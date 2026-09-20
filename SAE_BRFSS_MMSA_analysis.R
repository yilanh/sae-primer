# ============================================================
# BRFSS MMSA Small Area Estimation Example:
# Dental Care Utilization Across California MMSAs
#
# Outcome:
#   Dental visit within the past year
#
# Covariates:
#   - Sex
#   - Education (> high school vs. high school or less)
#   - Race/ethnicity (Hispanic vs. non-Hispanic)
#
# Data source:
#   Behavioral Risk Factor Surveillance System (BRFSS)
#
# Estimation methods:
#   - Design-based direct estimation
#   - Area-level Fay-Herriot model
#   - Unit-level generalized linear mixed model
#
# Purpose:
#   To compare direct estimates with model-based small area
#   estimates obtained from area-level and unit-level models.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(haven)
library(dplyr)
library(lme4)
library(MASS)
library(sae)


# ------------------------------------------------------------
# 2. Set working directory and import data
# ------------------------------------------------------------

data <- read_sas("mmsa2018.sas7bdat") %>%
  filter(grepl("CA", MMSANAME))


# ------------------------------------------------------------
# 3. Inspect MMSA variables and outcome coding
# ------------------------------------------------------------

xtabs(~ MMSANAME, data = data)
xtabs(~ `_MMSA`, data = data)
xtabs(~ `_DENVST3`, data = data)


# ------------------------------------------------------------
# 4. Prepare dental-visit outcome
# ------------------------------------------------------------

# _DENVST3:
#   1 = Yes
#   2 = No
# Recode to:
#   1 = Yes
#   0 = No

# Keep only valid Yes/No responses

data2 <- data %>%
  filter(`_DENVST3` %in% c(1, 2)) %>%
  mutate(
    DENVST = if_else(`_DENVST3` == 1, 1, 0)
  )


# ------------------------------------------------------------
# 5. Define MMSA domains and population sizes
# ------------------------------------------------------------

# Confirm that these MMSA codes and population sizes correspond
# to the intended five California MMSAs.

domain_info <- data.frame(
  MMSA = c(31080, 36084, 40140, 40900, 41940),
  dom = 1:5,
  popn = c(
    10476831,
    2207584,
    3423128,
    1808544,
    1556487
  )
)


# ------------------------------------------------------------
# 6. Simple weighted direct estimates by MMSA
# ------------------------------------------------------------

direct_summary <- data2 %>%
  mutate(MMSA = as.numeric(`_MMSA`)) %>%
  inner_join(domain_info %>% dplyr::select(MMSA, dom), by = "MMSA") %>%
  group_by(dom, MMSA) %>%
  summarize(
    direct_estimate = weighted.mean(
      DENVST,
      w = as.numeric(`_MMSAWT`),
      na.rm = TRUE
    ),
    sample_size = n(),
    .groups = "drop"
  )

print(direct_summary)


# ------------------------------------------------------------
# 7. Prepare data for sae::direct()
# ------------------------------------------------------------

sample_data <- data2 %>%
  transmute(
    y = as.numeric(DENVST),
    weight = as.numeric(`_MMSAWT`),
    MMSA = as.numeric(`_MMSA`)
  ) %>%
  inner_join(
    domain_info %>% dplyr::select(MMSA, dom),
    by = "MMSA"
  ) %>%
  filter(
    !is.na(y),
    !is.na(weight),
    weight > 0,
    !is.na(dom)
  ) %>%
  as.data.frame()

popnsize <- domain_info %>%
  dplyr::select(dom, popn) %>%
  as.data.frame()

# Check domain consistency
print(sort(unique(sample_data$dom)))
print(sort(unique(popnsize$dom)))
print(setdiff(unique(sample_data$dom), popnsize$dom))


# ------------------------------------------------------------
# 8. Design-based direct estimates and variances
# ------------------------------------------------------------

fit_direct <- sae::direct(
  y = y,
  dom = dom,
  sweight = weight,
  domsize = popnsize,
  replace = FALSE,
  data = sample_data
)

print(fit_direct)


# ------------------------------------------------------------
# 9. Area-level auxiliary covariates
# ------------------------------------------------------------

# Proportion male

gender <- c(
  5026532 / (5026532 + 5277993),
  631069  / (631069  + 668118),
  1643061 / (1643061 + 1688978),
  694527  / (694527  + 747093),
  748873  / (748873  + 740643)
)

# Proportion with more than high school education

education <- c(
  (2341617 + 3887387) / 10304525,
  (269308  + 634730)  / 1299187,
  (895462  + 887064)  / 3332039,
  (391315  + 563189)  / 1441620,
  (269073  + 813426)  / 1489516
)

# Proportion Hispanic or Latino

Latino <- c(
  5973798 / 13262234,
  369061  / 1643700,
  2282330 / 4518699,
  399523  / 1890100,
  495455  / 1922200
)


# ------------------------------------------------------------
# 10. Create Fay-Herriot area-level dataset
# ------------------------------------------------------------

fh_data <- fit_direct %>%
  as.data.frame() %>%
  mutate(
    var = SD^2,
    gender = gender,
    education = education,
    Latino = Latino
  )

print(fh_data)


# ------------------------------------------------------------
# 11. Fit Fay-Herriot model
# ------------------------------------------------------------

fit_FH <- mseFH(
  formula = Direct ~ gender + education + Latino,
  vardir = var,
  data = fh_data
)

FH_CV <- 100 * sqrt(fit_FH$mse) / fit_FH$est$eblup

fh_results <- fh_data %>%
  mutate(
    FH_estimate = fit_FH$est$eblup,
    FH_MSE = fit_FH$mse,
    FH_SE = sqrt(FH_MSE),
    FH_CV = FH_CV
  )

print(fh_results)


# ------------------------------------------------------------
# 12. Prepare unit-level data for logistic mixed model
# ------------------------------------------------------------

# SEX1:
#   1 = Male
#   2 = Female
#
# _EDUCAG:
#   1 = Did not graduate high school
#   2 = Graduated high school
#   3 = Attended college/technical school
#   4 = Graduated from college/technical school
#   9 = Don’t know/Not sure/ Missing
#
# Define education as more than high school:
#   0 = _EDUCAG 1 or 2
#   1 = _EDUCAG 3 or 4
#
# _RACE:
#   8 = Hispanic
#   9 = Don't know / refused category excluded here


data3 <- data2 %>%
  transmute(
    DENVST = as.numeric(DENVST),
    SEX1 = as.numeric(SEX1),
    EDUCAG = as.numeric(`_EDUCAG`),
    RACE = as.numeric(`_RACE`),
    MMSA = as.numeric(`_MMSA`),
    weight = as.numeric(`_MMSAWT`)
  ) %>%
  filter(RACE != 9) %>%
  mutate(
    gender = case_when(
      SEX1 == 1 ~ 1,
      SEX1 == 2 ~ 0,
      TRUE ~ NA_real_
    ),
    education = case_when(
      EDUCAG %in% c(1, 2) ~ 0,
      EDUCAG %in% c(3, 4) ~ 1,
      TRUE ~ NA_real_
    ),
    Latino = case_when(
      RACE == 8 ~ 1,
      !is.na(RACE) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  inner_join(
    domain_info %>% dplyr::select(MMSA, dom),
    by = "MMSA"
  ) %>%
  filter(
    !is.na(DENVST),
    !is.na(gender),
    !is.na(education),
    !is.na(Latino),
    !is.na(weight),
    weight > 0,
    !is.na(dom)
  ) %>%
  mutate(
    dom = factor(dom, levels = 1:5)
  ) %>%
  as.data.frame()

# Check sample size and weighted outcome by domain

glmm_sample_summary <- data3 %>%
  group_by(dom) %>%
  summarize(
    n = n(),
    dental_visit_percent = weighted.mean(
      DENVST,
      w = weight,
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  )

print(glmm_sample_summary)


# ------------------------------------------------------------
# 13. Prepare area-level covariates for GLMM prediction
# ------------------------------------------------------------

pred_dat <- fh_data %>%
  transmute(
    dom = as.integer(Domain),
    gender = gender,
    education = education,
    Latino = Latino
  ) %>%
  arrange(dom) %>%
  mutate(
    dom = factor(dom, levels = levels(data3$dom))
  ) %>%
  as.data.frame()

print(pred_dat)


# ------------------------------------------------------------
# 14. Fit weighted logistic mixed model
# ------------------------------------------------------------

fit_glmm <- glmer(
  DENVST ~ gender + education + Latino + (1 | dom),
  data = data3,
  family = binomial(link = "logit"),
  weights = weight,
  control = glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

print(summary(fit_glmm))

# Check convergence and singular fit
print(fit_glmm@optinfo$conv$lme4$messages)
print(isSingular(fit_glmm, tol = 1e-4))

# Note:
# Although lme4 reported a "Model is nearly unidentifiable" warning,
# the model converged successfully (convergence code = 0) and was not singular,
# so no further action was taken.


# ------------------------------------------------------------
# 15. Obtain GLMM area estimates
# ------------------------------------------------------------

pred_dat$eta <- predict(
  fit_glmm,
  newdata = pred_dat,
  re.form = NULL
)

pred_dat$glmm_estimate <- plogis(pred_dat$eta)

print(pred_dat[, c("dom", "glmm_estimate")])


# ------------------------------------------------------------
# 16. Approximate GLMM uncertainty by simulation
# ------------------------------------------------------------

# This simulation incorporates fixed-effect uncertainty while
# treating the estimated random effects as fixed plug-in values.
# It is therefore not a complete MSE estimator for unit-level SAE.

beta_hat <- fixef(fit_glmm)
V_beta <- as.matrix(vcov(fit_glmm))

u_hat <- ranef(fit_glmm)$dom[, 1]
names(u_hat) <- rownames(ranef(fit_glmm)$dom)

X_pred <- model.matrix(
  ~ gender + education + Latino,
  data = pred_dat
)

B <- 1000
set.seed(123)

sim_est <- matrix(
  NA_real_,
  nrow = nrow(pred_dat),
  ncol = B
)

for (b in seq_len(B)) {

  beta_b <- MASS::mvrnorm(
    n = 1,
    mu = beta_hat,
    Sigma = V_beta
  )

  u_b <- u_hat[as.character(pred_dat$dom)]

  eta_b <- as.vector(
    X_pred %*% beta_b + u_b
  )

  sim_est[, b] <- plogis(eta_b)
}

pred_dat <- pred_dat %>%
  mutate(
    glmm_mean = rowMeans(sim_est, na.rm = TRUE),
    glmm_se = apply(sim_est, 1, sd, na.rm = TRUE),
    glmm_cv = 100 * glmm_se / glmm_mean
  )

print(
  pred_dat[, c(
    "dom",
    "glmm_estimate",
    "glmm_mean",
    "glmm_se",
    "glmm_cv"
  )]
)


# ------------------------------------------------------------
# 17. Combine direct, FH, and GLMM estimates
# ------------------------------------------------------------

plot_data <- fit_direct %>%
  as.data.frame() %>%
  mutate(
    FH_estimate = fit_FH$est$eblup,
    FH_CV = 100 * sqrt(fit_FH$mse) / fit_FH$est$eblup
  ) %>%
  left_join(
    pred_dat %>%
      transmute(
        Domain = as.numeric(as.character(dom)),
        GLMM_estimate = glmm_mean,
        GLMM_CV = glmm_cv
      ),
    by = "Domain"
  )

print(plot_data)


# ------------------------------------------------------------
# 18. Prepare plotting variables
# ------------------------------------------------------------

area <- plot_data$Domain

direct_estimate <- plot_data$Direct * 100
fh_estimate <- plot_data$FH_estimate * 100
glmm_estimate <- plot_data$GLMM_estimate * 100

direct_cv <- plot_data$CV
fh_cv <- plot_data$FH_CV
glmm_cv <- plot_data$GLMM_CV


# ------------------------------------------------------------
# 19. Plot estimates and coefficients of variation
# ------------------------------------------------------------

# png(
#   "Figure_2_BRFSS_MMSA.png",
#   width = 1200,
#   height = 800,
#   res = 150
# )

par(mfrow = c(1, 2))

# Panel A: Estimates

estimate_range <- range(
  c(
    direct_estimate,
    fh_estimate,
    glmm_estimate
  ),
  na.rm = TRUE
)

plot(
  area,
  direct_estimate,
  type = "p",
  pch = 1,
  ylim = estimate_range + c(-5, 5),
  xlab = "MMSA",
  ylab = "Estimate (%)",
  main = "",
  col = "black"
)

points(
  area,
  fh_estimate,
  col = "blue",
  pch = 4
)

points(
  area,
  glmm_estimate,
  col = "green3",
  pch = 0
)

legend(
  "topleft",
  legend = c(
    "Direct",
    "Fay-Herriot",
    "Unit-level logistic"
  ),
  col = c(
    "black",
    "blue",
    "green3"
  ),
  pch = c(1, 4, 0),
  bty = "n"
)

mtext(
  "(a)",
  side = 3,
  line = 1,
  adj = 0,
  font = 2
)


# Panel B: Coefficients of variation

cv_range <- range(
  c(
    direct_cv,
    fh_cv,
    glmm_cv
  ),
  na.rm = TRUE
)

plot(
  area,
  direct_cv,
  type = "p",
  pch = 1,
  ylim = cv_range + c(-0.8, 0.8),
  xlab = "MMSA",
  ylab = "CV (%)",
  main = "",
  col = "black"
)

points(
  area,
  fh_cv,
  col = "blue",
  pch = 4
)

points(
  area,
  glmm_cv,
  col = "green3",
  pch = 0
)

legend(
  "topleft",
  legend = c(
    "Direct",
    "Fay-Herriot",
    "Unit-level logistic"
  ),
  col = c(
    "black",
    "blue",
    "green3"
  ),
  pch = c(1, 4, 0),
  bty = "n"
)

mtext(
  "(b)",
  side = 3,
  line = 1,
  adj = 0,
  font = 2
)

# dev.off()


# ============================================================
# End of script
# ============================================================
