# ============================================================
# Fay-Herriot Small Area Estimation Example:
# Life Expectancy at Birth Across European NUTS-2 Regions
#
# Outcome:
#   Life expectancy at birth, in years
#
# Area-level covariates:
#   - Gross domestic product
#   - Unemployment
#   - Education
#   - Population density
#
# Data source:
#   Eurostat regional databases
#
# Model:
#   Area-level Fay-Herriot model fitted using sae::eblupFH()
#
# Note:
#   Because sampling variances are not directly available in the
#   downloaded Eurostat table, a constant standard error of 0.20
#   years is used for illustration.
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(eurostat)
library(dplyr)
library(stringr)
library(tidyr)
library(sae)
library(sf)
library(ggplot2)
library(giscoR)


# ============================================================
# 2. Define helper functions
# ============================================================

# Identify NUTS-2 geographic codes.
is_nuts2 <- function(geo) {
  nchar(geo) == 4 &
    !str_detect(geo, "^(EU|EA)")
}

# Select the first preferred category available in a variable.
# If none of the preferred categories are available, select the
# most frequently occurring category.
pick_level <- function(x, preferred) {
  
  available_levels <- unique(x[!is.na(x)])
  
  for (level in preferred) {
    if (level %in% available_levels) {
      return(level)
    }
  }
  
  names(
    sort(
      table(x),
      decreasing = TRUE
    )
  )[1]
}


# ============================================================
# 3. Prepare the direct estimates
# ============================================================

# Download regional life expectancy data.
life_expectancy_raw <- get_eurostat(
  id = "demo_r_mlifexp",
  time_format = "num"
) %>%
  filter(is_nuts2(geo))

# Use the most recent year available in the life expectancy table.
year_use <- max(
  life_expectancy_raw$TIME_PERIOD,
  na.rm = TRUE
)

# Select dimension categories corresponding to:
#   - life expectancy at birth,
#   - total population,
#   - years.
#
# The helper function allows for minor differences in category
# naming across Eurostat table versions.
age_level <- pick_level(
  life_expectancy_raw$age,
  preferred = c("Y_LT0", "TOTAL", "Y0")
)

sex_level <- pick_level(
  life_expectancy_raw$sex,
  preferred = c("T", "TOTAL")
)

unit_level <- pick_level(
  life_expectancy_raw$unit,
  preferred = c("YR", "YEAR", "YEARS")
)

# Create the area-level direct-estimate dataset.
life_expectancy <- life_expectancy_raw %>%
  filter(
    TIME_PERIOD == year_use,
    age == age_level,
    sex == sex_level,
    unit == unit_level
  ) %>%
  transmute(
    geo,
    year = TIME_PERIOD,
    direct_estimate = values
  )


# ============================================================
# 4. Prepare the area-level covariates
# ============================================================

# ------------------------------------------------------------
# 4.1 Gross domestic product
# ------------------------------------------------------------

gdp_raw <- get_eurostat(
  id = "nama_10r_2gdp",
  time_format = "num"
) %>%
  filter(
    is_nuts2(geo),
    TIME_PERIOD == year_use
  )

gdp_unit <- pick_level(
  gdp_raw$unit,
  preferred = c(
    "EUR_HAB",
    "PPS_HAB",
    "MIO_EUR",
    "EUR"
  )
)

gdp <- gdp_raw %>%
  filter(unit == gdp_unit) %>%
  group_by(geo) %>%
  slice_max(
    order_by = values,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    geo,
    gdp = values
  )


# ------------------------------------------------------------
# 4.2 Unemployment
# ------------------------------------------------------------

unemployment_raw <- get_eurostat(
  id = "lfst_r_lfu3rt",
  time_format = "num"
) %>%
  filter(
    is_nuts2(geo),
    TIME_PERIOD == year_use
  )

unemployment <- unemployment_raw %>%
  group_by(geo) %>%
  slice_max(
    order_by = values,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    geo,
    unemployment = values
  )


# ------------------------------------------------------------
# 4.3 Education
# ------------------------------------------------------------

education_raw <- get_eurostat(
  id = "edat_lfse_04",
  time_format = "num"
) %>%
  filter(
    is_nuts2(geo),
    TIME_PERIOD == year_use
  )

education <- education_raw %>%
  group_by(geo) %>%
  slice_max(
    order_by = values,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    geo,
    education = values
  )


# ------------------------------------------------------------
# 4.4 Population density
# ------------------------------------------------------------

density_raw <- get_eurostat(
  id = "demo_r_d3dens",
  time_format = "num"
) %>%
  filter(
    is_nuts2(geo),
    TIME_PERIOD == year_use
  )

population_density <- density_raw %>%
  group_by(geo) %>%
  slice_max(
    order_by = values,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  transmute(
    geo,
    population_density = values
  )


# ============================================================
# 5. Construct the Fay-Herriot analysis dataset
# ============================================================

fh_data <- life_expectancy %>%
  left_join(
    gdp,
    by = "geo"
  ) %>%
  left_join(
    unemployment,
    by = "geo"
  ) %>%
  left_join(
    education,
    by = "geo"
  ) %>%
  left_join(
    population_density,
    by = "geo"
  ) %>%
  drop_na(
    direct_estimate,
    gdp,
    unemployment,
    education,
    population_density
  )


# ============================================================
# 6. Specify illustrative sampling variances
# ============================================================

# Working assumption:
# Each direct estimate has a standard error of 0.20 years.
direct_se <- 0.20

fh_data <- fh_data %>%
  mutate(
    sampling_variance = direct_se^2
  )

# sae::eblupFH() expects a standard data.frame.
fh_data <- as.data.frame(fh_data)


# ============================================================
# 7. Fit the Fay-Herriot model
# ============================================================

fit_fh <- eblupFH(
  formula =
    direct_estimate ~
    gdp +
    unemployment +
    education +
    population_density,
  vardir = sampling_variance,
  method = "REML",
  data = fh_data
)

# Store the empirical best linear unbiased predictors.
fh_data$fh_estimate <- as.numeric(
  fit_fh$eblup
)

# Estimated area-level random-effect variance.
random_effect_variance <- fit_fh$fit$refvar


# ============================================================
# 8. Display model information and estimates
# ============================================================

cat(
  "Outcome: Life expectancy at birth\n"
)

cat(
  "Year:", year_use, "\n"
)

cat(
  "Selected age category:", age_level, "\n"
)

cat(
  "Selected sex category:", sex_level, "\n"
)

cat(
  "Selected outcome unit:", unit_level, "\n"
)

cat(
  "Selected GDP unit:", gdp_unit, "\n"
)

cat(
  "Number of NUTS-2 regions:",
  nrow(fh_data),
  "\n"
)

cat(
  "Estimated random-effect variance:",
  random_effect_variance,
  "\n"
)

fh_results <- fh_data %>%
  transmute(
    geo,
    direct_estimate,
    fh_estimate
  ) %>%
  arrange(
    desc(fh_estimate)
  )

print(
  head(fh_results, 10)
)


# ============================================================
# 9. Download NUTS-2 geographic boundaries
# ============================================================

nuts2_boundaries <- gisco_get_nuts(
  year = "2021",
  epsg = "4326",
  nuts_level = "2",
  resolution = "20"
)


# ============================================================
# 10. Join the estimates to the geographic boundaries
# ============================================================

map_data <- nuts2_boundaries %>%
  left_join(
    fh_data,
    by = c("NUTS_ID" = "geo")
  ) %>%
  filter(
    !str_detect(NUTS_ID, "^FRY")
  )

# Note:
# French overseas NUTS regions beginning with "FRY" are excluded
# to improve the display of the mainland European map.


# ============================================================
# 11. Prepare the final plotting dataset
# ============================================================

plot_data <- map_data %>%
  dplyr::select(
    geo = NUTS_ID,
    geometry,
    fh_estimate,
    direct_estimate
  ) %>%
  pivot_longer(
    cols = c(
      fh_estimate,
      direct_estimate
    ),
    names_to = "estimate_type",
    values_to = "estimate"
  ) %>%
  mutate(
    estimate_type = factor(
      estimate_type,
      levels = c(
        "fh_estimate",
        "direct_estimate"
      ),
      labels = c(
        "(A) Fay-Herriot SAE",
        "(B) Direct estimate"
      )
    )
  )


# ============================================================
# 12. Create the comparison figure
# ============================================================

figure_s1 <- ggplot(plot_data) +
  geom_sf(
    aes(fill = estimate),
    color = "white",
    linewidth = 0.1
  ) +
  facet_wrap(
    ~ estimate_type,
    nrow = 1
  ) +
  scale_fill_gradient(
    name = "Life expectancy (years)",
    low = "#deebf7",
    high = "#08519c"
  ) +
  theme_minimal(
    base_size = 11
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

print(figure_s1)


# ============================================================
# 13. Optional: save the final figure
# ============================================================

# ggsave(
#   filename = "Figure_1_Eurostat_Fay_Herriot.png",
#   plot = figure_s1,
#   width = 12,
#   height = 6,
#   units = "in",
#   dpi = 600
# )


# ============================================================
# End of script
# ============================================================
