# ----------------------------------------------------------------------------------------------------------------------------------------------
# Code to generated time-series cross validated posterior predictions for each candidate model
# - Data from Jan 2000 - Jan 2009 is used for training, then for each week in the remaining dataset, we fit an INLA model using
#   data up to but not including the target week, and generate a posterior predictive distribution for cases in the target week, t.
# - If run in an early-warning set up, we fit an INLA model up to target week t minus the horizon h, and then use lagged predictor
#   variables to generate posterior predictive distributions for dengue cases with an h week lead time.
# - We also calculate the posterior predictive probability of an exceeding a pre-determined outbreak threshold
# - We use a seasonal moving 75th percentile as an outbreak threshold (calculated using the calculate_thresholds fn)
# - This script is currently set up so it can be run on an HPC (called by a bash script)
# - If running locally three c_args must be provided:
#    1. Forecast horizon in weeks: can be 0, 2, 4, 6 or 8
#    2. Which model formulae to run: can be 'all', 'covariate-only', 'sero_clim_form', 'clim_form', 'sero_form', 'baseline' or 'baseline_year'
#    3. Whether to estimate yearly RE or not: can be 'estimated' or 'na'
#    4. Whether to include lagged cases as a predictor in the model: can be 'no' or 'yes'

# Author: Emilie Finch
# ---------------------------------------------------------------------------------------------------------------------------------------------

source(here("R", "tscv-prediction_fn.R"))
source(here("R", "utils_fn.R"))
dir.create(here("output", paste0(Sys.Date())))

# c_args = commandArgs(trailingOnly = TRUE);
c_args <- list(0, "all", "estimated", "no") # To run all models for forecast horizon of 0 weeks without lagged cases as a predictor

# Prepare data ---------------------------------------------------------------------------------------------------------------------------------

if (!exists("df_model")) {
  source(here("R", "00_load-data.R"))
  source(here("R", "create-lagged-data_fn.R"))
  df_model <- lag_data(dengue_singapore)
}

df_eval <- df_model  |> 
  group_by(year, month) |>
  mutate(month_index = cur_group_id()) |>
  ungroup() |>
  mutate(date_index = row_number()) |>
  dplyr::select(
    date, date_index, year, year_index, month, month_index,
    eweek, cases, pop, 
    time_since_switch,
    time_since_switch_2,
    time_since_switch_4,
    time_since_switch_6,
    time_since_switch_8,
    max_t_scale_12_wk_avg_0,
    max_t_scale_10_wk_avg_2,
    max_t_scale_8_wk_avg_4,
    max_t_scale_6_wk_avg_6,
    max_t_scale_4_wk_avg_8,
    nino34_12_wk_avg_4,
    nino34_10_wk_avg_6,
    nino34_8_wk_avg_8,
    days_no_rain_12_wk_total_0,
    days_no_rain_10_wk_total_2,
    days_no_rain_8_wk_total_4,
    days_no_rain_6_wk_total_6,
    days_no_rain_4_wk_total_8,
    lag_cases,
    lag_cases_2,
    lag_cases_4,
    lag_cases_6,
    lag_cases_8
  ) |>
  # As days without rain is a cumulative variable, scale up the lagged versions to what would be expected over a 12 week period
  mutate(
    days_no_rain_10_wk_total_2 = days_no_rain_10_wk_total_2 * 12 / 10,
    days_no_rain_8_wk_total_4 = days_no_rain_8_wk_total_4 * 12 / 8,
    days_no_rain_6_wk_total_6 = days_no_rain_6_wk_total_6 * 12 / 6,
    days_no_rain_4_wk_total_8 = days_no_rain_4_wk_total_8 * 12 / 4
  ) |>
  mutate(
    time_since_switch_2 = time_since_switch_2 + 2,
    time_since_switch_4 = time_since_switch_4 + 4,
    time_since_switch_6 = time_since_switch_6 + 6,
    time_since_switch_8 = time_since_switch_8 + 8
  )

print("Data loaded")

# Add outbreak thresholds

year_month <- df_eval |>
  group_by(year_index, month) |>
  filter(year_index >= 10) |>
  summarise(.groups = "keep")
thresholds <- purrr::map2_df(year_month$month, year_month$year_index, possibly(calculate_thresholds), data_input = df_eval, quantile = 0.75, .progress = TRUE)

df_eval <- df_eval |> left_join(thresholds, by = c("year_index", "month"))

# Here we define an outbreak week where the number of cases is > seasonal moving 75th percentile
# threshold, and an outbreak year as having more than 8 outbreak weeks
df_eval <- df_eval |>
  mutate(outbreak_week = case_when(cases > threshold ~ 1, TRUE ~ 0)) |>
  group_by(year) |>
  mutate(outbreak_year = case_when(sum(outbreak_week) > 12 ~ 1, TRUE ~ 0)) |>
  mutate(outbreak_year = case_when(year == 2004 | year == 2005 | year == 2007 ~ 1, TRUE ~ outbreak_year)) |> # Add outbreak years pre 2010
  ungroup()

# Parameters to run evaluation
horizon <- c_args[[1]]
form_input <- c_args[[2]]
yearly_re <- c_args[[3]]
include_lag_cases <- c_args[[4]]

print("Generating time series cross-validated predictions.")

# Run TSCV predictions ----------------------------------------------------------------------------------------------------------------------------

tscv_predictions_weekly(
  data_input = df_eval,
  horizon = as.numeric(horizon), # One of: 0, 2, 4, 6, 8
  form_input = as.character(form_input), # One of: "all", "covariate-only", "sero_clim_form", "clim_form", "sero_form", "baseline", "baseline-year"
  yearly_re = as.character(yearly_re), # Either: "estimated", "na"
  include_lag_cases = as.character(include_lag_cases), # Either: "yes", "no"
  filename = "tscv-preds-weekly"
)

print(paste0("TSCV predictions finished for horizon ", horizon, "."))

