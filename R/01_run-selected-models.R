# -----------------------------------------------------------------------------------------------------------------------
# Code to fit INLA models to NEA dengue data for Singapore from 2000 to 2022
# Script is set up to run a: sero_climate model (with best set of covariates from model selection); a climate-only model
# and a serotype only model
# Author: Emilie Finch
# -----------------------------------------------------------------------------------------------------------------------

# Load data and fns for analysis

if (!exists("dengue_singapore")) {
  source(here("code", "00_read-data.R"))
}
source(here("R", "create-lagged-data_fn.R"))
source(here("R", "fit-inla_fn.R"))

# Wrangle data for model fitting

df_model <- lag_data(dengue_singapore) 

# Define formulae

sero_climate <- "+f(inla.group(time_since_switch, n = 18), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(max_t_scale_12_wk_avg_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(nino34_12_wk_avg_4, n = 12), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(days_no_rain_12_wk_total_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
sero_only <- "+f(inla.group(time_since_switch, n = 18), model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
climate_only <- "+f(inla.group(max_t_scale_12_wk_avg_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(nino34_12_wk_avg_4, n = 12), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(days_no_rain_12_wk_total_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"

forms <- c(sero_climate, climate_only, sero_only)

# Run models

mods_out <- fit_inla(df_model, forms)
saveRDS(mods_out, here("output", paste0("model-output_", Sys.Date(), ".rds")))
