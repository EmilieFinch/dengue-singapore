# -----------------------------------------------------------------------------------------------------------------------
# Code to conduct counterfactual analysis to estimate the potential impact of Project Wolbachia releases
# Author: Emilie Finch
# -----------------------------------------------------------------------------------------------------------------------

# Set up ----------------------------------------------------------------------------------------------------------------
if (!exists("dengue_wol")) {
  source(here("R", "00_load-data.R"))
}
source(here("R", "inla-counterfactual_fn.R"))

# Create output folder
dir.create(here("output", paste0(Sys.Date())))
output_folder <- here("output", paste0(Sys.Date()))

# Run counterfactual ----------------------------------------------------------------------------------------------------

# Define model formula

sero_climate <- "+f(inla.group(time_since_switch, n = 18), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(max_t_scale_12_wk_avg_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(nino34_12_wk_avg_4, n = 12), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(days_no_rain_12_wk_total_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
forms <- c(sero_climate)

# Run model 
## Fitting up to 2022-06-26 (start of Wolbachia RCT in Singapore)
## For sensitivity: 2020-05-10 (start of larger releases)
## And 2018-07-01 start of smaller releases

df_full <- lag_data(dengue_wol)
df_wol <- lag_data(dengue_wol) |> 
  filter(year <= 2023) |> 
  mutate(cases = case_when(date > as.Date("2022-06-26") ~ NA_integer_,
                           T ~ cases)) 

model_output <- run_scenario(df_wol, sero_climate)

wol_out <- cbind(df_full |> filter(year <= 2023) |> select(year,cases), model_output$exp_quantiles) |> 
  remove_rownames() |> 
  mutate(flag = case_when(date <= as.Date("2022-06-26") ~ "Tuning",
                          T ~ "Counterfactual"),
         training_end_date = "2022-06-26")

# Sensitivity analysis --------------------------------------------------------------------------------------------------------

sensitivity_table <- main_table_out
all_preds <- wol_out
training_end_dates <- c("2020-05-10", "2018-07-01")

for(end_date in training_end_dates){
  
  df_wol <- df_full |> 
    filter(year <= 2023) |> 
    mutate(cases = case_when(date > as.Date(end_date) ~ NA_integer_,
                             T ~ cases))  
  
  model_output <- run_scenario(df_wol, sero_climate)
  
  wol_out <- cbind(df_full |> filter(year <= 2023) |> select(year,cases), model_output$exp_quantiles) |> 
    remove_rownames() |> 
    mutate(flag = case_when(date <= as.Date(end_date) ~ "Tuning",
                            T ~ "Counterfactual")) |> 
    mutate(training_end_date = end_date)

  all_preds <- rbind(all_preds, wol_out)

}

saveRDS(all_preds, here(output_folder, paste0("counterfactual-predictions_", Sys.Date(), ".rds")))

