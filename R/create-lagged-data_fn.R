# -----------------------------------------------------------------------------------------------------------------------
# Function to create running averages and total for covariates, and add lags from 2 to 16 weeks
# Author: Emilie Finch
# Arguments:
# data (data frame)
# Output: data frame for model selection including lagged and aggregated climate and serotype variables
# - Variables created are: max_t (maximum temperature), mean_t (mean temperature), min_t (minimum temperature)
#   ab_hum (absolute humidity), prec (precipitation), days_no_rain (days without rain), nino34 (Niño 3.4 SSTA) and
#   time_since_switch_weekly (time since a switch in dominant serotype)
# - Running averages and totals are tagged in variable names with e.g. 4_wk_avg
# - Number of weeks the variable is lagged is tagged in the variable name with e.g. _4 for a four week lag
# - E.g. days_no_rain_4_wk_total_8 is a four week running total of days without rain, at an 8 week lag
# -----------------------------------------------------------------------------------------------------------------------

lag_data <- function(data) {
  
  df_model <- data |>
    rename(max_t = maximum_temperature, mean_t = mean_temperature, min_t = minimum_temperature, ab_hum = absolute_humidity, prec = rainfall, 
          pop = total_population, cases = dengue_cases) |>
    mutate(year_index = year - 1999, date_num = round(as.numeric(date)), max_t_scale = max_t - mean(max_t), mean_t_scale = mean_t - mean(mean_t), min_t_scale = min_t - mean(min_t),
           ab_hum_scale = ab_hum - mean(ab_hum), prec_scale = prec - mean(prec)) # Note need to round as otherwise can't match numeric date columns

  vars_avg <- c("max_t_scale", "min_t_scale", "mean_t_scale", "ab_hum_scale", "nino34")
  vars_sum <-c("prec_scale", "days_no_rain")
  
  df_model <- df_model |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_avg)), rollmean, k = 4, fill = NA, align = "right"),
                       c(paste0(c(vars_avg), "_4_wk_avg_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_avg)), rollmean, k = 6, fill = NA, align = "right"),
                       c(paste0(c(vars_avg), "_6_wk_avg_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_avg)), rollmean, k = 8, fill = NA, align = "right"),
                       c(paste0(c(vars_avg), "_8_wk_avg_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_avg)), rollmean, k = 10, fill = NA, align = "right"),
                       c(paste0(c(vars_avg), "_10_wk_avg_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_avg)), rollmean, k = 12, fill = NA, align = "right"),
                       c(paste0(c(vars_avg), "_12_wk_avg_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_sum)), rollsum, k = 4, fill = NA, align = "right"),
                       c(paste0(c(vars_sum), "_4_wk_total_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_sum)), rollsum, k = 6, fill = NA, align = "right"),
                       c(paste0(c(vars_sum), "_6_wk_total_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_sum)), rollsum, k = 8, fill = NA, align = "right"),
                       c(paste0(c(vars_sum), "_8_wk_total_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_sum)), rollsum, k = 10, fill = NA, align = "right"),
                       c(paste0(c(vars_sum), "_10_wk_total_0")))) |> 
    bind_cols(setNames(lapply(df_model  |> dplyr::select(all_of(vars_sum)), rollsum, k = 12, fill = NA, align = "right"),
                       c(paste0(c(vars_sum), "_12_wk_total_0"))))
  
  df_model <- df_model |> 
    bind_cols(setNames(shift(df_model$mean_t_scale,seq(2,16, by = 2)), c(paste0("mean_t_scale_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale,seq(2,16, by = 2)), c(paste0("min_t_scale_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale,seq(2,16, by = 2)), c(paste0("max_t_scale_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale,seq(2,16, by = 2)), c(paste0("prec_scale_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain,seq(2,16, by = 2)), c(paste0("days_no_rain_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale,seq(2,16, by = 2)), c(paste0("ab_hum_scale_", seq(2,16, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$nino34,seq(2,20, by = 2)), c(paste0("nino34_", seq(2,20, by = 2))))) |> 
    
    bind_cols(setNames(shift(df_model$mean_t_scale_4_wk_avg_0,seq(2,12, by = 2)), c(paste0("mean_t_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale_4_wk_avg_0,seq(2,12, by = 2)), c(paste0("min_t_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale_4_wk_avg_0,seq(2,12, by = 2)), c(paste0("max_t_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale_4_wk_total_0,seq(2,12, by = 2)), c(paste0("prec_scale_4_wk_total_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain_4_wk_total_0,seq(2,12, by = 2)), c(paste0("days_no_rain_4_wk_total_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale_4_wk_avg_0,seq(2,12, by = 2)), c(paste0("ab_hum_scale_4_wk_avg_", seq(2,12, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$nino34_4_wk_avg_0,seq(2,16, by = 2)), c(paste0("nino34_4_wk_avg_", seq(2,16, by = 2))))) |> 
    
    bind_cols(setNames(shift(df_model$mean_t_scale_6_wk_avg_0,seq(2,10, by = 2)), c(paste0("mean_t_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale_6_wk_avg_0,seq(2,10, by = 2)), c(paste0("min_t_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale_6_wk_avg_0,seq(2,10, by = 2)), c(paste0("max_t_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale_6_wk_total_0,seq(2,10, by = 2)), c(paste0("prec_scale_6_wk_total_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain_6_wk_total_0,seq(2,10, by = 2)), c(paste0("days_no_rain_6_wk_total_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale_6_wk_avg_0,seq(2,10, by = 2)), c(paste0("ab_hum_scale_6_wk_avg_", seq(2,10, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$nino34_6_wk_avg_0,seq(2,14, by = 2)), c(paste0("nino34_6_wk_avg_", seq(2,14, by = 2))))) |> 
    
    bind_cols(setNames(shift(df_model$mean_t_scale_8_wk_avg_0,seq(2,8, by = 2)), c(paste0("mean_t_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale_8_wk_avg_0,seq(2,8, by = 2)), c(paste0("min_t_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale_8_wk_avg_0,seq(2,8, by = 2)), c(paste0("max_t_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale_8_wk_total_0,seq(2,8, by = 2)), c(paste0("prec_scale_8_wk_total_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain_8_wk_total_0,seq(2,8, by = 2)), c(paste0("days_no_rain_8_wk_total_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale_8_wk_avg_0,seq(2,8, by = 2)), c(paste0("ab_hum_scale_8_wk_avg_", seq(2,8, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$nino34_8_wk_avg_0,seq(2,12, by = 2)), c(paste0("nino34_8_wk_avg_", seq(2,12, by = 2))))) |> 
    
    bind_cols(setNames(shift(df_model$mean_t_scale_10_wk_avg_0,seq(2,6, by = 2)), c(paste0("mean_t_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale_10_wk_avg_0,seq(2,6, by = 2)), c(paste0("min_t_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale_10_wk_avg_0,seq(2,6, by = 2)), c(paste0("max_t_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale_10_wk_total_0,seq(2,6, by = 2)), c(paste0("prec_scale_10_wk_total_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain_10_wk_total_0,seq(2,6, by = 2)), c(paste0("days_no_rain_10_wk_total_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale_10_wk_avg_0,seq(2,6, by = 2)), c(paste0("ab_hum_scale_10_wk_avg_", seq(2,6, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$nino34_10_wk_avg_0,seq(2,10, by = 2)), c(paste0("nino34_10_wk_avg_", seq(2,10, by = 2))))) |> 
    
    bind_cols(setNames(shift(df_model$mean_t_scale_12_wk_avg_0,seq(2,4, by = 2)), c(paste0("mean_t_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$min_t_scale_12_wk_avg_0,seq(2,4, by = 2)), c(paste0("min_t_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$max_t_scale_12_wk_avg_0,seq(2,4, by = 2)), c(paste0("max_t_scale_12_wk_avg_", seq(2,4, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$prec_scale_12_wk_total_0,seq(2,4, by = 2)), c(paste0("prec_scale_12_wk_total_", seq(2,4, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$days_no_rain_12_wk_total_0,seq(2,4, by = 2)), c(paste0("days_no_rain_12_wk_total_", seq(2,4, by = 2))))) |> 
    bind_cols(setNames(shift(df_model$ab_hum_scale_12_wk_avg_0,seq(2,4, by = 2)), c(paste0("ab_hum_scale_12_wk_avg_", seq(2,4, by = 2)))))  |> 
    bind_cols(setNames(shift(df_model$nino34_12_wk_avg_0,seq(2,8, by = 2)), c(paste0("nino34_12_wk_avg_", seq(2,8, by = 2)))))  
  
  # Add lagged version of time vars of interest

  df_model <- df_model |>
    do(data.frame(., setNames(shift(.$time_since_switch, c(2, 4, 6, 8, 12, 16)), c(paste0("time_since_switch_", c(2, 4, 6, 8, 12, 16))))))

  return(df_model)
}
