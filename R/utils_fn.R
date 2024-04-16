
# Helper functions to wrangle output ---------------------------------------------------------------------------------------------------------------------

# Extract posterior predictions from tscv output

get_pred_quantiles <- function(preds, lag, data_input) {
  pred_quantiles <- NULL
  for (p in 1:length(preds)) {
    add <- data.frame(
      median = matrixStats::rowQuantiles(preds[[p]]$post.samples, cols = c(1:1000), probs = 0.5),
      lower = matrixStats::rowQuantiles(preds[[p]]$post.samples, cols = c(1:1000), probs = 0.025),
      upper = matrixStats::rowQuantiles(preds[[p]]$post.samples, cols = c(1:1000), probs = 0.975),
      mod = preds[[p]]$mod,
      horizon = lag
    )
    add <- add |>
      cbind(date = data_input$date) |>
      cbind(true_cases = data_input$cases)
    pred_quantiles <- pred_quantiles |> bind_rows(add)
  }
  return(pred_quantiles)
}

# Extract outbreak probabilities from tscv output

get_outbreak_preds <- function(preds, data_input) {
  outbreak_preds <- NULL
  for (p in 1:length(preds)) {
    preds[[p]]$roc <- pROC::roc(preds[[p]]$outbreak.probs$outbreak,
                                preds[[p]]$outbreak.probs$prob,
                                direction = "<", ci = TRUE, quiet = TRUE
    )
    preds[[p]]$trigger <- pROC::coords(preds[[p]]$roc,
                                       x = "best",
                                       best.method = "closest.topleft"
    )
    
    add <- preds[[p]]$outbreak.probs |>
      rename(outbreak.occurred = outbreak) |>
      mutate(outbreak.predicted = case_when(
        prob >= preds[[p]]$trigger[[1]] ~ 1,
        TRUE ~ 0
      )) |>
      mutate(
        trigger = preds[[p]]$trigger[[1]],
        mod = preds[[p]]$mod
      )
    outbreak_preds <- outbreak_preds |> bind_rows(add)
  }
  outbreak_preds <- cbind(outbreak_preds, data_input[, !names(data_input) %in% c("cases")])
  
  return(outbreak_preds)
}

# Function to generate ROC curve coordinates

get_roc_coords <- function(preds) {
  thresh_out <- data.frame()
  for (p in unique(preds$mod)) {
    print(p)
    temp_outbreak <- preds |> filter(mod == p)
    roc_out <- roc(temp_outbreak$outbreak.occurred,
                   temp_outbreak$prob,
                   direction = "<", quiet = TRUE, ci = TRUE
    )
    roc_ci <- pROC::ci.se(roc_out, specificities = seq(0, 1, l = 50))
    dat_ci <- data.frame(
      specificity = as.numeric(rownames(roc_ci)),
      sensitivity_lower = roc_ci[, 1],
      sensitivity = roc_ci[, 2],
      sensitivity_upper = roc_ci[, 3]
    )
    rownames(dat_ci) <- NULL
    
    best_coords <- pROC::coords(roc_out,
                                x = "best", best.method = "youden", ci = TRUE
    )
    best_coords$best_flag <- 1
    temp_out <- bind_rows(dat_ci, best_coords |> dplyr::select(-threshold))
    temp_out$mod <- p
    temp_out$auc <- roc_out$auc
    temp_out$auc_lower <- min(ci.auc(roc_out))
    temp_out$auc_upper <- max(ci.auc(roc_out))
    thresh_out <- rbind(thresh_out, temp_out)
  }
  return(thresh_out)
}

# Function to extract hit rates and false alarm rates

get_hit_rates <- function(preds) {
  pred_rates <- preds |>
    group_by(mod, outbreak.occurred, trigger) |>
    mutate(total = n()) |>
    ungroup() |>
    group_by(mod, outbreak.occurred, outbreak.predicted, total, trigger) %>%
    summarise(events = n()) |>
    filter(outbreak.predicted == 1) |>
    mutate(rate = events / total) |>
    mutate(rate_type = case_when(
      outbreak.occurred == 1 & outbreak.predicted == 1 ~ "hit_rate",
      TRUE ~ "false_alarm"
    )) |>
    ungroup()
  return(pred_rates)
}

# Function to calculate CRPS

get_crps <- function(preds, data_input) {
  scores <- data.frame()
  for (i in 1:length(preds)) {
    preds_df <- as.data.frame(cbind(preds[[i]]$post.samples, true_value = data_input$cases, target_end_date = data_input$date))
    temp <- preds_df |>
      mutate(target_end_date = as.Date(target_end_date, origin = "1970-01-01")) |>
      pivot_longer(-c(true_value, target_end_date), names_to = "sample", values_to = "prediction") |>
      mutate(sample = as.numeric(substring(sample, 2))) |>
      transform_forecasts(append = TRUE, fun = log_shift, offset = 1) |>
      filter(!is.na(prediction)) |>
      score() |>
      summarise_scores(by = c("scale"))
    
    temp$mod <- preds[[i]]$model
    
    scores <- rbind(scores, temp)
  }
  return(scores)
}

# Function to calculate brier score

get_brier <- function(preds, data_input) {
  scores <- data.frame()
  for (i in 1:length(preds)) {
    preds_df <- as.data.frame(cbind(prediction = preds[[i]]$outbreak.probs$prob, true_value = preds[[i]]$outbreak.probs$outbreak, target_end_date = data_input$date))
    temp <- preds_df |>
      mutate(target_end_date = as.Date(target_end_date, origin = "1970-01-01")) |>
      filter(!is.na(prediction)) |>
      score(metrics = "brier_score") |>
      summarise_scores(by = "model") |>
      select(-model)
    
    temp$mod <- preds[[i]]$model
    
    scores <- rbind(scores, temp)
  }
  return(scores)
}

# Function to read latest model results

read_results <- function(date = NULL, type = "main-output", horizon = NULL){
  files <- list.files(here("output"))
  files <- files[str_detect(files, type)]
  dates <- as.Date(sub("^.*(202[2-9]-[0-9]*-[0-9]*).*$", "\\1", files))
  
  if (is.null(date)) {
    date <- max(dates, na.rm = TRUE)
  } 
  
  if(type == "tscv-preds"){
    date_files <- files[dates == date]
    horizons <- sub("^.*(horizon-[0-8]).*$", "\\1", date_files)
    horizons <- as.numeric(sub("^.*([0-9]+).*$", "\\1", horizons))
    final_file <- date_files[horizons == horizon]
    
  } else { 
    final_file <- max(files[dates == date])
  }
  results <- readRDS(here("output", final_file))
  print(paste0("Using ",final_file, " to generate results."))
  return(results)
}

# Wrapper function to wrangle scoring output

score_tscv <- function(tscv_output, data, horizon){
  
  tscv_crps <- get_crps(tscv_output, data)
  tscv_brier <- get_brier(tscv_output, data)
  tscv_outbreaks <- get_outbreak_preds(tscv_output, data)
  hit_rates <- get_hit_rates(tscv_outbreaks)
  roc_curves <- get_roc_coords(tscv_outbreaks)
  
  
  score_table <- tscv_crps |>
    select(crps, scale, mod) |>
    pivot_wider(names_from = scale, names_prefix = "crps_", values_from = crps) |>
    left_join(
      roc_curves |>
        select(mod, auc, auc_lower, auc_upper) |>
        unique()
    ) |>
    left_join(hit_rates |>
                select(mod, trigger, rate, rate_type) |>
                pivot_wider(names_from = rate_type, values_from = rate)) |>
    left_join(tscv_brier) |> 
    mutate(horizon = horizon)
  
  return(score_table)
}
  
  