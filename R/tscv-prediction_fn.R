# -----------------------------------------------------------------------------------------------------------------------
# Function to generate time series cross-validated predictions
# And helper function to calculate seasonal moving percentile threshold
# Arguments:
# data_input (data frame)
# horizon (numeric): number of weeks into the future to predict for each formula
# form_input (character): defines formula input, one of: "all", "covariate-only", "sero_clim_form", "sero_form", "clim_form", "baseline", "baseline_year"
# yearly_re (character): whether yearly RE to NA or estimated, one of: "na", "estimated"
# filename (character): output file name
# Outputs:
# List with:
# post.samples - posterior predictive samples for each week
# outbreak.probs - model outbreak probabilities, with cases and outbreak threshold for each week
# model - model formula name
# yearly_re - whether yearly RE estimated or not
# Author: Emilie Finch, adapted from code by Martin Lotto Batista and Rachel Lowe
# -----------------------------------------------------------------------------------------------------------------------

# Function to calculate outbreak thresholds for a given quantile

calculate_thresholds <- function(data_input, month_input, year_input, quantile) {
  # Using only data prior to the year/month of interest
  threshold <- data_input |>
    filter(year_index < as.numeric(year_input), month == month_input) |>
    mutate(threshold = quantile(cases, quantile)) |>
    pull(threshold)

  threshold_out <- data.frame(year_index = year_input, month = month_input, threshold = unique(threshold))

  return(threshold_out)
}

# Function to conduct expanding window TSCV at a weekly level -----------------------------------------------------------------------------------------------------------------------------------------------------------

tscv_predictions_weekly <- function(data_input, # Data frame
                                    horizon, # Numeric: number of weeks into the future to predict for each formula
                                    form_input, # Character string of formula input, one of: "all", "covariate-only", "sero_clim_form", "sero_form", "clim_form", "baseline", "baseline_year"
                                    yearly_re, # Yearly RE to NA or estimated, one of: "na", "estimated"
                                    filename) { # Character string of output file name
  start_time <- Sys.time()

  # Set up formulae for fitting and prediction
  sero_clim_form <- "+f(inla.group(time_since_switch, n = 18), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(max_t_scale_12_wk_avg_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(nino34_12_wk_avg_4, n = 12), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(days_no_rain_12_wk_total_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
  sero_form <- "+f(inla.group(time_since_switch, n = 18), model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
  clim_form <- "+f(inla.group(max_t_scale_12_wk_avg_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(nino34_12_wk_avg_4, n = 12), model = 'rw2', scale.model = TRUE, hyper = precision_prior)+f(inla.group(days_no_rain_12_wk_total_0),model = 'rw2', scale.model = TRUE, hyper = precision_prior)"
  baseline <- ""
  baseline_year <- ""

  if (form_input == "all") {
    forms <- list(sero_clim_form, sero_form, clim_form, baseline, baseline_year)
    form_names <- c("sero-climate", "sero-only", "climate-only", "seasonal-baseline", "seasonal-year-baseline")
  } else if (form_input == "covariate-only") {
    forms <- list(sero_clim_form, sero_form, clim_form)
    form_names <- c("sero-climate", "sero-only", "climate-only")
  } else if (form_input == "sero_clim_form") {
    forms <- get(form_input)
    form_names <- "sero-climate"
  } else if (form_input == "clim_form") {
    forms <- get(form_input)
    form_names <- "climate_only"
  } else if (form_input == "sero_form") {
    forms <- get(form_input)
    form_name <- "sero-only"
  } else if (form_input == "baseline") {
    forms <- get(form_input)
    form_names <- "seasonal-baseline"
  } else if (form_input == "baseline_year") {
    forms <- get(form_input)
    form_names <- "seasonal-year-baseline"
  } else {
    print("Inputted formula not regonised")
  }

  # Set variables to use for prediction
  if (horizon == 0) {
    predict_vars <- c("time_since_switch", "max_t_scale_12_wk_avg_0", "nino34_12_wk_avg_4", "days_no_rain_12_wk_total_0")
  } else if (horizon == 2) {
    predict_vars <- c("time_since_switch_2", "max_t_scale_10_wk_avg_2", "nino34_12_wk_avg_4", "days_no_rain_10_wk_total_2")
  } else if (horizon == 4) {
    predict_vars <- c("time_since_switch_4", "max_t_scale_8_wk_avg_4", "nino34_12_wk_avg_4", "days_no_rain_8_wk_total_4")
  } else if (horizon == 6) {
    predict_vars <- c("time_since_switch_6", "max_t_scale_6_wk_avg_6", "nino34_10_wk_avg_6", "days_no_rain_6_wk_total_6")
  } else if (horizon == 8) {
    predict_vars <- c("time_since_switch_8", "max_t_scale_4_wk_avg_8", "nino34_8_wk_avg_8", "days_no_rain_4_wk_total_8")
  }

  print(cat(c("Using these variables for prediction: ", predict_vars)))

  # Set up templates
  temp <- list(
    post.samples = NULL,
    outbreak.probs = NULL,
    trigger = NULL,
    roc = NULL,
    auc = NULL,
    model = NULL
  )

  temp_out <- list()

  # Set up parameters
  s <- 1000 # Number of samples
  y_pred_full <- matrix(NA, 1200, s) # Data frame to store posterior predictions for each week

  for (f in 1:length(forms)) {
    # Create formulae for random effects
    precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
    eweek_re <- " + f(eweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
    year_re <- " + f(year_index, model='iid', hyper = precision_prior)"

    cat(paste0("Running evaluation for model ", f, " with a lag of ", horizon, ".\n"))

    # Whether null models should include yearly random effect or seasonality only (note this will be an average yearly effect as year is set to NA below)
    if (form_names[[f]] != "seasonal-baseline") {
      res <- paste0(eweek_re, year_re)
    } else {
      res <- paste0(eweek_re)
    }

    cat(paste0("Random effects used are: ", res, "\n"))

    idx.all <- which(data_input$date_index >= 471)
    form <- forms[[f]]
    cat(paste0("Formula used is: ", form, "\n"))

    for (n in 471:1200) { # Should be 471:1200 (using Jan 2000 to Jan 2009 as training data)

      df_pred <- data_input |>
        filter(date_index <= n) |>
        mutate(
          cases = case_when(
            date_index > n - horizon ~ NA_integer_,
            TRUE ~ cases
          ),
          time_since_switch = case_when(
            date_index > n - horizon & date_index != n ~ NA_integer_,
            TRUE ~ time_since_switch
          ),
          max_t_scale_12_wk_avg_0 = case_when(
            date_index > n - horizon & date_index != n ~ NA_integer_,
            TRUE ~ max_t_scale_12_wk_avg_0
          ),
          nino34_12_wk_avg_4 = case_when(
            date_index > n - horizon & date_index != n ~ NA_integer_,
            TRUE ~ nino34_12_wk_avg_4
          ),
          days_no_rain_12_wk_total_0 = case_when(
            date_index > n - horizon & date_index != n ~ NA_integer_,
            TRUE ~ days_no_rain_12_wk_total_0
          )
        )

      # What to do with yearly RE for prediction
      if (yearly_re == "na" & form_names[[f]] != "seasonal-year-baseline") {
        df_pred <- df_pred |>
          mutate(year_index = case_when(date_index == n ~ NA, T ~ year_index))
        cat(paste0("Replacing year_index with NA \n"))
      }

      # If horizon != 0 then use variables from h weeks ago (lagged variables) to predict week n
      if (horizon != 0) {
        df_pred <- df_pred |>
          mutate(
            time_since_switch = case_when(date_index == n ~ get(predict_vars[1]), TRUE ~ time_since_switch),
            max_t_scale_12_wk_avg_0 = case_when(date_index == n ~ get(predict_vars[2]), TRUE ~ max_t_scale_12_wk_avg_0),
            nino34_12_wk_avg_4 = case_when(date_index == n ~ get(predict_vars[3]), TRUE ~ nino34_12_wk_avg_4),
            days_no_rain_12_wk_total_0 = case_when(date_index == n ~ get(predict_vars[4]), TRUE ~ days_no_rain_12_wk_total_0)
          )
      }

      # If horizon = 0 using data from this week to predict cases this week
      if (horizon == 0) {
        df_pred <- df_pred |>
          mutate(cases = case_when(date_index == n ~ NA_integer_, TRUE ~ cases))
      }

      # Run model
      mod <- inla(as.formula(paste0("cases~1", res, form)),
        family = "nbinomial",
        offset = log(pop / 100000),
        control.inla = list(strategy = "adaptive"),
        control.predictor = list(link = 1, compute = TRUE),
        control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
        control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
        num.threads = 4,
        verbose = FALSE,
        data = df_pred
      )

      # Extract samples
      xx <- inla.posterior.sample(s, mod)

      # Extract size parameter and predicted values
      xx_s <- inla.posterior.sample.eval(
        function(...) {
          c(
            theta[1], # This is the size parameter of the negative binomial distribution (overdispersion)
            Predictor
          )
        },
        xx
      )
      xx_s <- xx_s[c(1, n + 1), ] # Select size parameter and predictor of interest, note this is the n+1th row as the first row is the size parameter
      # Compute the posterior predictive distribution for each of the
      # weeks for which a prediction was made using the posterior mean
      # of the target week
      y_pred <- matrix(NA, 1, s)
      for (s_idx in 1:s) {
        xx_sample <- xx_s[, s_idx]
        y_pred[, s_idx] <- rnbinom(1,
          mu = exp(xx_sample[-1]),
          size = xx_sample[1]
        )

        if (is.na(y_pred[, s_idx])) {
          print(paste0("NA prediction generated with mu = ", exp(xx_sample[-1]), " and size = ", xx_sample[1]))
        }

        y_pred_full[n, ] <- y_pred
      }

      cat(paste0("Progress: ", n, "/", 1200, " for model ", f, "."))
    }

    # Calculate outbreak probabilities
    outbreak <- NULL
    for (o in 1:1200) {
      if (!is.na(sum(y_pred_full[o, ]))) {
        outbreak$prob[o] <- length(y_pred_full[o, ][y_pred_full[o, ] > data_input$threshold[o]]) / s
        print(outbreak$prob[o])
      } else {
        outbreak$prob[o] <- NA
      }

      outbreak$thresh[o] <- data_input$threshold[o]
      outbreak$cases[o] <- data_input$cases[o]
    }

    outbreak <- as.data.frame(outbreak)
    outbreak <- outbreak |>
      mutate(outbreak = case_when(
        cases > thresh ~ 1,
        TRUE ~ 0
      ))

    temp[["post.samples"]] <- y_pred_full
    temp[["outbreak.probs"]] <- outbreak
    temp[["model"]] <- form_names[f]
    temp[["yearly_re"]] <- yearly_re
    temp_out[[f]] <- temp

    output_dir <- here("output")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }

    saveRDS(temp_out, here("output", paste0(filename, "_lag-", horizon, "_", Sys.Date(), "_temp.rds")))

    # Remove temp file
    cat(paste0("Predictions completed for model ", f, " of ", length(forms), " at time ", Sys.time(), "\n"))
    cat(Sys.time() - start_time)
  }
  saveRDS(temp_out, here("output", paste0(filename, "_horizon-", horizon, "_", Sys.Date(), ".rds")))
  file.remove(here("output", paste0(filename, "_lag-", horizon, "_", Sys.Date(), "_temp.rds")))
}

