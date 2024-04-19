# -----------------------------------------------------------------------------------------------------------------------
# Function to run INLA models for a given list of formulae
# And helper function to calculate continuous ranked probility score using INLA model output
# Arguments:
# fit_data (data frame): data frame for model selection, must include variables specified in forms
# forms (list): list of strings of model formulae to test
# lims (numeric): upper limit for loop of model selection (if blank set to NULL)
# Outputs: List including:
# adeq_stats: data frame of adequacy statistics for all models
# fitted_vals: list of fitted values for each model
# fixed_effects: list of fixed effects for each model
# random_effects: list of random effects for each model
# Author: Emilie Finch, adapted from code by Martin Lotto Batista and Rachel Lowe
# -----------------------------------------------------------------------------------------------------------------------

# Helper function to calculate continuous ranked probability score

crps <- function(crps_mod, crps_data) {
  # Sample from the posterior
  s <- 10000 # Number of samples
  xx <- inla.posterior.sample(s, crps_mod)
  # Extract values of interest
  xx_s <- inla.posterior.sample.eval(
    function(...) c(theta[1], Predictor[1:1200]),
    xx
  )
  # Create posterior predictive sample
  y_pred <- matrix(NA, nrow(crps_data), s)
  for (s_idx in 1:s) {
    xx_sample <- xx_s[, s_idx]
    y_pred[, s_idx] <- rnbinom(nrow(crps_data),
      mu = exp(xx_sample[-1]), # Predicted means
      size = xx_sample[1]
    ) # Overdispersion parameter
  }

  # Calculate CRPS
  crps_val <- median(scoringutils::crps_sample(crps_data$cases, y_pred))

  return(crps_val)
}

# Helper function to wrangle model parameters

w_params <- function(x, n) {
  x %>%
    rownames_to_column(var = "var") %>%
    mutate(mod = n)
}

# Function to fit INLA models

fit_inla <- function(fit_data, forms, lims = FALSE) {

  start_time <- Sys.time()

  # Define objects for model output
  fits <- NULL
  params <- NULL
  random <- NULL
  adeq_stats <- NULL

  # Create formula for REs
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

  eweek_re <- " + f(eweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
  year_re <- " + f(year_index, model='iid', hyper = precision_prior)"
  res <- paste0(eweek_re, year_re)


  # Baseline model #
  # Run base model (REs only) and extract model adequacy statistics (DIC, WAIC, CPO, MAE and Rsq)
  null_mod <- inla(as.formula(paste0("cases~1", res)),
    family = "nbinomial",
    offset = log(pop / 100000),
    control.inla = list(strategy = "adaptive"),
    control.predictor = list(link = 1, compute = TRUE),
    control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
    control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
    verbose = FALSE,
    data = fit_data
  )

  fits[[1]] <- null_mod$summary.fitted.values %>% mutate(mod = 1)
  params[[1]] <- w_params(null_mod$summary.fixed, 1)
  random[[1]] <- lapply(null_mod$summary.random, mutate, mod = 1)

  # MAE base model
  mean_null_fit <- hydroGOF::mae(null_mod$summary.fitted.values$mean, fit_data$cases, na.rm = T)
  median_null_fit <- hydroGOF::mae(null_mod$summary.fitted.values$`0.5quant`, fit_data$cases, na.rm = T)

  # CRPS base model

  crps_null_fit <- crps(crps_mod = null_mod, crps_data = fit_data)

  add <- data.frame(
    mod = 1,
    form = "RE only",
    dic = null_mod$dic$dic,
    waic = null_mod$waic$waic,
    cpo = -mean(log(null_mod$cpo$cpo)),
    rsq = 0,
    mean_mae = mean_null_fit,
    median_mae = median_null_fit,
    mean_dif_mae = 0,
    median_dif_mae = 0,
    crps = crps_null_fit,
    crpss = 0
  )

  adeq_stats <- bind_rows(adeq_stats, add)

  print("Baseline model completed.")

  # Covariate models 
  
  if (lims == FALSE) {
    start_end <- 1:length(forms)
  } else {
    start_end <- lims
  }

  for (i in start_end) {
    tryCatch(
      {
        mod <- inla(as.formula(paste0("cases~1", res, forms[i])),
          family = "nbinomial",
          offset = log(pop / 100000),
          control.inla = list(strategy = "adaptive"),
          control.predictor = list(link = 1, compute = TRUE),
          control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
          control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
          verbose = FALSE,
          data = fit_data
        )

        fits[[i + 1]] <- mod$summary.fitted.values %>% mutate(mod = i + 1)
        params[[i + 1]] <- w_params(mod$summary.fixed, i + 1)
        random[[i + 1]] <- lapply(mod$summary.random, mutate, mod = i + 1)

        # Extract model adequacy criteria and save it in the table. Both MAE and Rsq are calculated relative to the base model (REs only)
        # R-square
        dev <- mod$dic$deviance.mean
        nulldev <- null_mod$dic$deviance.mean
        n <- nrow(mod$summary.fitted.values)
        x <- round(1 - exp((-2 / n) * ((dev / -2) - (nulldev / -2))), 3)

        # MAE
        mean_fit <- hydroGOF::mae(mod$summary.fitted.values$mean, fit_data$cases, na.rm = T)
        median_fit <- hydroGOF::mae(mod$summary.fitted.values$`0.5quant`, fit_data$cases, na.rm = T)

        mean_diff <- mean_null_fit - mean_fit
        median_diff <- median_null_fit - median_fit

        crps_mod <- crps(mod, fit_data)
        crpss_mod <- 1 - crps_mod / crps_null_fit

        # Create table
        add <- data.frame(
          mod = i + 1,
          form = forms[[i]],
          dic = mod$dic$dic,
          waic = mod$waic$waic,
          cpo = -mean(log(mod$cpo$cpo)),
          rsq = x,
          mean_mae = mean_fit,
          median_mae = median_fit,
          mean_dif_mae = mean_diff,
          median_dif_mae = median_diff,
          crps = crps_mod,
          crpss = crpss_mod
        )

        adeq_stats <- bind_rows(adeq_stats, add)

        print(paste0(i, "/", length(start_end)))
      },
      error = function(e) {
        cat(paste0("Error with model ", i, " and formula ", forms[i], ". Skipping to next model."))
      }
    )
  }
  print(Sys.time() - start_time)
  # FN outputs a list with i) model adequacy statistics, ii) fitted values, iii) fixed effects, iv) random effects
  return(list(
    adeq_stats = adeq_stats,
    fitted_vals = fits,
    fixed_effects = params,
    random_effects = random
  ))
}
