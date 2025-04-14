# -----------------------------------------------------------------------------------------------------------------------
# Function to run INLA models for counterfactual scenario analysis
# Arguments:
# data_input (data frame): data frame for model selection, must include variables specified in forms
# forms (list): string of model formula for INLA model
# Outputs: List including:
# samples_out: 1000 samples of predicted cases for tuning and counterfactual weeks
# exp_quantiles: median value and 95th percentile of predicted cases for tuning and counterfactual weeks
# Author: Emilie Finch, adapted from code by Martin Lotto Batista and Rachel Lowe
# -----------------------------------------------------------------------------------------------------------------------

run_scenario <- function(data_input, form) {
  # Arguments #
  #   data: a dataframe object
  #   form: formula for INLA models
  
  start_time <- Sys.time()
  
  # Define objects for model output
  fits <- NULL
  params <- NULL
  random <- NULL
  adeq_stats <- NULL
  linear_predictor <- NULL
  mods <- NULL
  
  # Create formula for REs
  precision_prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  eweek_re <- " + f(eweek, model='rw2', cyclic=TRUE, hyper = precision_prior)"
  year_re <- " + f(year_index, model='iid', hyper = precision_prior)"
  res <- paste0(eweek_re, year_re)
  
  set.seed(141)
  
 mod <- inla(as.formula(paste0("cases~1", res, form)),
              family = "nbinomial",
              offset = log(pop / 100000),
              control.inla = list(strategy = "adaptive"),
              control.predictor = list(link = 1, compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE), # which model assessment criteria to include
              control.fixed = list(correlation.matrix = TRUE, prec.intercept = 1, prec = 1),
              verbose = FALSE,
              data = data_input
        )
 
  # Extract samples
  xx <- inla.posterior.sample(1000, mod, seed = 141L, num.threads = "1:1")
 
  # Extract size parameter and predicted values
  xx_s <- inla.posterior.sample.eval(
    function(...) {
      c(
        Predictor
      )
    },
    xx
  )
  
  samples_out <- exp(xx_s) 
  
  # Wrangle to quantiles

      exp_quantiles <- NULL
      exp_quantiles <- data.frame(
        median = matrixStats::rowQuantiles(samples_out, cols = c(1:1000), probs = 0.5),
        lower = matrixStats::rowQuantiles(samples_out, cols = c(1:1000), probs = 0.025),
        upper = matrixStats::rowQuantiles(samples_out, cols = c(1:1000), probs = 0.975),
        mean = matrixStats::rowMeans2(samples_out, cols = c(1:1000)),
        date = data_input$date
      )
           
  print(Sys.time() - start_time)

  return(list(samples_out = samples_out,
         exp_quantiles = exp_quantiles))
}
