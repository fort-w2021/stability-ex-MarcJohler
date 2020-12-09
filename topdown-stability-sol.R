# b)
##### Pseudo-Code #####
# do not run #
##### Top-Level-function #####
compute_stability_paths(model,
                        data,
                        resampling_options) {
  
  ## loop over resampling_options$repetitions
  iteration_sample <- resample_data(data, resampling_options)
  # loop over different penalization parameters
  iteration_model <- refit_model(model, iteration_sample)
  selection_matrices[rightspot] <- save_selection(iteration_model)
  # loop end
  ## loop end
  
  # compute probability of feature selection 
  # by averaging over all selection matrices
  compute_proportion(selection_matrices)
}

#### Resampling ####
resample_data(data, 
              resampling_options) {
  
  # subsampling
  execute_subsampling # fraction required
  # OR
  execute_bootstrapping 
}

# resampling_options must contain fraction and strata
execute_subsampling(data,
                    fraction,
                    strata) {
  # loop over strata 
  resample_without_replacement(fraction)
  # loop end
}

execute_bootstrapping(data,
                    strata) {
  # loop over strata 
  resample_with_replacement
  # loop end
}

#### Refitting the model ####
refit_model(model,
            sample) {
  # get the model formula (or parameters) for chosen method
  modelspecifics <- get_modelspecifics(model)
  # apply that formula (or parameters) to the current sample
  apply_modelspecifics(modelspecifics, sample)
}

get_formula(model) {
  # What happens here depends on the chosen method 
  extract_parameters
  # OR
  extract_formula
}

apply_formula(formula, sample) {
  # execute for example:
  random_fores(modelspecifics)
  # OR
  lasso(modelspecifics)
  # OR
  randomized_lasso(modelspecifics)
  # OR ...
}

#### Save selection ####
save_selection(model) {
  # extract chosen features
  # vector of TRUEs and FALSEs
}


##### Implementation #####
# function to compute stability paths for 
# a specific model with a certain resampling technique
get_stability_paths <- function(model, data, reps = 100,
                                method = c("subsample", "bootstrap"),
                                strata = NULL, fraction = 0.5) {
  checkmate::assert_class(model, "regsubsets")
  checkmate::assert_data_frame(data)
  checkmate::assert_count(reps)
  method <- match.arg(method)
  checkmate::assert_vector(strata, any.missing = FALSE,
                           len = NROW(data), null.ok = TRUE)
  checkmate::assert_number(fraction, lower = 0, upper = 1)
  
  selected <- vector("list", reps)
  for (i in seq_len(reps)) {
    new_data <- resample(data, method = method, strata = strata,
                         fraction = fraction)
    new_model <- refit(model, new_data)
    selected[[i]] <- get_selected(new_model)
  }
  stability_paths <- make_paths(selected)
  stability_paths
}

############## resample ########################################################

resample <- function(data, method = c("subsample", "bootstrap"),
                     strata = NULL, fraction = 0.5) {
  nrows <- nrow(data)
  rows <- resample_rows(nrows, method, strata, fraction)
  data[rows, ]
}

resample_rows <- function(nrows, method, strata = NULL, fraction = 0.5) {
  switch(method,
         "bootstrap" = sample_with_replacement(nrows, strata),
         "subsample" = sample_without_replacement(nrows, strata,
                                                  fraction = fraction)
  )
}

sample_with_replacement <- function(nrows, strata = NULL) {
  if (is.null(strata)) {
    return(sample(nrows, replace = TRUE)) # --> early exit!
  }
  rows <- tapply(
    X = seq_len(nrows), INDEX = strata, FUN = sample, replace = TRUE
  )
  as.vector(rows)
}

############## refit ###########################################################

# redo subset selection <model> on <new_data>1
refit <- function(model, new_data) {
  # works by overwriting the data argument of the original model
  # and then re-doing the function call that produced the original model
  modelcall <- model$call
  modelcall$data <- new_data
  # use regsubsets-generic here instead of regsubsets.formula or other method as
  # these methods are not exported by {leaps}
  # (quote s.t. just the name of the function is handed over, not the
  # function code itself...)
  modelcall[[1]] <- quote(leaps::regsubsets)
  eval(modelcall)
}
