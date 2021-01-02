# b)
############## Pseudo-Code ##################
# do not run #
##### Top-Level-function #####
compute_stability_paths(model,
                        data,
                        resampling_options) {
  
  # loop over resampling_options$repetitions
  iteration_sample <- resample_data(data, resampling_options)
  # refit the model for different penelaization parameters
  iteration_model <- refit_model(model, iteration_sample)
  selection_matrices[rightspot] <- save_selection(iteration_model)
  # loop end
  
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

################### execute from here ##################

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
  checkmate::assert_vector(strata,
    any.missing = FALSE,
    len = nrow(data), null.ok = TRUE
  )
  checkmate::assert_number(fraction, lower = 0, upper = 1)

  selected <- vector("list", reps)
  for (i in seq_len(reps)) {
    new_data <- resample(data,
      method = method, strata = strata,
      fraction = fraction
    )
    new_model <- refit(model, new_data)
    selected[[i]] <- get_selected(data, new_model)
  }
  stability_paths <- make_paths(selected, reps)
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
      fraction = fraction
    )
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

sample_without_replacement <- function(nrows, fraction, strata = NULL) {
  # check if sampling is possible
  checkmate::assert(round(fraction * nrows) >= 1)

  abs_sample_size <- ceiling(fraction * nrows)

  if (is.null(strata)) {
    return(sample(nrows,
      size = abs_sample_size,
      replace = FALSE
    )) # --> early exit!
  }
  
  # for strata sample we need separate sizes
  abs_sample_sizes <- ceiling(table(strata) * fraction)
  row_sample <- c()
  
  for (i in names(abs_sample_sizes)) {
    row_indizes <- seq_len(nrows)[strata %in% c(i, as.numeric(i))]
    row_sample <- c(row_sample, sample(row_indizes,
                                      size = abs_sample_sizes[[i]], 
                                      replace = FALSE))
  }
  row_sample
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

############## get_selected ###########################################################

# extract information from model which features
# have been selected for which penalization parameter
get_selected <- function(data, model) {
  # currently the only implemented case
  if (class(model) == "regsubsets") {
    # extract select information
    selections_per_model <- summary(model)[["which"]]

    # generate output matrix
    selected_matrix <- matrix(FALSE,
      ncol = ncol(data),
      nrow = nrow(selections_per_model)
    )
    colnames(selected_matrix) <- colnames(data)

    # check for colnames of selected eatures
    selected_features <- colnames(selections_per_model)

    # check if model has intercept or not
    has_intercept <- model[["intercept"]]
    # if so, don't check for the first column
    for (i in seq(
      1 + has_intercept,
      ncol(selections_per_model)
    )) {
      # insert the information on variable i to 'selected_matrix'
      selected_matrix[, selected_features[i]] <-
        selections_per_model[, i]
    }
    # add FALSE-only row for Intercept-only model
    return(rbind(
      rep(FALSE, ncol(selected_matrix)),
      selected_matrix
    ))
  }
}

############## make_paths ###########################################################

# averages over a lists of selection matrices
# to compute the probability for each feature to be
# included in a model with specific penalization parameter
make_paths <- function(selected, reps) {
  # extract the first matrix to have a layout for the final output
  sum_matrix <- selected[[1]]
  # if there have been more than one repetition,
  # compute the averages
  if (reps > 1) {
    for (i in seq(2, reps)) {
      sum_matrix <- sum_matrix + selected[[i]]
    }
  }
  # this line of code can also be used as type conversion
  # for the case reps = 1
  paths <- sum_matrix / reps
  paths
}

############## plot_paths ###########################################################

## I assume this should look like this:

# transform stability_paths to plotable format
prepare_stability_path_toplot <- function(stability_paths) {
  plot_data <- data.frame(
    row.names =
      seq(1, ncol(stability_paths) * nrow(stability_paths))
  )
  plot_data$variable <- rep(colnames(stability_paths),
    each = nrow(stability_paths)
  )
  # if stability_paths contain information about the penalization parameter
  if (!is.null(rownames(stability_paths))) {
    plot_data$penalization <- rep(rownames(stability_paths),
      times = ncol(stability_paths)
    )
    # convert it to a numerical value for correct display on x-axis
    plot_data[["penalization"]] <- as.numeric(plot_data[["penalization"]])
  } else {
    # otherwise simply take the row number
    plot_data$penalization <- rep(seq(0, nrow(stability_paths) - 1),
      times = ncol(stability_paths)
    )
  }

  plot_data$probability <- NA
  # fill it with the values from stability_path
  for (i in seq_len(ncol(stability_paths))) {
    first_index <- (i - 1) * nrow(stability_paths) + 1
    last_index <- i * nrow(stability_paths)
    plot_data[first_index:last_index, "probability"] <- stability_paths[, i]
  }
  return(plot_data)
}


# plots stability paths
plot_stability_paths <- function(stability_paths) {
  # input checking - should be a numeric matrix with minimum of one column
  # and one row
  # columns should have names (to label the features)
  # no missing values are allowed
  checkmate::assert_matrix(stability_paths,
    min.cols = 1,
    min.rows = 1,
    mode = "numeric",
    any.missing = FALSE,
    col.names = "named"
  )

  # make stability_paths plotable
  plot_data <- prepare_stability_path_toplot(stability_paths)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = penalization,
      y = probability
    )
  ) +
    ggplot2::geom_line(ggplot2::aes(color = factor(variable))) +
    ggplot2::scale_x_continuous(breaks = seq(
      min(plot_data[["penalization"]]),
      max(plot_data[["penalization"]])
    ))
}

## Tests in topdown-stability-ex.Rmd scheinen auf
# nicht-vorhanden Dateien im Repo zugreifen zu wollen - wird deshalb ignoriert
# test-get-stability-paths.R funktioniert:
source("test-get-stability-paths.R")
stability_paths

# try to plot it
plot_stability_paths(stability_paths)
