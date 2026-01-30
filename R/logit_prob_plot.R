#' @title Plot Predicted Probabilities from a Logistic Model
#'
#' @description
#' This function visualizes the predicted probabilities from a Bayesian logistic regression model
#' fitted using `rstan`. It computes the posterior predicted probabilities over a grid defined by
#' two continuous predictor variables, and creates a plot showing how these probabilities vary
#' across their values. Color is used to represent the estimated probability, and the original
#' data points are overlaid for reference.
#'
#' @import ggplot2
#' @importFrom rstan extract
#' @importFrom dplyr rename
#'
#' @param data A data frame containing the original data used to fit the model.
#' @param model A fitted Stan model object of class `'stanfit'` (from `rstan`). Required if `ypredict` is not provided.
#' @param ypredict Optional. A matrix of posterior predictive simulations of the outcome variable
#'   (e.g., generated externally). If \code{NULL}, predictions will be simulated internally assuming a single-level logistic regression.
#'   The matrix should have dimensions \code{S x N}, where \code{S} is the number of posterior draws (rows) and \code{N} is the number of observations (columns).
#' @param parameters A named vector where the names are the predictor variable names (as in the data), and the values are the corresponding parameter names in the Stan model. Required if `ypredict` is not provided.
#' @param intercept Optional. A character string indicating the name of the intercept parameter
#'   in the Stan model (if present).
#' @param outcome A character string with the name of the binary outcome variable in the data.
#' @param predictors.plot A character vector of length 2 specifying which two predictor variables
#'   to use for the x and y axes of the plot.
#'
#' @return A `ggplot2` object showing the mean predicted probabilities across the grid of the
#'   two specified predictors, with the observed outcome overlaid as colored points.
#'
#' @export

logit_prob_plot <- function(data, ypredict = NULL, model = NULL, parameters = NULL, intercept = NULL, outcome, predictors.plot){
  # Check inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(outcome) || length(outcome) != 1) stop("`outcome` must be a single character string.")
  if (is.null(predictors.plot)) stop("`predictors.plot` must be a non-null vector.")
  if (length(predictors.plot) != 2) stop("`predictors.plot` must be a character vector of length 2.")

  ## Preparing the data for the probabilities
  x1 = predictors.plot[1]
  x2 = predictors.plot[2]

  grid <- expand.grid(
    X1 = seq(min(data[, x1]), max(data[, x1]), length.out = 100),
    X2 = seq(min(data[, x2]), max(data[, x2]), length.out = 100)
  )

  names_map <- c(X1 = x1, X2 = x2)
  grid <- dplyr::rename(grid, !!!setNames(names(names_map), names_map))

  if(is.null(ypredict)){
    if (!inherits(model, "stanfit")) stop("`model` must be an object of class 'stanfit' (from rstan).")
    if (is.null(parameters) || is.null(names(parameters)) || any(!nzchar(names(parameters)))) stop("`parameters` must be a named, non-null vector.")
    if (!all(names(parameters) %in% colnames(data))) stop("All names of `parameters` must be variables in the data.")
    ## Getting the parameters
    posterior_param <- as.matrix(model)
    N <- dim(grid)[1]
    S <- dim(posterior_param)[1]
    posterior_pred <- rstan::extract(model)
    var_names <- names(parameters)

    psim <- matrix(NA, nrow = S, ncol = N)

    ## Doing the matrix for the calculations
    #Checking if there is intercept
    if(!is.null(intercept)){
      if (!is.character(intercept) || length(intercept) != 1){
        stop("`intercept` must be a single character string.")
      }else{
        parameters <- c(parameters, intercept)
      }
    }

    if (!all(c(parameters) %in% names(posterior_pred))) stop("All values in `parameters` must be in the model.")

    matrix_model <- do.call(cbind, posterior_pred[parameters])
    data_model <- do.call(cbind, grid[var_names])

    # Calculating the probabilities
    for (s in 1:S) {
      if(is.null(intercept)){
        psim[s, ] <- plogis(matrix_model[s] %*% t(data_model))
      }else{
        psim[s, ] <- plogis(matrix_model[s, intercept] + matrix_model[s, -which(colnames(matrix_model) == intercept)] %*% t(data_model))
      }
    }
  }else{
    if (!is.matrix(ypredict)) stop("`ypredict` must be a matrix.")
    psim <- ypredict
  }

  observed_prob <- colMeans(psim)
  grid$prob <- observed_prob

  # Plot
  p <- ggplot() +
    geom_tile(data = grid, aes(x = .data[[x1]], y = .data[[x2]], fill = prob)) +  # adjust size!
    geom_point(data = data, aes(x = .data[[x1]], y = .data[[x2]], color = factor(.data[[outcome]])), size = 3) +
    scale_color_manual(values = c("lightgrey", "black"))+
    scale_fill_viridis_c(name = 'Mean of probabilities')  +
    labs(color = outcome) +
    theme_minimal()

  return(p)
}



