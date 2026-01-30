#' @title Summarize Bayesian Logistic Regression Model Results
#'
#' @description
#' This function provides a visual and numeric summary of a Bayesian logistic regression model
#' fitted with \code{rstan}. It displays posterior distributions of selected parameters and performs
#' posterior predictive checks based on either user-provided simulations or internal simulations.
#'
#' @import ggpubr
#' @import bayesplot
#' @importFrom rstan extract
#' @importFrom stats plogis
#'
#' @param model A fitted model object of class \code{stanfit}, from \code{rstan} or \code{stanreg}, from \code{rstanarm}, representing a single-level logistic regression.
#' @param ypredict Optional. A matrix of posterior predictive simulations of the outcome variable
#'   (e.g., generated externally). If \code{NULL}, predictions will be simulated internally assuming a single-level logistic regression.
#'   The matrix should have dimensions \code{S x N}, where \code{S} is the number of posterior draws (rows) and \code{N} is the number of observations (columns).
#' @param data A data frame containing the variables used in the model.
#' @param outcome A character string specifying the name of the binary outcome variable in `data`.
#' @param intercept Optional. A character string naming the intercept parameter in the model. Defaults to `NULL`.
#' @param var.param A named character vector mapping dataset variable names to model parameter names.
#' When `ypredict` is `NULL`, the variable names must exist in `data` for posterior prediction. If `ypredict` is provided, this check is skipped since predictions are supplied directly.
#' @param rounding An integer specifying the number of decimal places for printed parameter summaries. Must be a non-negative integer. Default is 2.
#' @param prob A numeric value between 0 and 1 specifying the width of the credible interval for posterior plots (e.g., 0.8 for 80% intervals). Default is 0.8.
#' @param point.est Character string, either `"mean"` or `"median"`, indicating the point estimate to use in posterior plots. Default is median.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{posterior_plot}{A `ggplot` object showing posterior distributions of the coefficients.}
#'   \item{ppc_mean}{A `ggplot` object showing the posterior predictive check for the point estimate.}
#'   \item{ppc_sd}{A `ggplot` object showing the posterior predictive check for the standard deviation.}
#' }
#' Also prints parameter summaries and displays plots.
#'
#'
#' @export

summarize_results <- function(model, ypredict =  NULL, data, outcome, intercept = NULL, var.param, rounding = 2, prob = 0.8, point.est = "median"){
  # Check inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(outcome) || length(outcome) != 1) stop("`outcome` must be a single character string.")
  if (!is.numeric(rounding) || rounding < 0) stop("`rounding` must be a non-negative integer.")
  if (!is.numeric(prob) || prob <= 0 || prob >= 1) stop("`prob` must be a number between 0 and 1.")
  if (!(point.est %in% c("mean", "median"))) stop("`point.est` must be either 'mean' or 'median'.")
  if (!inherits(model, c("stanfit", "stanreg"))) stop("`model` must be an object of class 'stanfit' (from rstan) or 'stanreg' (from rstanarm).")
  if (is.null(var.param) || is.null(names(var.param)) || any(!nzchar(names(var.param)))) stop("`var.param` must be a named, non-null vector.")


  var.names <- names(var.param)
  if(!is.null(intercept)){
    if (!is.character(intercept) || length(intercept) != 1){
      stop("`intercept` must be a single character string.")
    }else{
      var.param <- c(var.param, intercept)
    }
  }

  if (!all(c(outcome) %in% colnames(data))) stop("`outcome` must be a column in `data`.")

  # Summarize of the parameters
  if(inherits(model, "stanfit")){
    print(model, pars = var.param, digits = rounding)
  }else{
    print(rstan::summary(model), pars = var.param, digits = 6)
  }



  posterior <- as.matrix(model)
  p0 <- mcmc_areas(posterior, pars = var.param, prob = prob, point_est = point.est, area_method = "equal height") + ggtitle("Posterior distributions", paste0("with ", point.est, "s and ", prob * 100, "% intervals"))

  # Density plot of the probabilities
  ## Getting the posterior predictive
  if(is.null(ypredict)){
    if (!all(c(var.names) %in% colnames(data))) stop("All names in `var.param` (and `intercept`, if given) must be columns in `data`.")
    if (!all(c(var.param) %in% colnames(posterior))) stop("All values in `var.param` must be in the model.")
    S <- dim(posterior)[1]
    yrep <- predict_posterior(mod = model, dat = data, y = outcome, sim = S, inter = intercept, var_param = var.param, var_names = var.names)
  }else{
    if (!is.matrix(ypredict)) stop("`ypredict` must be a matrix.")
    if (ncol(ypredict) != nrow(data)) stop("Number of columns in `ypredict` must match the number of rows in `data`.")
    if (nrow(ypredict) < 2) warning("`ypredict` has fewer than 2 rows (posterior draws).")
    yrep <- ypredict
  }


  ## Plots
  color_scheme_set("mix-blue-pink")
  p1 <- ppc_stat(y = data[[outcome]], yrep = yrep, stat = 'mean') + ggtitle("Posterior Predictive Check:", "Simulated Means")
  p2 <- ppc_stat(y = data[[outcome]], yrep = yrep, stat = 'sd')+ ggtitle("Posterior Predictive Check:", "Simulated SD")

  print(p0)
  print(ggarrange(p1, p2, ncol = 2))

  invisible(list(posterior_plot = p0, ppc_pe = p1, ppc_sd = p2))
}
