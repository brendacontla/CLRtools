#' @title Generate MCMC Diagnostic Plots for a Bayesian Model
#'
#' @description
#' Produces a set of standard diagnostic plots to evaluate convergence and sampling efficiency
#' for a Bayesian model fitted using `rstan`.
#'
#' @importFrom bayesplot mcmc_rhat_hist mcmc_neff_hist mcmc_trace mcmc_rank_overlay color_scheme_set
#'
#' @param model A fitted Bayesian model object created with `rstan::stan()`.
#' @param var.param A character vector specifying the names of the model parameters (e.g., slopes or coefficients)
#' to include in the trace and rank plots.
#'
#' @return A list of four ggplot2 objects:
#' \describe{
#'   \item{plot_rhat}{A histogram of Rhat values (`mcmc_rhat_hist`).}
#'   \item{plot_neff}{A histogram of effective sample sizes (`mcmc_neff_hist`).}
#'   \item{plot_trace}{Trace plots of the MCMC chains for selected parameters (`mcmc_trace`).}
#'   \item{plot_trank}{Rank plots of the selected parameters across chains (`mcmc_rank_overlay`).}
#' }
#'
#' @details This function uses the `bayesplot` package to visualize diagnostics commonly used to assess
#' convergence and sampling performance in MCMC estimation.
#'
#'
#' @export

diagnostic_bayes <- function(model, var.param){
  #Check inputs
  if (!inherits(model, "stanfit")) stop("`model` must be an object of class 'stanfit' (from rstan).")
  if (any(!var.param %in% names(rstan::extract(model)))) stop("Some values in `var.param` do not match parameter names in the model.")

  rhat_values <- rstan::summary(model)$summary[, "Rhat"]
  neff_values <- rstan::summary(model)$summary[, "n_eff"]

  color_scheme_set("mix-blue-pink")
  p1 <- mcmc_rhat_hist(rhat_values) + ggtitle('Rhat statistics')
  p2 <- mcmc_neff_hist(neff_values) + ggtitle("Effective Sample Size (ESS) Histogram")
  p3 <- mcmc_trace(model,  pars = var.param)
  p4 <- mcmc_rank_overlay(model,  pars = var.param)

  print(p1)
  print(p2)
  print(p3)
  print(p4)

  invisible(list(plot_rhat = p1, plot_neff = p2, plot_trace = p3, plot_trank = p4))
}
