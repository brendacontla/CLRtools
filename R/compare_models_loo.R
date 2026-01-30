#' @title Compare Bayesian Models Using PSIS-LOO
#'
#' @description
#' This function compares multiple Bayesian models using PSIS-LOO (Pareto-smoothed importance sampling leave-one-out cross-validation)
#' from the `loo` package. It returns a comparison table and a plot of the estimated ELPD (expected log predictive density) with standard errors.
#'
#' @import ggplot2
#' @importFrom loo loo_compare
#' @importFrom rstan loo
#'
#' @param ... Two or more `stanfit` model objects to compare. Each model must include pointwise log-likelihood values (usually named `log_lik`) stored in the generated quantities or transformed parameters block.
#' @param k A numeric value specifying the Pareto-k diagnostic threshold. Default is 0.7.
#' @param name_log A character string specifying the name of the log-likelihood parameter in the model. Default is `"log_lik"`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{`p_loo`}{A ggplot object showing `elpd_loo` values and standard errors for each model.}
#'   \item{`comparison`}{A `loo_compare` table comparing the relative fit of the models.}
#'   \item{`models_loo`}{A named list of individual `loo` objects for each model.}
#' }
#'
#' @details This function performs PSIS-LOO diagnostics on each model, creates a visual summary, and ranks them using `loo_compare`.
#' Ensure that each model includes pointwise log-likelihood values named consistently (e.g., `"log_lik"`).
#'
#' @export

compare_models_loo <- function(..., k = 0.7, name_log = "log_lik"){
  models <- list(...)

  ## Check inputs
  if (length(models) < 2) stop("Please provide at least two fitted model objects.")
  if (!is.numeric(k) || length(k) != 1 || k < 0 || k > 1) stop("`k` must be a single numeric value between 0 and 1.")
  if (!is.character(name_log) || length(name_log) != 1) stop("`name_log` must be a single character string indicating the log-likelihood variable name.")

  models_name <- as.character(match.call(expand.dots = FALSE)$...)
  names(models) <- models_name
  models_loo <- list()

  for(m in models_name){
    if (!inherits(models[[m]], "stanfit")) stop("`model` must be an object of class 'stanfit' (from rstan).")
    if (!name_log %in% names(models[[m]]@par_dims)) stop("The model does not contain `name_log` in generated quantities.")
    l_model<-rstan::loo(models[[m]], pars = name_log, k_threshold = k)
    models_loo[[m]] <- l_model
    plot(l_model)
  }

  ## Plot
  loo_df <- do.call(rbind, lapply(models_loo, function(x) {
    est <- x$estimates["elpd_loo", "Estimate"]
    se <- x$estimates["elpd_loo", "SE"]
    data.frame(elpd_loo = est, se = se)
  }))

  # Add model names as a column
  loo_df$model <- rownames(loo_df)

  p1 <- ggplot(loo_df, aes(x = elpd_loo, y = model)) +
    geom_point(size = 4) +
    geom_errorbar(aes(xmin = elpd_loo - se, xmax = elpd_loo + se), width = 0.2) +
    labs(x = 'elpd (LOO)', y = "Model") +
    theme_bw()

  ## Compare the LOO
  comparison <- loo_compare(models_loo)

  return(list(p_loo = p1, comparison = comparison, models_loo = models_loo))
}

