#' @title Compare Bayesian Models by Predictor Using Posterior Predictive Simulations
#'
#' @description
#' This function compares posterior predictive distributions from several Bayesian models across levels of a selected predictor variable.
#' For numeric predictors, the variable is binned; for categorical predictors, the original factor levels are used directly.
#' The function visualizes the distribution of simulated means and standard deviations per predictor level alongside the observed values.
#'
#' @import ggplot2
#' @importFrom dplyr mutate group_by summarise bind_rows rename
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom ggpubr ggarrange
#' @importFrom stats sd
#'
#' @param data A data frame containing the original dataset.
#' @param models A named list of fitted `stanfit` model objects. Each name will be used as the model label.
#' @param parameters Optional. A **named list** mapping each model to a named character vector where each name is a variable in the data and the value is the name of the corresponding parameter/coefficient in the model.
#' Must have the same names as `models`. Required if `ypredict` is not provided.
#' @param var.plot A single character string. Name of the predictor variable in `data` used for binning and plotting. Must be in all the models.
#' @param intercept Optional. A named list with the intercept parameter names for each model. Each entry should be a character string or `NULL` if no intercept is used. Must have the same names as `models`.
#' @param ypredict Optional. A named list of posterior predictive matrices. Each matrix should have rows as posterior draws
#'   and columns as data points. If not provided, predictions are computed internally. Must have the same names as `models`.
#' @param outcome A character string. The name of the outcome variable in `data`.
#' @param mbreaks Number of bins if `var.plot` is numeric; ignored if it's a factor.
#'
#' @details
#' This function provides a visual diagnostic for comparing posterior predictive summaries across multiple Bayesian models.
#' The predictor variable can be either continuous or categorical. For continuous variables, the range is divided into bins using `cut()` and `mbreaks`; for categorical variables, no binning is applied.
#' Posterior predictive distributions are either precomputed via `ypredict` or generated internally using `parameters` and (optionally) `intercept`.
#'
#' @return A list containing:
#' \describe{
#'   \item{models_summary}{A data frame summarizing the posterior predictive means and standard deviations per draw, model, and predictor level (or bin).}
#'   \item{p_mean}{A `ggplot` showing the distribution of simulated outcome means for each bin and model, overlaid with the observed means and sample size per bin.}
#'   \item{p_sd}{A `ggplot` showing the distribution of simulated outcome standard deviations for each bin and model, overlaid with the observed standard deviations and sample size per bin.}
#' }
#' The function also prints the plots side by side using `ggarrange()`.
#'
#' @export
#'
compare_bayesm_by_predictor <- function(data, models, parameters = NULL, var.plot, intercept = NULL, ypredict = NULL, outcome, mbreaks){
  # Check inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.list(models) || is.null(names(models)) || any(!nzchar(names(models)))) stop("`models` must be a named list.")
  if(!is.null(ypredict)){
    if (!is.list(ypredict) || is.null(names(ypredict)) || any(!nzchar(names(ypredict)))) stop("`ypredict` must be a named list.")
    if (!all(names(ypredict) %in% names(models))) stop("All names in `ypredict` must be present in `models`.")
    models_ypredict <- names(ypredict)
  }else{
    if(!is.null(parameters)){
      if (!is.list(parameters) || is.null(names(parameters)) || any(!nzchar(names(parameters)))) stop("`parameters` must be a named list.")
      if (!all(names(parameters) %in% names(models))) stop("All names in `parameters` must be present in `models`.")
    }
    if(!is.null(intercept)){
      if (!is.list(intercept) || is.null(names(intercept)) || any(!nzchar(names(intercept)))) stop("`intercept` must be a named list.")
      if (!all(names(intercept) %in% names(models))) stop("All names in `intercept` must be present in `models`.")
    }
  }

  if (!is.character(outcome) || length(outcome) != 1) stop("`outcome` must be a single character string.")
  if (!is.character(var.plot) || length(var.plot) != 1) stop("`var.plot` must be a single character string.")
  if (!var.plot %in% names(data)) stop("`var.plot` must be present in `data`.")

  # Setting parameters
  model_names <- names(models)
  models_summary <- data.frame()

  # Getting the posterior predictive
  for(m in model_names){

    if(m %in% models_ypredict){

      m_predict <- ypredict[[m]]

    }else{
      ## Check for each model the parameters
      model <- models[[m]]
      if (!inherits(model, "stanfit")) stop("`model` must be an object of class 'stanfit' (from rstan).")

      ## Checking if the model has the parameters
      if (!m %in% names(parameters)) {
        stop(paste0("Model '", m, "' is not provided in `ypredict` and is also missing or unnamed in `var.param`. Cannot compute predictions."))
      }

      param_model <- parameters[[m]]
      if (is.null(param_model) || is.null(names(param_model)) || any(!nzchar(names(param_model)))) stop(paste0("`var.param[[", m, "]]` must be a named character vector."))

      ## Getting the names of parameters and checking intercept
      var.names <- names(param_model)
      if(!is.null(intercept)){
        intercept_model <- intercept[[m]]
        if (!is.character(intercept_model) || length(intercept_model) != 1){
          stop(paste0("`intercept[[", m, "]]` must be a single character string."))
        }else{
          param_model <- c(param_model, intercept_model)
        }
      }else{
        intercept_model <- NULL
      }

      ## Check that all variables exist in data
      posterior <- as.matrix(model)
      S <- dim(posterior)[1]

      missing_vars <- setdiff(c(outcome, var.names), colnames(data))
      if (length(missing_vars) > 0) {
        stop(paste0("The following variables are missing from `data`: ", paste(missing_vars, collapse = ", ")))
      }

      missing_params <- setdiff(param_model, colnames(posterior))
      if (length(missing_params) > 0) {
        stop(paste0("The following parameters are missing from the posterior samples of model '", m, "': ",
                    paste(missing_params, collapse = ", ")))
      }

      m_predict <- predict_posterior(mod = model, dat = data, sim = S, inter = intercept_model, var_param = param_model, var_names = var.names)
    }

    ## Getting the bins if the predictor variable is continuous
    if(inherits(data[, var.plot], 'numeric')){
      if (!is.numeric(mbreaks) || length(mbreaks) != 1) stop("`mbreaks` must be a numeric value.")
      df_total <- data.frame(data[, c(var.plot, outcome)], t(m_predict)) %>%
        mutate(E_bin = cut(!!sym(var.plot), breaks = mbreaks)) %>%
        pivot_longer(cols = starts_with("X"), names_to = "draw", values_to = "y_rep")
    }else{
      df_total <- data.frame(data[, c(var.plot, outcome)], t(m_predict)) %>%
        rename(E_bin = !!sym(var.plot)) %>%
        pivot_longer(cols = starts_with("X"), names_to = "draw", values_to = "y_rep")
    }



    # For total model
    summary_total <- df_total %>%
      group_by(E_bin, draw) %>%
      summarise(mean_r = mean(y_rep), sd_r = sd(y_rep), .groups = "drop_last") %>%
      mutate(model_name = m)

    models_summary <- bind_rows(models_summary, summary_total)
  }

  #Summary data per mean
  if(inherits(data[, var.plot], 'numeric')){
    observed_data_mean <- data %>%
      mutate(E_bin = cut(!!sym(var.plot), breaks = mbreaks)) %>%
      group_by(E_bin) %>%
      summarise(mean_observed = mean(!!sym(outcome)))

    bin_counts <- data %>%
      mutate(E_bin = cut(!!sym(var.plot), breaks = mbreaks)) %>%
      group_by(E_bin) %>%
      summarise(n = n())
  }else{
    observed_data_mean <- data %>%
      rename(E_bin = !!sym(var.plot)) %>%
      group_by(E_bin) %>%
      summarise(mean_observed = mean(!!sym(outcome)))

    bin_counts <- data %>%
      rename(E_bin = !!sym(var.plot)) %>%
      group_by(E_bin) %>%
      summarise(n = n())
  }

  # Plot for mean
  dodge_width <- position_dodge(width = 0.4)
  y_position_mean <- min(models_summary$mean_r)
  p_mean <- ggplot() +
    geom_boxplot(data=models_summary, aes(x = E_bin, y = mean_r, fill=model_name), position = dodge_width, alpha= 0.85) +
    geom_point(data=observed_data_mean,  aes(x = E_bin, y = mean_observed, color='Observed'), size = 3)+
    geom_text(data = bin_counts, aes(x = E_bin, y = y_position_mean, label = paste0("n=", n)), vjust = 2, size = 4) +
    scale_color_manual(values=c('red'))+
    scale_fill_brewer(palette = "BuPu")+
    labs(x = "Binned Predictor Variable", y = "Simulated Mean of Outcome per Bin", title = "Distribution of Posterior Predictive Means by Bin", color ='', fill='Model') +
    theme_bw()

  #Summary data per SD
  if(inherits(data[, var.plot], 'numeric')){
    observed_data_sd <- data %>%
      mutate(E_bin = cut(!!sym(var.plot), breaks = mbreaks)) %>%
      group_by(E_bin) %>%
      summarise(sd_observed = sd(!!sym(outcome)))
  }else{
    observed_data_sd <- data %>%
      rename(E_bin = !!sym(var.plot)) %>%
      group_by(E_bin) %>%
      summarise(sd_observed = sd(!!sym(outcome)))
  }


  # Plot for SD
  y_position_sd <- min(models_summary$sd_r)
  p_sd <- ggplot() +
    geom_boxplot(data=models_summary, aes(x = E_bin, y = sd_r, fill=model_name), position = dodge_width, alpha= 0.85) +
    geom_point(data=observed_data_sd,  aes(x = E_bin, y = sd_observed, color='Observed'), size = 3)+
    geom_text(data = bin_counts, aes(x = E_bin, y = y_position_sd, label = paste0("n=", n)), vjust = 2, size = 4) +
    scale_color_manual(values=c('red'))+
    scale_fill_brewer(palette = "BuPu")+
    labs(x = "Binned Predictor Variable", y = "Simulated SD of Outcome per Bin", title = "Distribution of Posterior Predictive SDs by Bin", color ='', fill='Model') +
    theme_bw()

  print(ggarrange(p_mean, p_sd, ncol = 2))

  invisible(return(list(models_summary = models_summary, p_mean = p_mean, p_sd = p_sd)))
}
