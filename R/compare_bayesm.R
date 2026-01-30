#' @title Posterior Predictive Check for Multiple Bayesian Models
#'
#' @description
#' This function performs a posterior predictive check for one or more Bayesian models by computing
#' the mean and standard deviations of simulated predictions, comparing them
#' to the observed outcome. It produces grouped histograms for each model.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom tidyr pivot_longer
#' @importFrom stats sd
#'
#' @param models A **named list** of fitted `rstan` models (objects of class `stanfit`).
#' @param ypredict Optional. A named list of posterior predictive matrices. Each matrix should have rows as posterior draws
#'   and columns as data points. If not provided, predictions are computed internally. Must have the same names as `models`.
#' @param data A data frame containing the predictor variables and the outcome used for model prediction.
#' @param outcome A character string naming the outcome variable in `data`.
#' @param intercept (Optional) A named list with the intercept parameter names for each model. Each entry should be a character string or `NULL` if no intercept is used. Must have the same names as `models`.
#' @param var.param A **named list** mapping each model to a named character vector where each name is a variable in the data and the value is the name of the corresponding parameter/coefficient in the model.
#' Must have the same names as `models`. Required if `ypredict` is not provided.
#'
#' @return Invisibly returns a list with two ggplot objects:
#' \describe{
#'   \item{ppc_pe}{A ggplot object showing the distribution of simulated mean of the outcome across models.}
#'   \item{ppc_sd}{A ggplot object showing the distribution of simulated standard deviations of the outcome across models.}
#' }
#'
#'
#' @export

compare_bayesm <- function(models, ypredict=NULL, data, outcome, intercept = NULL, var.param=NULL){
  # Check inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(outcome) || length(outcome) != 1) stop("`outcome` must be a single character string.")
  if (!outcome %in% colnames(data)) stop(paste0("The outcome variable '", outcome, "' is not in `data`."))
  if (!is.list(models) || is.null(names(models)) || any(!nzchar(names(models)))) stop("`models` must be a named list.")
  if(!is.null(ypredict)){
    if (!is.list(ypredict) || is.null(names(ypredict)) || any(!nzchar(names(ypredict)))) stop("`ypredict` must be a named list.")
    if (!all(names(ypredict) %in% names(models))) stop("All names in `ypredict` must be present in `models`.")
  }

  if(!is.null(var.param)){
    if (!is.list(var.param) || is.null(names(var.param)) || any(!nzchar(names(var.param)))) stop("`var.param` must be a named list.")
  }

  if (!is.null(intercept)) {
    if (!is.list(intercept) || is.null(names(intercept)) || any(!nzchar(names(intercept)))) stop("`intercept` must be a named list with non-empty names.")
  }

  model_names <- names(models)
  models_y <- names(ypredict)


  df_pe <- data.frame()
  df_sd <- data.frame()

  for(m in model_names){
    ## Calculating the posterior predictive
    if(m %in% models_y){

      m_predict <- ypredict[[m]]

    }else{
      ## Check for each model the parameters
      model <- models[[m]]
      if (!inherits(model, "stanfit")) stop("`model` must be an object of class 'stanfit' (from rstan).")

      ## Checking if the model has the parameters
      if (!m %in% names(var.param)) {
        stop(paste0("Model '", m, "' is not provided in `ypredict` and is also missing or unnamed in `var.param`. Cannot compute predictions."))
      }

      param_model <- var.param[[m]]
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

    ## Stat Mean or median
    ypre_pe <- apply(m_predict, 1, mean)
    ypre_pe <- data.frame(mod = m, t(ypre_pe))

    df_pe <- rbind(df_pe, ypre_pe)

    ## Stat SD
    ypre_sd <- apply(m_predict, 1, sd)
    ypre_sd <- data.frame(mod = m, t(ypre_sd))

    df_sd <- rbind(df_sd, ypre_sd)
  }

  ## Stat Mean or median
  df_long_pe <- pivot_longer(df_pe, cols = -mod, names_to = "data_point", values_to = "prediction")

  ## Stat SD
  df_long_sd <- pivot_longer(df_sd, cols = -mod, names_to = "data_point", values_to = "prediction")

  ########### Plots
  obs_mean <- mean(data[[outcome]])
  obs_sd <- sd(data[[outcome]])

  p1 <- ggplot(df_long_pe, aes(x = prediction, fill = mod)) +
    geom_histogram(alpha = 0.5, position = "identity", color="#b97c9b") +
    geom_vline(xintercept = obs_mean, color = "darkred", linewidth = 1) +
    labs(x = NULL, y=NULL, fill = "Model")+
    scale_fill_brewer(palette = "Blues")+
    ggtitle("Posterior Predictive Check of the Outcome:", "Simulated Means")+
    theme_classic()

  print(p1)

  p2 <- ggplot(df_long_sd, aes(x = prediction, fill = mod)) +
    geom_histogram(alpha = 0.5, position = "identity", color="#b97c9b") +
    geom_vline(xintercept = obs_sd, color = "darkred", linewidth = 1) +
    labs(x = NULL, y=NULL, fill = "Model")+
    scale_fill_brewer(palette = "Blues")+
    ggtitle("Posterior Predictive Check of the Outcome:", "Simulated SD")+
    theme_classic()

  print(p2)

  invisible(list(ppc_pe = p1, ppc_sd = p2))

}




