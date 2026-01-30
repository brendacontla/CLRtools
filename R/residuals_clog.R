#' @title Model Diagnostic for Conditional Logistic Regression
#'
#' @description
#' This function computes various diagnostic measures for a conditional logistic regression model,
#' including Pearson residuals, leverage, standardized Pearson residuals, delta chi, and delta beta values.
#' Additionally, it generates plots to visualize these diagnostics, helping to assess the model's fit and
#' identify potential influential data points.
#' This function assumes that the \code{strata} and \code{y} vectors are provided in the same order as they were used
#' when fitting the model with the \code{clogit} function.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @param model A fitted conditional logistic regression model of class \code{clogit}.
#' @param strata A vector of the strata (matching variables) that defines the pairs in the matched case-control design.
#' @param y A vector of the binary outcome variable (case-control indicator, where 1 represents the case and 0 represents the control).
#' @param plot A logical value indicating whether diagnostic plots should be generated. Defaults to TRUE.
#'
#' @details The diagnostic measures are calculated based on the formulas provided by Hosmer et al. (2013).
#' The following measures are computed:
#' \itemize{
#'   \item Pearson residuals
#'   \item Standardized Pearson residuals
#'   \item Leverage
#'   \item Delta statistics: Delta beta (change in Cook's distance) and Delta chi (change in Pearson chi-squared)
#' }
#'
#' Additionally, the function generates the following plots:
#' \itemize{
#'   \item Leverage vs. Estimated Probability
#'   \item Change in Pearson chi-squared vs. Estimated Probability
#'   \item Cook's Distance vs. Estimated Probability
#'   \item Change in Pearson chi-squared and Change in Cook's Distance vs Pair ID
#' }
#'
#' @return A list containing:
#' \item{data.results}{A data frame with the calculated residuals, leverage, and other diagnostic measures for each observation.}
#' \item{leverage}{A ggplot object displaying leverage vs. estimated probability.}
#' \item{change.Pearsonchi}{A ggplot object displaying the change in Pearson chi-squared vs. estimated probability.}
#' \item{cooks}{A ggplot object displaying Cook's distance vs. estimated probability.}
#' \item{m.Pearsonchi}{A ggplot object displaying the total pairwise change in Pearson chi-squared.}
#' \item{mcooks}{A ggplot object displaying the total pairwise change in Cook's distance.}
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 7
#'
#' # Convert 'fracture' to binary (0 = No, 1 = Yes)
#' glow11m$fracture <- ifelse(glow11m$fracture == "Yes", 1, 0)
#'
#' # Load required package
#' library(survival)
#'
#' # Fit a conditional logistic regression model
#' mod7.3 <- clogit(
#'   fracture ~ weight + bmi + priorfrac + momfrac + armassist + strata(pair),
#'   data = glow11m)
#'
#' # Run diagnostics for the conditional logistic model
#' residuals_clog(mod7.3, strata = glow11m$pair, y = glow11m$fracture)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc.
#'
#' @export

residuals_clog<-function(model,strata, y, plot=TRUE){
  #Calculate estimated N8 (estimated probability)
  x.design<-model.matrix(model)
  coeff<-model$coefficients

  ## Checking inputs
  if (!inherits(model, "clogit")) stop("The input model must be a 'clogit' object from the survival package.")
  if (length(strata) != nrow(x.design)) stop("Length of 'strata' must match the number of rows in model matrix.")
  if (length(y) != nrow(x.design)) stop("Length of response variable 'y' must match the number of rows in model matrix.")
  if (any(is.na(y))) stop("The outcome 'y' contains NA values. Please remove or impute them before proceeding.")
  if (is.factor(y)) y <- as.numeric(as.character(y))
  if (is.character(y)) y <- as.numeric(y)

  ##Checking the outcome is binary (0/1)
  unique_vals <- unique(y)
  if (!all(unique_vals %in% c(0, 1))) {
    stop("The outcome 'y' must be binary (only 0 and 1). Found values: ", paste(unique_vals, collapse = ", "))
  }

  ## Removing possible NA's in the coefficients
  if (anyNA(coeff)) {
    dropped_vars <- names(coeff[is.na(coeff)])
    warning("The following variables were dropped from the model (likely due to collinearity) and will be excluded from calculations: ",
            paste(dropped_vars, collapse = ", "))
  }

  non_na_vars <- names(coeff[!is.na(coeff)])
  x.design <- x.design[, non_na_vars, drop = FALSE]
  coeff <- coeff[non_na_vars]

  data.results <- data.frame(x.design, y = y, pair = strata, check.names = FALSE) #Strata needs to be in the same order than the datapoints wthat were used in the model

  ##Calculating the estimated probability that the subject is the case among the m subjects
  num<-exp(x.design %*% coeff)
  denom<-ave(num, strata, FUN = sum)
  data.results$theta<-num/denom


  ######Pearson residuals
  data.results$pearson<-(data.results$y-data.results$theta)/sqrt(data.results$theta)

  ######Calculating leverage h
  variables<-names(coeff)
  variables<-variables[variables != 'strata']
  x_centered <- data.results %>%
    group_by(pair) %>%
    mutate(across(all_of(variables), ~ . - sum(. * theta), .names = "centered_{.col}")) %>%
    ungroup() %>%
    dplyr::select(starts_with("centered_"))
  x_centered<-as.matrix(x_centered)

  u_j<-c(data.results$theta)
  U<-diag(u_j)

  inv<-solve(t(x_centered) %*% U %*% x_centered)

  data.results$leverage <- data.results$theta * rowSums((x_centered %*% inv) * x_centered)

  #Standardized Pearson residuals
  data.results$spearson <- (data.results$y - data.results$theta) / sqrt(data.results$theta * (1 - data.results$leverage))

  ################## Delta chi calculation
  data.results$deltachi<-data.results$spearson^2

  ################## Delta BETA calculation
  data.results$deltabeta<-data.results$deltachi*(data.results$leverage/(1-data.results$leverage))

  ###PLOTS
  if(plot==TRUE){
    leverage=ggplot()+
      geom_point(data=data.results, aes(x=theta, y=leverage, shape = factor(y), color = factor(y)), size=3)+
      scale_color_manual(values = c('black', 'purple'), labels= c('control', 'case'))+
      scale_shape_manual(values = c(1,2), labels= c('control', 'case'))+
      theme_classic()+
      labs(x='Estimated probability', y='Leverage')+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))+
      theme(legend.title=element_blank(), legend.position = 'top')

    change.Pearsonchi=ggplot()+
      geom_point(data=data.results, aes(x=theta, y=deltachi, shape = factor(y), color = factor(y)), size=3)+
      scale_color_manual(values = c('black', 'purple'), labels= c('control', 'case'))+
      scale_shape_manual(values = c(1,2), labels= c('control', 'case'))+
      theme_classic()+
      labs(x='Estimated probability', y='Change in Pearson chi-squared')+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))+
      theme(legend.title=element_blank(), legend.position = 'top')

    cooks=ggplot()+
      geom_point(data=data.results, aes(x=theta, y=deltabeta, shape = factor(y), color = factor(y)), size=3)+
      scale_color_manual(values = c('black', 'purple'), labels= c('control', 'case'))+
      scale_shape_manual(values = c(1,2), labels= c('control', 'case'))+
      theme_classic()+
      labs(x='Estimated probability', y='Cooks distance')+
      theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))+
      theme(legend.title=element_blank(), legend.position = 'top')

    ## For matched data
    df.match<-data.results %>% group_by(pair) %>% dplyr::summarise(sum_deltachi=sum(deltachi), sum_deltabeta=sum(deltabeta))

    m.Pearsonchi=ggplot()+
      geom_point(data=df.match, aes(x=pair, y=sum_deltachi))+
      theme_classic()+
      labs(x='Pair ID', y='Pairwise total change in Pearson chi-squared')

    mcooks=ggplot()+
      geom_point(data=df.match, aes(x=pair, y=sum_deltabeta))+
      theme_classic()+
      labs(x='Pair ID', y='Pairwise total change in Cooks distance')

    return(list(data.results=data.results,leverage=leverage,change.Pearsonchi=change.Pearsonchi,cooks=cooks,m.Pearsonchi=m.Pearsonchi,mcooks=mcooks))

  }else{
    return(data.results)
  }

}

