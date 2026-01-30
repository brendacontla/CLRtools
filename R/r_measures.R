#' @title Model Fit Evaluation: R^2-like Measures for Logistic Regression Models with J=n
#'
#' @description
#' This function computes several R^2-like measures for logistic regression models, including the squared Pearson correlation
#' coefficient, pseudo R-squared, adjusted R-squared, shrinkage-adjusted R-squared, and log-likelihood-based R-squared.
#' These metrics, adapted from Hosmer et al. (2013), are specifically designed for datasets where the number of covariate patterns (J)
#' is equal to the number of data points (n), offering a comprehensive assessment of model fit in such contexts.
#'
#' @import lmtest
#'
#' @param model A fitted logistic regression model of class \code{glm}.
#'
#' @return A list containing the following measures:
#' \item{squared.Pearson.cc}{The squared Pearson correlation coefficient.}
#' \item{R_ss}{The linear regression-like measure or pseudo R-squared value.}
#' \item{ajusted_R}{The adjusted R-squared value, adjusted for the number of covariates and sample size.}
#' \item{adjust_shrink_R}{The shrinkage-adjusted R-squared value, which accounts for potential shrinkage in model estimates.}
#' \item{R_L}{The log-likelihood-based R-squared value.}
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Section 5.2.5
#'
#' # Recode 'raterisk' into a binary variable 'raterisk_cat'
#' glow500 <- dplyr::mutate(
#'   glow500,
#'   raterisk_cat = dplyr::case_when(
#'     raterisk %in% c("Less", "Same") ~ "C1",
#'     raterisk == "Greater" ~ "C2"
#'   )
#' )
#'
#' # Fit a multiple logistic regression model with interactions
#' model.int <- glm(
#'   fracture ~ age + height + priorfrac + momfrac + armassist +
#'     raterisk_cat + age * priorfrac + momfrac * armassist,
#'   family = binomial,
#'   data = glow500
#' )
#'
#' # Compute R-squared-like model fit statistics assuming J = n
#' r_measures(model.int)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc.
#'
#' @export

r_measures<-function(model){
  ## Checking inputs
  if (!inherits(model, "glm")) stop("The input 'model' must be a glm object.")
  if (family(model)$family != "binomial") stop("This function is intended for binomial logistic regression models.")
  if (is.null(model$y)) stop("The model must be fitted with y=TRUE.")

  ######################################################### Using data
  y_j<-model$y
  est.prob<-model[["fitted.values"]]

  y_bar <- mean(y_j)
  pi_bar <- mean(est.prob)

  ################################# Calculating squared Pearson correlation coefficient
  num <- sum((y_j - y_bar) * (est.prob - pi_bar))^2
  denom <- sum((y_j - y_bar)^2) * sum((est.prob - pi_bar)^2)
  r_c <- num / denom

  ###################################### Calculating linear regression like sum-of-squares measure or pseudo R_SQUARED
  num <- sum((y_j - est.prob)^2)
  denom <- sum((y_j - pi_bar)^2)
  R_SS <- 1-(num / denom)

  ###################################### Calculating R_SS_adj when sample size is small and the model contains many covariates
  n<-length(y_j)
  #Identify if there is intercept
  has_intercept <- attr(terms(model), "intercept") == 1
  p <- length(coef(model)) - if (has_intercept) 1 else 0
  if (n <= (p + 1)) warning("Sample size too small for adjusted R-squared calculation. Result may be unstable.")

  num <- sum((y_j - est.prob)^2)/(n-p-1)
  denom <- sum((y_j - pi_bar)^2)/(n-1)
  R_SS_adj<-1-(num/denom)

  ###################################### Calculating R_SS_shr for shrinkage in the estimates, a condition that is typically present when the model is fit with a small sample
  G<-lmtest::lrtest(model)$Chisq[2]
  if (G == 0){
    shr <- 0
  }else{
    shr <- (G - p) / G
  }

  R_SS_SHR<-R_SS*shr

  ###################################### Calculating  the log-likelihood-based R^2
  model.0<-glm(y_j ~ 1, family = binomial)
  L_0<-as.numeric(logLik(model.0))
  L_P<-as.numeric(logLik(model))

  R_L<-1-(L_P/L_0)


  return(list(squared.Pearson.cc=r_c, R_ss=R_SS, adjusted_R=R_SS_adj, adjust_shrink_R=R_SS_SHR, R_L=R_L))

}








