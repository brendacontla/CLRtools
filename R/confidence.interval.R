#' @title Compute Wald-Based Confidence Intervals for Logit and Predicted Probability
#'
#' @description
#' Calculates Wald-based confidence intervals for the predicted logit values and estimated
#' probabilities for new data points, based on a fitted logistic regression model.
#'
#' @param model A fitted logistic regression model from \code{glm()} with \code{family = binomial}.
#' @param data A data frame containing the new observations for which predictions and confidence intervals are needed.
#' @param confidence.level Confidence level for the interval estimates. Defaults to 0.95.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{logit}}{A data frame with predicted logit values and their lower and upper confidence limits.}
#'   \item{\code{estimated.prob}}{A data frame with predicted probabilities and their lower and upper confidence limits.}
#' }
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 2
#'
#' # Fit logistic regression model with selected predictors, Table 2.3
#' mod2.3 <- glm(
#'   fracture ~ age + priorfrac + raterisk,
#'   family = binomial,
#'   data = glow500
#' )
#'
#' # Create a new data point for prediction
#' new.x <- data.frame(
#'   age = 65,
#'   priorfrac = "Yes",
#'   raterisk = "Same"
#' )
#'
#' # Compute the 95% confidence interval for the predicted probability
#' confidence.interval(mod2.3, data = new.x, confidence.level = 0.95)
#'
#' @export

confidence.interval<-function(model, data, confidence.level=0.95){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!inherits(model, "glm")) stop("`model` must be a fitted `glm` object.")
  if (!is.numeric(confidence.level) || confidence.level <= 0 || confidence.level >= 1) stop("`confidence.level` must be a number between 0 and 1.")

  level<-(1-confidence.level)/2
  z.val<-qnorm(1-level)

  #Obtaining logit prediction and SE
  prediction<-predict(model, newdata = data, se.fit=TRUE, type = 'link')[["fit"]]
  se.fit<-predict(model, newdata = data, se.fit=TRUE, type = 'link')[["se.fit"]]

  #Calculating the confidence interval of logit
  up.limit<-as.numeric(prediction+z.val*se.fit)
  lo.limit<-as.numeric(prediction-z.val*se.fit)
  confidence.logit<-data.frame(prediction=prediction,lower.limit=lo.limit,upper.limit=up.limit)

  #Calculating the estimated probability and confidence intervals
  est.prob<-exp(prediction)/(1+exp(prediction))
  up.limit.prob<-exp(up.limit)/(1+exp(up.limit))
  lo.limit.prob<-exp(lo.limit)/(1+exp(lo.limit))
  confidence.prob<-data.frame(est.probability=est.prob, lower.limit=lo.limit.prob, upper.limit=up.limit.prob)

  results<-list(logit=confidence.logit,estimated.prob=confidence.prob)

  return(results)
}
