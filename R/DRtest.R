#' @title Deviance Residuals Test (HL Test)
#'
#' @description
#' The function performs the Hosmer-Lemeshow (HL) goodness-of-fit test, which is used to evaluate
#' how well a logistic regression model fits the data by comparing observed and expected frequencies
#' across different groups based on predicted probabilities. This function calculates the test using
#' either the model's predictions or user-supplied predictions. It divides the data into `g` groups,
#' based on the predicted values, and calculates the chi-squared statistic to assess the fit.
#'
#' @param model A fitted logistic regression model from \code{glm()} with \code{family = binomial}. If model is supplied, the yvar and yhatvar arguments will not be used.
#' @param yvar A vector of observed response values. Required if \code{model} is \code{NULL}.
#' @param yhatvar A vector of predicted probabilities. Required if \code{model} is \code{NULL}.
#' @param g The number of groups (quantiles) to divide the data into for the Hosmer-Lemeshow test. Default is 10.
#'
#' @return A list containing the following components:
#' \item{results}{A data frame with observed and expected values, and total counts for each group. The data frame includes columns for the observed and expected counts for both the outcome (1 and 0) and total counts.}
#' \item{chisq}{The chi-squared statistic used to assess the fit of the model.}
#' \item{df}{The degrees of freedom for the chi-squared test.}
#' \item{p.value}{The p-value of the chi-squared test, which assesses the overall goodness of fit.}
#' \item{groups}{The number of groups (quantiles) used in the test.}
#'
#' @details
#' The Hosmer-Lemeshow test compares the observed and expected frequencies of events across different quantiles
#' of predicted probabilities. The test statistic follows a chi-squared distribution.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Table 5.2
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
#' # Perform Hosmer-Lemeshow test with default 10 groups
#' DRtest(model.int)
#'
#' @export
#'
#'
DRtest <- function(model=NULL,yvar=NULL, yhatvar=NULL, g=10) {
  # This function is based on the original function `HLtest` from the package vcdExtra (GPL >= 2)
  # Original author: Michael Friendly.
  # Modifications: Extended to include additional features

  ## Checking inputs
  if (!is.null(model)) {
    if (is.null(model$y)){
      stop("The model must be fitted with y=TRUE.")
    }else{
      y<-model$y
      estimated.y<-predict(model, type='response')
    }
  }else{
    if (is.null(yvar) || is.null(yhatvar)) {
      stop("If 'model' is not provided, both 'yvar' and 'yhatvar' must be supplied.")
    }else{
      y<-yvar
      estimated.y<-yhatvar
    }
    }

  if (!all(y %in% c(0, 1))) stop("y must contain only 0s and 1s.")
  if (!is.numeric(estimated.y)) stop("estimated.y must be numeric.")

  esty.cut <- cut(estimated.y, breaks = quantile(estimated.y, probs=seq(0, 1, 1/g)), include.lowest=TRUE)
  obs <- xtabs(cbind(1 - y, y) ~ esty.cut)
  exp <- xtabs(cbind(1 - estimated.y, estimated.y) ~ esty.cut)
  n<- as.numeric(apply(obs, 1, sum))

  chi <- (obs - exp)/sqrt(exp)
  results <- data.frame(cut=dimnames(obs)$esty.cut,
                        obs.1=as.numeric(as.character(n-obs[,1])),
                        exp.1=round(as.numeric(as.character(n-exp[,1])),2),
                        obs.0=as.numeric(as.character(obs[,1])),
                        exp.0=round(as.numeric(as.character(exp[,1])),2),
                        total= n
  )

  rownames(results) <- 1:g
  chisq = sum(chi^2)
  p = 1 - pchisq(chisq, g-2 )
  result <- list(results=results, chisq=chisq, df=g, p.value=p, groups=g)
  return(result)
}

