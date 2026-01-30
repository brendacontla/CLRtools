#' @title Delta-beta hat percentage: Change in Coefficients when Adding a Variable
#'
#' @description
#' This function checks the percentage change in the coefficient of a variable when
#' it is added to a model. The comparison is made between two models: model1 (with
#' the new variable included) and model0 (without the new variable).
#'
#' @param model1 A fitted model (e.g., `glm`, `lm`) with the new variable included.
#' @param model0 A fitted model (e.g., `glm`, `lm`) without the new variable.
#' @param variable A vector of variable names to check. If `NULL` (default), checks all coefficients.
#'
#' @return A vector of percentage changes in the coefficients between model1 and model0.
#'
#' @details This function compares the coefficients of variables between two logistic regression models or conditional logistic regressions:
#' one that includes a variable of interest (\code{model1}) and one that excludes it (\code{model0}).
#' It calculates how much the coefficients of other variables change, expressed as a percentage.
#' The change is calculated as ((beta0 - beta1) / beta1) * 100, where beta0 is the coefficient
#' from model0 and beta1 is the coefficient from model1.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 3, Table 3.10
#'
#' mod3.1 <- glm(fracture ~ priorfrac, family = binomial, data = glow500)
#' mod3.2 <- glm(fracture ~ priorfrac + height, family = binomial, data = glow500)
#'
#' # delta-beta-coefficient
#' delta.coefficient(model1 = mod3.2, model0 = mod3.1, variable = 'priorfrac')
#'
#' @export

delta.coefficient<-function(model1, model0, variable=NULL){
  ## Checking inputs
  if (!inherits(model0, c("glm", "clogit")) || !inherits(model1, c("glm", "clogit"))) stop("Both `model0` and `model1` must be `glm` or `clogit` objects.")

  beta0<-c()
  beta1<-c()
  # Get coefficients
  if (is.null(variable)) {
    coefficients <- names(coef(model0))
  } else {
    ncoefficients <- names(coef(model0))
    data <- model.frame(model0)
    coefficients <- match_variables(data, variable, ncoefficients)
  }

  if (length(coefficients) == 0) stop("No matching coefficients found for the specified variables.")

  ## Obtaining the coefficient values
  for(i in coefficients){
    beta0[i]<-data.frame(model0$coefficients)[i,]
    beta1[i]<-data.frame(model1$coefficients)[i,]
  }

  dbc<-((beta0-beta1)/beta1)*100

  return(dbc)

}


