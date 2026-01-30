#' @title Compute Odds Ratios for Logistic and Conditional Logistic Regression Models
#'
#' @description
#' Calculates odds ratios and confidence intervals for one or more coefficients in a fitted
#' logistic regression model (from \code{glm}) or conditional logistic regression model
#' (from \code{clogit}). Also supports scaling coefficients to meaningful units (e.g., per 5 years, 5 kg, or 5 cm)
#' for clearer interpretation.
#'
#' @param model A fitted model object from \code{glm()} (with \code{family = binomial}) or \code{survival::clogit()}.
#' @param variable Optional character vector specifying the variable(s) for which to compute odds ratios.  If `NULL` (default), checks all coefficients
#' @param c A numeric value specifying the increment(s) for the variables. Either a single (unnamed) numeric value
#'   applied to all selected variables, or a named numeric vector specifying a separate increment for each. The names of `c` must match
#'   the coefficient names in the model. Defaults to 1.
#' @param confidence.level Confidence level for the confidence interval. Defaults to 0.95.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{odds.ratio}}{Exponentiated coefficients (odds ratios).}
#'   \item{\code{lower.limit}}{Lower bounds of the confidence intervals.}
#'   \item{\code{upper.limit}}{Upper bounds of the confidence intervals.}
#' }
#'
#' @details Supports both logistic regression models fit with \code{glm()} and conditional
#' logistic regression models fit with \code{clogit()} from the \code{survival} package.
#' Confidence intervals are calculated on the log-odds scale and transformed back using \code{exp()}.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Table 5.18
#'
#' # Recode 'raterisk' into a binary categorical variable 'raterisk_cat'
#' glow500<-dplyr::mutate(
#'   glow500,
#'   raterisk_cat = dplyr::case_when(
#'     raterisk %in% c('Less', 'Same') ~ 'C1',
#'     raterisk == 'Greater' ~ 'C2'))
#'
#'  model.int <- glm(
#'    fracture ~ age + height + priorfrac + momfrac +
#'      armassist + raterisk_cat + age*priorfrac + momfrac*armassist,
#'   family = binomial, data = glow500)
#'
#'  # Specify variables and interpretation units for OR computation
#'  var.or<-c('raterisk_catC2','height')
#'  units.var<-c('raterisk_catC2'=1,'height'=5)
#'
#'  # Calculate and interpret adjusted odds ratios
#'  coeff.OR(model.int,variable = var.or, c = units.var)
#'
#' @export


coeff.OR<-function(model, variable=NULL, c=1, confidence.level=0.95){ #C increment of units in a numerical variable
  ## Checking inputs
  if (!inherits(model, c("glm", "clogit"))) stop("`model` must be a `glm` or `clogit` object.")
  if (!(length(c) == 1 || (!is.null(names(c)) && all(names(c) != "")))) stop("`c` must be either a single numeric value or a *named* numeric vector.")
  if (!is.numeric(confidence.level) || confidence.level <= 0 || confidence.level >= 1) stop("`confidence.level` must be a number between 0 and 1.")

  level<-(1-confidence.level)/2
  z.val<-qnorm(1-level)
  beta<-c()
  se.beta<-c()

  ## Getting the variables
  if(is.null(variable)){
    coefficients <- names(coef(model))
  }else{
    ncoefficients <- names(coef(model))
    data <- model.frame(model)
    coefficients<- match_variables(data, variable, ncoefficients)
  }

  if (length(coefficients) == 0) stop("No matching coefficients found for the specified variables.")

  ## Obtaining beta and SE from the model
  if(inherits(model, "clogit")){

    for(i in coefficients){
      if(length(c)==1){
        multiplier=c
      }else{
        multiplier=c[i]
      }
      beta[i]<-as.numeric(model$coefficients[i])*multiplier
      se.beta[i]<-summary(model)$coefficients[i, "se(coef)"]*multiplier
    }

  }else{

    for(i in coefficients){
      if(length(c)==1){
        multiplier=c
      }else{
        multiplier=c[i]
      }
      beta[i]<-coef(model)[i]*multiplier
      se.beta[i]<-summary(model)$coefficients[i,'Std. Error']*multiplier
    }

  }

  or<-exp(beta) # odds ratio
  lo.limit<-exp(beta-z.val*se.beta)
  up.limit<-exp(beta+z.val*se.beta)

  results<-list(odds.ratio=or, lower.limit=lo.limit, upper.limit=up.limit)

  return(results)

}
