#' @title Univariable Conditional Logistic Regression Models
#'
#' @description
#' This function fits separate univariable conditional logistic regression models for a specified outcome and a
#' list of explanatory variables, using the `clogit` function from the `survival` package. It returns a summary
#' table with coefficients, standard errors, p-values, and optionally odds ratios with confidence intervals.
#' The strata vector must be in the same order as the observations in the \code{data} frame.
#'
#' @import survival
#' @import dplyr
#'
#' @details
#' Each predictor in `xval` is fit in a separate model with the outcome and stratification variable.
#' When `OR = TRUE`, the user must supply `inc.or`, a vector of increments matching the number of coefficients
#' across all univariable models (i.e., including levels of factors if applicable). The function assumes that
#' the outcome is binary and numeric (0/1); it will attempt to coerce it if not.
#'
#' @param data A data frame containing the outcome, predictors, and strata variables
#' @param yval A string indicating the name of the binary outcome variable.
#' @param xval A character vector of predictor variable names for which univariable models will be fit.
#' @param strata A character string specifying the name of the stratification variable (e.g., matching ID).
#' @param OR Logical; if TRUE, odds ratios and their confidence intervals are returned. Default to FALSE.
#' @param inc.or A numeric vector specifying the unit increment to apply when calculating odds ratios for each coefficient.
#'   Required if `OR = TRUE`. The length of `inc.or` must match the number of coefficients produced across all univariable models
#'   (including dummy variables for factor levels, if applicable). It must align with the order of `xval`: for each variable,
#'   include one value per resulting coefficient.
#' @param confidence.level The confidence level to use for interval estimation of odds ratios. Defaults to 0.95.
#'
#' @return A data frame with coefficients from the univariable models, standard errors, p-values, and, if requested, odds ratios with lower and upper confidence limits.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 7, Table 7.1
#'
#'  # Convert 'fracture' to binary (0 = No, 1 = Yes)
#' glow11m$fracture <- ifelse(glow11m$fracture == "Yes", 1, 0)
#'
#' # Define variables to evaluate
#' unvariables <- c(
#'   "height", "weight", "bmi", "priorfrac", "premeno", "momfrac",
#'   "armassist", "smoke", "raterisk")
#'
#' # Define value ranges used to interpret odds ratios (Optional)
#' val.pe <- c(10, 10, 3, 1, 1, 1, 1, 1, 1, 1)
#'
#' # Run univariable conditional logistic regressions
#' univariable.clogmodels(glow11m, yval = 'fracture', xval = unvariables,
#'                         strata = 'pair', OR = TRUE, inc.or = val.pe)
#'
#' @export

univariable.clogmodels<-function(data, yval, xval, strata, OR=FALSE, inc.or=NULL, confidence.level=0.95){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!all(c(yval, xval, strata) %in% colnames(data))) stop("yval and xval must exist in 'data'.")
  if (OR && (is.null(inc.or) )) {
    stop("When OR = TRUE, 'inc.or' must be provided.")
  }

  level<-(1-confidence.level)/2
  z.val<-qnorm(1-level)
  n<-length(xval)
  coeff<-c()
  p.val<-c()
  se<-c()
  estim.OR<-c()
  low.limit<-c()
  up.limit<-c()
  xval_expanded <- c()

  for(var.fit in xval){
    data.model<-data %>% dplyr::select(all_of(yval),all_of(var.fit)) %>% dplyr::rename(y = all_of(yval))
    data.model$strata<-data[[strata]]

    if(!inherits(data.model$y, "numeric")){
      data.model$y<-as.numeric(data.model$y)
    }

    #Model
    model.constant<-clogit(y ~ 1+strata(strata), data = data.model)
    model.univariate <- clogit(y ~ .+strata(strata), data = data.model)

    var.model<-names(coef(model.univariate))
    var.model<-var.model[var.model != 'strata']
    xval_expanded <- c(xval_expanded, var.model)

    for(j in var.model){
      coeff[j]<-data.frame(model.univariate$coefficients)[j,]
      se[j]<-data.frame(summary(model.univariate)$coefficients[,3])[j,]
      p.val[j]<-data.frame(summary(model.univariate)$coefficients[,5])[j,]
    }

  }

  #OR
  i=1
  if(OR==FALSE){
    table<-round(data.frame(coeff,se,p.val),3)
  }else{
    ##Checking vector
    if (length(inc.or) != length(xval_expanded)) {
      stop(paste0(
        "Length of 'inc.or' (", length(inc.or), ") must match the number of variables including all levels in 'xval_expanded' (", length(xval_expanded), "). Please provide one multiplier per coefficient." ))
    }

    for(var.fit in xval_expanded){
      estim.OR[var.fit]<-exp(inc.or[i]*coeff[var.fit])
      low.limit[var.fit]<-exp(inc.or[i]*coeff[var.fit]-z.val*inc.or[i]*se[var.fit])
      up.limit[var.fit]<-exp(inc.or[i]*coeff[var.fit]+z.val*inc.or[i]*se[var.fit])
      i=i+1
    }
    table<-data.frame(coeff,se,p.val,estim.OR,low.limit,up.limit)
  }


  return(table)
}
