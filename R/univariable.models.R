#' @title Univariable Logistic Regression Summary Table
#'
#' @description
#' Fits univariable logistic regression models `glm` function for a set of predictors and summarizes coefficients,
#' standard errors, likelihood ratio test statistics, and optionally odds ratios with confidence intervals.
#'
#' @import dplyr
#' @import lmtest
#'
#' @param data A data frame containing the outcome and predictor variables.
#' @param yval A string indicating the name of the binary outcome variable.
#' @param xval A character vector with the names of the explanatory variables to be evaluated univariately.
#' @param OR Logical; if TRUE, odds ratios and their confidence intervals are returned. Default to FALSE.
#' @param inc.or A numeric vector of scaling factors to be applied to each coefficient to obtain odds ratios.
#' @param confidence.level The confidence level to use for interval estimation of odds ratios. Defaults to 0.95.
#'
#' @return A data frame with coefficients from the univariable models, standard errors, p-values, and, if requested, odds ratios with lower and upper confidence limits.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 4, Table 4.7
#'
#' # Define variables to evaluate
#' unvariables <- c(
#'   'age','weight','height', 'bmi', 'priorfrac', 'premeno', 'momfrac',
#'   'armassist','smoke', 'raterisk')
#'
#' # Define value ranges used to interpret odds ratios (Optional)
#' val.pe <- c(5, 5, 10, 5, 1, 1, 1, 1, 1, 1, 1)
#'
#' # Run univariable conditional logistic regressions
#' univariable.models(glow500, yval = 'fracture', xval = unvariables, OR = TRUE, inc.or = val.pe)
#'
#' @export

univariable.models<-function(data, yval, xval, OR=FALSE, inc.or=NULL, confidence.level=0.95){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (OR && (is.null(inc.or) )) stop("When OR = TRUE, 'inc.or' must be provided.")
  missing_vars <- setdiff(c(yval, xval), colnames(data))
  if (length(missing_vars) > 0) stop(paste("These variables are missing in 'data':", paste(missing_vars, collapse = ", ")))

  level<-(1-confidence.level)/2
  z.val<-qnorm(1-level)
  n<-length(xval)
  coeff<-c()
  se<-c()
  G<-c()
  p_val<-c()
  estim.OR<-c()
  low.limit<-c()
  up.limit<-c()
  xval_expanded <- c()

  for(var.fit in xval){
    data.model<-data %>% dplyr::select(all_of(yval),all_of(var.fit)) %>% dplyr::rename(y = all_of(yval))
    #Model
    model.constant<-glm(y ~ 1, family = binomial, data = data.model)
    model.univariate <- glm(y ~ ., family = binomial, data = data.model)

    var.model<-names(coef(model.univariate))[-1]
    xval_expanded <- c(xval_expanded, var.model)

    for(j in var.model){
      coeff[j]<-data.frame(model.univariate$coefficients)[j,]
      se[j]<-data.frame(summary(model.univariate)$coefficients[,2])[j,]
      G[j]<-lmtest::lrtest(model.univariate,model.constant)$Chisq[2]
      p_val[j]<-lmtest::lrtest(model.univariate,model.constant)$`Pr(>Chisq)`[2]
    }

  }

  #OR
  i=1
  if(OR==FALSE){
    table<-round(data.frame(coeff,se,G,p_val),3)
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
    table<-data.frame(coeff,se,estim.OR,low.limit,up.limit,G,p_val)
  }


  return(table)
}





