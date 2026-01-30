#' @title Check Significance of Excluded Variables
#'
#' @description
#' Obtain summary statistics, including Wald test z-values and p-values, for coefficients
#' of variables added one at a time to an existing logistic regression model.
#' Supports both standard and conditional logistic regression.
#'
#' @import dplyr
#' @import survival
#'
#' @param data A data frame containing the outcome, predictors, and optionally a stratification variable.
#' @param yval A string naming the binary outcome variable.
#' @param xpre A character vector of variables included in the preliminary model.
#' @param xcheck A character vector of variables to be tested for significance when added individually to the preliminary model.
#' @param strata (Optional) A string naming the stratification variable. If provided, conditional logistic regression is used via \code{clogit()}.
#'
#' @return A matrix containing the coefficient estimates, standard errors, z-values,
#' and p-values for each variable in \code{xcheck} when added to the preliminary model.
#'
#' @details For each variable in \code{xcheck}, this function fits a model that includes the variable
#' alongside \code{xpre}. It extracts and returns the coefficient summary for each added variable
#' to assess statistical significance. When \code{strata} is provided, a conditional logistic regression
#' model is used instead of standard logistic regression.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 4
#'
#' # Variables selected for the full model
#' preliminar <- c('age', 'height', 'priorfrac', 'momfrac', 'armassist', 'raterisk')
#'
#' # Variables to test
#' excluded <- c('weight', 'bmi', 'premeno', 'smoke')
#'
#' # Assess whether any excluded variables become significant when added to the preliminary model
#' check_coef_significant(data = glow500, yval = 'fracture', xpre = preliminar, xcheck = excluded)
#'
#' @export

check_coef_significant<-function(data, yval, xpre, xcheck, strata = NULL){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.character(yval) || length(yval) != 1) stop("`yval` must be a character string naming the outcome variable.")
  if (!is.character(xpre) || length(xpre) == 0) stop("`xpre` must be a non-empty character vector.")
  if (!is.character(xcheck) || length(xcheck) == 0) stop("`xcheck` must be a non-empty character vector.")

  required_vars <- c(yval, xpre, xcheck)
  if (!all(required_vars %in% names(data))) {
    missing_vars <- required_vars[!required_vars %in% names(data)]
    stop("The following variables are missing from `data`: ", paste(missing_vars, collapse = ", "))
  }

  # Warn if any xcheck vars already in xpre
  overlap_vars <- intersect(xpre, xcheck)
  if (length(overlap_vars) > 0) {
    warning("Variables present in both `xpre` and `xcheck`: ", paste(overlap_vars, collapse = ", "))
  }

  if (!is.null(strata)) {
    if (!(strata %in% names(data))) {
      stop(paste0("Stratification variable '", strata, "' not found in `data`."))
    }
  }

  ## Calculating the significance in coefficients
  check_coeff<-NULL
  coeff<-c()

  for (var in xcheck) {
    xpc <- c(xpre, var)

    if (is.null(strata)) {
      # Standard logistic regression
      data.pc<-data %>% dplyr::select(all_of(yval),all_of(xpc)) %>% dplyr::rename(y = all_of(yval))
      model.pc<-glm(y ~ ., family = binomial, data = data.pc)
      coefficients <- names(coef(model.pc))
      coefficients<- coefficients[grepl(var, coefficients)]

      for(i in coefficients){
        check_coeff<-cbind(check_coeff,summary(model.pc)[["coefficients"]][i,])
        coeff[i]<-i
      }
    } else {
      # Conditional logistic regression
      data.pc <- data %>% dplyr::select(all_of(yval), all_of(xpc), all_of(strata)) %>% dplyr::rename(y = all_of(yval))
      f_pc <- as.formula(paste("as.numeric(y) ~", paste(xpc, collapse = " + "), "+ strata(", strata, ")"))
      model.pc <- clogit(f_pc, data = data.pc)
      coefficients <- names(coef(model.pc))
      coefficients <- coefficients[grepl(var, coefficients)]

      for (i in coefficients) {
        check_coeff <- cbind(check_coeff, summary(model.pc)[["coefficients"]][i, ])
        coeff[i] <- i
      }
    }
  }
  colnames(check_coeff)<-coeff

  return(check_coeff)
}

