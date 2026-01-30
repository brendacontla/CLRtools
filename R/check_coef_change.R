#' @title Assess Coefficient Change After Variable Removal
#'
#' @description
#' Computes the percentage change in logistic regression coefficients (\eqn{\Delta \hat{\beta}\%}) as each additional variable is introduced
#' to the model one at a time.
#' Supports both standard logistic regression (via \code{glm}) and conditional logistic regression (via \code{clogit}) when a stratification variable is provided.
#'
#' @import dplyr
#' @import survival
#'
#' @param data A data frame containing the outcome, predictors, and optional stratification variable.
#' @param yval A string naming the binary outcome variable.
#' @param xpre A character vector of variable names that are already selected for the model.
#' @param xcheck A character vector of variable names to be added one-by-one for comparison.
#' @param strata Optional; a string specifying the name of the stratification variable for conditional logistic regression.
#'
#' @return A data frame showing how the coefficients change when each variable in \code{xcheck}
#' is added to the model containing \code{xpre}.
#'
#' @details This function fits a logistic regression model using variables in \code{xpre},
#' and then adds each variable in \code{xcheck} one at a time to assess how the coefficients
#' of the model change with the delta beta hat percentages. Useful for evaluating confounding or additional variable contribution.
#' When \code{strata} is \code{NULL}, the function uses standard logistic regression (\code{glm} with binomial family).
#' When \code{strata} is specified, conditional logistic regression is used instead via \code{survival::clogit}.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 4
#'
#' # Variables selected to evaluate
#' preliminar <- c('age', 'height', 'priorfrac', 'momfrac', 'armassist')
#'
#' # Variable to evaluate for potential confounding
#' excluded <- c('raterisk')
#'
#' # Assess coefficient change after adding 'raterisk'
#' check_coef_change(data = glow500, yval = 'fracture', xpre = preliminar, xcheck = excluded)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc. The formulas for calculating residuals and diagnostics are adapted from this source.
#'
#' @seealso \code{\link{delta.coefficient}}
#' @export
#'


check_coef_change<-function(data, yval, xpre, xcheck, strata = NULL){
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

  ## Calculating the changes in coefficients
  if (is.null(strata)) {
    ##Model with all significant variables
    data.pre<-data %>% dplyr::select(all_of(yval),all_of(xpre)) %>% dplyr::rename(y = all_of(yval))
    model.pre<-glm(y ~ ., family = binomial, data = data.pre)

    coefficients <- names(coef(model.pre))
    check_coeff<-setNames(data.frame(matrix(ncol = length(xcheck), nrow = length(coefficients))), xcheck)
    rownames(check_coeff)<-coefficients

    for(var in xcheck){
      xpc<-c(xpre,var)
      data.pc<-data %>% dplyr::select(all_of(yval),all_of(xpc)) %>% dplyr::rename(y = all_of(yval))
      model.pc<-glm(y ~ ., family = binomial, data = data.pc)
      check_coeff[,var]<-delta.coefficient(model.pc,model.pre)
    }
  } else {
    ##Model with all significant variables
    data.pre<-data %>% dplyr::select(all_of(yval),all_of(xpre),all_of(strata)) %>% dplyr::rename(y = all_of(yval))
    f <- as.formula(paste("as.numeric(y) ~", paste(xpre, collapse = " + "), "+ strata(", strata, ")"))
    model.pre <- clogit(f, data = data.pre)

    coefficients <- names(coef(model.pre))
    check_coeff<-setNames(data.frame(matrix(ncol = length(xcheck), nrow = length(coefficients))), xcheck)
    rownames(check_coeff)<-coefficients

    for(var in xcheck){
      xpc<-c(xpre,var)
      data.pc<-data %>% dplyr::select(all_of(yval),all_of(xpc),all_of(strata)) %>% dplyr::rename(y = all_of(yval))
      f_pc <- as.formula(paste("as.numeric(y) ~", paste(xpc, collapse = " + "), "+ strata(", strata, ")"))
      model.pc <- clogit(f_pc, data = data.pc)
      check_coeff[,var]<-delta.coefficient(model.pc,model.pre)
    }
  }


  return(check_coeff)
}

