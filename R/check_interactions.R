#' @title Check Pairwise Interactions in Logistic Regression
#'
#' @description
#' Evaluates all pairwise interactions between provided predictor variables in a logistic regression model.
#' It compares models that include individual interaction terms to a main-effects model using likelihood ratio tests.
#' Supports both standard logistic regression and conditional logistic regression (CLR) when stratification is provided.
#'
#' @import dplyr
#' @import lmtest
#' @import survival
#'
#' @param data A data frame containing the outcome variable and predictor variables.
#' @param variables A character vector of predictor variables to test for pairwise interactions.
#' @param yval A string naming the binary outcome variable.
#' @param rounding Optional integer specifying the number of decimal places to round the output. Defaults to \code{NULL} (no rounding).
#' @param strata A string naming the stratification variable for conditional logistic regression (CLR). Defaults to \code{NULL}, meaning standard logistic regression is used.
#'
#' @return A data frame showing the log-likelihood, likelihood ratio test statistic (G), and p-value
#' for each pairwise interaction tested. The first row corresponds to the main-effects model.
#'
#' @details This function fits a logistic regression model using the variables in \code{variables},
#' then iteratively fits additional models that include each pairwise interaction. The models are
#' compared using likelihood ratio tests via \code{lmtest::lrtest}. If \code{strata} is provided,
#' conditional logistic regression (CLR) is used instead of standard logistic regression.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 4, Table 4.14
#'
#' # Evaluate potential interaction terms among predictors.
#' # Variables included in the final main effects model to evaluate potential interactions
#' var.names <- c('age', 'height', 'priorfrac', 'momfrac', 'armassist', 'raterisk_cat')
#'
#' # Recode 'raterisk' into a binary categorical variable 'raterisk_cat'
#' glow500<-dplyr::mutate(
#'   glow500,
#'   raterisk_cat = dplyr::case_when(
#'     raterisk %in% c('Less', 'Same') ~ 'C1',
#'     raterisk == 'Greater' ~ 'C2'))
#'
#' # Run the interaction-checking procedure
#' check_interactions(data = glow500, variables = var.names, yval = 'fracture', rounding = 4)
#'
#' @export

check_interactions<-function(data, variables, yval,rounding=NULL, strata=NULL){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.character(yval) || length(yval) != 1) stop("`yval` must be a single character string.")
  if (!is.character(variables) || length(variables) < 2)  stop("`variables` must be a character vector of at least two variable names.")

  required_vars <- c(yval, variables)
  if (!all(required_vars %in% names(data))) {
    missing_vars <- required_vars[!required_vars %in% names(data)]
    stop("The following variables are missing from `data`: ", paste(missing_vars, collapse = ", "))
  }

  if (!is.null(strata)) {
    if (!(strata %in% names(data))) {
      stop(paste0("Stratification variable '", strata, "' not found in `data`."))
    }
  }

  ##Get the combinations of interactions
  combinations <- combn(variables, 2, simplify = FALSE)

  log.likh<-c()
  G<-c()
  p_value<-c()

  if (is.null(strata)) {
    # Standard logistic regression
    data.model<-data %>% dplyr::select(all_of(yval),all_of(variables)) %>% dplyr::rename(y = all_of(yval))
    model.multivariate <- glm(y ~ ., family = binomial, data = data.model)

    # Checking each interaction
    for (pair in combinations) {
      pair_str <- paste(pair[1], pair[2], sep = "*")
      model.multi.int <- glm(paste('y ~ .+',pair_str), family = binomial, data = data.model)

      log.likh[pair_str]<-logLik(model.multi.int)
      G[pair_str]<-lmtest::lrtest(model.multi.int,model.multivariate)[2,4]
      p_value[pair_str]<-lmtest::lrtest(model.multi.int,model.multivariate)[2,5]
    }
    base_loglik <- logLik(model.multivariate)

  } else {
    # Conditional logistic regression
    data.model <- data %>% dplyr::select(all_of(yval), all_of(variables), all_of(strata)) %>% dplyr::rename(y = all_of(yval))

    f_base <- as.formula(paste("as.numeric(y) ~", paste(variables, collapse = " + "), "+ strata(", strata, ")"))
    model.multivariate <- clogit(f_base, data = data.model)

    for (pair in combinations) {
      pair_str <- paste(pair[1], pair[2], sep = "*")
      interaction_formula <- as.formula(paste("as.numeric(y) ~", paste(c(variables, pair_str), collapse = " + "), "+ strata(", strata, ")"))
      model.multi.int <- clogit(interaction_formula, data = data.model)

      log.likh[pair_str] <- logLik(model.multi.int)
      G[pair_str] <- lmtest::lrtest(model.multi.int, model.multivariate)[2, 4]
      p_value[pair_str] <- lmtest::lrtest(model.multi.int, model.multivariate)[2, 5]
    }
    base_loglik <- logLik(model.multivariate)
  }

  if(!is.null(rounding)){
    df<-round(data.frame(log.likh,G,p_value),rounding)
    df<-rbind(c(round(logLik(model.multivariate), rounding),NA_real_,NA_real_),df)
  }else{
    df<-data.frame(log.likh,G,p_value)
    df<-rbind(c(logLik(model.multivariate),NA_real_,NA_real_),df)
  }

  return(df)

}




