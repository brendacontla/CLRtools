#' @title Extract Unique Covariate Patterns from a Logistic Regression Model
#'
#' @description
#' Returns a summary of unique covariate patterns from a fitted logistic regression model, including
#' the number of observations per pattern and the estimated probability.
#'
#' @import dplyr
#'
#' @param model A fitted logistic regression model object, typically from \code{glm()} with \code{family = binomial}.
#'
#' @return A data frame where each row corresponds to a unique covariate pattern. The output includes:
#' \describe{
#'   \item{\code{y_j}}{The number of observed events (cases) for each pattern.}
#'   \item{\code{m}}{The number of observations sharing that pattern.}
#'   \item{\code{est.prob}}{The estimated probability from the fitted model for that pattern.}
#'   \item{\code{Covariate columns}}{The covariates used in the model, showing the pattern structure.}
#' }
#'
#' @details This function summarizes the unique covariate patterns in the data, capturing the combinations of
#' predictor values across observations. It calculates the frequency of each pattern, the number of events (cases),
#' and the estimated probability for each combination of covariates.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5
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
#' # Examine covariate patterns from the fitted model
#' X.cv <- cov.patterns(model.int)
#' head(X.cv, n=10)
#'
#' @export

cov.patterns<-function(model){
  ## Checking inputs
  if (!inherits(model, "glm")) stop("`model` must be a fitted `glm` object.")
  if (is.null(model$y)) stop("The model must be fitted with y=TRUE.")

  X<-model.matrix(model)
  covariates <- data.frame(y=model$y, est.prob=model[["fitted.values"]], X)
  covariates$id <- seq_len(nrow(covariates))
  X.names <- colnames(covariates)[-c(1:2, ncol(covariates))]

  # Create a data frame with unique combinations of covariates
  covariate_patterns <- covariates %>%
    group_by(across(all_of(X.names))) %>%
    dplyr::summarise(id= first(id),y_j = sum(y), m = n(), est.prob=first(est.prob), .groups = "drop")%>%
    arrange(id)  %>%
    dplyr::select(-id)

  return(covariate_patterns)
}

