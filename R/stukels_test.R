#' @title Stukel’s Test for Logistic Regression Model Fit
#'
#' @description
#' This function performs Stukel’s test to assess the goodness-of-fit of a logistic regression model,
#' as adapted from Hosmer et al. (2013) by using the likelihood ratio test .
#'
#' @import lmtest
#'
#' @param model A fitted logistic regression model of class \code{glm}.
#'
#' @return A list with the following components:
#' \item{G}{Likelihood ratio test statistic for added z1 and z2.}
#' \item{p_value}{P-value for the likelihood ratio test.}
#' \item{model_summary}{Summary of the logistic regression model including z1 and z2.}
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Section 5.2.2
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
#' # Apply the Stukels test for goodness-of-fit
#' stukels_test(model.int)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression* (3rd ed.).
#' John Wiley & Sons, Inc.
#'
#' @export

stukels_test<-function(model){
  ## Checking inputs
  if (!inherits(model, "glm") || family(model)$family != "binomial") stop("Input must be a logistic regression model of class 'glm' with binomial family.")

  # Step 1. Save the fitted values
  pi_hat <- fitted(model)
  if (any(pi_hat <= 0 | pi_hat >= 1)) {
    stop("Fitted probabilities contain 0 or 1, leading to infinite logits. Stukel's test cannot be performed.")
  }

  # Step 2
  g_hat <- log(pi_hat / (1 - pi_hat))

  # Step 3. Compute two new covariates z
  z1 <- 0.5 * g_hat^2 * as.numeric(pi_hat >= 0.5)
  z2 <- -0.5 * g_hat^2 * as.numeric(pi_hat < 0.5)

  # Step 4. Perform the Score test for the addition of z1 and/or z2 to the model
  X.design<-model.frame(model)
  X.design$z1<- z1
  X.design$z2<- z2

  model_original<- glm(formula(model), data = X.design, family = binomial)
  model_z1z2<- update(model_original, . ~ . + z1 + z2)

  G<-lmtest::lrtest(model, model_z1z2)$Chisq[2]
  p_value<-lmtest::lrtest(model, model_z1z2)$`Pr(>Chisq)`[2]

  results<-list(G=G,p_value=p_value, model_summary = summary(model_z1z2))

  return(results)
}

