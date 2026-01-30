#' @title Osius and Rojek Goodness-of-Fit Test for Logistic Regression
#'
#' @description
#' Applies the Osius and Rojek normal approximation to assess the overall goodness-of-fit
#' of a logistic regression model. This includes:
#' - A normal approximation to the distribution of the Pearson chi-square statistic.
#' - A normal approximation to the distribution of the sum-of-squares statistic.
#'
#' @import dplyr
#'
#' @param model A logistic regression model fitted with \code{glm()}.
#'
#' @return A list containing the following measures:
#'  \item{z_chisq}{The standardized Z-statistic based on the Pearson chi-square approximation.}
#'  \item{p_value}{The two-sided p-value associated with \code{z_chisq}.}
#'  \item{z_s}{The standardized Z-statistic based on the sum-of-squares approximation.}
#'  \item{p_value_S}{The two-sided p-value associated with \code{z_s}.}
#' @details This function implements the Osius and Rojek test as described in Hosmer et al. (2013)
#' for assessing the goodness-of-fit of logistic regression models.
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
#' # Apply the Osius and Rojek test for goodness-of-fit
#' osius_rojek(model.int)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc.
#'
#' @seealso \code{\link{cov.patterns}}
#' @export


osius_rojek<-function(model){
  ## Checking inputs
  if (!inherits(model, "glm") || family(model)$family != "binomial") stop("Input must be a logistic regression model of class 'glm' with binomial family.")

  #Covariate patterns dataset
  X.cv<-cov.patterns(model)

  intercept_name <- grep("Intercept", colnames(X.cv), value = TRUE)
  X.design <- X.cv %>% dplyr::select(-all_of(c('y_j', 'm', 'est.prob', intercept_name)))
  n=as.numeric(dim(X.design)[1])
  p=as.numeric(dim(X.design)[2])

  ########################## Normal approximation to the distribution of the Pearson chi-square statistic
  ##Step 2
  v_j<- X.cv$m*X.cv$est.prob*(1-X.cv$est.prob)
  if (any(v_j == 0)) stop("Estimated probabilities of 0 or 1 lead to division by zero in variance.")

  #Step 3
  c_j=(1-2*X.cv$est.prob)/v_j

  #Step 4
  chisq<-sum((X.cv$y_j-X.cv$m*X.cv$est.prob)^2/v_j)

  #Step 5
  model.rss <- lm(c_j ~ .,data=X.design, weights = v_j)
  RSS <- sum(resid(model.rss)^2*v_j)

  #Step 6
  A=2*(n-sum(1/X.cv$m))

  #Step 7
  z_chisq=(chisq-(n-p-1))/sqrt(A+RSS)
  p_value <- 2 * (1 - pnorm(abs(z_chisq)))

  ########################## Normal approximation to the distribution of S (sum-of-squares)
  d_j=(1-2*X.cv$est.prob)
  S=sum((X.cv$y_j-X.cv$m*X.cv$est.prob)^2)

  modelS.rss <- lm(d_j ~ .,data=X.design, weights = v_j)
  RSS.s <- sum(resid(modelS.rss)^2*v_j)

  z_s=(S-sum(v_j))/sqrt(A+RSS.s)
  p_value_S <- 2 * (1 - pnorm(abs(z_s)))

  results<-list(z_chisq=z_chisq,p_value=p_value,z_s=z_s,p_value_S=p_value_S)

  return(results)
}

