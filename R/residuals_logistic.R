#' @title Model Diagnostic for Logistic Regression Models
#'
#' @description
#' This function computes various residuals and diagnostic measures for logistic regression models,
#' including Pearson residuals, standardized Pearson residuals, deviance residuals, leverage, and delta statistics
#' (change in Pearson chi-squared, change in deviance, and change in Cook's distance).
#' The function also generates diagnostic plots for evaluating the model fit and identifying influential data points.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @param model A fitted logistic regression model of class \code{glm}.
#'
#' @details The residuals and diagnostic measures are computed using standard formulas for logistic regression
#' (adapted from Hosmer et al., 2013, *Applied Logistic Regression*, 3rd ed.). The following measures are computed:
#' \itemize{
#'   \item Pearson residuals
#'   \item Standardized Pearson residuals
#'   \item Deviance residuals
#'   \item Leverage
#'   \item Delta statistics: Delta beta (change in Cook's distance), Delta chi (change in Pearson chi-squared) and Delta deviance
#' }
#'
#' Additionally, the function generates the following plots:
#' \itemize{
#'   \item Leverage vs. Estimated Probability
#'   \item Change in Pearson chi-squared vs. Estimated Probability
#'   \item Change in Deviance vs. Estimated Probability
#'   \item Cook's Distance vs. Estimated Probability
#'   \item Change in Pearson chi-squared vs. Estimated Probability
#'   \item Change in Pearson chi-squared vs. Estimated Probability with size determinate by Cook's distance.
#' }
#'
#' @return A list containing:
#' \item{x.cv}{A data frame with the computed residuals and diagnostic measures for each observation, and their respective covariates.}
#' \item{leverage}{A ggplot object displaying leverage vs. estimated probability.}
#' \item{change.Pearsonchi}{A ggplot object displaying the change in Pearson chi-squared vs. estimated probability.}
#' \item{change.deviance}{A ggplot object displaying the change in deviance vs. estimated probability.}
#' \item{cooks}{A ggplot object displaying Cook's distance vs. estimated probability.}
#' \item{change.Pb}{A ggplot object displaying the change in Pearson chi-squared vs. Estimated Probability with size determinate by Cook's distance.}
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Section 5.3
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
#' # Run diagnostics for the logistic model
#' residuals_logistic(model.int)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc. The formulas for calculating residuals and diagnostics are adapted from this source.
#'
#' @seealso \code{\link{cov.patterns}}
#' @export


residuals_logistic<-function(model){
  ## Checking inputs
  if (!inherits(model, "glm")) stop("The model must be a 'glm' object.")
  if (family(model)$family != "binomial") stop("Only binomial glm models are supported.")

  #Covariate patterns dataset
  X.cv<-cov.patterns(model)

  X.design<-X.cv %>% dplyr::select(-all_of(c('y_j',"m", 'est.prob')))
  X.design<-as.matrix(X.design)
  n=as.numeric(dim(X.design)[1])

  ############ Calculating, v_j, h_j.
  v_j<- X.cv$m*X.cv$est.prob*(1-X.cv$est.prob)
  if (any(v_j == 0)) stop("Estimated probabilities of 0 or 1 lead to division by zero in variance.")

  V<-diag(v_j)
  inv<-solve(t(X.design) %*% V %*% X.design)
  b_j<-c()
  i=1
  while(i<=n){
    x_j<-as.matrix(X.design[i,])
    b_j[i]<-t(x_j) %*% inv %*% x_j
    i=i+1
  }
  h_j<-v_j*b_j  #Leverage

  ### Residuals
  pearson<-(X.cv$y_j-X.cv$m*X.cv$est.prob)/sqrt(v_j)

  spearson<-pearson/(sqrt(1-h_j))

  #Calculating the deviance residuals
  deviance<-c()
  for (i in 1:length(X.cv$y_j)) {
    yj <- X.cv$y_j[i]
    mj <- X.cv$m[i]
    pi_j_hat <- X.cv$est.prob[i]

    if (yj == 0) {
      # Deviance residual for yj = 0
      deviance[i] <- -sqrt(2 * mj * abs(log(1 - pi_j_hat)))
    } else if (yj == mj) {
      # Deviance residual for yj = mj
      deviance[i] <- sqrt(2 * mj * abs(log(pi_j_hat)))
    } else {
      # General case: deviance residual for 0 < yj < mj
      sign_term <- ifelse((yj - mj * pi_j_hat) >= 0, 1, -1)
      term1 <- yj * log(yj / (mj * pi_j_hat))
      term2 <- (mj - yj) * log((mj - yj) / (mj * (1 - pi_j_hat)))
      deviance[i] <- sign_term * sqrt(2 * (term1 + term2))
    }
  }

  ##Delta
  delta.beta<-(spearson^2*h_j)/(1-h_j)
  delta.chi<-spearson^2
  delta.deviance<-deviance^2/(1-h_j)

  X.cv<-data.frame(X.cv, leverage=h_j, pearson=pearson, spearson=spearson, deviance=deviance, delta.beta=delta.beta, delta.chi=delta.chi, delta.deviance=delta.deviance)


  ###PLOTS
    leverage<-ggplot()+
      geom_point( aes(x=X.cv$est.prob, y=X.cv$leverage))+
      theme_classic()+
      labs(x='Estimated probability', y='Leverage')+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))

    change.Pearsonchi<-ggplot()+
      geom_point( aes(x=X.cv$est.prob, y=X.cv$delta.chi))+
      theme_classic()+
      geom_hline(aes(yintercept = 4))+
      labs(x='Estimated probability', y='Change in Pearson chi-squared')+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))


    change.deviance<-ggplot()+
      geom_point( aes(x=X.cv$est.prob, y=X.cv$delta.deviance))+
      theme_classic()+
      labs(x='Estimated probability', y='Change in deviance')+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))

    cooks<-ggplot()+
      geom_point( aes(x=X.cv$est.prob, y=X.cv$delta.beta))+
      theme_classic()+
      labs(x='Estimated probability', y='Cooks distance')+
      theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))

    change.Pb<-ggplot()+
      geom_point( aes(x=X.cv$est.prob, y=X.cv$delta.chi, size=X.cv$delta.beta), shape=1)+
      scale_size_area(max_size = 23)+
      theme_classic()+
      labs(x='Estimated probability', y='Change in Pearson chi-squared')+
      theme(legend.position="none")+
      scale_x_continuous(breaks=seq(0,1,by=0.1), limits = c(0,1))

    return(list(x.cv=X.cv,leverage=leverage,change.Pearsonchi=change.Pearsonchi,change.deviance=change.deviance,cooks=cooks,change.Pb=change.Pb))

}







