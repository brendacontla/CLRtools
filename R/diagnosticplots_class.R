#' @title Diagnostic Plots for Model Discrimination
#'
#' @description
#' Generates four diagnostic plots to evaluate the discriminatory ability of a logistic regression model.
#' These plots follow the style of those presented by Hosmer et al. (2013) and are helpful for visually
#' assessing model performance.
#'
#' @import ggplot2
#' @import ggpubr
#'
#' @param model A fitted logistic regression model object of class \code{glm}, with a binary outcome.
#'
#' @return This function displays a 2x2 grid of diagnostic plots.
#'
#' @details This function creates the following diagnostic plots:
#' \enumerate{
#'   \item A jitter plot showing the distribution of the estimated probabilities against the observed outcome.
#'   \item An ROC curve displaying sensitivity versus 1-specificity.
#'   \item A histogram of estimated probabilities for observations with outcome = 0.
#'   \item A histogram of estimated probabilities for observations with outcome = 1.
#' }
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Section 5.2.4
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
#' # Generate classification diagnostic plots
#' diagnosticplots_class(model.int)
#'
#' @references
#' Hosmer, D. W., Lemeshow, S., & Sturdivant, R. X. (2013). *Applied Logistic Regression*
#' (3rd ed.). John Wiley & Sons, Inc.
#'
#' @seealso \code{\link{cutpoints}}
#' @export

diagnosticplots_class<-function(model){
  ## Checking inputs
  if (!inherits(model, c("glm", "clogit"))) stop("The function only supports binomial glm or clogit models.")
  if (is.null(model$y)) stop("The model must be fitted with y=TRUE.")

  y<-factor(model$y)
  est.prob<-predict(model, type='response')

  df<-data.frame(y,est.prob)
  table_cutoffs<-cutpoints(model,0,1,0.05,FALSE)

  p1<-ggplot()+
    geom_jitter(data=df, aes(x=est.prob, y=y),height = 0.1)+
    theme_classic()+
    labs(x='Estimated probability', y='Jittered outcome')

  p2<-ggplot(data = table_cutoffs, aes(x=Specificity1, y=Sensitivity))+
    geom_line()+
    geom_abline()+
    theme_classic()+
    labs(x='1-Specificity', y='Sensitivity')

  p3<-ggplot(data = subset(df, y %in% c('0')), aes(x=est.prob))+
    geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray", color = "black", alpha = 0.5) +
    theme_classic()+
    labs(y='Density', x='Estimated probability')+
    ggtitle('Outcome = 0') +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks=seq(0,1,0.1), limits=c(0,1))

  p4<-ggplot(data = subset(df, y %in% c('1')), aes(x=est.prob))+
    geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray", color = "black", alpha = 0.5) +
    theme_classic()+
    labs(y='Density', x='Estimated probability')+
    ggtitle('Outcome = 1') +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(breaks=seq(0,1,0.1), limits=c(0,1))

  print(ggarrange(p1,p3,p2,p4, ncol = 2, nrow = 2))

}

