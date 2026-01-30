#' @title Table with Sensitivity and Specificity at Different Cutpoints
#'
#' @description
#' This function computes the sensitivity and specificity at various cutpoints for a given logistic regression model.
#' It generates a table summarizing the performance metrics (sensitivity, specificity) at different probability cutoffs
#' and optionally plots these metrics and the distribution of probabilities for each class.
#' This is useful for selecting an optimal threshold for classification.
#'
#' @import ggplot2
#' @import patchwork
#' @importFrom caret sensitivity specificity
#'
#' @param model A fitted logistic regression model (either `glm` or `clogit`).
#' @param cmin The minimum cutoff value for the predicted probabilities. Defaults to 0.
#' @param cmax The maximum cutoff value for the predicted probabilities. Defaults to 1.
#' @param byval The increment for cutpoints. Defaults to 0.05.
#' @param plot Logical value indicating whether to generate plots. Defaults to TRUE.
#'
#' @return A data frame containing cutpoints, sensitivity, specificity, and specificity complement for each cutoff. If
#' \code{plot = TRUE}, a ggplot2-based visualization is also printed,
#' showing sensitivity and specificity curves and the distribution of predicted probabilities
#' by outcome class, with the optimal cutoff (where sensitivity and specificity are closest)
#' indicated on the histogram.
#'
#' @details The function calculates sensitivity and specificity for a range of cutpoints from `cmin` to `cmax` with
#' a step size of `byval`. It then plots the relationship between sensitivity and specificity, as well as histograms
#' of estimated probabilities. The cutpoint with the smallest difference between sensitivity and
#' specificity is also marked on the histogram plots. This can aid in finding an optimal classification threshold.
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 5, Table 5.8
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
#' # Compute sensitivity and specificity at multiple cutpoints
#' cutpoints(model.int, cmin = 0.05, cmax = 0.75, byval = 0.05, plot = FALSE)
#'
#' @export
#'
cutpoints<-function(model, cmin=0, cmax=1, byval=0.05, plot=TRUE){
  ## Checking inputs
  if (!inherits(model, "glm")) stop("`model` must be a `glm` object.")
  if (family(model)$family != "binomial") stop("`model` must be a logistic regression (binomial family).")
  if (!is.numeric(cmin) || !is.numeric(cmax) || !is.numeric(byval)) stop("`cmin`, `cmax`, and `byval` must be numeric.")
  if (cmin < 0 || cmax > 1 || cmin >= cmax) stop("`cmin` and `cmax` must be between 0 and 1, with `cmin` < `cmax`.")
  if (is.null(model$y)) stop("The model must be fitted with y=TRUE.")


  cp<-seq(cmin,cmax, by=byval)
  n<-length(cp)
  summ<-data.frame(Cutpoint=cp, Sensitivity=replicate(0,n=n), Specificity=replicate(0,n=n), Specificity1=replicate(0,n=n))
  y<-factor(model$y)
  est.prob<-predict(model, type='response')

  i=1
  while(i<=n){
    cutpoint.i<-summ[i,1]
    y.class<-factor(ifelse(est.prob<=cutpoint.i,'0','1'))
    summ[i,2]<-caret::sensitivity(y.class,y, positive=levels(y)[2])*100
    summ[i,3]<-caret::specificity(y.class,y, negative=levels(y)[1])*100
    summ[i,4]<-100-caret::specificity(y.class,y, negative=levels(y)[1])*100
    i=i+1
  }

  if(plot==TRUE){
    df<-data.frame(y,est.prob)
    plot.data <- data.frame(summ,difference=abs(summ$Sensitivity - summ$Specificity))
    crossvalue<-plot.data [which.min(plot.data $difference), ]
    crossvalue<-crossvalue$Cutpoint

    p1<-ggplot()+
      geom_line(data=summ, aes(y=Sensitivity, x=Cutpoint, col='purple'))+
      geom_line(data=summ, aes(y=Specificity, x=Cutpoint, col='black'))+
      labs(y='Specificity/Sensitivity', x='Probability cutoff')+
      scale_color_manual(name=NULL,values=c('purple','black'),labels=c('Sensitivity','Specificity'))+
      theme_classic()+
      theme(legend.position = 'top')

    p2<-ggplot(data = subset(df, y %in% c('0')), aes(x=est.prob))+
      geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray", color = "black", alpha = 0.5) +
      geom_vline(xintercept = crossvalue, color='red', linewidth=1)+
      theme_classic()+
      labs(y='Density', x='Estimated probability')+
      ggtitle('Outcome = 0') +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks=seq(0,1,0.1), limits=c(0,1))

    p3<-ggplot(data =subset(df, y %in% c('1')), aes(x=est.prob))+
      geom_histogram(aes(y = after_stat(density)), bins = 25, fill = "gray", color = "black", alpha = 0.5) +
      geom_vline(xintercept = crossvalue, color='red', linewidth=1)+
      theme_classic()+
      labs(y='Density', x='Estimated probability')+
      ggtitle('Outcome = 1') +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks=seq(0,1,0.1), limits=c(0,1))

    layout <- p1 | (p2 / p3)
    print(layout)
  }


  return(summ)
}
