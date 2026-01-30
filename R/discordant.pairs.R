#' @title Count Discordant Pairs in Matched Case-Control Data
#'
#' @description
#' This function identifies and counts discordant/non-informative pairs for categorical variables within matched case-control data.
#' It supports flexible matching designs, including 1:1 and 1:m, by evaluating each stratum individually.
#' Only strata with both case and control observations are included in the evaluation.
#'
#' @param data A data frame containing the matched dataset.
#' @param outcome A character string giving the name of the binary outcome variable.
#' @param strata A character string giving the name of the variable that identifies matched strata or sets.
#' @param variables Optional. A character vector of variable names to evaluate for discordant pairs. If `NULL`, all categorical (factor or character) variables in `data` are used.
#'
#'
#' @return A named numeric vector:
#' \itemize{
#'   \item Each element corresponds to the number of discordant strata for one variable.
#'   \item The final element, `total.valid.pairs`, gives the total number of strata that were eligible for evaluation.
#' }
#'
#' @examples
#' # Example from Hosmer et al., 2013
#' # Applied Logistic Regression (3rd ed.), Chapter 7, Table 7.1
#'
#' discordant.pairs(glow11m, outcome = 'fracture', strata = 'pair')
#'
#' @export

discordant.pairs<-function(data, outcome, strata, variables=NULL){
  ## Checking inputs
  if (!is.data.frame(data)) stop("`data` must be a data frame.")
  if (!is.character(outcome) || length(outcome) != 1) stop("`outcome` must be a single character string.")
  if (!is.character(strata) || length(strata) != 1) stop("`strata` must be a single character string.")
  if (!is.null(variables) && !all(variables %in% names(data))) stop("Some variables provided are not in the dataset.")
  if (!all(outcome %in% names(data))) stop("Outcome is not in the dataset.")
  if (!all(strata %in% names(data))) stop("Strata is not in the dataset.")

  n.strata<-unique(data[,strata])

  ## Optional for the user to give the variables ro review either a vector or only a name.
  if(is.null(variables)){
    cat.names<-names(data)[sapply(data, function(x) is.factor(x) || is.character(x))]
    if (length(cat.names) == 0) {
      stop("No categorical variables found in the dataset.")
    }
  }else{
    cat.names<-variables
  }

  result<-c()
  count.total <- 0

  ##Looping for each variable in a dataset
  for(variable in cat.names){
    ## Checking the strata for one variable
    count.pairs<-0
    for(pair in n.strata){
      pair.subset<-data[ data[[strata]] == pair, ]

      check.pairs <- length(unique(pair.subset[[variable]]))
      check.outcome <- length(unique(pair.subset[[outcome]]))

      ## Checking that the outcome has control and case per stratum
      if(check.outcome==1){
        warning(paste0("Stratum ", pair, " is missing either a case or a control and will be excluded from the analysis."))
      }else if(nrow(pair.subset) < 2){
        warning(paste0("Stratum ", pair, " has fewer than two observations and will be skipped."))
      }else{
        ##Getting the total number of valid pairs
        if (variable == cat.names[1]) {
          count.total <- count.total + 1
        }

        ## Getting the total number of discordant pairs
        if(check.pairs>1){
          count.pairs<-count.pairs+1
        }

      }
    }
    result[variable]<-count.pairs
  }

  result<- result[names(result) != outcome]
  result <- c(result, total.valid.pairs = count.total)

  return(result)
}

