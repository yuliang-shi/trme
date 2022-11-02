#' @title Summary() Method for Objects of Class 'trme'
#' @param object Objects of Class 'trme'
#' @return a data frame containing summarized estimates of parameters.
#' @export
#'
summary.trme <- function(object){

  ##print summarized data frame
  print(object$results)

}


#' @title print() Method for Objects of Class 'trme'
#' @param object Objects of Class 'trme'.
#' @return print the summarized results.
#' @export
print.trme <- function(object){

  ##self defined print for this class or list
  summary(object)
}


