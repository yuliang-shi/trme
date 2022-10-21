#' @title summary() Method for Objects of Class 'trmd'
#' @param object Objects of Class 'trmd'
#' @return a list containing estimates of parameters.
#' @export
#'
summary.trmd <- function(object){

  ##print summarized data frame
  print(object$results)

}


#' @title print() Method for Objects of Class 'trmd'
#' @param object Objects of Class 'trmd'
#' @return invisible() with printing
#' @export
print.trmd <- function(object){

  ##self defined print for this class or list
  summary(object)
}


