#' @title plot() Method for Objects of Class 'trme'
#' @param object Objects of Class 'trme'
#' @param ... additional arguments
#' @return a histogram plots.
#' @export
#'
#'
plot.trme<- function(object,...)
{
    ##self defined plot function for "trme" class

    ##remove NA
    data_naomit=na.omit(object$data)

    ##draw hist plots
    ##cannot store the label
    par(mfrow=c(1,2))
    p1 <- hist(object$fit_ps[data_naomit$A==0],main="Fitted PS Values",
                           xlab="non-exposure group",ylab="density",freq=F)

    p2 <- hist(object$fit_ps[data_naomit$A==1],main="Fitted PS Values",
                       xlab="exposure group",ylab="density",freq=F)

    par(mfrow=c(1,1))

    ##output
    out=list(p1,p2)
    invisible(out)

}


