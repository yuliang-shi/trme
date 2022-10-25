#' @title plot() Method for Objects of Class 'trmd'
#' @param object Objects of Class 'trmd'
#' @param ... additional arguments
#' @return a histogram plots.
#' @export
#'
plot.trmd<- function(object,...)
{
    ##self defined plot function for "trmd" class

    ##remove NA
    df_mar_naomit=na.omit(object$df_mar)

    ##draw hist plots
    ##cannot store the label
    par(mfrow=c(1,2))
    p1 <- hist(object$fit_ps_all[df_mar_naomit$A==0],main="Fitted PS Values",
                           xlab="non-exposure group",ylab="density",freq=F,...)

    p2 <- hist(object$fit_ps_all[df_mar_naomit$A==1],main="Fitted PS Values",
                       xlab="exposure group",ylab="density",freq=F,...)

    par(mfrow=c(1,1))

    ##output
    out=list(p1,p2)
    invisible(out)

}


