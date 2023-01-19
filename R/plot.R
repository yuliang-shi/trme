#' @title plot() Method for Objects of Class 'trme'
#' @param object Objects of Class 'trme'
#' @param ... additional arguments
#' @return histogram plots for fitted PS values and density plot for bootstrap estimate.
#' @export
#'
#'
plot.trme= function(object,...)
{
    ##self defined plot function for "trme" class
    ##... other auguments

    ##remove NA
    data_naomit=na.omit(object$data)

    ##draw hist plots for fitted PS values in two groups
    ##cannot store the label
    par(mfrow=c(1,2))
    p1 = hist(object$fit_ps[data_naomit$A==0],main="Fitted PS Values",
                           xlab="non-exposure group",ylab="density",freq=F)

    p2 = hist(object$fit_ps[data_naomit$A==1],main="Fitted PS Values",
                       xlab="exposure group",ylab="density",freq=F)


    ##density plot for bootstrap point estimate
    par(mfrow=c(1,1))
    p3=hist(object$boot_est,prob=T,xlab="bootstrap estimated values",ylab="density",
            main=paste0("Density Plot for Bootstrap Using ", object$method))
    lines(density(object$boot_est,kernel = "gaussian"),col="black",lwd=2,lty=2,...)

    ##output as invisible list
    out=list(p1,p2,p3)
    invisible(out)

}




