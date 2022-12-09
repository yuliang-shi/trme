###IPW for missing exposure core function ######
#' @title IPW Estimation for Missing Exposure
#' @author Yuliang Shi
#' @description Estimate the causal effect of the exposure on the outcome when the exposure is MAR. Adjust for both missing and confounding issues via IPW methods. Provide estimated standard errors for inference purposes.
#'
#'
#' @param covs \code{\link{character}}  the required vector name of all confounders which
#'   will affect both exposure and outcome.  Any interaction term or non-linear
#'   form can be added after creating those variables based on the data.
#' @param Y \code{\link{character}} the required name of the outcome variable. The trme
#'   function currently only works on the binary outcome.
#' @param A \code{\link{character}} the required name of binary treatment or exposure
#'   variable.
#' @param data a required \code{\link{data.frame}} contains the variables in the models.
#' @param shrink_rate \code{\link{numeric}} shrinkage rate. By default, no shrinkage applies.
#'   In some cases, we can shrink weights larger than 99\% quantile into exact 99\% quantile to avoid
#'   extreme weights.
#' @param method \code{\link{character}} the method to be used. either "IPW-IPW" for traditional IPW; "IPW-DR" for IPW with DR estimator; or "IPW WEE" for IPW estimator using weighted estimating equations to avoid extreme weights.
#'
#' @param ci_alpha \code{\link{numeric}} 0\%-100\% percentage of confidence interval. By default,  95\% CI is present.
#' @param bootstrap a \code{\link{logical}} value either \code{TRUE} or \code{FALSE} for bootstrap method. If \code{TRUE}, provide statistical reference. If \code{FALSE}, only return point estimate.
#' @param B a \code{\link{numeric}} value for replication times of bootstrap.
#'
#'
#' @details The methods currently only work for the missing
#'   exposure with binary outcome and exposure variables. If either the
#'   covariates or outcome is also missing, please impute missing values via \code{\link{mice}}.
#'
#'   The basic assumption is that the exposure variable is
#'   missing at random (MAR), i.e. given all observed covariates and
#'   outcomes, the missingness (or missing indicator) is independent of the
#'   missing value itself. If the exposure is missing not at random (MNAR), try
#'   another method instead.
#'
#'   \code{method="IPW-IPW"} is for the traditional IPW estimator, which requires **both missingness and treatment models to be correct**.
#'
#'   \code{method="IPW-DR"} is a IPW with DR estimator, which requires **missingness model is correct and either treatment or outcome model is correct**.
#'
#'   \code{method="IPW-WEE"} is a IPW using WEE estimator (recommended), which keeps **IPW-DR properties and avoids some effects of extreme weights** in the finite samples, but it still keeps the same TR properties as the complex form. For more details, please review the reference paper.
#'
#' All IPW estimators utilize Bootstrap approach to estimate standard errors which can protect against
#'   the misspecification of model based on the simulation studies. By default, parallel computing will be applied to speed up computing process based on the operating system.
#'
#' @return Use \code{summary()} function to print out a data frame including summarized results.
#' \itemize{
#' \item{\code{Estimate: }}{estimated causal effect  (odds ratio) of exposure on the outcome. }
#' \item{\code{95\% CI: }}{95\% two-sided confidence interval.}
#' \item{\code{p.value: }}{p values for two-sided Wald test.}}
#'
#' In addition, other fitted values are also saved in the list.
#' \itemize{
#' \item{\code{fit_ps: }}{fitted propensity scores for all subjects, which are used to adjust for the confounding issue.}
#' \item{\code{miss_weights: }}{fitted inverse weights of missingness used to adjust for the missing issue.}
#' \item{\code{hist_ps_control, hist_ps_trt: } use \code{plot()} function to draw density plots for fitted propensity score between observed control and treatment groups.}
#' }
#'
#' @keywords regression, robust.
#'
#' @note For more details, please review \href{https://github.com/yuliang-shi/trme}{Yuliang's Github}.
#' For citation, please cite the package as **Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Triple Robust Estimator}. Statistics in Medicine.**
#'
#' @seealso \code{\link{summary.trme}}, \code{\link{print.trme}} for summarized result; \code{\link{plot.trme}} for histograms of fitted propensity score; \code{\link{covid19}} for description of real data set.
#'
#' @references Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Triple Robust Estimator}. Statistics in Medicine. Submitted (12/2022).
#'
#' Zhang, Z., Liu, W., Zhang, B., Tang, L., and Zhang, J. (2016). \emph{Causal inference with missing exposure information: Methods and applications to an obstetric study}. Statistical Methods in Medical Research 25, 2053–2066.
#'
#' Williamson, E. J., Forbes, A., and Wolfe, R. (2012). \emph{Doubly robust estimators of causal exposure effects with missing data in the outcome, exposure or a confounder}. Statistics in medicine 31, 4382–4400.
#'
#'
#' @examples
#' ########The first example for simulated data##########
#' require("trme")
#' set.seed(2000)
#' n = 1000 #sample size
#'
#' #####generate some continuous covariates
#' id = seq(1, n, by = 1)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = rnorm(n)
#'
#' #generate binary exposure from the PS model. x1,x2,x3 as confounders
#' a = -0.2 + 0.8 * x1 + 0.9 * x2 + 1 * x3
#' prob_a = 1 / (1 + exp(-a))
#' A = rbinom(n, 1, prob_a)
#'
#' #generate binary outcome from OR model
#' z = 0.9 + 1 *A + 0.9 * x1 + 0.6 * x2 + 0.5 * x3   #-1.5 is ok
#' Y = rbinom(n, 1, 1 / (1 + exp(-z)))      #' bernoulli response variable
#'
#' ##df: before remove missing data
#' df = data.frame(
#'   id = id,
#'   Y = Y,
#'   A = factor(A),
#'   x1 = x1,
#'   x2 = x2,
#'   x3 = x3,
#' prob_a = prob_a
#' )
#'
#' ##control the miss rate as 20\%
#' r_sim = -1 +0.4 * x1 +0.6 * x2 +0.8 * x3 - 0.9 * Y #add random error
#' r = rbinom(n, 1, 1 / (1 + exp(-r_sim))) #miss rate
#'
#' ##data: after include missing data
#' data = cbind(df, "r" = r)
#' data$A[which(data$r == 1)] = NA
#'
#' ##test in the simulated data
#' ##true causal (odds ratio) is 2.247 in this case
#' ##use recommended IPW-WEE as estimate
#' ipw_wee = ipwe(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   data = data,
#'   shrink_rate = 1,
#'   ci_alpha=0.95,
#'   method = "IPW-WEE",
#'   bootstrap=T,
#'   B=200
#' )
#'
#' ##print out results and plots
#' summary(ipw_wee)
#' plot(ipw_wee)
#'
#' ##use IPW-DR method
#' ipw_dr = ipwe(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   data = data,
#'   shrink_rate = 1,
#'   ci_alpha=0.95,
#'   method = "IPW-DR",
#'   bootstrap=T,
#'   B=200
#' )
#'
#' ##print out results
#' summary(ipw_dr)
#'
#'##use IPW-IPW method
#'ipw_ipw = ipwe(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   data = data,
#'   shrink_rate = 1,
#'   ci_alpha=0.95,
#'   method = "IPW-IPW",
#'   bootstrap=T,
#'   B=200
#' )
#'
#' ##print out results
#' summary(ipw_ipw)
#'
#'
#' ########The second example for real data##########
#' require("trme")
#' data(covid19)
#'
#' ##use IPW WEE estimator
#' ipw_wee=ipwe(covs = c("age","sex","diabetes"),Y="Y",A="CVD", data=covid19,
#'             shrink_rate = 1,ci_alpha=0.95,
#'             method="IPW-WEE",bootstrap=T,B=200)
#'
#' ##obtain estimate of causal effect and CI
#' summary(ipw_wee)
#'
#' ##use IPW DR method
#' ipw_dr=ipwe(covs = c("age","sex","diabetes"),Y="Y",A="CVD", data=covid19,
#'            shrink_rate = 1,ci_alpha=0.95,
#'            method="IPW-DR",bootstrap=T,B=200)
#'
#' ##use TR AIPW method
#' summary(ipw_dr)
#'
#' ##use traditional IPW IPW method
#' ipw_ipw=ipwe(covs = c("age","sex","diabetes"),Y="Y",A="CVD", data=covid19,
#'            shrink_rate = 1,ci_alpha=0.95,
#'            method="IPW-IPW",bootstrap=T,B=200)
#'
#' ##use TR AIPW method
#' summary(ipw_ipw)

#' @importFrom rootSolve multiroot
#' @importFrom survey svydesign svyglm
#' @importFrom parallel mclapply detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export


ipwe=function(covs,Y,A,data,shrink_rate=1,ci_alpha=0.95,
              method=c("IPW-WEE","IPW-DR","IPW-IPW"),bootstrap=T,B=200)
{
  ##match the method
  method=match.arg(method)

  ###start the main functions
  ##rename variables in data
  names(data)[names(data)==Y]="Y"
  names(data)[names(data)==A]="A"


  ###check input
  if((nrow(data)<1)|(ncol(data)<1))
    stop("Input dataset has either zero row or column.")

  if(length(levels(factor(data$Y)))>2)
    stop("not a binary outcome.")

  if(length(levels(factor(data$A)))>2)
    stop("not a binary treatment.")

  if((method!="IPW-IPW")&(method!="IPW-DR")&(method!="IPW-WEE"))
    stop("Specify a wrong method.")


  ##test missing values
  data_covomit=na.omit(data[,covs])
  data_outomit=na.omit(data[,"Y"])

  if(nrow(data_covomit)<nrow(data))
    stop("covariates have missing values. Impute missing covariates first.")

  if(length(data_outomit)<nrow(data))
    stop("The outcome has missing values. Impute missing outcome first.")


  #####ipwe key core function starts####

  ipwe_core=function(covs=covs,Y=Y,A=A,data=data,
                     shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method)
  {
    ##convert data without factor
    df_type=lapply(data[,c("A","Y",covs)], FUN=class)
    fac_index=c(which(df_type=="character"),which(df_type=="factor"))

    if(length(fac_index)==1)
    {
      data[,c("A","Y",covs)][,fac_index]=as.numeric(as.factor(data[,c("A","Y",covs)][,fac_index]))-1

    }else{

      data[,c("A","Y",covs)][,fac_index]=lapply(data[,c("A","Y",covs)][,fac_index],
                                                FUN=function(x)
                                                {x=as.factor(x)
                                                y=as.numeric(x)-1
                                                return(y)})
    }


    ##create miss ind
    data$r=ifelse(is.na(data$A)==T,1,0)

    ##create baseline
    data$x0=1

    ##fit miss model
    glm_mar=glm(as.formula(paste("r~Y+",paste(covs,collapse = "+"))),data=data,family = binomial)
    data$fit_mar=glm_mar$fitted.values

    ##inverse weights of missingessness
    data$weight_miss=1/(1-glm_mar$fitted.values)

    #shrinkage
    data$weight_miss[which(data$weight_miss>=quantile(data$weight_miss,shrink_rate,na.rm = TRUE))]=quantile(data$weight_miss,shrink_rate,na.rm = TRUE)

    ##remove all NA
    data_naomit=na.omit(data)


    #####step1: fit alpha ipw constructed from Miss+PS model#######
    ##Miss model must be correct, so alpha ee becomes WLA IPW
    ipw_design=svydesign(ids=~1, weights=~weight_miss, data=data_naomit)
    alpha_ipwglm= svyglm(as.formula(paste("A~",paste(covs,collapse = "+"))),design=ipw_design,family=quasibinomial(link="logit"),data=data_naomit)

    ##alpha ipw
    alpha_ipw=coefficients(alpha_ipwglm)

    ##fitted values
    fit_ps_ipw=as.vector(predict(alpha_ipwglm,newdata =data,type = "response"))


    #####step2: fit beta IPW #####
    ##miss model must be correct, so beta ee becomes beta IPW
    beta_glmipw=svyglm(as.formula(paste("Y~A+",paste(covs,collapse = "+"))),design=ipw_design,family=quasibinomial(link="logit"),data=data_naomit)

    ##beta ipw
    beta_ipw=coefficients(beta_glmipw)


    ##fitted response values in outcome model
    data_trt=data
    data_trt$A=1

    data_con=data
    data_con$A=0

    m1_ipw=as.vector(predict(beta_glmipw,newdata =data_trt,type = "response"))
    m0_ipw=as.vector(predict(beta_glmipw,newdata =data_con,type = "response"))

    ##X matrix as input variable values
    ##without the treatment
    X_all=as.matrix(data[,c("x0",covs)])

    ##with trt
    XA_all=as.matrix(data[,c("x0","A",covs)])

    ##set diff of F all-Fobs
    F_diff=rep(0,ncol(X_all))

    ##reset NA value
    data$A[is.na(data$A)] = -100

    ####Traditional IPW estimator####
    ##must require miss+ps are correct
    if(method=="IPW-IPW"){

      tau1_ipw=mean((1-data$r)*data$weight_miss*data$A*data$Y/fit_ps_ipw)
      tau0_ipw=mean((1-data$r)*data$weight_miss*(1-data$A)*data$Y/(1-fit_ps_ipw))

      ipw_est=(tau1_ipw/(1-tau1_ipw))/(tau0_ipw/(1-tau0_ipw))

    }

    if(method=="IPW-DR"){

      ####IPW-DR Estimator####
      ##Miss must be correct+ either PS/ OR is correct
      ##trt group dr estimator
      fit_dr1=data$A*data$Y/fit_ps_ipw-(data$A-fit_ps_ipw)/fit_ps_ipw*m1_ipw
      tau1_ipw_dr=mean((1-data$r)*data$weight_miss*fit_dr1)

      ##control group dr estimator
      fit_dr0=(1-data$A)*data$Y/(1-fit_ps_ipw)-(fit_ps_ipw-data$A)/(1-fit_ps_ipw)*m0_ipw
      tau0_ipw_dr=mean((1-data$r)*data$weight_miss*fit_dr0)

      ipw_est=(tau1_ipw_dr/(1-tau1_ipw_dr))/(tau0_ipw_dr/(1-tau0_ipw_dr))

    }


    if(method=="IPW-WEE"){

      ##for trt group
      score_beta_wee1=function(beta)
      {
        ##unknown beta
        link=X_all%*%beta
        exbit_beta=exp(link)/(1+exp(link))

        #five beta_ee_miss score functions
        score_ij=function(x_all)
        {
          F_all=sum((1-data$r)*data$weight_miss*data$A/fit_ps_ipw*x_all*(data$Y-exbit_beta))

          return(F_all)
        }

        for (j in 1:ncol(X_all)) {

          F_diff[j]=score_ij(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee1_sol=multiroot(f = score_beta_wee1, start =c(beta_ipw[1]+beta_ipw[2],beta_ipw[-c(1,2)]))


      ##for control group
      score_beta_wee0=function(beta)
      {
        ##unknown beta
        link=X_all%*%beta
        exbit_beta=exp(link)/(1+exp(link))

        #five beta_ee_miss score functions
        score_ij=function(x_all)
        {
          F_all=sum((1-data$r)*data$weight_miss*(1-data$A)/(1-fit_ps_ipw)*x_all*(data$Y-exbit_beta))

          return(F_all)
        }

        for (j in 1:ncol(X_all)) {

          F_diff[j]=score_ij(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee0_sol=multiroot(f = score_beta_wee0, start =beta_ipw[-2])

      ##estimated beta WEE
      beta_ipw_wee1=beta_wee1_sol$root
      beta_ipw_wee1_converge=beta_wee1_sol$estim.precis

      beta_ipw_wee0=beta_wee0_sol$root
      beta_ipw_wee0_converge=beta_wee0_sol$estim.precis


      ##first step: estimate fitted response from beta WEE
      link_wee1=exp(X_all%*%beta_ipw_wee1)/(1+exp(X_all%*%beta_ipw_wee1))
      link_wee0=exp(X_all%*%beta_ipw_wee0)/(1+exp(X_all%*%beta_ipw_wee0))

      #####final IPW DR WEE estimate####
      if((is.na(beta_ipw_wee1_converge)==T)|(is.na(beta_ipw_wee0_converge)==T))
      {
        ##if estimate tau is too large, remove TR WEE
        tau_ipw_wee=NA
      }else{

        ##estimate IPW WEE only when beta_ipw_wee is converged
        tau1_ipw_wee=mean((1-data$r)*data$weight_miss*link_wee1)
        tau0_ipw_wee=mean((1-data$r)*data$weight_miss*link_wee0)

        ipw_est=(tau1_ipw_wee/(1-tau1_ipw_wee))/(tau0_ipw_wee/(1-tau0_ipw_wee))
      }

    }

    ##return NA
    data$A[which(data$A==-100)] =NA


    ##point est
    out=list("Estimate"=ipw_est, "fit_ps"=as.vector(fit_ps_ipw),"miss_weights"=data$weight_miss,data=data)
    return(out)

  }


  ##point estimate
  ipw_est=ipwe_core(covs=covs,Y=Y,A=A,data=data,
                    shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method)
  point_est=ipw_est$Estimate


  ###########Bootstrap Function ##############


  if(bootstrap==T){

    ##variance
    boot_fun=function(j) {

      set.seed(j)

      ##size of data
      n=nrow(data)

      ##step1: Function for Bootstrap j loop: B=200, sample size=1000, with replace
      data_boot=data[sample(1:n,size=n,replace = T),]

      ##step2: for each data set, just get 200 times point estimation, 4x200 matrix
      beta_boot=ipwe_core(covs=covs,Y=Y,A=A,data=data_boot,
                          shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method)$Estimate
      return(beta_boot) #return a vector 3x1
    }

    ##detect system
    sys=(Sys.info()[['sysname']])

    ##detectCores
    n_cores=detectCores()-1

    if(sys=="Windows"){

      ##run foreach for windows.run parallel
      cl = makeCluster(n_cores)
      registerDoParallel(cl)

      boot_mat=foreach(i=1:B,.packages="rootSolve",.combine=cbind) %dopar% boot_fun(i)
      boot_vec=as.vector(boot_mat)

      stopCluster(cl)

    }else{

      ##run mclapply for mac and linux. run parallel
      boot_ls=mclapply(1:B, FUN=boot_fun,mc.cores=n_cores)

      ##cbind for each list, to get matrix
      boot_vec=as.vector(do.call(cbind,boot_ls))

    }

    ##remove NA in bootstrap
    boot_na_col=which(is.na(boot_vec))

    if(length(boot_na_col) != 0)
    {
      boot_vec=boot_vec[-boot_na_col]
    }

    ##step3: bootstrap standard error 4x1 vector
    # boot_se=sd(boot_vec)
    #
    # ##summary dataframe
    # df_sum=data.frame(
    #   "Estimate"=point_est,
    #   "BSE"=boot_se
    #   # ,"na_boot"=length(boot_na_col)
    #   )
    #
    # df_sum=round(df_sum,3)
    # row.names(df_sum)=paste0("TR ",method)
    # df_sum
    #
    # ##95% CI. paste0 without dropping 0
    # tr_ci_low=round(point_est-qnorm(0.5+0.5*ci_alpha,0,1)*boot_se,3)
    # tr_ci_up=round(point_est+qnorm(0.5+0.5*ci_alpha,0,1)*boot_se,3)
    # tr_ci=paste0("(",format(tr_ci_low,drop0Trailing = F),",",format(tr_ci_up,drop0Trailing = F),")")
    # df_sum=cbind(df_sum,"ci"=tr_ci)
    # colnames(df_sum)[which(colnames(df_sum)=="ci")]=paste0(100*ci_alpha,"% CI")

    ##add p value as 3 digit.
    #change character <0.001 when smaller than 0.001
    # df_sum$p.value=2*(1-pnorm(q=abs(point_est/boot_se-1),mean=0,sd=1))
    # df_sum$p.value=round(df_sum$p.value,3)

    ####bootstrap percentile CI ####
    ci_low_hybrid=round(quantile(boot_vec,na.rm = T,probs=0.025,type=1),3)
    ci_up_hybrid=round(quantile(boot_vec,na.rm = T,probs=0.975,type=1),3)
    tr_ci=paste0("(",format(ci_low_hybrid,drop0Trailing = F),",",format(ci_up_hybrid,drop0Trailing = F),")")

    ####percentile pvalue
    p.value=sum(boot_vec>=point_est)/B

    if(p.value<0.001){
      p.value="<0.001"

    }else{
      p.value=round(p.value,3)
    }

    ###summary df
    df_sum=data.frame(round(point_est,3),tr_ci,p.value)
    rownames(df_sum)=paste0(method,": ",A)
    colnames(df_sum)=c("Estimate",paste0(100*ci_alpha,"% CI"),"p.value")
  }

  if(bootstrap==F)
  {

    ##if not use bootstrap, only return point est
    df_sum=data.frame("Estimate"=point_est)
    df_sum=round(df_sum,3)
    rownames(df_sum)=paste0(method,": ",A)

  }

  ##return final list
  final=list("results"=df_sum,"fit_ps"=ipw_est$fit_ps,"miss_weights"=ipw_est$miss_weights,data=data)

  structure(final,class="trme")

}
