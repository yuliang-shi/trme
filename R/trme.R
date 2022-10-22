###triple robust for missing exposure (TRME) core function ######
#' @title Triple Robust Estimation for Missing Exposure
#' @author Yuliang Shi
#' @description Estimate the causal effect of the exposure on the outcome when the exposure is MAR. Adjust for both missing and confounding issues via simplified estimating equations with triple robust (TR) properties. Provide robust standard errors for inference purposes.
#'
#'
#' @param covs \code{\link{character}}  the vector name of all confounders which
#'   will affect both exposure and outcome.  Any interaction term or non-linear
#'   form can be added after creating those variables based on the data.
#' @param Y \code{\link{character}} the name of the outcome variable. The trme
#'   function currently only works on the binary outcome.
#' @param A \code{\link{character}} the name of binary treatment or exposure
#'   variable.
#' @param df_mar a data frame contains the variables in the models.
#' @param ps_model,imp_model a logical value either \code{TRUE} or \code{FALSE} ps_model, imp_model refers to the
#'   treatment and imputation model, respectively. If the model=\code{FALSE}, the Bayes
#'   rule will be applied to estimate it.
#' @param quan_value \code{\link{numeric}} shrinkage value. By default, we
#'   shrink weights larger than 99\% quantile into exact 99\% quantile to avoid
#'   extreme weights.
#' @param method the method to be used. either proposed "new" TR estimator or "ee" complex TR estimator.
#'
#'
#' @details The methods currently only work for the situation of missing
#'   exposure and require binary outcome and exposure variables. If either the
#'   covariates or outcome is also missing, first, impute missing values based on
#'   \code{\link{mice}}.
#'
#'   The basic assumption is the exposure variable is
#'   missing at random (MAR), which means that given all observed covariates and
#'   outcomes, the missingness (or missing indicator) is independent of the
#'   missing value itself. If the exposure is missing not at random (MNAR), try
#'   another method instead.
#'
#'   \code{method="new"} is a simplified estimating
#'   equation for the triple robust (TR) estimator (recommended), which speeds
#'   up the computation process and avoids some effects of extreme weights in
#'   the finite samples, but still keeps the same TR properties as the complex
#'   form. For more details, please review the reference paper.
#'
#'   \code{method="ee"} is a complex estimating equation for TR estimator, which
#'   may be influenced by extreme weights.
#'
#'   Both two TR estimators have the same
#'   asymptotical consistency when the sample size is large. To achieve
#'   consistency, both TR estimators require at least **one of two models**
#'   condition, which means if the missingness model is correct, we require
#'   either the treatment or outcome model to be correct; If the missingness
#'   model is wrong, but the outcome model is correct, we require either the
#'   imputation or the treatment model to be correct.
#'
#'   The asymptotical standard
#'   errors may not be the same. However, both TR estimators utilize robust
#'   standard error (RSE) based on the sandwich formula, which is quite robust
#'   to the model's misspecification.
#'
#' @return Use \code{summary()} function to print out a data frame including summarized results.
#' \itemize{
#' \item{\code{Estimate: }}{estimated causal effect (log odds ratio) of exposure on the outcome. }
#' \item{\code{Robust.SE: }}{estimated robust standard errors used for inference.}
#' \item{\code{p.value: }}{p values for two-sided Wald test.}
#' \item{\code{95\% CI: }}{95\% two-sided confidence interval.}}
#'
#' In addition, other fitted values are stored in the full list.
#' \itemize{
#' \item{\code{vcov: }}{variance-covariance matrix among exposure and covariates.}
#' \item{\code{fit_ps_all: }}{predicted propensity score values for all subjects used to adjust for the confounding issue.}
#' \item{\code{fit_weightmiss: }}{fitted inverse weights of missingness used to adjust for the missing issue.}}
#'
#' Plots are also provided by calling \code{plot()} function. \code{hist_ps_control, hist_ps_trt: } histogram for predicted propensity score between control and treatment groups.
#'
#'
#' @keywords regression, robust.
#'
#' @note For more details, please review my GitHub website: https://github.com/yuliang-shi/trmd.
#' For citation, please cite the package as \dQuote{Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Triple Robust Estimator}. Statistics in Medicine.}
#'
#' @seealso \code{\link{summary.trmd}} for summarized result and \code{\link{plot.trmd}} for drawing histograms of fitted propensity score.
#' Other useful functions include \code{\link{svyglm}} for inverse-probability weighting or double robust methods and \code{\link{mice}} for multiple imputation on missing data.
#'
#' @references Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Triple Robust Estimator}. Statistics in Medicine. Submitted (11/2022).
#'
#' Zhang, Z., Liu, W., Zhang, B., Tang, L., and Zhang, J. (2016). \emph{Causal inference with missing exposure information: Methods and applications to an obstetric study}. Statistical Methods in Medical Research 25, 2053–2066.
#'
#'
#' @examples
#' set.seed(2000)
#' n = 2000 #sample size
#'
#' #####generate some continuous covariates
#' id = seq(1, n, by = 1)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = rnorm(n)
#'
#' #generate binary exposure from the PS model. x1,x2,x3 as confounders
#' a = -0.5 + 0.3 * x1 + 0.6 * x2 + 0.6 * x3
#' prob_a = 1 / (1 + exp(-a))
#' A = rbinom(n, 1, prob_a)
#'
#' #generate binary outcome from OR model
#' z = 0.9 + 1 *A + 0.1 * x1 + 0.3 * x2 + 0.2 * x3   #-1.5 is ok
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
#'
#' ##control the miss rate as 60\%
#' r_sim = 1.1 - 0.2 * x1 - 0.3 * x2 - 0.6 * x3 - 0.9 * Y #add random error
#' r = rbinom(n, 1, 1 / (1 + exp(-r_sim))) #miss rate
#'
#' ##df_mar: after include missing data
#' df_mar = cbind(df, "r" = r)
#' df_mar$A[which(df_mar$r == 1)] = NA
#'
#' ##test in the simulated data
#' ##use simplified method
#' tr_est = trme(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   df_mar = df_mar,
#'   ps_model = T,
#'   imp_model = T,
#'   quan_value = 0.99,
#'   method = "new"
#' )
#'
#' ##print out results and plots
#' summary(tr_est)
#' plot(tr_est)
#'
#' ##use complex method
#' tr_est = trme(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   df_mar = df_mar,
#'   ps_model = T,
#'   imp_model = T,
#'   quan_value = 0.99,
#'   method = "ee"
#' )
#'
#' ##print out results and plots
#' summary(tr_est)
#' plot(tr_est)
#'
#' @importFrom rootSolve multiroot
#' @importFrom matrixcalc is.singular.matrix
#' @export


trme=function(covs,Y,A,df_mar,ps_model=T,imp_model=T,quan_value=0.99,method=c("new","ee"))
{

  ##match the method
  method=match.arg(method)

  ##rename variables in df_mar
  names(df_mar)[names(df_mar)==Y]="Y"
  names(df_mar)[names(df_mar)==A]="A"


  ###check input
  if((nrow(df_mar)<1)|(ncol(df_mar)<1))
    stop("Input dataset has either zero row or column.")

  if(length(levels(factor(df_mar$Y)))>2)
    stop("not a binary outcome.")

  if(length(levels(factor(df_mar$A)))>2)
    stop("not a binary treatment.")

  if((method!="new")&(method!="ee"))
    stop("Specify a wrong method.")


  ##test missing values
  df_mar_covomit=na.omit(df_mar[,covs])
  df_mar_outomit=na.omit(df_mar[,"Y"])

  if(nrow(df_mar_covomit)<nrow(df_mar))
    stop("covariates have missing values. Impute missing covariates first.")

  if(length(df_mar_outomit)<nrow(df_mar))
    stop("The outcome has missing values. Impute missing outcome first.")


  ##TRME Core function starts
  ##convert data without factor
  df_type=lapply(df_mar[,c("A","Y",covs)], FUN=class)
  fac_index=c(which(df_type=="character"),which(df_type=="factor"))

  if(length(fac_index)==1)
  {
    df_mar[,c("A","Y",covs)][,fac_index]=as.numeric(as.factor(df_mar[,c("A","Y",covs)][,fac_index]))-1

  }else{

    df_mar[,c("A","Y",covs)][,fac_index]=lapply(df_mar[,c("A","Y",covs)][,fac_index],
                                                FUN=function(x)
                                                {x=as.factor(x)
                                                y=as.numeric(x)-1
                                                return(y)})
  }



  ##create miss ind
  df_mar$r=ifelse(is.na(df_mar$A)==T,1,0)

  ##create baseline
  df_mar$x0=1

  ##fit miss model
  glm_mar=glm(as.formula(paste("r~Y+",paste(covs,collapse = "+"))),data=df_mar,family = binomial)
  glm_mar


  ##inverse weights of missingessness
  df_mar$weight_miss=1/(1-glm_mar$fitted.values)

  #shrinkage to 99%
  df_mar$weight_miss[which(df_mar$weight_miss>=quantile(df_mar$weight_miss,quan_value,na.rm = TRUE))]=quantile(df_mar$weight_miss,quan_value,na.rm = TRUE)


  ##fit naive ps model without adjusting for missingness
  glm_ps=glm(as.formula(paste("A~",paste(covs,collapse = "+"))),data=df_mar,family = binomial(link="logit"))
  glm_ps
  alpha_ps=coefficients(glm_ps)

  ##fit naive or model without adjusting for missingness
  glm_or=glm(as.formula(paste("Y~A+",paste(covs,collapse = "+"))),data=df_mar,family = binomial)
  glm_or


  ##fit imp model
  ##imp model A|X,Y by MAR, not affected by missing data
  if(imp_model==T){
    ##else: imp model is right.
    glm_ps_xy=glm(as.formula(paste("A~Y+",paste(covs,collapse = "+"))),data=df_mar,family = binomial(link="logit")) #for A|X,Y

    ###Fit on complete data. Then, predict A|X,Y for all subjects
    df_mar$pred_ps_xy=predict(glm_ps_xy,newdata = df_mar,type="response")

  }else{

    ##if imp model is wrong, but ps and or models are correct.
    ##Use Bayes to predict A|X,Y

    ##note: for users, as long as they use imp_model=F, we will apply Bayes rule.
    ##in fact, we assume or model and ps models are correct, so we can use Bayes rule to transfer.
    ##if imp is wrong, or and PS models are also wrong, using Bayes rule does not help. Always biased P(A|X,Y)

    df_mar$pred_ps_xy=rep(0,nrow(df_mar))

    ##or model: subset on trt & observed
    glm_or_trt=glm(as.formula(paste("Y~",paste(covs,collapse = "+"))),data=subset(df_mar,A==1),family=binomial(link="logit"))
    beta_trt=coefficients(glm_or_trt)

    ##or model: subset on control & observed
    glm_or_con=glm(as.formula(paste("Y~",paste(covs,collapse = "+"))),data=subset(df_mar,A==0),family=binomial(link="logit"))
    beta_con=coefficients(glm_or_con)

    ##trt model: subset on observed alpha_ps

    ##Bayes transfer A|X,Y. predict for all subjects
    ##design X mat
    X_all=as.matrix(df_mar[,c("x0",covs)])
    link_glm_or_trt=exp(X_all%*%beta_trt)/(1+exp(X_all%*%beta_trt))
    link_glm_or_con=exp(X_all%*%beta_con)/(1+exp(X_all%*% beta_con))
    link_glm_ps=exp(X_all%*%alpha_ps)/(1+exp(X_all%*%alpha_ps))


    ##predict for all subjects in two different y=1 and y=0 model
    df_mar$pred_ps_xy_y1=(link_glm_or_trt*link_glm_ps)/(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps))
    df_mar$pred_ps_xy_y0=((1-link_glm_or_trt)*link_glm_ps)/(1-(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps)))

    ##match the fitted value with the correct position with observed y=1 or y=0
    df_mar$pred_ps_xy[df_mar$Y==1]=df_mar$pred_ps_xy_y1[df_mar$Y==1]
    df_mar$pred_ps_xy[df_mar$Y==0]=df_mar$pred_ps_xy_y0[df_mar$Y==0]

  }



  ##remove all NA
  df_mar_naomit=na.omit(df_mar)


  ##X matrix as input variable values
  ##without the treatment
  X_all=as.matrix(df_mar[,c("x0",covs)])
  X_obs=as.matrix(df_mar_naomit[,c("x0",covs)]) #design matrix nxp

  ##with trt
  XA_all=as.matrix(df_mar[,c("x0","A",covs)])
  XA_obs=as.matrix(df_mar_naomit[,c("x0","A",covs)])

  #input matrix for the subset on trt and control
  XA_con_all=XA_all
  XA_con_all[,"A"]=0
  XA_trt_all=XA_all
  XA_trt_all[,"A"]=1

  ##set diff of F all-Fobs
  F_diff=rep(0,ncol(X_all))


  #######Fit alpha EE ######

  alpha=coefficients(glm_ps)

  ##score functions for alpha ee
  score_alpha_ee=function(alpha)
  {

    #dimnesion problem must use observed values to input the data
    link=X_obs%*%alpha
    link_all=X_all%*%alpha #for all subjects

    #five beta_ee_miss score functions
    score_alphaj=function(x_obs,x_all)
    {
      #only observed subjects
      u=(df_mar_naomit$A*x_obs-(1-df_mar_naomit$A)*x_obs*exp(link))/(1+exp(link))
      F_obs=sum((1-df_mar_naomit$r)*df_mar_naomit$weight_miss*u)


      ##all subjects
      u_fit=(df_mar$pred_ps_xy*x_all-(1-df_mar$pred_ps_xy)*x_all*exp(link_all))/(1+exp(link_all))
      F_all=sum((glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*u_fit)

      return(F_obs-F_all)

    }

    ##loop all columns for all covs and observed covs
    ##return F diff=F_obs-F_all to multiroot
    ##Then try the possible values s.t F diff is small enough to 0 and not change much
    for (j in 1:ncol(X_all)) {

      F_diff[j]=score_alphaj(x_obs=X_obs[,j],x_all=X_all[,j])

    }

    return(F_diff)
  }


  #solve alpha
  alpha_ee_sol=multiroot(f = score_alpha_ee, start =coefficients(glm_ps))
  alpha_ee=alpha_ee_sol$root
  alpha_ee




  #######Fit beta EE ######

  ##score functions for beta miss
  score_beta_miss=function(beta)
  {


    #dimenesion problem must use observed values to input the data
    link=XA_obs%*%beta
    link_con=XA_con_all%*%beta
    link_trt=XA_trt_all%*%beta

    #five beta_ee_miss score functions
    score_ij=function(x_obs,x_all,type)
    {
      F_obs=sum((1-df_mar_naomit$r)*df_mar_naomit$weight_miss*(df_mar_naomit$Y*x_obs-(1-df_mar_naomit$Y)*x_obs*exp(link))/(1+exp(link)))

      ##F_all for all subjects n=1000
      ##v control group
      v_con=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_con))/(1+exp(link_con))
      prob_con=1-df_mar$pred_ps_xy

      ##v trt group
      v_trt=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_trt))/(1+exp(link_trt))
      prob_trt=df_mar$pred_ps_xy
      trt=v_trt*prob_trt

      if(type=="cov"){
        control=v_con*prob_con
      } else
      {  control=0 }

      F_all=sum((glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*(control+trt))
      return(F_obs-F_all)

    }


    ##loop all columns for all covs and observed covs
    ##return F diff=F_obs-F_all to multiroot
    ##Then try the possible values s.t F diff is small enough to 0 and not change much
    for (j in 1:ncol(X_all)) {

      F_diff[j]=score_ij(x_obs=X_obs[,j],x_all=X_all[,j],type = "cov")

    }

    ## especially, solve the equation for exposure
    F_trt=score_ij(x_obs=df_mar_naomit$A,x_all=1,type="expose")

    ##combine results with baseline, exposure, and covariates
    return(c(F_diff[1], F_trt, F_diff[-1]))
  }

  #beta
  beta_ee_miss_sol=multiroot(f = score_beta_miss, start =coefficients(glm_or))
  beta_ee_miss=beta_ee_miss_sol$root
  beta_ee_miss_coverage=beta_ee_miss_sol$estim.precis




  ###fitted values for PS
  if(ps_model==T){

    ###if PS mdoel is correct.
    linear_dr_alpha=X_obs%*%alpha_ee
    fit_ps_dr=exp(linear_dr_alpha)/(1+exp(linear_dr_alpha))

    link_dr_alpha_all=X_all%*%alpha_ee
    fit_ps_dr_all=exp(link_dr_alpha_all)/(1+exp(link_dr_alpha_all))

  }else{
    ##Similarly, as long as users use ps_mdoel=F, we will apply the Bayes rule to transfer A|X
    ##However, it needs both OR and Imp models are correct to use Bayes rule.
    ##if PS model is wrong, both OR and Imp models are also wrong, using the Bayes rule doesn't help. Always a biased estimate.

    ##fit two imp models on a different subset

    glm_ps_xy_trt=glm(as.formula(paste("A~",paste(covs,collapse = "+"))),data=subset(df_mar,Y==1),family = binomial)
    glm_ps_xy_con=glm(as.formula(paste("A~",paste(covs,collapse = "+"))),data=subset(df_mar,Y==0),family = binomial)

    ##predict for all subjects due to MAR
    pred_glm_ps_xy_trt=predict(glm_ps_xy_trt,newdata=df_mar,type="response")
    pred_glm_ps_xy_con=predict(glm_ps_xy_con,newdata=df_mar,type="response")


    ##use beta ee miss to predict for all subjects
    pred_beta_ee_trt=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))
    pred_beta_ee_con=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))

    ##Bayes predict ps from correct or imp models
    ##need stronger assumption that P(Y=1|A=0,X) P(Y=1|A=1,X) has same format
    fit_ps_dr_all=-(pred_glm_ps_xy_trt*pred_beta_ee_con)/(pred_glm_ps_xy_trt*pred_beta_ee_trt-pred_glm_ps_xy_trt*pred_beta_ee_con-pred_beta_ee_trt)
    fit_ps_dr=fit_ps_dr_all[df_mar$r==0]

  }


  #inverse weights of ps_dr values and shrink
  ipw_dr=rep(0,dim(df_mar_naomit)[1])
  ipw_dr[which(df_mar_naomit$A==1)]=1/fit_ps_dr[which(df_mar_naomit$A==1)]
  ipw_dr[which(df_mar_naomit$A==0)]=1/(1-fit_ps_dr[which(df_mar_naomit$A==0)])

  ipw_dr_aug=rep(0,dim(df_mar_naomit)[1])
  ipw_dr_aug[which(df_mar_naomit$A==1)]=(1-fit_ps_dr[which(df_mar_naomit$A==1)])/fit_ps_dr[which(df_mar_naomit$A==1)]
  ipw_dr_aug[which(df_mar_naomit$A==0)]=(fit_ps_dr[which(df_mar_naomit$A==0)])/(1-fit_ps_dr[which(df_mar_naomit$A==0)])

  ####inverse fit ps for all subjects in control and treatment
  inv_ps_con=1/(1-fit_ps_dr_all)
  inv_ps_trt=1/fit_ps_dr_all


  ##for observed subjects 1/fit_ps_dr
  # inv_ps_dr_con=1/(1-fit_ps_dr_all)
  # inv_ps_dr_trt=1/fit_ps_dr_all

  #augmented term
  inv_ps_con_aug=fit_ps_dr_all/(1-fit_ps_dr_all)
  inv_ps_trt_aug=(1-fit_ps_dr_all)/fit_ps_dr_all

  #shrink
  ipw_dr[which(ipw_dr>=quantile(ipw_dr,quan_value,na.rm = TRUE))]=quantile(ipw_dr,quan_value,na.rm = TRUE)
  ipw_dr_aug[which(ipw_dr_aug>=quantile(ipw_dr_aug,quan_value,na.rm = TRUE))]=quantile(ipw_dr_aug,quan_value,na.rm = TRUE)
  inv_ps_con[which(inv_ps_con>=quantile(inv_ps_con,quan_value,na.rm = TRUE))]=quantile(inv_ps_con,quan_value,na.rm = TRUE)
  inv_ps_trt[which(inv_ps_trt>=quantile(inv_ps_trt,quan_value,na.rm = TRUE))]=quantile(inv_ps_trt,quan_value,na.rm = TRUE)
  inv_ps_con_aug[which(inv_ps_con_aug>=quantile(inv_ps_con_aug,quan_value,na.rm = TRUE))]=quantile(inv_ps_con_aug,quan_value,na.rm = TRUE)
  inv_ps_trt_aug[which(inv_ps_trt_aug>=quantile(inv_ps_trt_aug,quan_value,na.rm = TRUE))]=quantile(inv_ps_trt_aug,quan_value,na.rm = TRUE)



  #####New TR Core function ######
  if(method=="new")
  {
    score_beta_tr_new=function(beta)
    {
      #linear fit
      # link_fit=exp(glm_or$linear.predictors)

      #dimenesion problem must use observed values to input the data
      link=XA_obs%*%beta

      link_con=XA_con_all%*%beta_ee_miss

      link_trt=XA_trt_all%*%beta_ee_miss

      #five beta_ee_miss score functions
      score_tr_new=function(x_obs,x_all,type)
      {
        F_obs=sum((1-df_mar_naomit$r)*df_mar_naomit$weight_miss*ipw_dr*(df_mar_naomit$Y*x_obs-(1-df_mar_naomit$Y)*x_obs*exp(link))/(1+exp(link)))

        ##F_all for all subjects n=1000
        ##v control group
        v_con=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_con))/(1+exp(link_con))
        prob_con=1-df_mar$pred_ps_xy

        ##v trt group
        v_trt=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_trt))/(1+exp(link_trt))
        prob_trt=df_mar$pred_ps_xy
        trt=inv_ps_trt*v_trt*prob_trt

        if(type=="cov"){
          control=inv_ps_con*v_con*prob_con
        } else
        {  control=inv_ps_con*0*prob_con}

        F_all=sum((glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*(control+trt))
        return(F_obs-F_all)

      }


      ##loop all columns for all covs and observed covs
      ##return F diff=F_obs-F_all to multiroot
      ##Then try the possible values s.t F diff is small enough to 0 and not change much
      for (j in 1:ncol(X_all)) {

        F_diff[j]=score_tr_new(x_obs=X_obs[,j],x_all=X_all[,j],type = "cov")

      }

      ## especially, solve equations for exposure
      F_trt=score_tr_new(x_obs=df_mar_naomit$A,x_all=1,type="expose")

      ##combine results with baseline, exposure, and covariates
      return(c(F_diff[1], F_trt, F_diff[-1]))

    }

    #beta TR New
    beta_tr_sol=multiroot(f = score_beta_tr_new,start =beta_ee_miss)
    beta_tr=beta_tr_sol$root
    beta_tr_converge=beta_tr_sol$estim.precis

  }





  #####TR EE complex form######
  if(method=="ee")
  {

    score_beta_tr_ee=function(beta)
    {
      #link
      link=XA_obs%*%beta

      #fitted link value
      link_fit=XA_obs%*%beta_ee_miss
      y_fit_dr=exp(link_fit)/(1+exp(link_fit))

      #fit in the control group
      link_fit_allcon=XA_con_all%*%beta_ee_miss

      #fit in treat group
      link_fit_alltrt=XA_trt_all%*%beta_ee_miss

      #equation x_j
      score_j=function(x_obs,x_all,type)
      {
        ##x obs is observed data, x_all is all data

        #score function
        v=(df_mar_naomit$Y*x_obs-(1-df_mar_naomit$Y)*x_obs*exp(link))/(1+exp(link))
        v_fit=(y_fit_dr*x_obs-(1-y_fit_dr)*x_obs*exp(link_fit))/(1+exp(link_fit))

        #observed part 1
        F_obs=sum(df_mar_naomit$weight_miss*(ipw_dr*v-ipw_dr_aug*v_fit))

        ###all subjects part 2
        #a=0 control group
        v_fit_con=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_allcon))/(1+exp(link_fit_allcon))
        v_fit_con_aug=(-x_all)/(1+exp(link_fit_allcon))^2
        prob_con=1-df_mar$pred_ps_xy

        if(type=="cov")
          control=(inv_ps_con*v_fit_con-inv_ps_con_aug*v_fit_con_aug)*prob_con

        else
          control=0
        #when xij=aij exposure #ai=0 no values in control. only have values in trt

        ##a=1 in treatment
        v_fit_trt=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_alltrt))/(1+exp(link_fit_alltrt))
        v_fit_trt_aug=(-x_all)/(1+exp(link_fit_alltrt))^2
        prob_trt=df_mar$pred_ps_xy
        trt=(inv_ps_trt*v_fit_trt-inv_ps_trt_aug*v_fit_trt_aug)*prob_trt

        F_all=sum((glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*(control+trt))

        return(F_obs-F_all)
      }

      ##loop all columns for all covs and observed covs
      ##return F diff=F_obs-F_all to multiroot
      ##Then try the possible values s.t F diff is small enough to 0 and not change much
      for (j in 1:ncol(X_all)) {

        F_diff[j]=score_j(x_obs=X_obs[,j],x_all=X_all[,j],type = "cov")

      }

      ## especially, solve equation for exposure
      F_trt=score_j(x_obs=df_mar_naomit$A,x_all=1,type="expose")

      ##combine results with baseline, exposure, and covariates
      return(c(F_diff[1], F_trt, F_diff[-1]))

    }

    #solve
    beta_tr_sol=multiroot(f = score_beta_tr_ee, start=coefficients(glm_or))
    beta_tr=beta_tr_sol$root
    beta_tr_converge=beta_tr_sol$estim.precis
  }



  ##check for converge
  if((is.na(beta_tr_converge)==T)|(is.na(beta_ee_miss_coverage)==T))
  {
    beta_tr_new=NA
    se_beta_tr_new=NA

    stop("Results are not converged. Enlarge dataset or specify the correct models.")
  }


  ###########Variance Function ##############

  ##setup
  n_obs=dim(df_mar)[1]
  l=length(coefficients(glm_or)) #length of total parameters including exposure
  I=matrix(0,nrow=l,ncol=l) #information matrix pxp
  # S_mat=matrix(0,nrow=n_obs,ncol=l) #observed score nxp
  F_obs=rep(0,dim(df_mar)[1]) #store values for observed subjects. others are 0.
  F_diff_mat=matrix(0,nrow=n_obs,ncol=ncol(X_all)) ##difference between F_obs and F_all. Include baseline and covs
  colnames(F_diff_mat)=c("x0",covs)


  Var_fun=function(beta_tr,method="new")
  {
    #beta_fit is fitted beta value
    #information matrix is all same
    ###robust SE=sqrt(diag(B%*%M%*%B))
    ##B is the inverse information matrix

    #fitted link value
    link_fit=XA_obs%*%beta_tr
    y_fit_dr=exp(link_fit)/(1+exp(link_fit))

    #fit in control group
    link_fit_allcon=XA_con_all%*%beta_tr

    #fit in treat group
    link_fit_alltrt=XA_trt_all%*%beta_tr



    ##bread matrix pxp only includes observed data,
    for (j in 1:l) {
      for (k in 1:l) {
        I[j,k]=sum((1-df_mar_naomit$r)*df_mar_naomit$weight_miss*ipw_dr*XA_obs[,j]*XA_obs[,k]*exp(link_fit)/(1+exp(link_fit))^2)
      }
    }


    ######Score fun for TR EE New#########
    ##Note: remove all sum for F obs and F all. We just need nxp fitted score
    ##in TSE, link should replace by link_fit
    ##all link_fit is just a linear combination of x*beta without exp. so add exp()
    ##F obs only has observed data. others are 0.
    #fitted score equation x_j
    if(method=="new")
    {
      score_j=function(x_obs,x_all,type)
      {
        ##x obs is observed data, x_all is all data
        F_obs[which(df_mar$r==0)]=(1-df_mar_naomit$r)*df_mar_naomit$weight_miss*ipw_dr*(df_mar_naomit$Y*x_obs-(1-df_mar_naomit$Y)*x_obs*exp(link_fit))/(1+exp(link_fit))

        ##F_all for all subjects n=1000
        ##v control group
        v_con=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_allcon))/(1+exp(link_fit_allcon))
        prob_con=1-df_mar$pred_ps_xy

        ##v trt group
        v_trt=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_alltrt))/(1+exp(link_fit_alltrt))
        prob_trt=df_mar$pred_ps_xy
        trt=inv_ps_trt*v_trt*prob_trt

        if(type=="cov"){
          control=inv_ps_con*v_con*prob_con
        } else
        {  control=inv_ps_con*0*prob_con}

        F_all=(glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*(control+trt)

        return(F_obs-F_all)

      }
    }




    if(method=="ee")
    {

      #fitted score equation x_j
      score_j=function(x_obs,x_all,type)
      {
        ##x obs is observed data, x_all is all data

        #fitted score function
        v_hat=(df_mar_naomit$Y*x_obs-(1-df_mar_naomit$Y)*x_obs*exp(link_fit))/(1+exp(link_fit))
        v_fit=(y_fit_dr*x_obs-(1-y_fit_dr)*x_obs*exp(link_fit))/(1+exp(link_fit))

        #observed part 1
        F_obs[which(df_mar$r==0)]=df_mar_naomit$weight_miss*(ipw_dr*v_hat-ipw_dr_aug*v_fit)

        ###all subjects part 2
        #a=0 control group
        v_fit_con=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_allcon))/(1+exp(link_fit_allcon))
        v_fit_con_aug=(-x_all)/(1+exp(link_fit_allcon))^2
        prob_con=1-df_mar$pred_ps_xy

        if(type=="cov")
          control=(inv_ps_con*v_fit_con-inv_ps_con_aug*v_fit_con_aug)*prob_con

        #when xij=aij exposure #ai=0 no values in control. only have values in trt
        else
          control=0

        ##a=1 in treatment
        v_fit_trt=(df_mar$Y*x_all-(1-df_mar$Y)*x_all*exp(link_fit_alltrt))/(1+exp(link_fit_alltrt))
        v_fit_trt_aug=(-x_all)/(1+exp(link_fit_alltrt))^2
        prob_trt=df_mar$pred_ps_xy
        trt=(inv_ps_trt*v_fit_trt-inv_ps_trt_aug*v_fit_trt_aug)*prob_trt

        F_all=(glm_mar$fitted.values-df_mar$r)*df_mar$weight_miss*(control+trt)

        return(F_obs-F_all)

      }
    }


    ##loop all columns for all covs and observed covs
    ##return F diff=F_obs-F_all to multiroot
    ##Then try the possible values s.t F diff is small enough to 0 and not change much
    for (j in 1:ncol(X_all)) {

      F_diff_mat[,j]=score_j(x_obs=X_obs[,j],x_all=X_all[,j],type = "cov")

    }

    ## especially, solve equations for exposure
    F_trt=score_j(x_obs=df_mar_naomit$A,x_all=1,type="expose")

    ##combine results with baseline, exposure, and covariates
    #observed score nx(p+2)
    S_mat=cbind(F_diff_mat[,1],F_trt,F_diff_mat[,-1])


    #second condition for TSE: check whether it is singular matrix.
    #if it is singular. cannot solve it. set as NA
    if(is.singular.matrix(I)==T)
    {
      se_beta_tr=NA
      p_tr=NA
      cov_beta=NA
    }

    #if it is NOT singular. can solve it.
    if(is.singular.matrix(I)==F)
    {
      #inverse of information matrix pxp
      B=solve(I, tol = 1e-20)

      #meat matrix pxn x nxp=pxp
      M=t(S_mat)%*%S_mat


      ##covariance matrix using RSE
      cov_beta= B%*%M%*%B

      rownames(cov_beta)=c("x0","A",covs)
      colnames(cov_beta)=c("x0","A",covs)


      #covariance  matrix
      se_beta_tr=sqrt(diag(cov_beta))
      p_tr=2*(1-pnorm(q=beta_tr/se_beta_tr,mean=0,sd=1))

    }

    ##output
    out=list("se_beta_tr"=se_beta_tr,"p_tr"=p_tr,"cov_beta"=cov_beta)
    return(out)
  }


  ##variance and covariance matrix
  var_beta_tr=Var_fun(beta_tr=beta_tr,method = method)

  #robust se output
  # se_beta_tr_new= sqrt(diag(var_beta_tr_new))
  # se_beta_tr_ee= sqrt(diag(var_beta_tr_ee))
  #
  # ##approx normal. Wald test.
  # p_tr_new=2*(1-pnorm(q=beta_tr_new/se_beta_tr_new,mean=0,sd=1))
  # p_tr_ee=2*(1-pnorm(q=beta_tr_ee/se_beta_tr_ee,mean=0,sd=1))



  ######final output######
  df_sum=data.frame("Estimate"=beta_tr,"Robust SE"=var_beta_tr$se_beta_tr,
                    "p value"=var_beta_tr$p_tr)
  df_sum=round(df_sum,3)
  row.names(df_sum)=c("(Intercept)",A,covs)

  ##95% CI
  tr_ci_low=round(df_sum$Estimate-1.96*df_sum$Robust.SE,3)
  tr_ci_up=round(df_sum$Estimate+1.96*df_sum$Robust.SE,3)
  tr_ci=paste0("(",tr_ci_low,",",tr_ci_up,")")
  df_sum=cbind(df_sum,"95% CI"=tr_ci)


  ###return list with a summary table, covariance matrix, fitted ps values, fitted miss weights
  final=list("results"=df_sum,"vcov"=var_beta_tr$cov_beta,
             "fit_ps_all"=fit_ps_dr_all,"fit_weightmiss"=df_mar$weight_miss,
             "df_mar"=df_mar)

  ##return
  structure(final, class = "trmd") # S3 class
  # return(final)

}


