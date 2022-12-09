###triple robust for missing exposure (TRME) core function ######
#' @title Triple Robust Estimation for Missing Exposure
#' @author Yuliang Shi
#' @description Estimate the causal effect of the exposure on the outcome when the exposure is MAR. Adjust for both missing and confounding issues via simplified estimating equations with triple robust (TR) properties. Provide estimated standard errors for inference purposes.
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
#' @param imp_model a \code{\link{logical}} value either \code{TRUE} or \code{FALSE} for correct or wrong imputation model. If the model=\code{FALSE}, the Bayes rule will be applied to estimate it.
#' @param shrink_rate \code{\link{numeric}} shrinkage rate. By default, no shrinkage applies.
#'   In some cases, we can shrink weights larger than 99\% quantile into exact 99\% quantile to avoid
#'   extreme weights.
#' @param method \code{\link{character}} the method to be used. either "AIPW" for TR estimator using augmented inverse-probability weighting or "WEE" for TR estimator using weighted estimating equations.
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
#'   \code{method="AIPW"} is for the TR AIPW estimator, which requires "two correct models from three groups (missing/imputation model group, treatment model group, outcome model group)".
#'
#'   \code{method="WEE"} is a weighted estimating equation for TR estimator (recommended), which avoids some effects of extreme weights in the finite samples, but it still keeps the same TR properties as the complex form. For more details, please review the reference paper.
#'
#'   Both two TR estimators have the same
#'   asymptotic consistency when the sample size is large. To achieve
#'   consistency, both TR estimators require at least **two correct models from three groups**
#'   condition, which means if the missingness model is correct, we require
#'   either the treatment or outcome model to be correct; If the missingness
#'   model is wrong, but the imputation model is correct, we require either the treatment or outcome model to be correct. If both imputation and missingness models are wrong, to acheive consistency, both treatment and outcome models should be correct, in order to apply Bayes approach.
#'
#' Both TR estimators utilize Bootstrap approach to estimate standard errors which can protect against
#'   the misspecification of model based on the simulation studies. By default, parallel computing will be applied to speed up computing process based on the operating system.
#'
#' @return Use \code{summary()} function to print out a data frame including summarized results.
#' \itemize{
#' \item{\code{Estimate: }}{estimated causal effect (odds ratio) of exposure on the outcome. }
#' \item{\code{95\% CI: }}{95\% two-sided confidence interval.}
#' \item{\code{p.value: }}{p values for two-sided Wald test.}}
#'
#' In addition, other fitted values are also saved in the list.
#' \itemize{
#' \item{\code{fit_ps: }}{fitted propensity scores for all subjects, which are used to adjust for the confounding issue.}
#' \item{\code{miss_weights: }}{fitted inverse weights of missingness used to adjust for the missing issue.}
#' \item{\code{hist_ps_control, hist_ps_trt: } use \code{plot()} function to draw density plots for fitted propensity score between control and treatment groups.}
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
#' ##use WEE as estimate
#' tr_wee = trme(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   data = data,
#'   imp_model = T,
#'   shrink_rate = 1,
#'   ci_alpha=0.95,
#'   method = "WEE",
#'   bootstrap=T,
#'   B=200
#' )
#'
#' ##print out results and plots
#' summary(tr_wee)
#'
#' ##use AIPW method
#' tr_aipw = trme(
#'   covs = c("x1", "x2", "x3"),
#'   Y = "Y",
#'   A = "A",
#'   data = data,
#'   imp_model = T,
#'   shrink_rate = 1,
#'   ci_alpha=0.95,
#'   method = "AIPW",
#'   bootstrap=T,
#'   B=200
#' )
#'
#' ##print out results and PS histogram plots
#' summary(tr_aipw)
#' plot(tr_aipw)
#'
#'
#' ########The second example for real data##########
#' require("trme")
#' data(covid19)
#'
#' ##use new TR estimator
#' tr_wee=trme(covs = c("age","sex","diabetes"),Y="Y",A="CVD", data=covid19,
#'             imp_model=T,shrink_rate = 1,ci_alpha=0.95,
#'             method="WEE",bootstrap=T,B=200)
#'
#' ##obtain estimate of causal effect and CI
#' summary(tr_wee)
#'
#' ##use TR WEE method
#' tr_aipw=trme(covs = c("age","sex","diabetes"),Y="Y",A="CVD", data=covid19,
#'            imp_model=T,shrink_rate = 1,ci_alpha=0.95,
#'            method="AIPW",bootstrap=T,B=200)
#'
#' ##use TR AIPW method
#' summary(tr_aipw)

#' @importFrom rootSolve multiroot
#' @importFrom parallel mclapply detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @export


trme=function(covs,Y,A,data,imp_model=T,shrink_rate=1,ci_alpha=0.95,
              method=c("AIPW","WEE"),bootstrap=T,B=200)
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

  if((method!="WEE")&(method!="AIPW"))
    stop("Specify a wrong method.")


  ##test missing values
  data_covomit=na.omit(data[,covs])
  data_outomit=na.omit(data[,"Y"])

  if(nrow(data_covomit)<nrow(data))
    stop("covariates have missing values. Impute missing covariates first.")

  if(length(data_outomit)<nrow(data))
    stop("The outcome has missing values. Impute missing outcome first.")


  #####trme core function starts####

  trme_core=function(covs=covs,Y=Y,A=A,data=data,imp_model=imp_model,
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


    ##fit naive ps model without adjusting for missingness
    glm_ps=glm(as.formula(paste("A~",paste(covs,collapse = "+"))),data=data,family = binomial(link="logit"))
    glm_ps
    alpha_ps=coefficients(glm_ps)

    ##fit naive or model without adjusting for missingness
    glm_or=glm(as.formula(paste("Y~A+",paste(covs,collapse = "+"))),data=data,family = binomial)
    glm_or_beta=coefficients(glm_or)


    ##fit imp model
    ##imp model A|X,Y by MAR, not affected by missing data
    if(imp_model==T){
      ##else: imp model is right.
      glm_ps_xy=glm(as.formula(paste("A~Y+",paste(covs,collapse = "+"))),data=data,family = binomial(link="logit")) #for A|X,Y

      ###Fit on complete data. Then, predict A|X,Y for all subjects
      data$pred_ps_xy=predict(glm_ps_xy,newdata = data,type="response")

    }else{

      ##note: for users, as long as they use imp_model=F, we will apply Bayes rule.
      ##if imp model wrong, but ps and or models are correct.
      ##Use Bayes to predict A|X,Y
      data$pred_ps_xy=rep(0,nrow(data))

      ##or model: subset on trt & observed
      ##trt model: subset on observed alpha_ps

      ##Bayes transfer A|X,Y. predict for all subjects
      link_glm_or_trt=exp(X_all%*%tr_estt)/(1+exp(X_all%*%tr_estt))
      link_glm_or_con=exp(X_all%*%beta_con)/(1+exp(X_all%*% beta_con))
      link_glm_ps=exp(X_all%*%alpha_ps)/(1+exp(X_all%*%alpha_ps))

      ##predict for all subjects in two different y=1 and y=0 model
      pred_ps_xy_y1=(link_glm_or_trt*link_glm_ps)/(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps))
      pred_ps_xy_y0=((1-link_glm_or_trt)*link_glm_ps)/(1-(link_glm_or_trt*link_glm_ps+link_glm_or_con*(1-link_glm_ps)))

      ##match the fitted value with correct position with observed y=1 or y=0
      data$pred_ps_xy[data$Y==1]=pred_ps_xy_y1[data$Y==1]
      data$pred_ps_xy[data$Y==0]=pred_ps_xy_y0[data$Y==0]


      ##for extreme case, when sample is very small, and glm does not converge,
      ##so pred_ps_xy is na or nan. bayes rule can't help. we use the wrong imp model
      if(sum(is.nan(data$pred_ps_xy))>=1|sum(is.na(data$pred_ps_xy))>=1)
      {
        ##warning
        warning("Bayes method does not converged. Apply wrong imputation model.")

        ##fit wrong imp models without Y
        glm_ps_xy=glm(A~x1+x2+x3,data=data,family = binomial(link="logit")) #for A|X,Y

        ###Fit on complete data. Then, predict A|X,Y for all subjects
        data$pred_ps_xy=predict(glm_ps_xy,newdata = data,type="response")

      }
    }


    ##remove all NA
    data_naomit=na.omit(data)


    ##X matrix as input variable values
    ##without the treatment
    X_all=as.matrix(data[,c("x0",covs)])
    X_obs=as.matrix(data_naomit[,c("x0",covs)]) #design matrix nxp

    ##with trt
    XA_all=as.matrix(data[,c("x0","A",covs)])
    XA_obs=as.matrix(data_naomit[,c("x0","A",covs)])

    #input matrix for the subset on trt and control
    XA_con_all=XA_all
    XA_con_all[,"A"]=0
    XA_trt_all=XA_all
    XA_trt_all[,"A"]=1

    ##set diff of F all-Fobs
    F_diff=rep(0,ncol(X_all))


    #######Fit alpha EE ######

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
        u=(data_naomit$A*x_obs-(1-data_naomit$A)*x_obs*exp(link))/(1+exp(link))
        F_obs=sum((1-data_naomit$r)*data_naomit$weight_miss*u)


        ##all subjects
        u_fit=(data$pred_ps_xy*x_all-(1-data$pred_ps_xy)*x_all*exp(link_all))/(1+exp(link_all))
        F_all=sum((glm_mar$fitted.values-data$r)*data$weight_miss*u_fit)

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
    alpha_ee_converge=alpha_ee_sol$estim.precis
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
        F_obs=sum((1-data_naomit$r)*data_naomit$weight_miss*(data_naomit$Y*x_obs-(1-data_naomit$Y)*x_obs*exp(link))/(1+exp(link)))

        ##F_all for all subjects n=1000
        ##v control group
        v_con=(data$Y*x_all-(1-data$Y)*x_all*exp(link_con))/(1+exp(link_con))
        prob_con=1-data$pred_ps_xy

        ##v trt group
        v_trt=(data$Y*x_all-(1-data$Y)*x_all*exp(link_trt))/(1+exp(link_trt))
        prob_trt=data$pred_ps_xy
        trt=v_trt*prob_trt

        if(type=="cov"){
          control=v_con*prob_con
        } else
        {  control=0 }

        F_all=sum((glm_mar$fitted.values-data$r)*data$weight_miss*(control+trt))
        return(F_obs-F_all)

      }


      ##loop all columns for all covs and observed covs
      ##return F diff=F_obs-F_all to multiroot
      ##Then try the possible values s.t F diff is small enough to 0 and not change much
      for (j in 1:ncol(X_all)) {

        F_diff[j]=score_ij(x_obs=X_obs[,j],x_all=X_all[,j],type = "cov")

      }

      ## especially, solve the equation for exposure
      F_trt=score_ij(x_obs=data_naomit$A,x_all=1,type="expose")

      ##combine results with baseline, exposure, and covariates
      return(c(F_diff[1], F_trt, F_diff[-1]))
    }

    #beta
    beta_ee_miss_sol=multiroot(f = score_beta_miss, start =coefficients(glm_or))
    beta_ee_miss_coverage=beta_ee_miss_sol$estim.precis
    beta_ee_miss=beta_ee_miss_sol$root




    ###fitted values for PS
    ###if PS model is correct. we can use original PS model to get fitted PS
    ##if PS model is wrong, it is biased. So we require Miss+OR or OR+imp is correct.
    linear_dr_alpha=X_obs%*%alpha_ee
    fit_ps_dr=exp(linear_dr_alpha)/(1+exp(linear_dr_alpha))

    link_dr_alpha_all=X_all%*%alpha_ee
    fit_ps_dr_all=exp(link_dr_alpha_all)/(1+exp(link_dr_alpha_all))


    ####inverse fit ps for all subjects in control and treatment
    inv_ps_dr=1/fit_ps_dr
    inv_ps_dr_con=1/(1-fit_ps_dr)

    inv_ps_dr_all=1/fit_ps_dr_all
    inv_ps_dr_all_con=1/(1-fit_ps_dr_all)

    #shrink
    inv_ps_dr[which(inv_ps_dr>=quantile(inv_ps_dr,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr,shrink_rate,na.rm = TRUE)
    inv_ps_dr_con[which(inv_ps_dr_con>=quantile(inv_ps_dr_con,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_con,shrink_rate,na.rm = TRUE)
    inv_ps_dr_all[which(inv_ps_dr_all>=quantile(inv_ps_dr_all,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_all,shrink_rate,na.rm = TRUE)
    inv_ps_dr_all_con[which(inv_ps_dr_all_con>=quantile(inv_ps_dr_all_con,shrink_rate,na.rm = TRUE))]=quantile(inv_ps_dr_all_con,shrink_rate,na.rm = TRUE)


    #####TR AIPW Method######


    ##reset A=NA as -100
    data$A[is.na(data$A)] = -100

    if(method=="AIPW"){


      ###predicted response for trt and control groups
      data$m1=exp(XA_trt_all%*%beta_ee_miss)/(1+exp(XA_trt_all%*%beta_ee_miss))
      data$m0=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))

      if((is.na(alpha_ee_converge)==F)&(is.na(beta_ee_miss_coverage)==F))
      {
        ##for trt group
        dr_tau1=data$A*data$Y*inv_ps_dr_all-(data$A-fit_ps_dr_all)*inv_ps_dr_all*data$m1
        dr_aug_tau1=data$pred_ps_xy*data$Y*inv_ps_dr_all-(data$pred_ps_xy-fit_ps_dr_all)*inv_ps_dr_all*data$m1

        ##for control group
        dr_tau0=(1-data$A)*data$Y*inv_ps_dr_all_con-(fit_ps_dr_all-data$A)*inv_ps_dr_all_con*data$m0
        dr_aug_tau0=(1-data$pred_ps_xy)*data$Y*inv_ps_dr_all_con-(fit_ps_dr_all-data$pred_ps_xy)*inv_ps_dr_all_con*data$m0

        ###TR AIPW for trt group
        tr_est_tau1=mean((1-data$r)*data$weight_miss*dr_tau1)-
          mean((data$fit_mar-data$r)*data$weight_miss*dr_aug_tau1)

        ###TR AIPW for control group
        tr_est_tau0=mean((1-data$r)*data$weight_miss*dr_tau0)-
          mean((data$fit_mar-data$r)*data$weight_miss*dr_aug_tau0)


        ##TR AIPW
        tr_est= (tr_est_tau1/(1-tr_est_tau1))/(tr_est_tau0/(1-tr_est_tau0))
        tr_est

      }else{

        ##if alpha ee or beta ee is not converged
        tr_est=NA
        warning("Algorithm does not converged! Either enlarge sample size or specify the correct model.")
      }

    }



    #####TR WEE Method######
    if(method=="WEE"){

      ##predicted response for trt or control
      exbit_trt_fit=exp(XA_trt_all%*%beta_ee_miss)/(1+exp(XA_trt_all%*%beta_ee_miss))
      exbit_con_fit=exp(XA_con_all%*%beta_ee_miss)/(1+exp(XA_con_all%*%beta_ee_miss))

      ##score functions for beta miss
      score_beta_wee1=function(beta)
      {
        #dimenesion problem must use observed values to input the data
        link_beta=X_all%*%beta
        exbit_beta=exp(link_beta)/(1+exp(link_beta))

        #five beta_ee_miss score functions
        beta_wee1=function(x_all)
        {
          F_trt=data$A*(data$Y-exbit_beta)-data$pred_ps_xy*(data$Y-exbit_trt_fit)

          F_all=sum((1-data$r)*data$weight_miss*inv_ps_dr_all*x_all*F_trt)

          return(F_all)
        }

        ##loop all columns for all covs and observed covs
        ##return F diff=F_obs-F_all to multiroot
        ##Then try the possible values s.t F diff is small enough to 0 and not changes much
        for (j in 1:ncol(X_all)) {
          F_diff[j]=beta_wee1(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee1_sol=multiroot(f = score_beta_wee1, start =c(beta_ee_miss[1]+beta_ee_miss[2],beta_ee_miss[-c(1,2)]))


      ###score function for beta wee0
      # beta=beta_ee_miss[-2]
      # x_all=X_all[,1]

      score_beta_wee0=function(beta){

        exbit_beta=exp(X_all%*%beta)/(1+exp(X_all%*%beta))

        beta_wee0=function(x_all){

          F_con=(1-data$A)*(data$Y-exbit_beta)-(1-data$pred_ps_xy)*(data$Y-exbit_con_fit)

          F_all=sum((1-data$r)*data$weight_miss*inv_ps_dr_all_con*x_all*F_con)

          return(F_all)
        }

        for (j in 1:ncol(X_all)) {
          F_diff[j]=beta_wee0(x_all=X_all[,j])
        }

        ##combine results with baseline, expousure, and covariates
        return(F_diff)
      }

      #beta
      beta_wee0_sol=multiroot(f = score_beta_wee0, start =beta_ee_miss[-2])

      ##estimated beta WEE
      tr_est_wee1=beta_wee1_sol$root
      tr_est_wee1_converge=beta_wee1_sol$estim.precis

      tr_est_wee0=beta_wee0_sol$root
      tr_est_wee0_converge=beta_wee0_sol$estim.precis

      ##pred response use WEE
      link_wee1=exp(X_all%*%tr_est_wee1)/(1+exp(X_all%*%tr_est_wee1))
      link_wee0=exp(X_all%*%tr_est_wee0)/(1+exp(X_all%*%tr_est_wee0))

      ##beta WEE estimate1
      tr_tau_wee1=mean((1-data$r)*data$weight_miss*(link_wee1-exbit_trt_fit))+
        mean(data$pred_ps_xy*data$Y*inv_ps_dr_all-(data$pred_ps_xy-fit_ps_dr_all)*inv_ps_dr_all*exbit_trt_fit)

      tr_tau_wee0=mean((1-data$r)*data$weight_miss*(link_wee0-exbit_con_fit))+
        mean((1-data$pred_ps_xy)*data$Y*inv_ps_dr_all_con-(fit_ps_dr_all-data$pred_ps_xy)*inv_ps_dr_all_con*exbit_con_fit)


      #####final TR WEE estimate####
      if((is.na(tr_est_wee1_converge)==T)|(is.na(tr_est_wee0_converge)==T))
      {
        ##if estimate tau is too large, remove TR WEE
        tr_est=NA
        warning("Algorithm does not converged! Either enlarge sample size or specify the correct model.")

      }else{

        ##estimate TR WEE only when alpha ee, beta ee, beta wee are converged
        ##and tr_tau_wee1,0 are not larger than 1
        tr_est=(tr_tau_wee1/(1-tr_tau_wee1))/(tr_tau_wee0/(1-tr_tau_wee0))
      }

    }

    data$A[which(data$A==-100)] =NA


    ##point est
    out=list("Estimate"=tr_est, "fit_ps"=as.vector(fit_ps_dr_all),"miss_weights"=data$weight_miss,data=data)
    return(out)

  }


  ##point estimate
  tr_est=trme_core(covs=covs,Y=Y,A=A,data=data,imp_model=imp_model,
                   shrink_rate=shrink_rate,ci_alpha=ci_alpha,method=method)
  point_est=tr_est$Estimate

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
      beta_boot=trme_core(covs=covs,Y=Y,A=A,data=data_boot,imp_model=imp_model,
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
    ci_low_per=round(quantile(boot_vec,na.rm = T,probs=0.025,type=1),3)
    ci_up_per=round(quantile(boot_vec,na.rm = T,probs=0.975,type=1),3)
    tr_ci=paste0("(",format(ci_low_per,drop0Trailing = F),",",format(ci_up_per,drop0Trailing = F),")")

    ####percentile pvalue
    p.value=sum(boot_vec>=point_est)/B

    if(p.value<0.001){
      p.value="<0.001"

    }else{
      p.value=round(p.value,3)
    }

    ###summary df
    df_sum=data.frame(round(point_est,3),tr_ci,p.value)
    rownames(df_sum)=paste0("TR ",method,": ",A)
    colnames(df_sum)=c("Estimate",paste0(100*ci_alpha,"% CI"),"p.value")
  }

  if(bootstrap==F)
  {

    ##if not use bootstrap, only return point est
    df_sum=data.frame("Estimate"=point_est)
    df_sum=round(df_sum,3)
    rownames(df_sum)=paste0("TR ",method,": ",A)

  }

  ##return final list
  final=list("results"=df_sum,"fit_ps"=tr_est$fit_ps,"miss_weights"=tr_est$miss_weights,data=data)

  structure(final,class="trme")

}


