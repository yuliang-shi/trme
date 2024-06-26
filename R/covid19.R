#' @title COVID-19 Data for Missing Exposure
#'
#' @description
#' Data is collected from an observational study in Brazil during the second wave of pandemic (from 2020-09-20 to 2020-12-01).
#' In our study, the binary outcome is the mortality and the exposure is cardiovascular disease (CVD).
#' Other clinical variables are collected including age, sex, and diabetes, which may be considered as confounders to adjust for.
#' However, around 18\% patients are missing on CVD, which also leads to missing data problem.
#'
#'
#' @docType data
#'
#' @usage \code{data(covid19)}
#'
#' @format a \code{\link{data.frome}} including:
#' \itemize{
#'   \item \code{age: } \code{\link{numeric}} values.
#'   \item \code{sex: } \code{\link{character}} binary covariates, either male or female.
#'   \item \code{CVD: } \code{\link{numeric}}  binary exposure/treatment variable, either 0=non-CVD or 1=CVD.
#'   \item \code{diabetes: } \code{\link{character}} binary covariates, either No or Yes.
#'   \item \code{Y: } \code{\link{numeric}} the binary outcome for the mortality status, either 0=alive or 1=death}
#'
#' @keywords datasets
#'
#' @author Yuliang Shi
#'
#' @source View the original raw data in \href{https://integrasus.saude.ce.gov.br/}{IntegraSUS website}
#'
#' @references For more details, please review and cite the reference paper: Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Robust Estimator}. Biometrical Journal. Submitted.
#'
#' @examples
#' ##attach the data
#' require("trme")
#' data(covid19)
#'
#' ##show the structure of dataset
#' head(covid19)
#' str(covid19)
#'
"covid19"





