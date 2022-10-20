# trme: Triple Robust Estimation for Missing Exposure
Goal: Estimate the causal effect of the exposure on the outcome when the exposure is MAR. 

We develop a simplified estimating equations to adjust for both missing and confounding issues. The new triple robust estimators contains TR robust properties: it can achieve consistency if  **one of two models** condition satisfies, which means if the missingness model is correct, we require either the treatment or outcome model to be correct; If the missingness model is wrong, but the outcome model is correct, we require either the imputation or the treatment model to be correct.

Compared with previous TR estimators from complex estimating equation, the new TR estimamators is a simplified  speeds up the computation process and avoids some effects of extreme weights in the finite samples, but still keeps the same TR properties as the complex form. 


*To install packages, you need to write down the commands in R:*
require("devtools")
devtools::install_github("yuliang-shi/trme" ,ref="main" ,auth_token = "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")



For more details, please review the reference paper: Yuliang Shi, Yeying Zhu, Joel Dubin. \emph{Causal Inference on Missing Exposure via Triple Robust Estimator}. Statistics in Medicine.
