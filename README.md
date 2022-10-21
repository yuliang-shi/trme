# trme: Triple Robust Estimation for Missing Exposure
Goal: Estimate the causal effect of the exposure on the outcome when the exposure is missing at random (MAR). 

Method: We develop simplified estimating equations to adjust for both missing and confounding issues. The new triple robust estimators contain TR robust properties: it can achieve consistency if the ``**one of two models**" condition satisfies, which means if the missingness model is correct, we require either the treatment or outcome model to be correct; If the missingness model is wrong, but the outcome model is correct, we require either the imputation or the treatment model to be correct.

Both two TR estimators have the same asymptotical consistency when the sample size is large. The asymptotical standard errors may not be the same. However, both TR estimators utilize robust standard error (RSE) based on the sandwich formula, which is quite robust to the misspecification of the model. Therefore, RSE will be useful to make statistical inferences.

Compared with previous TR estimators from the complex estimating equation, the new TR estimator has a simpler form, so it speeds up the computation process and avoids some effects of extreme weights in the finite samples, but still keeps the same TR properties as the complex form. 


**To install the ``trme" package from the GitHub, you need to write down the commands in R or Rstudio:**
```
require("devtools")
devtools::install_github("yuliang-shi/trme" ,ref="main" ,auth_token = "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")
```

Note: the current version of the package only works for the case when the outcome and exposure variables are binary. In addition, the current package only supports the case when only the exposure variable is missing. If either the covariates or the outcome is also missing, we suggest imputing missing values based on the data. In the future, the package may be updated to deal with other cases such as the continuous outcome or the outcome is MAR.

For more details, please cite the reference paper: *Yuliang Shi, Yeying Zhu, Joel Dubin. Causal Inference on Missing Exposure via Triple Robust Estimator. Statistics in Medicine.*
