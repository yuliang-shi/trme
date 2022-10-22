# trmd: Triple Robust Estimation for Missing Data
**Goal: Estimate the causal effect of the exposure on the outcome when the exposure is missing at random (MAR). **

Method: We develop simplified estimating equations to adjust for both missing and confounding issues. The new triple robust (TR) estimators contain TR robust properties: it can achieve consistency if the ``**one of two models**" condition satisfies, which means *if the missingness model is correct, we require either the treatment or outcome model to be correct; or if the missingness model is wrong, but the outcome model is correct, we require either the imputation or the treatment model to be correct.*

Both two TR estimators have the same asymptotical consistency when the sample size is large. Their asymptotical variance may not be the same. However, both TR estimators utilize robust standard error (RSE) based on the sandwich formula, which is quite robust to the misspecification of the model. Therefore, RSE will be reliable to conduct statistical inferences.

Advantages: compared with previous TR estimators from the complex estimating equation, the new TR estimator has a simpler form, so it **avoids large computational cost and some effects of extreme weights in the finite samples.** Its general form is easy to be built after we simplify the estimation process, but new TR estimator can keep the same properties and perform even better than the complex form. 


**Use: To install the ``trmd" package from the GitHub, you need to write down the following commands in R or Rstudio:**
```
install.packages("devtools")
library("devtools")
devtools::install_github("yuliang-shi/trmd" ,ref="main" ,auth_token = "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")
```

Note 1: within the "trmd" package, the main function is **"trme" only working for missing exposure**. The main code should run well on R or Rstudio under either Windows, MacOS, or LINUX operating system. If the ``devtools" does not work, 1) check and update your R version, 2) check your path, and 3) restart your computer or R. If it still doesn’t work, try the following commands:
```
devtools::build github devtools()
require("devtools")
devtools::install_github("yuliang-shi/trmd" ,ref="main" ,auth_token = "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")
```

Note 2: **the current version of the "trmd" package only supports the case when only the exposure variable is missing, i.e. using "trme" main function. In addition, the current package only works for the case when the outcome and exposure variables are binary.** If either the covariates or the outcome is also missing, we suggest firstly imputing missing values based on the data. In the future, the package may be updated to deal with other cases such as the continuous outcome or the outcome is MAR.

Note 3: for more technical details and supporting our work, please review and cite our paper: **Yuliang Shi, Yeying Zhu, Joel Dubin. Causal Inference on Missing Exposure via Triple Robust Estimator. Statistics in Medicine.**

To contact the author, please visit my personal website: https://uwaterloo.ca/scholar/y323shi/home or send me an email: yuliang.shi@uwaterloo.ca.
