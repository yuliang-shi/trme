# trmd: Triple Robust Estimation for Missing Data
**Goal: Estimate the causal effect of the exposure on the outcome when the exposure is missing at random (MAR). **

Method: we develop simplified triple robust (TR) estimator to adjust for both missing and confounding issues in the observational studies. The new estimator contains TR robust properties, i.e. it can achieve consistency if the ``**one of two models**" condition satisfies, which means *if the missingness model is correct, we require either the treatment or outcome model to be correct; or if the missingness model is wrong, but the outcome model is correct, we require either the imputation or the treatment model to be correct.*

Both two TR estimators have the same asymptotical consistency when the sample size is large. Their asymptotical variance may not be the same. However, both TR estimators utilize robust standard error (RSE) based on the sandwich formula, which is quite robust to the misspecification of the model. Therefore, RSE will be reliable to conduct statistical inferences.

Advantages: compared with previous TR estimator from the complex estimating equation, the new TR estimator removes redundant terms, so it **avoids large computational cost and some effects of extreme weights in the finite samples.** Its simpler form helps us largely simplify the estimation process, but new TR estimator can keep the same properties and perform even better than the previous estimator in the complex form. 


**Use: to install the ``trmd" package from the GitHub, you need to run the following commands in R or Rstudio:**
```
install.packages("devtools")
library("devtools")
devtools::install_github("yuliang-shi/trmd" ,ref="main" ,auth_token = "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")
```

Within the "trmd" package, the main function is **"trme" only working for missing exposure**. After installation, to understand how to use the "trme" main function, please run the following codes. **We strongly suggest users to read the detailed instruction and try the example codes carefully.**
```
library("trmd")
help(trme) 

```
If the ``devtools" does not work, 1) check and update your R version, 2) check your path, and 3) restart your computer or R. If it still doesn’t work, try the following commands:
```
devtools::build github devtools()
require("devtools")
devtools::install_github("yuliang-shi/trmd" ,ref="main" ,auth_token = "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")
```

Note 1: The main code should run well on R or Rstudio under either Windows, MacOS, or LINUX operating system. However, sometimes the algorithm may not be converged because the sample size is small or the missing rates are very large. We would suggest to try different form of models, enlarge the sample size, impute missing data or reduce the missing rates, if possible.

Note 2: **the current version of the "trmd" package only supports the case when only the exposure variable is missing, i.e. using "trme" main function. In addition, the current package only works for the case when the outcome and exposure variables are binary.** If either the covariates or the outcome is also missing, we suggest firstly imputing missing values based on the data. In the future, the package may be updated to deal with other cases such as the continuous outcome or the outcome is MAR.

Note 3: for more technical details, please review and cite our paper: **Yuliang Shi, Yeying Zhu, Joel Dubin. Causal Inference on Missing Exposure via Triple Robust Estimator. Statistics in Medicine.** To support our work, please cite the paper correctly.

To contact the author, please visit my personal website: https://uwaterloo.ca/scholar/y323shi/home or send me an email on workday: yuliang.shi@uwaterloo.ca.
