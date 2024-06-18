# trme: Triple Robust Estimation for Missing Exposure

**Goal: Estimate the causal effect of the exposure on the outcome when the exposure is missing at random (MAR).**

Method: we develop simplified triple robust (TR) estimator to adjust for both missing and confounding issues in the observational studies. The new estimator contains TR robust properties, i.e. it can achieve consistency and asymptotic normality if the "**two out of three group models**" condition satisfies, which means *if either the missingness model or the imputation is correct, we require either the treatment or outcome model to be correct; or if both the missingness and imputation models are wrong, we require the treatment and outcome model are correct*.

Advantages: compared with previous TR estimator from the complex estimating equation, the new IPW and TR estimator **avoids large computational cost and some effects of extreme weights in the finite samples.** All estimators utilize bootstrap approach to conduct statistical inferences, which is reliable from the simulation studies.

**Use: to install the "trme" package from the GitHub, you need to run the following commands in R or Rstudio:**
```
install.packages("devtools")
devtools::install_github("yuliang-shi/trme", ref="main" ,auth_token = "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")
```

If the "devtools" does not work, 1) check and update your R version, 2) check the path of directory, and 3) restart R or your computer. If it still doesn’t work, please try the following commands:
```
devtools::build github devtools()
devtools::install_github("yuliang-shi/trme" ,ref="main" ,auth_token = "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")
```

Within the "trme" package, the main function is **"ipwe"  and "trme" only working for missing exposure**. After installation, to understand how to use the "trme" main function, please run the following codes. **We strongly suggest users to read the detailed instruction and run the examples carefully.**
```
##attach the package and read the instructions of trme function
library("trme")
help(ipwe)
help(trme) 
```

Note 1: The main code should run well on R or Rstudio under either Windows, MacOS, or LINUX operating system. However, sometimes the algorithm may not be converged because the sample size is small, or the missing rates are large, or extreme propensity scores occur. We would suggest to try different models, enlarge the sample size, or impute missing covariates to reduce the missing rates.

Note 2: **the current version of the "trme" package only supports the case when the exposure variable is missing. In addition, the current package only works for the case when the outcome is binary.** If either the covariates or the outcome is also missing, we suggest firstly imputing missing values based on the data. In the future, the package may be updated to deal with other cases such as the continuous outcome.

Note 3: for more technical details, please review and properly cite the paper: **Yuliang Shi, Yeying Zhu, Joel Dubin. Causal Inference on Missing Exposure via Robust Estimation.** published at arxiv: https://arxiv.org/abs/2406.08668.

To contact the author, please visit [Yuliang's Website](https://uwaterloo.ca/scholar/y323shi/home) or send me an email: yuliang.shi@uwaterloo.ca.
