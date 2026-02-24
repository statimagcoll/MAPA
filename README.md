**The package has been renamed from semicor to MAPA**

This R package allows users to compute the one-step semiparametric estimator
for maximum achievable predictive accuracy (MAPA), introduced in [Jones et al. (2025)](https://www.biorxiv.org/content/10.1101/2025.11.26.690778v1). The functions
assume the user has already run their ML models (using cross-fitting).

Correlation (continuous phenotypes) and AUC (binary phenotypes) have been implemented thus far.

The package and documentation is under active development, but the main functionality is available.
Here is a simple example:

```r
# Install package
devtools::install_github("statimagcoll/semicor")
library(semicor)
# Simple example (detailed package documentation in development)
# 10 sample splits, n = 100, 5-fold cross-fitting
set.seed(1)

# Fold indices
# Each column is a different sample split
test_fold = matrix(replicate(10, sample(rep(1:5, each = 20))), ncol = 10)
# Observed phenotype
y = runif(100)
# Predicted values resulting from ML model 1
muhat = matrix(runif(1000), ncol = 10)
# Predicted values resulting from ML model 2
muhat2 = matrix(runif(1000), ncol = 10)


# For only summarizing one ML model
split_res = sampleSplitOS(y, muhat, test_fold)
split_res
# For final estimates, take median
apply(split_res$result[,1:3], 2, median)

# For summarizing and comparing two ML models
split_res2 = sampleSplitOS(y, muhat, test_fold, muhat2)
split_res2 
# For final estimates, take median
apply(split_res2$result[split_res2$result$metric == "yhat1",1:3], 2, median)
apply(split_res2$result[split_res2$result$metric == "yhat2",1:3], 2, median)
apply(split_res2$result[split_res2$result$metric == "Diff",1:3], 2, median)
apply(split_res2$result[split_res2$result$metric == "Ratio",1:3], 2, median)
```
