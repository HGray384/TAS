# TAS
Target-Averaged linear Shrinkage R-package (dev version)

Implementation of TAS from Gray et al., (2018) for high-dimensional covariance matrix estimation. Conceptually,
TAS performs linear shrinkage by weighting multiple target matrices by their log-marginal likelihood values in
a Gaussian conjugate (inverse Wishart prior) model. This is a data-driven way of either performing regularisation
through a restricted target matrix, or using a strongly informative prior through a data-derived target matrix 
for potentially greater estimation accuracy. The model weighting procedure eliminates the need to choose a single 
target matrix, which is the case in standard linear shrinkage.

See documentation files for information on the usage of TAS.

This version can be directly installed from the R console using:

```
# install.packages("devtools")
devtools::install_github("HGray384/TAS")
```
