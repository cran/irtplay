
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtplay\_1.2.0 (2019-10-16)

o Fixed a problem of the function ‘est\_score’, which returned an error
message when only one dichotomous item or one polytomous item was
included in the item meta data set.

o Fixed a problem of the function ‘est\_item’, which returned an error
message when the inverse of hessian matrix is not obtainable.

o Included the ‘maximum likelihood estimation with fences’ scoring
method (Han, 2016) in the function ‘est\_score’.

o Included the ‘inverse test characteristic curve (TCC)’ scoring method
(e.g., Stocking, 1996) in the function ‘est\_score’.

o Included the function ‘llike\_item’ to compute the loglikelihood
values of items.

# irtplay\_1.1.0 (2019-09-15)

o For the function ‘est\_item’, default parameters of a-parameter prior
distribution were revised

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_score’ to estimate an ability in a brute
force way when MLE or MAP fails to find the solution.

o Updated the function ‘irtfit’ to compute the likelihood ratio
chi-square fit statistic (G2; Mckinley & Mills, 1985).

# irtplay\_1.0.0 (2019-08-21)

o initial release on CRAN
