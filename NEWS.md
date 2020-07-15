
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtplay\_1.6.0 (2020-07-14)

The package has been updated significantly in this verstion. In this
version, I have:

o Updated ‘est\_score’ function to estimate ability parameters much
faster than the previous version of the function.

o Updated ‘est\_irt’ and ‘est\_item’ functions to estimate item
parameters much faster than the previous version of the functions.

o Updated ‘test.info’ function to compute items infomation and test
information much faster than the previous version of the function.

o Added an option to use a prior distribution of the item difficulty (or
threshold) parameters in ‘est\_irt’, ‘est\_item’, and ‘llike\_item’
functions.

o Solved unstable item parameter estimation of ‘est\_irt’ and
‘est\_item’ functions which occured when the scaling factor of ‘D’
is other than 1.0 and ‘use.aprior = TRUE’.

o Fixed an error which occured in the function ‘est\_irt’ when the data
set contains missing values and ‘fix.a.1pl = FALSE’.

# irtplay\_1.5.1 (2020-06-16)

o Included ‘summary’ method to summarize the IRT calibration results
from ‘est\_irt’ or ‘est\_item’ objects.

o Included a new function of ‘getirt’ to extract various estimates
results from ‘est\_irt’ or ‘est\_item’ objects.

o Fixed an error which happens when “DRM” is specified in the model name
in the function ‘est\_irt’.

o Included total computation time in the function ‘est\_irt’.

# irtplay\_1.5.0 (2020-04-12)

o Changed the title of ‘irtplay’ package to “Unidimensional Item
Response Theory Modeling”.

o Included a new function of ‘est\_irt’ to fit unidimensional IRT models
to mixture of dichotomous and polytomous item data using the marginal
maximum likelihood estimation with expectation-maximization (MMLE-EM;
Bock & Aitkin, 1981) algorithm.

o Included the fixed item parameter calibration (FIPC; Kim, 2006)
approach, which is one of useful online calibration methods, in the
function ‘est\_irt’.

o Updated the documentation to explain how to implement the new function
‘est\_irt’.

o Included well-known LSAT6 dichotomous response data set from Thissen
(1982).

o Fixed a problem of inaccurate item parameter estimation in the
function ‘est\_item’ when a prior distribution of the slope parameter is
used with a scaling factor other than D = 1.

o Updated the function ‘bring.flexmirt’ to read the item parameters of
the generalized partial credit model when the number of score categories
are two.

o Updated the function ‘est\_score’ to find a smart starting value when
MLE is used. More specifically, the smart starting value is a theta
value where the log-likelihood is the maximum at the highest peak.

# irtplay\_1.4.1 (2020-02-21)

o Included the function ‘run\_flexmirt’ to implement flexMIRT software
(Cai, 2017) through R.

o Applied a prior distribution to the slope parameters of the IRT 1PL
model when the slope parameters are constrained to be equal in the
function of ‘est\_item’.

o Fixed a problem of using staring values to estimate item parameters in
the function of ‘est\_item’.

# irtplay\_1.4.0 (2020-01-23)

o Fixed a non-convergence problem of the maximum likelihood estimation
with fences (MLEF) in the function of ‘est\_score’.

o Updated the description and introduction of the package.

o Updated the documentation to explain how to implement the function
“est\_item” in more detail.

o Updated the README.md file to explain how to implement the function
“est\_item” in more detail.

# irtplay\_1.3.0 (2019-11-17)

o Included the function ‘llike\_score’ to compute the loglikelihood
function of ability for an examinee.

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_item’ to exclude items that contains no
item response data during the item parameter estimation.

o Updated the function ‘est\_item’ to count the number of item responses
for each item used to estimate the item parameters.

o Updated the function ‘est\_score’ to find better starting values when
MLE is used.

o Updated the function ‘est\_score’ to address NaNs of gradient values
and NaNs of hessian values when MLE, MLEF, or MAP is used.

o Fixed a problem of the function ‘est\_score’, which returned an error
message when a vector of an examinee’s response data was used in the
argument of ‘x’.

# irtplay\_1.2.0 (2019-10-16)

o Fixed a problem of the function ‘est\_score’, which returned an error
message when only one dichotomous item or one polytomous item was
included in the item meta data set.

o Fixed a problem of the function ‘est\_item’, which returned an error
message when the inverse of hessian matrix is not obtainable.

o Included the ‘maximum likelihood estimation with fences scoring method
(Han, 2016) in the function ’est\_score’.

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
