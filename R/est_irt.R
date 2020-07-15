#' Item parameter estimation using MMLE-EM algorithm
#'
#' @description This function fits unidimensional item response (IRT) models to mixture of dichotomous and polytomous data using
#' marginal maximum likelihood estimation with expectation-maximization (MMLE-EM) algorithm (Bock & Aitkin, 1981). This function also
#' implements the fixed item parameter calibration (FIPC; Kim, 2006). As Method A (Stocking, 1988), FIPC is one of useful online item
#' calibration methods for computerized adaptive testing (CAT) to put the parameter estimates of pretest items on the same scale of
#' operational item parameter estimates (Ban, Hanson, Wang, Yi, & Harris, 2001). For dichotomous items, IRT one-, two-, and three-parameter
#' logistic models are available. For polytomous items, the graded response model (GRM) and the (generalized) partial credit model (GPCM)
#' are available.
#'
#' @param x A data.frame containing the item meta data. This meta data is necessary to obtain the information of
#' each item (i.e., number of score categories and IRT model) to be calibrated. You can easily create an empty
#' item meta data using the function \code{\link{shape_df}}. When \code{use.startval = TRUE}, the item parameters
#' specified in the item meta data are used as the starting values for the item parameter estimation.
#' If \code{x = NULL}, the arguments of \code{model} and \code{cats} must be specified. Note that when \code{fipc = TRUE}
#' to implement the FIPC method, the item meta data of a test form must be provided in the argument \code{x}.
#' See below for details.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param model A vector of character strings indicating what IRT model is used to calibrate each item. Available IRT models are
#' "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function.
#' If a single character of the IRT model is specified, that model will be recycled across all items. This information is only required
#' when \code{x = NULL} and \code{fipc = FALSE}.
#' @param cats A numeric vector specifying the number of score categories for each item. For example, a dichotomous
#' item has two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items
#' are binary items (i.e., dichotomous items), it assumes that all items have two score categories. This information is only required
#' when \code{x = NULL} and \code{fipc = FALSE}.
#' @param fix.a.1pl A logical value. If TRUE, the slope parameters of the 1PLM items are fixed to a specific value specified in the argument
#' \code{a.val.1pl}. Otherwise, the slope parameters of all 1PLM items are constrained to be equal and estimated. Default is FALSE.
#' @param fix.a.gpcm A logical value. If TRUE, the GPCM items are calibrated with the partial credit model and the slope parameters of
#' the GPCM items are fixed to a specific value specified in the argument \code{a.val.gpcm}. Otherwise, the slope parameter of each GPCM item
#' is estimated. Default is FALSE.
#' @param fix.g A logical value. If TRUE, the guessing parameters of the 3PLM items are fixed to a specific value specified in the argument
#' \code{g.val}. Otherwise, the guessing parameter of each 3PLM item is estimated. Default is FALSE.
#' @param a.val.1pl A numeric value. This value is used to fixed the slope parameters of the 1PLM items.
#' @param a.val.gpcm A numeric value. This value is used to fixed the slope parameters of the GPCM items.
#' @param g.val A numeric value. This value is used to fixed the guessing parameters of the 3PLM items.
#' @param use.aprior A logical value. If TRUE, a prior distribution for the slope parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.bprior A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.gprior A logical value. If TRUE, a prior distribution for the guessing parameters is used for the parameter calibration
#' across all 3PLM items. Default is TRUE.
#' @param aprior A list containing the information of the prior distribution for item slope parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param bprior A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param gprior A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param Quadrature A numeric vector of two components specifying the number of quadrature points (in the first component) and
#' the symmetric minimum and maximum values of these points (in the second component). For example, a vector of c(49, 6) indicates 49 rectangular
#' quadrature points over -6 and 6. The quadrature points are used in the E step of the EM algorithm. Default is c(49, 6).
#' @param weights A two-column matrix or data.frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL, a normal prior density is used based on the information provided in the arguments
#' of \code{Quadrature}, \code{group.mean}, and \code{group.var}). Default is NULL.
#' @param group.mean A numeric value to set the mean of latent variable prior distribution. Default is 0. This value is fixed to remove
#' the indeterminancy of item parameter scale when calibrating items. However, the scale of prior distribution is updated when FIPC is implemented.
#' @param group.var A positive numeric value to set the variance of latent variable prior distribution. Default is 1. This value is fixed to remove
#' the indeterminancy of item parameter scale when calibrating items. However, the scale of prior distribution is updated when FIPC is implemented.
#' @param EmpHist A logical value. If TRUE, the empirical histogram of the latent variable prior distribution is simultaneously estimated with
#' the item parameters using Woods's (2007) approach. The item parameters are calibrated against the estimated empirical histogram prior distribution.
#' See below for details.
#' @param use.startval A logical value. If TRUE, the item parameters provided in the item meta data (i.e., the argument \code{x}) are used as
#' the starting values for the item parameter estimation. Otherwise, internal starting values of this function are used. Default is FALSE.
#' @param Etol A positive numeric value. This value sets the convergence criterion for E steps of the EM algorithm. Default is 1e-4.
#' @param MaxE A positive integer value. This value determines the maximum number of the E steps in the EM algorithm. Default is 500.
#' @param control A list of control parameters to be passed to the optimization function of \code{nlminb()} in the \pkg{stats} package. The control parameters
#' set the conditions of M steps of the EM algorithm. For example, the maximum number of iterations in each of the iterative M steps can
#' be set by \code{control = list(iter.max=200)}. Default maximum number of iterations in each M step is 200. See \code{nlminb()} in the \pkg{stats} package
#' for other control parameters.
#' @param fipc A logical value. If TRUE, FIPC is implemented for item parameter estimation. See below for details.
#' @param fipc.method A character string specifying the FIPC method. Available methods include "OEM" for "No Prior Weights Updating and One EM Cycle
#' (NWU-OEM; Wainer & Mislevy, 1990)" and "MEM" for "Multiple Prior Weights Updating and Multiple EM Cycles (MWU-MEM; Kim, 2006)."
#' When \code{fipc.method = "OEM"}, the maximum number of the E steps of the EM algorithm is set to 1 no matter what number is specified
#' in the argument \code{MaxE}.
#' @param fix.loc A vector of positive integer values specifying the location of the items to be fixed in the item meta data (i.e., \code{x})
#' when the FIPC is implemented. For example, suppose that five items located in the 1st, 2nd, 4th, 7th, and 9th rows of the item meta data \code{x}
#' should be fixed. Then \code{fix.loc = c(1, 2, 4, 7, 9)}.
#'
#' @details A specific form of a data.frame should be used for the argument \code{x}. The first column should have item IDs,
#' the second column should contain the number of score categories of the items, and the third column should include IRT models.
#' The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items.
#' Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. From the fourth column, item parameters should be included.
#' For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
#' item guessing parameters, respectively. When "1PLM" or "2PLM" is specified for any items in the third column, NAs should be inserted
#' for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be contained in the
#' fourth column and the item threshold (or step) parameters should be included from the fifth to the last columns.
#' When the number of categories differs between items, the empty cells of item parameters should be filled with NAs.
#' In this package, item step parameters should be used for the (generalized) partial credit model. The item step parameter is the
#' overall item difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for each category.
#' Thus, the number of step parameters for an item with m categories is m-1 because a step parameter for the first category does not
#' affect the category probabilities. For example, if an item has five categories under the (generalized) partial credit model,
#' four step parameters should be specified. An example of a data.frame with a single-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#' And an example of a data.frame for a mixed-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
#'   ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
#'   ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
#'   ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
#'   ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
#'   ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
#' }
#' For more details about the parameterization of the (generalized) partial credit model, see \code{IRT Models} section in
#' the page of \code{\link{irtplay-package}} for more details about the IRT models. An easier way to create a data.frame for
#' the argument \code{x} is by using the function \code{\link{shape_df}}.
#'
#' To fit the IRT models to data, the IRT model and the number of score category information for the estimated items must be provided as well as
#' the item response data. There are two way to provide the IRT model and score category information. The first way is to provide the item meta data
#' to the argument \code{x}. As explained above, the item meta data can be easily created by the function \code{\link{shape_df}}. The second way is
#' specify the IRT models and the score category information into the arguments of \code{model} and \code{cats}. Thus, if \code{x=NULL}, the specified
#' information in \code{model} and \code{cats} are used.
#'
#' To implement FIPC, however, the item meta data must be provided in the argument \code{x}. This is because the item parameters of the fixed items
#' in the item meta data are used to estimate the characteristic of the underlying latent variable prior distribution when calibrating the rest of freely estimated items.
#' More specifically, the underlying latent variable prior distribution of the fixed items is estimated during the calibration of the freely estimated items
#' to put the item parameters of the freely estimated items on the scale of the fixed item parameters (Kim, 2006).
#'
#' In terms of approaches for FIPC, Kim (2006) described five different methods. Among them, two methods are available in the
#' function \code{\link{est_irt}}. The first method is "NWU-OEM" where uses just one E step in the EM algorithm, involving data from only the fixed items, and
#' just one M step, involving data from only non-fixed items. This method is suggested by Wainer and Mislevy (1990) in the context of online calibration. This method
#' can be implemented by setting \code{fipc.method = "OEM"}. The second method is "MWU-MEM" which iteratively updates the latent variable prior distribution and
#' finds the parameter estimates of the non-fixed items. In this method, the first EM cycle is conducted in the same way as in the NWU-OEM method. From the second
#' EM cycle, both the parameters of all items (fixed and non-fixed) and the weights of the prior distribution are concurrently updated. See Kim (2006) for more details.
#'
#' When \code{EmpHist = TRUE}, the emprical histogram of latent variable prior distribution is simultaneously estimated with the item parameters. If \code{fipc = TRUE}
#' given \code{EmpHist = TRUE}, the scale parameters (e.g., mean and variance) of the empirical prior distribution are estimated as well. If \code{fipc = FALSE} given
#' \code{EmpHist = TRUE}, the scale parameters of the empirical prior distribution are fixed to the values specified in the arguments of \code{group.mean} and \code{group.var}.
#' When \code{EmpHist = FALSE}, the normal prior distribution is used during the item parameter estimation. If \code{fipc = TRUE} given \code{EmpHist = FALSE},
#' the scale parameters of the normal prior distribution are estimated as well as the item parameters. If \code{fipc = FALSE} given \code{EmpHist = FALSE},
#' the scale parameters of the normal prior distribution are fixed to the values specified in the arguments of \code{group.mean} and \code{group.var}.
#'
#' @return This function returns an object of class \code{\link{est_irt}}. Within this object, several internal objects are contained such as:
#' \item{estimates}{A data.frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
#' \item{par.est}{A data.frame containing the item parameter estimates.}
#' \item{se.est}{A data.frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
#' observed information functions. The standard errors are estimated using the cross-production approximation method (Meilijson, 1989).}
#' \item{pos.par}{A data.frame containing the position number of item parameters being estimated. The position information is useful
#' when interpreting the variance-covariance matrix of item parameter estimates.}
#' \item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
#' \item{loglikelihood}{A sum of the log-likelihood values of the observed data set (marginal log-likelihood) across all estimated items.}
#' \item{group.par}{A data.frame containing the mean, variance, and standard deviation of latent variable prior distribution.}
#' \item{weights}{A two-column matrix or data.frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the (updated) latent variable prior distribution.}
#' \item{posterior.dist}{A matrix of normalized posterior densities for all the response patterns at each of the quadrature points.
#' The row and column indicate the response pattern and the quadrature point, respectively.}
#' \item{data}{A data.frame of the examinees' response data set.}
#' \item{scale.D}{A scaling factor in IRT models.}
#' \item{ncase}{A total number of response patterns.}
#' \item{nitem}{A total number of items included in the response data.}
#' \item{Etol}{A convergence criteria for E steps of the EM algorithm.}
#' \item{MaxE}{The maximum number of E steps in the EM algorithm.}
#' \item{aprior}{A list containing the information of the prior distribution for item slope parameters.}
#' \item{gprior}{A list containing the information of the prior distribution for item guessing parameters.}
#' \item{npar.est}{A total number of the estimated item parameters.}
#' \item{niter}{The number of EM cycles completed.}
#' \item{maxpar.diff}{A maximum item parameter change when the EM cycles were completed.}
#' \item{EMtime}{Time (in seconds) spent for the EM cycles.}
#' \item{SEtime}{Time (in seconds) spent for computing the standard errors of the item parameter estimates.}
#' \item{TotalTime}{Time (in seconds) spent for total compuatation.}
#' \item{test.1}{Status of the first-order test to report if the gradients has vanished sufficiently for the solution to be stable.}
#' \item{test.2}{Status of the second-order test to report if the information matrix is positive definite, which is a prerequisite
#' for the solution to be a possible maximum.}
#' \item{var.note}{A note to report if the variance-covariance matrix of item parameter estimates is obtainable from the information matrix.}
#' \item{fipc}{A logical value to indicate if FIPC was used.}
#' \item{fipc.method}{A method used for the FIPC.}
#' \item{fix.loc}{A vector of integer values specifying the location of the fixed items when the FIPC was implemented.}
#'
#' The internal objects can be easily extracted using the function \code{\link{getirt}}.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{est_item}}, \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{sx2_fit}},
#' \code{\link{traceline.est_item}}, \code{\link{getirt}}
#'
#' @references
#' Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
#' in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm.
#' \emph{Psychometrika, 46}, 443-459.
#'
#' Kim, S. (2006). A comparative study of IRT fixed parameter calibration methods.
#' \emph{Journal of Educational Measurement, 43}(4), 355-381.
#'
#' Meilijson, I. (1989). A fast improvement to the EM algorithm on its own terms.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological), 51}, 127-138.
#'
#' Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.
#'
#' Wainer, H., & Mislevy, R. J. (1990). Item response theory, item calibration, and proficiency estimation. In H. Wainer (Ed.),
#' \emph{Computer adaptive testing: A primer} (Chap. 4, pp.65-102). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Woods, C. M. (2007). Empirical histograms in item response theory with ordinal data. \emph{Educational and Psychological Measurement, 67}(1), 73-87.
#'
#' @examples
#' \donttest{
#'
#' ##------------------------------------------------------------------------------
#' # 1. item parameter estimation for the dichotomous item data (LSAT6)
#' ##------------------------------------------------------------------------------
#' # fit the 1PL model to LSAT6 data and constrain the slope parameters to be equal
#' (mod.1pl.c <- est_irt(data=LSAT6, D=1, model="1PLM", cats=2, fix.a.1pl=FALSE))
#'
#' # summary of the estimation
#' summary(mod.1pl.c)
#'
#' # extract the item parameter estimates
#' getirt(mod.1pl.c, what="par.est")
#'
#' # extract the standard error estimates
#' getirt(mod.1pl.c, what="se.est")
#'
#' # fit the 1PL model to LSAT6 data and fix the slope parameters to 1.0
#' (mod.1pl.f <- est_irt(data=LSAT6, D=1, model="1PLM", cats=2, fix.a.1pl=TRUE, a.val.1pl=1))
#'
#' # summary of the estimation
#' summary(mod.1pl.f)
#'
#' # fit the 2PL model to LSAT6 data
#' (mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2))
#'
#' # summary of the estimation
#' summary(mod.2pl)
#'
#' # assess the fit of the 2PL model to the LSAT5 data using S-X2 fit statistic
#' (sx2fit.2pl <- sx2_fit(x=mod.2pl))
#'
#' # compute the item and test information at several theta points
#' theta <- seq(-4, 4, 0.1)
#' (info.2pl <- test.info(x=mod.2pl, theta=theta))
#'
#' # draw the test characteristic curve plot
#' (trace.2pl <- traceline(x=mod.2pl, theta=theta))
#' plot(trace.2pl)
#'
#' # draw the item characteristic curve for the 1st item
#' plot(trace.2pl, item.loc=1)
#'
#' # fit the 2PL model to LSAT6 data and
#' # estimate the empirical histogram of latent variable prior distribution
#' # also use a less stringent convergence criterion for E-step
#' (mod.2pl.hist <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2, EmpHist=TRUE, Etol=0.001))
#' (emphist <- getirt(mod.2pl.hist, what="weights"))
#' plot(emphist$weight ~ emphist$theta)
#'
#' # fit the 3PL model to LSAT6 data and use the Beta prior distribution for
#' # the guessing parameters
#' (mod.3pl <- est_irt(data=LSAT6, D=1, model="3PLM", cats=2, use.gprior=TRUE,
#'                     gprior=list(dist="beta", params=c(5, 16))))
#'
#' # summary of the estimation
#' summary(mod.3pl)
#'
#' # fit the 3PL model to LSAT6 data, but fix the guessing parameters to be 0.2
#' (mod.3pl.f <- est_irt(data=LSAT6, D=1, model="3PLM", cats=2, fix.g=TRUE, g.val=0.2))
#'
#' # summary of the estimation
#' summary(mod.3pl.f)
#'
#' # fit the differnt dichotomous models to each item of LSAT6 data
#' # fit the constrained 1PL model to the 1st, 2nd, and 3rd items, fit the 2PL model to
#' # the 4th item, and fit the 3PL model to the 5th item with the Beta prior of
#' # the guessing parameter
#' (mod.drm.mix <- est_irt(data=LSAT6, D=1, model=c("1PLM", "1PLM", "1PLM", "2PLM", "3PLM"),
#'                         cats=2, fix.a.1pl=FALSE, use.gprior=TRUE,
#'                         gprior=list(dist="beta", params=c(5, 16))))
#' # summary of the estimation
#' summary(mod.drm.mix)
#'
#' ##------------------------------------------------------------------------------
#' # 2. item parameter estimation for the mixed-item format data (simulation data)
#' ##------------------------------------------------------------------------------
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # modify the item meta data so that the 39th and 40th items follow GPCM
#' x[39:40, 3] <- "GPCM"
#'
#' # generate 1,000 examinees' latent abilities from N(0, 1)
#' set.seed(37)
#' score1 <- rnorm(1000, mean=0, sd=1)
#'
#' # simulate the response data
#' sim.dat1 <- simdat(x=x, theta=score1, D=1)
#'
#' # fit the 3PL model to all dichotmous items, fit the GPCM model to 39th and 40th items,
#' # and fit the GRM model to the 53th, 54th, 55th items.
#' # use the beta prior distribution for the guessing parameters, use the log-normal
#' # prior distribution for the slope parameters, and use the normal prior distribution
#' # for the difficulty (or thereshold) parameters.
#' # also, specify the argument 'x' to provide the IRT model and score category information
#' # for items
#' item.meta <- shape_df(item.id=x$id, cats=x$cats, model=x$model, empty.par=TRUE)
#' (mod.mix1 <- est_irt(x=item.meta, data=sim.dat1, D=1, use.aprior=TRUE, use.bprior=TRUE,
#'                      use.gprior=TRUE,
#'                      aprior=list(dist="lnorm", params=c(0.0, 0.5)),
#'                      bprior=list(dist="norm", params=c(0.0, 2.0)),
#'                      gprior=list(dist="beta", params=c(5, 16))))
#'
#' # summary of the estimation
#' summary(mod.mix1)
#'
#' # estimate examinees' latent scores given the item parameter estimates using the MLE
#' (score.mle <- est_score(x=mod.mix1, method = "MLE", range = c(-4, 4), ncore=2))
#'
#' # compute the traditional fit statistics
#' (fit.mix1 <- irtfit(x=mod.mix1, score=score.mle$est.theta, group.method="equal.width",
#'                     n.width=10, loc.theta="middle"))
#'
#' # residual plots for the first item (dichotomous item)
#' plot(x=fit.mix1, item.loc=1, type = "both", ci.method = "wald",
#'      show.table=TRUE, ylim.sr.adjust=TRUE)
#'
#' # residual plots for the last item (polytomous item)
#' plot(x=fit.mix1, item.loc=55, type = "both", ci.method = "wald",
#'      show.table=FALSE, ylim.sr.adjust=TRUE)
#'
#' # fit the 2PL model to all dichotmous items, fit the GPCM model to 39th and 40th items,
#' # and fit the GRM model to the 53th, 54th, 55th items.
#' # also, specify the arguments of 'model' and 'cats' to provide the IRT model and
#' # score category information for items
#' (mod.mix2 <- est_irt(data=sim.dat1, D=1,
#'                      model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'                      cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3))))
#'
#' # summary of the estimation
#' summary(mod.mix2)
#'
#' # fit the 2PL model to all dichotmous items, fit the GPCM model to 39th and 40th items,
#' # fit the GRM model to the 53th, 54th, 55th items, and estimate the empirical histogram
#' # of latent variable prior distribution.
#' # also, specify the arguments of 'model' and 'cats' to provide the IRT model and
#' # score category information for items
#' (mod.mix3 <- est_irt(data=sim.dat1, D=1,
#'                      model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'                      cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)), EmpHist=TRUE))
#' (emphist <- getirt(mod.mix3, what="weights"))
#' plot(emphist$weight ~ emphist$theta)
#'
#' # fit the 2PL model to all dichotmous items,
#' # fit the PCM model to 39th and 40th items by fixing the slope parameters to 1,
#' # and fit the GRM model to the 53th, 54th, 55th items.
#' # also, specify the arguments of 'model' and 'cats' to provide the IRT model and
#' # score category information for items
#' (mod.mix4 <- est_irt(data=sim.dat1, D=1,
#'                      model=c(rep("2PLM", 38), rep("GPCM", 2), rep("2PLM", 12), rep("GRM", 3)),
#'                      cats=c(rep(2, 38), rep(5, 2), rep(2, 12), rep(5, 3)),
#'                      fix.a.gpcm=TRUE, a.val.gpcm=1))
#'
#' # summary of the estimation
#' summary(mod.mix4)
#'
#' ##------------------------------------------------------------------------------
#' # 3. fixed item parameter calibration (FIPC) for the mixed-item format data
#' #    (simulation data)
#' ##------------------------------------------------------------------------------
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # generate 1,000 examinees' latent abilities from N(0.4, 1.3)
#' set.seed(20)
#' score2 <- rnorm(1000, mean=0.4, sd=1.3)
#'
#' # simulate the response data
#' sim.dat2 <- simdat(x=x, theta=score2, D=1)
#'
#' # fit the 3PL model to all dichotmous items, fit the GRM model to all polytomous data,
#' # fix the five 3PL items (1st - 5th items) and three GRM items (53th to 55th items)
#' # also, estimate the empirical histogram of latent variable
#' # use the MEM method.
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix1 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
#'                      gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE,
#'                      Etol=1e-3, fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
#' (prior.par <- mod.fix1$group.par)
#' (emphist <- getirt(mod.fix1, what="weights"))
#' plot(emphist$weight ~ emphist$theta)
#'
#' # summary of the estimation
#' summary(mod.fix1)
#'
#' # fit the 3PL model to all dichotmous items, fit the GRM model to all polytomous data,
#' # fix the five 3PL items (1st - 5th items) and three GRM items (53th to 55th items)
#' # at this moment, do estimate the empirical histogram of latent variable.
#' # instead, estimate the scale of normal prior distribution of latent variable
#' # use the MEM method.
#' fix.loc <- c(1:5, 53:55)
#' (mod.fix2 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
#'                      gprior=list(dist="beta", params=c(5, 16)), EmpHist=FALSE,
#'                      Etol=1e-3, fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc))
#' (prior.par <- mod.fix2$group.par)
#' (emphist <- getirt(mod.fix2, what="weights"))
#' plot(emphist$weight ~ emphist$theta)
#'
#' # fit the 3PL model to all dichotmous items, fit the GRM model to all polytomous data,
#' # at this moment fix only the five 3PL items (1st - 5th items)
#' # and estimate the empirical histogram of latent variable.
#' # use the OEM method. Thus, only 1 EM cycle is used.
#' fix.loc <- c(1:5)
#' (mod.fix3 <- est_irt(x=x, data=sim.dat2, D=1, use.gprior=TRUE,
#'                      gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE,
#'                      Etol=1e-3, fipc=TRUE, fipc.method="OEM", fix.loc=fix.loc))
#' (prior.par <- mod.fix3$group.par)
#' (emphist <- getirt(mod.fix3, what="weights"))
#' plot(emphist$weight ~ emphist$theta)
#'
#' # summary of the estimation
#' summary(mod.fix3)
#'
#' }
#'
#' @import purrr
#' @import dplyr
#'
#' @export
#'
est_irt <- function(x=NULL, data, D=1, model=NULL, cats=NULL, fix.a.1pl=FALSE, fix.a.gpcm=FALSE, fix.g=FALSE,
                    a.val.1pl=1, a.val.gpcm=1, g.val=.2, use.aprior=FALSE, use.bprior=FALSE, use.gprior=TRUE,
                    aprior=list(dist="lnorm", params=c(0.0, 0.5)), bprior=list(dist="norm", params=c(0.0, 1.0)),
                    gprior=list(dist="beta", params=c(5, 16)), missing=NA, Quadrature=c(49, 6.0), weights=NULL,
                    group.mean=0.0, group.var=1.0, EmpHist=FALSE,
                    use.startval=FALSE, Etol=1e-04, MaxE=500, control=list(iter.max=200),
                    fipc=FALSE, fipc.method="MEM", fix.loc=NULL) {

  # match.call
  cl <- match.call()

  # item parameter estimation
  if(!fipc) {
    # item parameter estimation using MMLE-EM algorithm
    est_par <- est_irt_em(x=x, data=data, D=D, model=model, cats=cats, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g,
                          a.val.1pl=a.val.1pl, a.val.gpcm=a.val.gpcm, g.val=g.val, use.aprior=use.aprior, use.bprior=use.bprior,
                          use.gprior=use.gprior, aprior=aprior, bprior=bprior, gprior=gprior, missing=missing, Quadrature=Quadrature,
                          weights=weights, group.mean=group.mean, group.var=group.var, EmpHist=EmpHist, use.startval=use.startval,
                          Etol=Etol, MaxE=MaxE, control=control)

  } else {
    # implement FIPC method
    est_par <- est_irt_fipc(x=x, data=data, D=D, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g,
                            a.val.1pl=a.val.1pl, a.val.gpcm=a.val.gpcm, g.val=g.val, use.aprior=use.aprior, use.bprior=use.bprior,
                            use.gprior=use.gprior, aprior=aprior, bprior=bprior, gprior=gprior, missing=missing, Quadrature=Quadrature,
                            weights=weights, group.mean=group.mean, group.var=group.var, EmpHist=EmpHist, use.startval=use.startval,
                            Etol=Etol, MaxE=MaxE, control=control, fipc=TRUE, fipc.method=fipc.method, fix.loc=fix.loc)

  }

  # return the estimation results
  class(est_par) <- "est_irt"
  est_par$call <- cl
  est_par

}


# This function estimates item parameters using MMLE-EM algorithm
est_irt_em <- function(x=NULL, data, D=1, model=NULL, cats=NULL, fix.a.1pl=FALSE, fix.a.gpcm=FALSE, fix.g=FALSE,
                       a.val.1pl=1, a.val.gpcm=1, g.val=.2, use.aprior=FALSE, use.bprior=FALSE, use.gprior=TRUE,
                       aprior=list(dist="lnorm", params=c(0.0, 0.5)), bprior=list(dist="norm", params=c(0.0, 1.0)),
                       gprior=list(dist="beta", params=c(5, 16)), missing=NA, Quadrature=c(49, 6.0), weights=NULL,
                       group.mean=0, group.var=1, EmpHist=FALSE, use.startval=FALSE, Etol=1e-04, MaxE=500,
                       control=list(eval.max=200, iter.max=200)) {

  # check start time
  start.time <- Sys.time()

  ##---------------------------------------------------------------------
  # prepare the item parameter estimation
  ##---------------------------------------------------------------------
  # check if the starting values are available
  if(use.startval & is.null(x)) {
    stop("To use starting values for item parameter estimation, the item meta data must be specified in the argument 'x'.", call.=FALSE)
  }

  # extract information about the number of score cetegories and models
  cat("Parsing input...", '\n')
  if(!is.null(x)) {
    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

    # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
    if(ncol(x[, -c(1, 2, 3)]) == 2) {
      x <- data.frame(x, par.3=NA)
    }

    # clear the item meta data set
    x <- back2df(metalist2(x))
    id <- x[, 1]
    cats <- x[, 2]
    model <-
      as.character(x[, 3]) %>%
      toupper()
    if(!use.startval) {
      x <- startval_df(cats=cats, model=model)
      x$id <- id
    }
  } else {
    model <- toupper(model)
    if(length(model) == 1) {
      model <- rep(model, ncol(data))
    }
    if(is.null(cats)) {
      if(all(model %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
        cats <- rep(2, ncol(data))
      } else {
        stop("The number of score categories for the items should be specified in the argument 'cats'.", call.=FALSE)
      }
    }
    if(length(cats) == 1) {
      cats <- rep(cats, ncol(data))
    }
    x <- startval_df(cats=cats, model=model)
    id <- x$id
  }

  # check whether included data are correct
  if(nrow(x) != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)

  # consider DRM as 3PLM
  if("DRM" %in% model) {
    model[model == "DRM"] <- "3PLM"
    x$model <- model
    memo <- "All 'DRM' items were considered as '3PLM' items during the item parameter estimation. \n"
    warning(memo, call.=FALSE)
  }

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if(length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call.=TRUE)
  }

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
  if("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x)
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(x=x, loc_1p_const=loc_1p_const, loc_else=loc_else,
                      fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g)

  ##---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ##---------------------------------------------------------------------
  # check the total number of examinees
  nstd <- nrow(data)

  # create initial weights of prior ability distribution when it is not specified
  if(is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out=Quadrature[1])

    # create the data.frame containing the quadurature points and weights
    weights <- gen.weight(dist="norm", mu=group.mean, sigma=sqrt(group.var), theta=quadpt)

  } else {
    quadpt <- weights[, 1]
  }

  # factorize the response values and create a list of item responses across items
  resp <- purrr::map2(.x=data.frame(data), .y=cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))

  # calculate the score categories for each examinee
  std.id <- 1:nstd
  freq.cat <- purrr::map(.x=resp, .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))

  # transform the score category data.frame to a matrix
  freq.cat <- purrr::map(.x=freq.cat,
                         .f=function(k) data.matrix(k) %>%
                           unname())

  # delete 'resp' object
  rm(resp, envir=environment(), inherits = FALSE)

  # divide the data set for the mixed-item format
  datlist <- divide_data(data=data, cats=cats, freq.cat=freq.cat)
  data1_drm <- datlist$data1_drm
  data2_drm <- datlist$data2_drm
  data_plm <- datlist$data_plm

  # delete 'datlist' object
  rm(datlist, envir=environment(), inherits = FALSE)

  # estimation
  cat("Estimating item parameters...", '\n')

  # implement EM algorithm
  time1 <- Sys.time()
  par.history <- list()
  for(r in 1:MaxE) {

    # implement E-step
    estep <- Estep(x=x, data1_drm=data1_drm, data2_drm=data2_drm, data_plm=data_plm, freq.cat=freq.cat, weights=weights, D=D)

    # implement M-step
    mstep <- Mstep(estep=estep, id=id, cats=cats, model=model, quadpt=quadpt, D=D, loc_1p_const=loc_1p_const, loc_else=loc_else,
                   EmpHist=EmpHist, weights=weights, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g, a.val.1pl=a.val.1pl,
                   a.val.gpcm=a.val.gpcm, g.val=g.val, use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                   aprior=aprior, bprior=bprior, gprior=gprior,  group.mean=group.mean, group.var=group.var, nstd=nstd,
                   Quadrature=Quadrature, use.startval=TRUE, control=control, iter=r, fipc=FALSE,
                   reloc.par=param_loc$reloc.par, info.mstep=FALSE)

    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$par_df[, -c(1, 2, 3)] - x[, -c(1, 2, 3)]
    max.diff <- abs(max(diff_par, na.rm=TRUE))

    # loglikehood value
    llike <- mstep$loglike

    # print
    cat("\r", paste0("EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall=4), ", Max-Change: ", format(round(max.diff, 6), nsmall=5)))

    # check the convergence of EM algorithm
    converge <- max.diff <= Etol

    # end the EM step if the convergence criterion is satisfied
    if(converge | r == MaxE) {
      # extract the final item parameter estimates
      par_df <- mstep$par_df
      par.history[[r]] <- mstep$par_vec

      # extract the final quadrature points and the corresponding weights of the prior population density
      weights <- mstep$weights
      break

    } else {
      # extract the interim item parameter estimates
      x <- mstep$par_df
      par.history[[r]] <- mstep$par_vec

      # extract the interim quadrature points and the corresponding weights of the prior population density
      weights <- mstep$weights
    }

  }
  cat("", "\n")
  time2 <- Sys.time()

  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # the first order test: check convergence-criteria test
  test_1st <- all(c(all(mstep$convergence == 0L), r < MaxE))
  if(test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3,  " \n"), call.=FALSE)
  }

  # conduct one more E-step to update the posterior distribution using the final item parameter estimates
  estep <- Estep(x=x, data1_drm=data1_drm, data2_drm=data2_drm, data_plm=data_plm, freq.cat=freq.cat, weights=weights, D=D)

  # compute the final log of marginal likelihood
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2])))

  ##---------------------------------------------------------------------
  # estimates the information matrix and standard errors of item parameter estimates
  ##---------------------------------------------------------------------
  # listrize the item parameter estimate data.frame
  meta <- metalist2(par_df)

  # extract the finalized posterior density
  post_dist <- estep$post_dist

  # compute the infomation matrix of item parameter estimates using the cross-product method
  cat("Computing item parameter var-covariance matrix...", '\n')
  time1 <- Sys.time()
  info.data <- info_xpd(meta=meta, freq.cat=freq.cat, post_dist=post_dist, cats=cats, model=model, quadpt=quadpt,
                        D=D, loc_1p_const=loc_1p_const, loc_else=loc_else, nstd=nstd,
                        fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g, a.val.1pl=a.val.1pl,
                        a.val.gpcm=a.val.gpcm, g.val=g.val, reloc.par=param_loc$reloc.par)

  # compute the infomation matrix of item parameter priors
  info.prior <- info_prior(meta=meta, cats=cats, model=model, D=D, loc_1p_const=loc_1p_const,
                           loc_else=loc_else, nstd=nstd, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g,
                           a.val.1pl=a.val.1pl, a.val.gpcm=a.val.gpcm, g.val=g.val, aprior=aprior, bprior=bprior,
                           gprior=gprior, use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                           reloc.par=param_loc$reloc.par)

  # sum of two information matrixs
  info.mat <- info.data + info.prior

  # the second-order test: check if the information matrix is positive definite
  test_2nd <- all(eigen(info.mat, only.values=TRUE)$values > 1e-20)
  if(test_2nd) {
    if(test_1st) {
      memo4 <- "Solution is a possible local maximum."
    } else {
      memo4 <- "Information matrix of item parameter estimates is positive definite."
    }
  } else {
    memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
    warning(paste0(memo4, " \n"), call.=FALSE)
  }

  # compute the variance-covariance matrix of the item parameter estimates, and
  # check if the hessian matrix can be inversed
  cov_mat <- suppressWarnings(tryCatch({solve(info.mat, tol=1e-200)}, error = function(e) {NULL}))

  # compute the standard errors of item parameter estimates
  if(is.null(cov_mat)) {
    se_par <- rep(99999, length(diag(info.mat)))
    memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
    warning(paste0(memo5, " \n"), call.=FALSE)
  } else {
    se_par <- suppressWarnings(sqrt(diag(cov_mat)))
    memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
  }

  # prevent showing NaN values of standard errors
  if(any(is.nan(se_par))) {
    se_par[is.nan(se_par)] <- 99999
  }

  # set an upper bound of standard error
  se_par <- ifelse(se_par > 99999, 99999, se_par)
  time2 <- Sys.time()

  # record the standard error computationtime
  est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # deploy the standard errors on the location of matrix as the item parameter estimates
  se_df <- loc.par <- param_loc$loc.par
  for(i in 1:nrow(loc.par)) {
    num.loc <- which(!is.na(loc.par[i, ]))
    se.loc <- loc.par[i, ][num.loc]
    se_df[i, num.loc] <- se_par[se.loc]
  }

  # create a full data.frame for the standard error estimates
  se_df <- data.frame(par_df[, 1:3], se_df)
  colnames(se_df) <- c("id", "cats", "model", paste0("se.", 1:ncol(loc.par)))

  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  loc_df <- data.frame(par_df[, 1:3], loc.par)
  colnames(loc_df) <- c("id", "cats", "model", paste0("par.", 1:ncol(loc.par)))

  ##---------------------------------------------------------------------
  # summarize the estimation results
  ##---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow=nrow(loc.par), ncol=2*ncol(loc.par)))
  all_df[, seq(1, 2*ncol(loc.par), 2)] <- par_df[, -c(1:3)]
  all_df[, seq(2, 2*ncol(loc.par), 2)] <- se_df[, -c(1:3)]
  col.names <- rep(NA, 2*ncol(loc.par))
  col.names[seq(1, 2*ncol(loc.par), 2)] <- paste0("par.", 1:ncol(loc.par))
  col.names[seq(2, 2*ncol(loc.par), 2)] <- paste0("se.", 1:ncol(loc.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x[, 1:3], all_df)

  # population density parameters
  moments <- c(mu=group.mean, sigma2=group.var, sigma=sqrt(group.var))
  moments.se <- rep(NA, 3)

  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")

  # prior information
  if(use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if(use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if(use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  ##---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- list(estimates=full_all_df, par.est=par_df, se.est=se_df, pos.par=loc_df, covariance=cov_mat, loglikelihood=llike, group.par=group.par,
              weights=weights, posterior.dist=post_dist, data=data, scale.D=D, ncase=nstd, nitem=nrow(par_df), Etol=Etol, MaxE=MaxE,
              aprior=aprior.dist, bprior=bprior.dist, gprior=gprior.dist, npar.est=length(param_loc$reloc.par), niter=r, maxpar.diff=max.diff,
              EMtime=est_time1, SEtime=est_time2, TotalTime=est_time3, test.1=memo3, test.2=memo4, var.note=memo5, fipc=FALSE,
              fipc.method=NULL, fix.loc=NULL)

  cat("Estimation is finished in", est_time3, "seconds.",'\n')
  return(rst)

}


# This function implement the fixed item parameter calibration (FIPC) using MMLE-EM algorithm
est_irt_fipc <- function(x=NULL, data, D=1, fix.a.1pl=FALSE, fix.a.gpcm=FALSE, fix.g=FALSE,
                         a.val.1pl=1, a.val.gpcm=1, g.val=.2, use.aprior=FALSE, use.bprior=FALSE, use.gprior=TRUE,
                         aprior=list(dist="lnorm", params=c(0.0, 0.5)), bprior=list(dist="norm", params=c(0.0, 1.0)),
                         gprior=list(dist="beta", params=c(5, 16)), missing=NA, Quadrature=c(49, 6.0), weights=NULL,
                         group.mean=0.0, group.var=1.0, EmpHist=FALSE, use.startval=FALSE, Etol=1e-04, MaxE=500,
                         control=list(eval.max=200, iter.max=200), fipc=TRUE, fipc.method="MEM", fix.loc=NULL) {

  # check start time
  start.time <- Sys.time()

  ##---------------------------------------------------------------------
  # prepare the item parameter estimation
  ##---------------------------------------------------------------------
  # check if item meta data argument of 'x' is provided
  if(is.null(x)) {
    stop("To implement the fixed item parameter calibration, the item meta data must be specified in the argument 'x'.", call.=FALSE)
  }

  # extract information about the number of score cetegories and models
  cat("Parsing input...", '\n')

  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }

  # clear the item meta data set
  x <- back2df(metalist2(x))

  # transform the model variables into character type
  x$model <-
    as.character(x$model) %>%
    toupper()

  # consider DRM as 3PLM
  if("DRM" %in% x$model) {
    x$model[x$model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items were considered as '3PLM' items during the item parameter estimation. \n"
    warning(memo, call.=FALSE)
  }

  # count the number of items
  nitem <- nrow(x)

  # check the location of items whose item parameters are estimated
  nofix.loc <- c(1:nitem)[!c(1:nitem) %in% fix.loc]

  # divide the meta data.frame into two groups: fixed (x_fix) and non-fixed (x_new)
  x_fix <- x[fix.loc, ]
  x_new <- x[nofix.loc, ]

  # clear the two item meta data sets
  x_fix <- back2df(metalist2(x_fix))
  x_new <- back2df(metalist2(x_new))

  # record the score categories and model information of the new items to be estimated
  id <- x_new$id
  cats <- x_new$cats
  model <- x_new$model

  # generate the empty meta data with starting values
  if(!use.startval) {
    x_new <- startval_df(cats=cats, model=model)
    x_new$id <- id
  }

  # create the total item meta data to be used in the further estimation process
  x_all <- startval_df(cats=x$cats, model=x$model)
  x_all[fix.loc, 1:ncol(x_fix)] <- x_fix
  x_all[nofix.loc, 1:ncol(x_new)] <- x_new

  # check whether included data are correct
  if(nrow(x_all) != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to matrix
  data <- data.matrix(data)

  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if(length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call.=TRUE)
  }

  # save the item response data into the different object
  data_all <- data
  rm(data, envir=environment(), inherits = FALSE)

  # divide the item response data into two groups: fixed (x_fix) and non-fixed (x)
  data_fix <- data_all[, fix.loc]
  data_new <- data_all[, nofix.loc]

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
  if("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")

  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x_new)
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(x=x_new, loc_1p_const=loc_1p_const, loc_else=loc_else,
                      fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g)

  ##---------------------------------------------------------------------
  # conduct item parameter estimation using MMLE-EM algorithm
  ##---------------------------------------------------------------------
  # check the total number of examinees
  nstd <- nrow(data_all)

  # create initial weights of prior ability distribution when it is not specified
  if(is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out=Quadrature[1])

    # create the data.frame containing the quadurature points and weights
    weights <- gen.weight(dist="norm", mu=group.mean, sigma=sqrt(group.var), theta=quadpt)

  } else {
    quadpt <- weights[, 1]
  }

  # factorize the response values and create a list of item responses across items
  resp_new <- purrr::map2(.x=data.frame(data_new), .y=cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))
  resp_fix <- purrr::map2(.x=data.frame(data_fix), .y=x_fix$cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))
  resp_all <- purrr::map2(.x=data.frame(data_all), .y=x_all$cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))

  # calculate the score categories for each examinee
  std.id <- 1:nstd
  freq_new.cat <- purrr::map(.x=resp_new, .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))
  freq_fix.cat <- purrr::map(.x=resp_fix, .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))
  freq_all.cat <- purrr::map(.x=resp_all, .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))

  # transform the score category data.frame to a matrix
  freq_new.cat <- purrr::map(.x=freq_new.cat,
                             .f=function(k) data.matrix(k) %>%
                               unname())
  freq_fix.cat <- purrr::map(.x=freq_fix.cat,
                             .f=function(k) data.matrix(k) %>%
                               unname())
  freq_all.cat <- purrr::map(.x=freq_all.cat,
                             .f=function(k) data.matrix(k) %>%
                               unname())

  # delete 'resp' object
  rm(resp_new, envir=environment(), inherits = FALSE)
  rm(resp_fix, envir=environment(), inherits = FALSE)
  rm(resp_all, envir=environment(), inherits = FALSE)

  # divide the item response data set into the dichotomous (correct & incorrect) and polytomous parts
  datlist_new <- divide_data(data=data_new, cats=cats, freq.cat=freq_new.cat)
  datlist_fix <- divide_data(data=data_fix, cats=x_fix$cats, freq.cat=freq_fix.cat)
  datlist_all <- divide_data(data=data_all, cats=x_all$cats, freq.cat=freq_all.cat)

  # estimation
  cat("Estimating item parameters...", '\n')

  # set the number of EM iteration to one when OEM (Wainer & Mislevy, 1990) method is used
  if(fipc.method=="OEM") MaxE <- 1


  # implement EM algorithm
  time1 <- Sys.time()
  par.history <- list()
  for(r in 1:MaxE) {

    # implement E-step
    if(r == 1L) {
      estep <- Estep_fipc(x1=x_new, x2=x_fix,
                          data1_drm1=datlist_new$data1_drm, data2_drm1=datlist_new$data2_drm, data_plm1=datlist_new$data_plm,
                          data1_drm2=datlist_fix$data1_drm, data2_drm2=datlist_fix$data2_drm, data_plm2=datlist_fix$data_plm,
                          freq.cat=freq_new.cat, weights=weights, D=D)
    } else {
      estep <- Estep_fipc(x1=x_new, x2=x_all,
                          data1_drm1=datlist_new$data1_drm, data2_drm1=datlist_new$data2_drm, data_plm1=datlist_new$data_plm,
                          data1_drm2=datlist_all$data1_drm, data2_drm2=datlist_all$data2_drm, data_plm2=datlist_all$data_plm,
                          freq.cat=freq_new.cat, weights=weights, D=D)
    }

    # implement M-step
    mstep <- Mstep(estep=estep, id=id, cats=cats, model=model, quadpt=quadpt, D=D, loc_1p_const=loc_1p_const, loc_else=loc_else,
                   EmpHist=EmpHist, weights=weights, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g, a.val.1pl=a.val.1pl,
                   a.val.gpcm=a.val.gpcm, g.val=g.val, use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                   aprior=aprior, bprior=bprior, gprior=gprior, group.mean=group.mean, group.var=group.var, nstd=nstd,
                   Quadrature=Quadrature, use.startval=TRUE, control=control, iter=r, fipc=TRUE,
                   reloc.par=param_loc$reloc.par, info.mstep=FALSE)

    # compute the difference between previous and updated item parameter estimates
    diff_par <- mstep$par_df[, -c(1, 2, 3)] - x_new[, -c(1, 2, 3)]
    max.diff <- abs(max(diff_par, na.rm=TRUE))

    # loglikehood value
    llike <- mstep$loglike

    # print
    cat("\r", paste0("EM iteration: ", r, ", Loglike: ", format(round(mstep$loglike, 4), nsmall=4), ", Max-Change: ", format(round(max.diff, 5), nsmall=4)))

    # check the convergence of EM algorithm
    converge <- max.diff <= Etol

    # end the EM step if the convergence criterion is satisfied
    if(converge | r == MaxE) {
      # extract the final item parameter estimates
      par_df <- mstep$par_df
      x_all[nofix.loc, 1:ncol(par_df)] <- par_df
      par.history[[r]] <- mstep$par_vec

      # extract the final quadrature points and the corresponding weights of the prior population density
      weights <- mstep$weights
      break

    } else {
      # extract the interim item parameter estimates
      x_new <- mstep$par_df
      x_all[nofix.loc, 1:ncol(x_new)] <- x_new
      par.history[[r]] <- mstep$par_vec

      # extract the interim quadrature points and the corresponding weights of the prior population density
      weights <- mstep$weights
    }

  }
  cat("", "\n")
  time2 <- Sys.time()

  # record the item parameter estimation time
  est_time1 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # the first order test: check convergence-criteria test
  test_1st <- all(c(all(mstep$convergence == 0L), r < MaxE))
  if(test_1st) {
    memo3 <- "Convergence criteria are satisfied."
  } else {
    memo3 <- "Convergence criteria are not satisfied."
    warning(paste0(memo3,  " \n"), call.=FALSE)
  }

  # conduct one more E-step to update the posterior distribution using the final item parameter estimates
  estep <- Estep_fipc(x1=par_df, x2=x_all,
                      data1_drm1=datlist_new$data1_drm, data2_drm1=datlist_new$data2_drm, data_plm1=datlist_new$data_plm,
                      data1_drm2=datlist_all$data1_drm, data2_drm2=datlist_all$data2_drm, data_plm2=datlist_all$data_plm,
                      freq.cat=freq_new.cat, weights=weights, D=D)

  # compute the final log of marginal likelihood
  llike <- sum(log(estep$likehd %*% matrix(weights[, 2])))

  # compute the mean and variance of the estimated density distribution
  pop_moments <- cal_moment(node=quadpt, weight=weights[, 2])

  ##---------------------------------------------------------------------
  # estimates the standard errors of item parameter estimates
  ##---------------------------------------------------------------------
  # listrize the item parameter estimate data.frame
  meta <- metalist2(par_df)

  # extract the finalized posterior density
  post_dist <- estep$post_dist

  # compute the infomation matrix of item parameter estimates using the cross-product method
  cat("Computing item parameter var-covariance matrix...", '\n')
  time1 <- Sys.time()
  info.data <- info_xpd(meta=meta, freq.cat=freq_new.cat, post_dist=post_dist, cats=cats, model=model, quadpt=quadpt,
                        D=D, loc_1p_const=loc_1p_const, loc_else=loc_else, nstd=nstd,
                        fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g, a.val.1pl=a.val.1pl,
                        a.val.gpcm=a.val.gpcm, g.val=g.val, reloc.par=param_loc$reloc.par)

  # compute the infomation matrix of item parameter priors
  info.prior <- info_prior(meta=meta, cats=cats, model=model, D=D, loc_1p_const=loc_1p_const,
                           loc_else=loc_else, nstd=nstd, fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g,
                           a.val.1pl=a.val.1pl, a.val.gpcm=a.val.gpcm, g.val=g.val, aprior=aprior, bprior=bprior,
                           gprior=gprior, use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                           reloc.par=param_loc$reloc.par)

  # sum of two information matrixs
  info.mat <- info.data + info.prior

  # the second-order test: check if the information matrix is positive definite
  test_2nd <- all(eigen(info.mat, only.values=TRUE)$values > 1e-8)
  if(test_2nd) {
    if(test_1st) {
      memo4 <- "Solution is a possible local maximum."
    } else {
      memo4 <- "Information matrix of item parameter estimates is positive definite."
    }
  } else {
    memo4 <- "Information matrix of item parameter estimates is not positive definite; unstable solution."
    warning(paste0(memo4, " \n"), call.=FALSE)
  }

  # compute the variance-covariance matrix of the item parameter estimates, and
  # check if the hessian matrix can be inversed
  cov_mat <- suppressWarnings(tryCatch({solve(info.mat, tol=1e-200)}, error = function(e) {NULL}))

  # compute the standard errors of item parameter estimates
  if(is.null(cov_mat)) {
    se_par <- rep(99999, length(diag(info.mat)))
    memo5 <- "Variance-covariance matrix of item parameter estimates is not obtainable; unstable solution."
    warning(paste0(memo5, " \n"), call.=FALSE)
  } else {
    se_par <- suppressWarnings(sqrt(diag(cov_mat)))
    memo5 <- "Variance-covariance matrix of item parameter estimates is obtainable."
  }

  # prevent showing NaN values of standard errors
  if(any(is.nan(se_par))) {
    se_par[is.nan(se_par)] <- 99999
  }

  # set an upper bound of standard error
  se_par <- ifelse(se_par > 99999, 99999, se_par)
  time2 <- Sys.time()

  # record the standard error computationtime
  est_time2 <- round(as.numeric(difftime(time2, time1, units = "secs")), 2)

  # deploy the standard errors on the location of matrix as the item parameter estimates
  # 1) for the only new items
  se_df <- loc.par <- param_loc$loc.par
  for(i in 1:nrow(loc.par)) {
    num.loc <- which(!is.na(loc.par[i, ]))
    se.loc <- loc.par[i, ][num.loc]
    se_df[i, num.loc] <- se_par[se.loc]
  }

  # 2) for the a total test form
  se_all_df <- loc_all.par <- matrix(NA, nrow=nitem, ncol=max(x_all$cats))
  if(ncol(se_df) < ncol(se_all_df)) {
    n2add <- ncol(se_all_df) - ncol(se_df)
    se_df <- cbind(se_df, matrix(NA, nrow=nrow(se_df), ncol=n2add))
  }
  se_all_df[nofix.loc, 1:ncol(se_df)] <- se_df
  loc_all.par[nofix.loc, 1:ncol(loc.par)] <- loc.par

  # create a full data.frame for the standard error estimates
  se_all_df <- data.frame(x_all[, 1:3], se_all_df)
  colnames(se_all_df) <- c("id", "cats", "model", paste0("se.", 1:ncol(loc_all.par)))

  # create a full data.frame containing the position of item parameter estimates
  # this data.frame is useful when interpreting the variance-covariance matrix of item parameter estimates
  loc_all_df <- data.frame(x_all[, 1:3], loc_all.par)
  colnames(loc_all_df) <- c("id", "cats", "model", paste0("par.", 1:ncol(loc_all.par)))

  ##---------------------------------------------------------------------
  # summarize the estimation results
  ##---------------------------------------------------------------------
  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow=nrow(loc_all.par), ncol=2*ncol(loc_all.par)))
  all_df[, seq(1, 2*ncol(loc_all.par), 2)] <- x_all[, -c(1:3)]
  all_df[, seq(2, 2*ncol(loc_all.par), 2)] <- se_all_df[, -c(1:3)]
  col.names <- rep(NA, 2*ncol(loc_all.par))
  col.names[seq(1, 2*ncol(loc_all.par), 2)] <- paste0("par.", 1:ncol(loc_all.par))
  col.names[seq(2, 2*ncol(loc_all.par), 2)] <- paste0("se.", 1:ncol(loc_all.par))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x_all[, 1:3], all_df)

  # population density parameters
  mu <- pop_moments[1]
  sigma2 <- pop_moments[2]
  sigma <- sqrt(pop_moments[2])
  moments.est <- c(mu, sigma2, sigma)

  # compute the standard errors of population density parameters
  se.mu <- sigma / sqrt(nstd)
  se.sigma2 <- sigma2 * sqrt(2 / (nstd - 1))
  se.sigma <- (1 / (2 * sigma)) * se.sigma2 # using Delta method
  moments.se <- c(se.mu, se.sigma2, se.sigma)

  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(moments.est, moments.se))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  rownames(group.par) <- c("estimates", "se")

  # prior information
  if(use.aprior) aprior.dist <- aprior else aprior.dist <- NULL
  if(use.bprior) bprior.dist <- bprior else bprior.dist <- NULL
  if(use.gprior) gprior.dist <- gprior else gprior.dist <- NULL

  ##---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time3 <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- list(estimates=full_all_df, par.est=x_all, se.est=se_all_df, pos.par=loc_all_df, covariance=cov_mat, loglikelihood=llike, group.par=group.par,
              weights=weights, posterior.dist=post_dist, data=data_all, scale.D=D, ncase=nstd, nitem=nitem, Etol=Etol, MaxE=MaxE,
              aprior=aprior.dist, bprior=bprior.dist, gprior=gprior.dist, npar.est=length(param_loc$reloc.par), niter=r, maxpar.diff=max.diff,
              EMtime=est_time1, SEtime=est_time2, TotalTime=est_time3, test.1=memo3, test.2=memo4, var.note=memo5, fipc=TRUE,
              fipc.method=fipc.method, fix.loc=fix.loc)

  cat("Estimation is finished in", est_time3, "seconds.",'\n')
  return(rst)

}

