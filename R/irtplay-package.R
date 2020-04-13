#' irtplay: Unidimensional Item Response Theory Modeling
#'
#' @description
#' Fit unidimensional item response theory (IRT) models to mixture of dichotomous and polytomous data,
#' calibrate online item parameters (i.e., pretest and operational items), estimate examinees abilities, and examine the IRT model-data
#' fit on item-level in different ways as well as provide useful functions related to unidimensional IRT.
#'
#' For the item parameter estimation, the marginal maximum likelihood estimation with expectation-maximization (MMLE-EM) algorithm
#' (Bock & Aitkin, 1981) is used. For the online calibration, Stocking's Method A (Ban, Hanson, Wang, Yi, & Harris, 2001; stocking, 1988)
#' and the fixed item parameter calibration (FIPC) method (Kim, 2006) are provided. For the ability estimation, several popular
#' scoring methods (e.g., MLE, EAP, and MAP) are implemented. In terms of assessing the IRT model-data fit, one of distinguished features
#' of this package is that it gives not only item fit statistics (e.g., \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
#' likelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley & Mills, 1985), infit and outfit statistics (Ames et al., 2015),
#' and \eqn{S-X^{2}} (Orlando & Thissen, 2000, 2003)) but also graphical displays to look at residuals between the observed data and
#' model-based predictions (Hambleton, Swaminathan, & Rogers, 1991).
#'
#' In addition, there are many useful functions such as computing asymptotic variance-covariance matrices of item parameter estimates,
#' importing item and/or ability parameters from popular IRT software, running flexMIRT (Cai, 2017) through R, generating simulated data,
#' computing the conditional distribution of observed scores using the Lord-Wingersky recursion formula, computing the loglikelihood of
#' individual items, computing the loglikelihood of abilities, computing item and test information functions, computing item and test
#' characteristic curve functions, and plotting item and test characteristic curves and item and test information functions.
#'
#' \tabular{ll}{ Package: \tab irtplay\cr Version: \tab 1.5.0\cr Date: \tab
#' 2020-04-12\cr Depends: \tab R (>= 3.6)\cr License: \tab GPL (>= 2)\cr }
#'
#' @details
#' Following five sections describe a) how to implement the online item calibration using FIPC, a) how to implement the online item
#' calibration using Method A, b) the process of evaluating the IRT model-data fit, c) two examples for the online calibration and
#' evaluating the IRT model-data fit, and d) IRT Models used in \pkg{irtplay} package.
#'
#'
#' @section Online item calibration with the fixed item parameter calibration (FIPC) method (e.g., Kim, 2006):
#'
#' The fixed item parameter calibration (FIPC) is one of useful online item calibration methods for computerized adaptive testing (CAT)
#' to put the parameter estimates of pretest items on the same scale of operational item parameter estimates without post hoc
#' linking/scaling (Ban, Hanson, Wang, Yi, & Harris, 2001; Chen & Wang, 2016). In FIPC, the operational item parameters are fixed to
#' estimate the characteristic of the underlying latent variable prior distribution when calibrating the pretest items. More specifically,
#' the underlying latent variable prior distribution of the operational items is estimated during the calibration of the pretest
#' items to put the item parameters of the pretest items on the scale of the operational item parameters (Kim, 2006). In \pkg{irtplay}
#' package, FIPC is implemented with two main steps:
#'
#' \enumerate{
#'   \item Prepare a response data set and the item meta data of the fixed (or operational) items.
#'   \item Implement FIPC to estimate the item parameters of pretest items using the \code{\link{est_irt}} function.
#' }
#'
#' \describe{
#'   \item{1. Preparing a data set}{
#'   To run the \code{\link{est_irt}} function, it requires two data sets:
#'
#'     \enumerate{
#'       \item Item meta data set (i.e., model, score category, and item parameters. see the desciption of the argument \code{x} in the function \code{\link{est_irt}}).
#'       \item Examinees' response data set for the items. It should be a matrix format where a row and column indicate the examinees and the items, respectively.
#'       The order of the columns in the response data set must be exactly the same as the order of rows of the item meta data.
#'     }
#'   }
#'
#'   \item{2. Estimating the pretest item parameters}{
#'   When FIPC is implemented in \code{\link{est_irt}} function, the pretest item parameters are estimated by fixing the operational item parameters. To estimate the item
#'   parameters, you need to provide the item meta data in the argument \code{x} and the response data in the argument \code{data}.
#'
#'   It is worthwhile to explain about how to prepare the item meta data set in the argument \code{x}. A specific form of a data.frame should be used for
#'   the argument \code{x}. The first column should have item IDs, the second column should contain the number of score categories of the items, and the third
#'   column should include IRT models. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous
#'   items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded response model and
#'   (generalized) partial credit model, respectively. From the fourth column, item parameters should be included. For dichotomous items, the fourth, fifth,
#'   and sixth columns represent the item discrimination (or slope), item difficulty, and item guessing parameters, respectively. When "1PLM" or "2PLM" is
#'   specified for any items in the third column, NAs should be inserted for the item guessing parameters. For polytomous items, the item discrimination (or slope)
#'   parameters should be contained in the fourth column and the item threshold (or step) parameters should be included from the fifth to the last columns.
#'   When the number of categories differs between items, the empty cells of item parameters should be filled with NAs. See `est_irt` for more details about
#'   the item meta data.
#'
#'   Also, you should specify in the argument \code{fipc = TRUE} and a specific FIPC method in the argument \code{fipc.method}. Finally, you should provide
#'   a vector of the location of the items to be fixed in the argument \code{fix.loc}. For more details about implementing FIPC, see the
#'   description of the function \code{\link{est_irt}}.
#'
#'   When implementing FIPC, you can estimate both the emprical histogram and the scale of latent variable prior distribution by setting \code{EmpHist = TRUE}.
#'   If \code{EmpHist = FALSE}, the normal prior distribution is used during the item parameter estimation and the scale of the normal prior distribution is
#'   updated during the EM cycle.
#'
#'   The \code{\link{est_item}} function requires a vector of the number of score categories for the items in the argument \code{cats}. For example, a dichotomous item has
#'   two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items are binary items
#'   (i.e., dichotomous items), it assumes that all items have two score categories.
#'
#'   If necessary, you need to specify whether prior distributions of item slope and guessing parameters (only for the IRT 3PL model) are used in the arguments of
#'   \code{use.aprior} and \code{use.gprior}, respectively. If you decide to use the prior distributions, you should specify what distributions will be used for the prior
#'   distributions in the arguments of \code{aprior} and \code{gprior}, respectively. Currently three probability distributions of Beta, Log-normal, and Normal
#'   distributions are available.
#'
#'   In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.
#'
#'   Once the \code{\link{est_irt}} function has been implemented, you'll get a list of several internal objects such as the item parameter estimates,
#'   standard error of the parameter estimates.
#'   }
#' }
#'
#'
#' @section Online item calibration with Method A (Stocking, 1988):
#' In CAT, Method A is the relatively simplest and most straightforward online calibration method,
#' which is the maximum likelihood estimation of the item parameters given the proficiency estimates. In CAT, Method A can be used
#' to put the parameter estimates of pretest items on the same scale of operational item parameter estimates and recalibrate
#' the operational items to evaluate the parameter drifts of the operational items (Chen & Wang, 2016; Stocking, 1988).
#' Also, Method A is known to result in accurate, unbiased item parameters calibration when items are randomly rather than
#' adaptively administered to examinees, which occurs most commonly with pretest items (Ban, Hanson, Wang, Yi, & Harris, 2001; Chen & Wang, 2016).
#' Using \pkg{irtplay} package, Method A is implemented to calibrate the items with two main steps:
#'
#' \enumerate{
#'   \item Prepare a data set for the calibration of item parameters (i.e., item response data and ability estimates).
#'   \item Implement Method A to estimate the item parameters using the \code{\link{est_item}} function.
#' }
#'
#' \describe{
#'   \item{1. Preparing a data set}{
#'   To run the \code{\link{est_item}} function, it requires two data sets:
#'
#'     \enumerate{
#'       \item Examinees' ability (or proficiency) estimates. It should be in the format of a numeric vector.
#'       \item response data set for the items. It should be in the format of matrix where a row and column indicate
#'     the examinees and the items, respectively. The order of the examinees in the response data set must be exactly the same as that of the examinees' ability estimates.
#'     }
#'   }
#'
#'   \item{2. Estimating the pretest item parameters}{
#'   The \code{\link{est_item}} function estimates the pretest item parameters given the proficiency estimates. To estimate the item parameters,
#'   you need to provide the response data in the argument \code{data} and the ability estimates in the argument \code{score}.
#'
#'   Also, you should provide a string vector of the IRT models in the argument \code{model} to indicate what IRT model is used to calibrate each item.
#'   Available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent
#'   the graded response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function. If a single
#'   character of the IRT model is specified, that model will be recycled across all items.
#'
#'   The \code{\link{est_item}} function requires a vector of the number of score categories for the items in the argument \code{cats}. For example, a dichotomous item has
#'   two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items are binary items
#'   (i.e., dichotomous items), it assumes that all items have two score categories.
#'
#'   If necessary, you need to specify whether prior distributions of item slope and guessing parameters (only for the IRT 3PL model) are used in the arguments of
#'   \code{use.aprior} and \code{use.gprior}, respectively. If you decide to use the prior distributions, you should specify what distributions will be used for the prior
#'   distributions in the arguments of \code{aprior} and \code{gprior}, respectively. Currently three probability distributions of Beta, Log-normal, and Normal
#'   distributions are available.
#'
#'   In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.
#'
#'   Once the \code{\link{est_item}} function has been implemented, you'll get a list of several internal objects such as the item parameter estimates,
#'   standard error of the parameter estimates.
#'   }
#' }
#'
#'
#' @section The process of evaluating the IRT model-data fit:
#' One way to assess goodness of IRT model-data fit is through an item fit analysis by examining the traditional item fit statistics
#' and looking at the discrepancy between the observed data and model-based predictions. Using \pkg{irtplay} package, the traditional approach
#' of evaluating the IRT model-data fit on item-level can be implemented with three main steps:
#'
#' \enumerate{
#'   \item Prepare a data set for the IRT item fit analysis (i.e., item meta data, ability estimates, and response data).
#'   \item Obtain the IRT fit statistics such as \eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit statistics using the function \code{\link{irtfit}}.
#'   \item Based on the results of IRT model fit analysis (i.e., an object of class \code{\link{irtfit}}) obtained in step 2,
#' draw the IRT residual plots (i.e., raw residual and standardized residual plots) using the function \code{\link{plot.irtfit}}.
#' }
#'
#' \describe{
#'   \item{1. Preparing a data set}{
#'   Before conducting the IRT model fit analysis, it is necessary to prepare a data set. To run the function \code{\link{irtfit}}, it requires
#'   three data sets:
#'
#'    \enumerate{
#'    \item Item meta data including the item ID, number of score categories, IRT models, and item parameters. The item meta data should be in the format of
#'    data.frame. You can prepare the data either by using the function \code{\link{shape_df}} or by creating a data.frame of the item meta data by yourself.
#'    If you have output files of item parameter estimates obtained from one of the IRT software such as BILOG-MG 3, PARSCALE 4, flexMIRT, and mirt (R package),
#'    the item meta data can be easily obtained using the functions of \code{\link{bring.bilog}}, \code{\link{bring.parscale}}, \code{\link{bring.flexmirt}},
#'    and \code{\link{bring.mirt}}. See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data format.
#'    \item Examinees' ability (or proficiency) estimates. It should be in the format of a numeric vector.
#'    \item Examinees' response data set for the items. It should be in the format of matrix where a row and column indicate the examinees and the items,
#'    respectively. The order of the examinees in the response data set must be exactly the same as that of the examinees' ability estimates. The order of the items
#'    in the response data set must be exactly the same as that of the items in the item meta data.
#'    }
#'
#'  }
#'
#'   \item{2. Computing the IRT model-data fit statistics}{
#'   The function \code{\link{irtfit}} computes the traditional IRT item fit statistics such as \eqn{\chi^{2}}, \eqn{G^{2}}, infit, and outfit statistics.
#'   To calculate the \eqn{\chi^{2}} and \eqn{G^{2}} statistics, two methods are available to divide the ability scale into several groups. The two methods are "equal.width"
#'   for dividing the scale by an equal length of the interval and "equal.freq" for dividing the scale by an equal frequency of examinees. Also, you need to
#'   specify the location of ability point at each group (or interval) where the expected probabilities of score categories are calculated from the IRT models.
#'   Available locations are "average" for computing the expected probability at the average point of examinees' ability estimates in each group and "middle" for
#'   computing the expected probability at the midpoint of each group.
#'
#'   To use the function \code{\link{irtfit}}, you need to insert the item meta data in the argument \code{x}, the ability estimates in the argument \code{score},
#'   and the response data in the argument \code{data}. If you want to divide the ability scale into other than ten groups, you need to specify the number of groups
#'   in the argument \code{n.width}. In addition, if the response data include missing values, you must indicate the missing value in argument \code{missing}.
#'
#'   Once the function \code{\link{irtfit}} has been implemented, you'll get the fit statistic results and the contingency tables for every item used
#'   to calculate the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
#'  }
#'
#'   \item{3. Drawing the IRT residual plots}{
#'   Using the saved object of class \code{\link{irtfit}}, you can use the \code{\link{plot}} method to evaluate the IRT raw residual and standardized residual plots.
#'
#'   Because the \code{\link{plot}} method can draw the residual plots for an item at a time, you have to indicate which item will be examined. For this,
#'   you can specify an integer value, which is the location of the studied item, in the argument \code{item.loc}.
#'
#'   In terms of the raw residual plot, the argument \code{ci.method} is used to select a method to estimate the confidence intervals among four methods.
#'   Those methods are "wald" for the Wald interval, which is based on the normal approximation (Laplace, 1812), "cp" for Clopper-Pearson interval
#'   (Clopper & Pearson, 1934), "wilson" for Wilson score interval (Wilson, 1927), and "wilson.cr" for Wilson score interval with continuity correction
#'   (Newcombe, 1998).
#'  }
#' }
#'
#'
#' @section Three examples of R script:
#'
#' The example code below shows how to implement the online calibration and how to evaluate the IRT model-data fit:\preformatted{
#' ##---------------------------------------------------------------
#' # Attach the packages
#' library(irtplay)
#'
#' ##----------------------------------------------------------------------------
#' # 1. The example code below shows how to prepare the data sets and how to
#' #    implement the fixed item parameter calibration (FIPC):
#' ##----------------------------------------------------------------------------
#'
#' ## Step 1: prepare a data set
#' ## In this example, we generated examinees' true proficiency parameters and simulated
#' ## the item response data using the function "simdat".
#'
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # generate 1,000 examinees' latent abilities from N(0.4, 1.3)
#' set.seed(20)
#' score <- rnorm(1000, mean=0.4, sd=1.3)
#'
#' # simulate the response data
#' sim.dat <- simdat(x=x, theta=score, D=1)
#'
#' ## Step 2: Estimate the item parameters
#' # fit the 3PL model to all dichotmous items, fit the GRM model to all polytomous data,
#' # fix the five 3PL items (1st - 5th items) and three GRM items (53th to 55th items)
#' # also, estimate the empirical histogram of latent variable
#' fix.loc <- c(1:5, 53:55)
#' mod.fix1 <- est_irt(x=x, data=sim.dat, D=1, use.gprior=TRUE,
#'                     gprior=list(dist="beta", params=c(5, 16)), EmpHist=TRUE, Etol=1e-3,
#'                     fipc=TRUE, fipc.method="MEM", fix.loc=fix.loc)
#' print(mod.fix1)
#'
#' # plot the estimated empirical histogram of latent variable prior distribution
#' emphist <- mod.fix1$weights
#' plot(emphist$weight ~ emphist$theta, xlab="Theta", ylab="Density")
#'
#'
#' ##----------------------------------------------------------------------------
#' # 2. The example code below shows how to prepare the data sets and how to estimate
#' #    the item parameters using Method A:
#' ##----------------------------------------------------------------------------
#'
#' ## Step 1: prepare a data set
#' ## In this example, we generated examinees' true proficiency parameters and simulated
#' ## the item response data using the function "simdat". Because, the true
#' ## proficiency parameters are not known in reality, however, the true proficiencies
#' ## would be replaced with the proficiency estimates for the calibration.
#'
#' # import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # modify the item meta data so that some items follow 1PLM, 2PLM and GPCM
#' x[c(1:3, 5), 3] <- "1PLM"
#' x[c(1:3, 5), 4] <- 1
#' x[c(1:3, 5), 6] <- 0
#' x[c(4, 8:12), 3] <- "2PLM"
#' x[c(4, 8:12), 6] <- 0
#' x[54:55, 3] <- "GPCM"
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' ## Step 2: Estimate the item parameters
#' # 1) item parameter estimation: constrain the slope parameters of the 1PLM to be equal
#' mod1 <- est_item(x, data, score, D=1, fix.a.1pl=FALSE, use.gprior=TRUE,
#'                  gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE)
#' print(mod1)
#'
#' # 2) item parameter estimation: fix the slope parameters of the 1PLM to 1
#' mod2 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, a.val.1pl=1, use.gprior=TRUE,
#'                  gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE)
#' print(mod2)
#'
#' # 3) item parameter estimation: fix the guessing parameters of the 3PLM to 0.2
#' mod3 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, fix.g=TRUE, a.val.1pl=1, g.val=.2,
#'                  use.startval=FALSE)
#' print(mod3)
#'
#'
#' ##----------------------------------------------------------------------------
#' # 3. The example code below shows how to prepare the data sets and how to conduct
#' #    the IRT model-data fit analysis:
#' ##----------------------------------------------------------------------------
#'
#' ## Step 1: prepare a data set for IRT
#' ## In this example, we use the simulated mixed-item format of CAT Data
#' ## But, only items that have examinees' responses more than 1,000 are assessed.
#'
#' # find the location of items that have more than 1,000 item responses
#' over1000 <- which(colSums(simCAT_MX$res.dat, na.rm=TRUE) > 1000)
#'
#' # (1) item meta data
#' x <- simCAT_MX$item.prm[over1000, ]
#'
#' # (2) examinee's ability estimates
#' score <- simCAT_MX$score
#'
#' # (3) response data
#' data <- simCAT_MX$res.dat[, over1000]
#'
#' ## Step 2: Compute the IRT mode-data fit statistics
#' # (1) the use of "equal.width"
#' fit1 <- irtfit(x=x, score=score, data=data, group.method="equal.width",
#'                n.width=10, loc.theta="average", range.score=NULL, D=1,
#'                alpha=0.05, missing=NA)
#'
#' # what kinds of internal objects does the results have?
#' names(fit1)
#'
#' # show the results of the fit statistics
#' fit1$fit_stat[1:10, ]
#'
#' # show the contingency tables for the first item (dichotomous item)
#' fit1$contingency.fitstat[[1]]
#'
#' # (2) the use of "equal.freq"
#' fit2 <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
#'                n.width=10, loc.theta="average", range.score=NULL, D=1,
#'                alpha=0.05, missing=NA)
#'
#' # show the results of the fit statistics
#' fit2$fit_stat[1:10, ]
#'
#' # show the contingency table for the fourth item (polytomous item)
#' fit2$contingency.fitstat[[4]]
#'
#' ## Step 3: Draw the IRT residual plots
#' # 1. for the dichotomous item
#' # (1) both raw and standardized residual plots using the object "fit1"
#' plot(x=fit1, item.loc=1, type = "both", ci.method = "wald",
#'      ylim.sr.adjust=TRUE)
#'
#' # (2) the raw residual plots using the object "fit1"
#' plot(x=fit1, item.loc=1, type = "icc", ci.method = "wald",
#'      ylim.sr.adjust=TRUE)
#'
#' # (3) the standardized residual plots using the object "fit1"
#' plot(x=fit1, item.loc=113, type = "sr", ci.method = "wald",
#'      ylim.sr.adjust=TRUE)
#'
#' # 2. for the polytomous item
#' # (1) both raw and standardized residual plots using the object "fit1"
#' plot(x=fit1, item.loc=113, type = "both", ci.method = "wald",
#'      ylim.sr.adjust=TRUE)
#'
#' # (2) the raw residual plots using the object "fit1"
#' plot(x=fit1, item.loc=113, type = "icc", ci.method = "wald",
#'      layout.col=2, ylim.sr.adjust=TRUE)
#'
#' # (3) the standardized residual plots using the object "fit1"
#' plot(x=fit1, item.loc=113, type = "sr", ci.method = "wald",
#'      layout.col=4, ylim.sr.adjust=TRUE)
#' }
#'
#'
#' @section IRT Models:
#'
#' In \pkg{irtplay} package, both dichotomous and polytomous IRT models are available.
#' For dichotomous items, IRT one-, two-, and three-parameter logistic models are used. For polytomous items,
#' the graded response model (GRM) and the (generalized) partial credit model (GPCM) are used.
#' Note that the item discrimination (or slope) parameters should be fixed to 1 when the partial credit model is used.
#' The IRT model parameter labels use the conventions shown below. Note that \eqn{M} is the total number of categories
#' for each of polytomous item.
#'
#' \describe{
#'   \item{IRT 1-3PL models}{
#'     For the IRT 1-3PL models, the probability that an examinee with \eqn{\theta} provides a correct answer for an item is given by,
#'      \deqn{P(u = 1|\theta) = g + \frac{(1 - g)}{1 + exp(-Da(\theta - b))},}
#'     where \eqn{D} is a scaling factor. When the IRT 1PL model is used, \eqn{a} is fixed to a constant and \eqn{g = 0}
#'     for all items. When the IRT 2PL model is used, \eqn{g = 0}.
#'   }
#'   \item{GRM}{
#'     For the graded response model, the probability that an examinee with \eqn{\theta} scores \eqn{k} on an item with \eqn{M} categories
#'     is a given by,
#'     \deqn{P(u = k | \theta) = P(x \ge k | \theta) - P(x \ge k + 1 | \theta),}
#'     \deqn{P(x \ge k | \theta) = \frac{1}{1 + exp(-Da(\theta - b_{k}))}, }
#'     \deqn{P(x \ge k + 1 | \theta) = \frac{1}{1 + exp(-Da(\theta - b_{k+1}))}, }
#'     where \eqn{b_{k}} is the \eqn{k}th threshold parameter for an item.
#'   }
#'   \item{GPCM}{
#'     For the (generalized) partial credit model, the probability that an examinee with \eqn{\theta} scores \eqn{k} on an item with \eqn{M}
#'     categories is a given by,
#'      \deqn{P(u = k | \theta) = \frac{exp(\sum_{u=0}^{k}{Da(\theta - b_{u})})}{\sum_{j=0}^{M-1}{exp(\sum_{u=0}^{j}{Da(\theta - b_{u})})}},}
#'     where \eqn{b_{u}} is the \eqn{u}th step parameter for an item. In other contexts, the step parameter \eqn{b_{u}} can
#'     also be parameterized as \eqn{b_{u} = \beta - \tau_{u}}, where \eqn{\beta} denotes the location (or overall difficulty) parameter
#'     and \eqn{\tau_{u}} represents a difficulty (or threshold) parameter for each category. In this package, only the step parameter of \eqn{b_{u}}
#'     is only used. Also, note that M-1 step parameters are necessary when an item has M categories because a step parameter for the first category
#'     does not affect the category probabilities. When a partial credit model is used, \eqn{a = 1}.
#'    }
#'
#'}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional Module on Item-Fit Statistics for Item Response Theory Models.
#' \emph{Educational Measurement: Issues and Practice, 34}(3), 39-48.
#'
#' Baker, F. B., & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques.} CRC Press.
#'
#' Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
#' in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
#' \emph{Statistical theories of mental test scores} (pp. 397-479). Reading, MA: Addison-Wesley.
#'
#' Bock, R.D. (1960), \emph{Methods and applications of optimal scaling}. Chapel Hill, NC: L.L. Thurstone Psychometric Laboratory.
#'
#' Bock, R. D., & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters: Application of an EM algorithm.
#' \emph{Psychometrika, 46}, 443-459.
#'
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Psychometrika, 35}, 179-198.
#'
#' Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
#' Chapel Hill, NC: Vector Psychometric Group.
#'
#' Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment.
#' \emph{Journal of Statistical Software, 48}(6), 1-29.
#'
#' Chen, P., & Wang, C. (2016). A new online calibration method for multidimensional computerized adaptive testing.
#' \emph{Psychometrika, 81}(3), 674-701.
#'
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial.
#' \emph{Biometrika, 26}(4), 404-413.
#'
#' Hambleton, R. K., & Swaminathan, H., & Rogers, H. J. (1991) \emph{Fundamentals of item response theory}.
#' Newbury Park, CA: Sage.
#'
#' Han, K. T. (2016). Maximum likelihood score estimation method with fences for short-length tests and computerized adaptive tests.
#' \emph{Applied psychological measurement, 40}(4), 289-301.
#'
#' Kang, T., & Chen, T. T. (2008). Performance of the generalized S-X2 item fit index for polytomous IRT models.
#' \emph{Journal of Educational Measurement, 45}(4), 391-406.
#'
#' Kim, S. (2006). A comparative study of IRT fixed parameter calibration methods.
#' \emph{Journal of Educational Measurement, 43}(4), 355-381.
#'
#' Kolen, M. J. & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
#' Springer.
#'
#' Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
#' \emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.
#'
#' Laplace, P. S. (1820).\emph{Theorie analytique des probabilites} (in French). Courcier.
#'
#' Li, Y. & Lissitz, R. (2004). Applications of the analytically derived asymptotic standard errors of item response theory
#' item parameter estimates. \emph{Journal of educational measurement, 41}(2), 85-117.
#'
#' Lord, F. & Wingersky, M. (1984). Comparison of IRT true score and equipercentile observed score equatings.
#' \emph{Applied Psychological Measurement, 8}(4), 453-461.
#'
#' McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
#' \emph{Applied Psychological Measurement, 9}, 49-57.
#'
#' Meilijson, I. (1989). A fast improvement to the EM algorithm on its own terms.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological), 51}, 127-138.
#'
#' Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
#' scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com
#'
#' Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of seven methods.
#' \emph{Statistics in medicine, 17}(8), 857-872.
#'
#' Orlando, M., & Thissen, D. (2000). Likelihood-based item-fit indices for dichotomous item response theory models.
#' \emph{Applied Psychological Measurement, 24}(1), 50-64.
#'
#' Orlando, M., & Thissen, D. (2003). Further investigation of the performance of S-X2: An item fit index for use with
#' dichotomous item response theory models. \emph{Applied Psychological Measurement, 27}(4), 289-298.
#'
#' Pritikin, J. (2018). \emph{rpf: Response Probability Functions}. R package version 0.59.
#' https://CRAN.R-project.org/package=rpf.
#'
#' Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
#' \emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.
#'
#' Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.
#'
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' \emph{Psychometrika, 47}, 175-186.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
#' for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
#' Measurement, 19}(1), 39-49.
#'
#' Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
#' \emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Wainer, H., & Mislevy, R. J. (1990). Item response theory, item calibration, and proficiency estimation. In H. Wainer (Ed.),
#' \emph{Computer adaptive testing: A primer} (Chap. 4, pp.65-102). Hillsdale, NJ: Lawrence Erlbaum.
#'
#' Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests Using IRT-Based Methods.
#' \emph{Journal of Statistical Software, 35}(12), 1-33. URL http://www.jstatsoft.org/v35/i12/.
#'
#' Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric procedure for assessing goodness-of-fit in
#' item response theory. \emph{Applied Measurement in Education, 21}(1), 22-40.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
#' \emph{Journal of the American Statistical Association, 22}(158), 209-212.
#'
#' Woods, C. M. (2007). Empirical histograms in item response theory with ordinal data. \emph{Educational and Psychological Measurement, 67}(1), 73-87.
#'
#' Yen, W. M. (1981). Using simulation results to choose a latent trait model. \emph{Applied Psychological Measurement, 5}, 245-262.
#'
#' Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). BILOG-MG 3: Multiple-group
#' IRT analysis and test maintenance for binary items [Computer Program]. Chicago, IL: Scientific
#' Software International. URL http://www.ssicentral.com
#'
#'
#' @docType package
#' @name irtplay-package
#' @keywords package
NULL
