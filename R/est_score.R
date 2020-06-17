#' Estimate examinees' ability (proficiency) parameters
#'
#' @description This function estimates examinees' latent ability parameters. Available scoring methods are maximum likelihood estimation (MLE),
#' maximum likelihood estimation with fences (MLEF; Han, 2016), maximum a posteriori estimation (MAP; Hambleton et al., 1991),
#' expected a posteriori estimation (EAP; Bock & Mislevy, 1982), EAP summed scoring (Thissen et al., 1995; Thissen & Orlando, 2001),
#' and inverse test characteristic curve (TCC) scoring (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...) or an object of
#' class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}. See \code{\link{irtfit}}, \code{\link{test.info}},
#' or \code{\link{simdat}} for more details about the item meta data. This data.frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param data A matrix or vector containing examinees' response data for the items in the argument \code{x}. When a matrix is used, a row and column indicate
#' the examinees and items, respectively. When a vector is used, it should contains the item response data for an examinee.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param method A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation,
#' "MLEF" for the maximum likelihood estimation with fences, "MAP" for the maximum a posteriori estimation,
#' "EAP" for the expected a posteriori estimation, "EAP.SUM" for the expected a posteriori summed scoring, and "INV.TCC" for the inverse TCC scoring.
#' Default method is "MLE".
#' @param range A numeric vector of two components to restrict the range of ability scale for the MLE. Default is c(-4, 4).
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "MLE", "MLEF", or "INV.TCC".
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' Ignored if \code{method} is "MLE", "MLEF", "MAP", or "INV.TCC".
#' @param weights A two-column matrix or data.frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL and \code{method} is "EAP" or "EAP.SUM", default values are used (see the arguments
#' of \code{norm.prior} and \code{nquad}). Ignored if \code{method} is "MLE", "MLEF", "MAP", or "INV.TCC".
#' @param fence.a A numeric value specifying the item slope parameter (i.e., \emph{a}-parameter) for the two imaginary items in MLEF. See below for details.
#' Default is 3.0.
#' @param fence.b A numeric vector of two components specifying the lower and upper fences of item difficulty parameters (i.e., \emph{b}-parameters)
#' for the two imaginary items, respectively, in MLEF. When \code{fence.b = NULL}, the lower and upper fences of item difficulty parameters were
#' automatically set. See below for details. Default is NULL.
#' @param se A logical value. If TRUE, the standard errors of ability estimates are computed. However, if \code{method} is "EAP.SUM" or "INV.TCC", the standard
#' errors are always returned.
#' @param constant A numeric value used to adjust zero and perfect raw sum scores, or the raw sum score equal to the sum of item guessing parameters,
#' if necessary, to find estimable solutions for those raw sum scores when \code{method = "INV.TCC"}. The zero raw score is forced to become the score of "zero raw score + constant"
#' and the perfect raw score is forced to become the score of "perfect raw score - constant". If the 3PLM items are included in the item meta data,
#' the raw sum score equal to the sum of item guessing parameters is forced to become the score of "the raw sum score + constant". Default is .1.
#' @param constraint A logical value indicating whether the ability estimates will be restricted within a specific ability range
#' specified in the argument \code{range.tcc} when \code{method = "INV.TCC"}. If \code{constraint = TRUE}, all ability estimates less than the first value in the vector specified in
#' the argument \code{range.tcc} are transformed to the first value and all ability estimates greater than the second value in the vector specified in
#' the argument \code{range.tcc} are transformed to the second value. Also, when \code{constraint = TRUE} and the 3PLM items are contained
#' in the item meta data, linear interpolation method is used to find the ability estimates for the raw sum scores less than the sum of item guessing
#' parameters. When \code{constraint = FALSE} and the 3PLM items are contained in the item meta data, linear extrapolation method is used to find
#' the ability estimates for the raw sum scores less than the sum of item guessing parameters. See below for details. Default is FALSE.
#' @param range.tcc A numeric vector of two components to be used as the lower and upper bounds of ability estimates when \code{method = "INV.TCC"} and
#' \code{constraint = TRUE}. Default is c(-7, 7).
#' @param missing A value indicating missing values in the response data set. Default is NA. See below for details.
#' @param ncore The number of logical CPU cores to use. Default is 1. See below for details.
#' @param ... additional arguments to pass to \code{parallel::makeCluster}.
#'
#' @details For MAP scoring method, only the normal prior distribution is available for the population distribution.
#'
#' When there are missing data in the response data set, the missing value must be specified in \code{missing}. The missing data are taken into account
#' when either of MLE, MLEF, MAP, and EAP is used. However, there must be no missing data in the response data set when "EAP.SUM" or "INV.TCC" is used.
#' One of possible ways to use "EAP.SUM" or "INV.TCC" method when missing values exist is to remove rows with any missing values.
#'
#' In the maximum likelihood estimation with fences (MLEF; Han, 2016), two 2PLM imaginary items are necessary. The first imaginary item serves as the lower
#' fence and its difficulty parameter (i.e., \emph{b}-parameters) should be lower than any difficulty parameter values in the test form. Likewise, the second
#' imaginary item serves as the upper fence and its difficulty parameter should be greater than any difficulty parameter values in the test form. Also, the two
#' imaginary items should have a very high item slope parameter (i.e., \emph{a}-parameter) value. See Han (2016) for more details.
#'
#' When \code{fence.b = NULL} in MLEF, the function automatically sets the lower and upper fences of item difficulty parameters using two steps. More specifically,
#' in the first step, the lower fence of the item difficulty parameter is set to the greatest integer value less than the minimum of item difficulty parameters
#' in the item meta data and the upper fence of the item difficulty parameter is set to the smallest integer value greater than the maximum of item difficulty
#' parameters in the item meta data. Then, in the second step, if the lower fence set in the first step is greater than -3.5, the lower fence is constrained to -3.5
#' and if the upper fence set in the first step is less than 3.5, the upper fence is constrained to 3.5. Otherwise, the fence values of item difficulty parameters
#' set in the first step are used.
#'
#' When "INV.TCC" method is used employing the IRT 3-parameter logistic model (3PLM) in a test, ability estimates for the raw sum scores less than the sum of item
#' guessing parameters are not attainable. In this case, either of linear interpolation and linear extrapolation can be applied. Note that
#' if \code{constraint = TRUE}, linear interpolation method is used. Otherwise, linear extrapolation method is used. Let \eqn{\theta_{min}} and
#' \eqn{\theta_{max}} be the minimum and maximum ability estimates and \eqn{\theta_{X}} be the ability estimate for the smallest raw score, X, greater than or equal
#' to the sum of item guessing parameters. When linear interpolation method is used, a linear line is constructed between two points of (x=\eqn{\theta_{min}}, y=0) and
#' (x=\eqn{\theta_{X}}, y=X). Because \code{constraint = TRUE}, \eqn{\theta_{min}} is the first value in the argument \code{range.tcc}.
#' When linear extrapolation method is used, a linear line is constructed using two points of (x=\eqn{\theta_{X}}, y=X) and
#' (x=\eqn{\theta_{max}}, y=maximum raw score). Then, ability estimates for the raw sum scores between zero and the smallest raw score greater than or equal
#' to the sum of item guessing parameters are found using the constructed linear line.
#'
#' To speed up the ability estimation for MLE, MLEF, MAP, and EAP methods, this function applies a parallel process using multiple logical CPU cores.
#' You can set the number of logical CPU cores by specifying a positive integer value in the argument \code{ncore}. Default value is 1.
#'
#' Note that the standard errors of ability estimates are computed using observed information functions for MLE, MLEF, and MAP methods.
#'
#' @return A list including a vector of the ability estimates and a vector of the standard errors of ability estimates. When \code{method} is
#' "EAP.SUM" or "INV.TCC", raw sum scores of examinees and a table with the possible raw sum scores and corresponding ability estimates are returned as well.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
#'
#' @references
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a microcomputer environment. \emph{Psychometrika, 35}, 179-198.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.
#'
#' Han, K. T. (2016). Maximum likelihood score estimation method with fences for short-length tests and computerized adaptive tests.
#' \emph{Applied psychological measurement, 40}(4), 289-301.
#'
#' Kolen, M. J. & Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
#' Springer
#'
#' Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
#' \emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.
#'
#' Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
#' \emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.
#'
#' Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
#' \emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
#' for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
#' Measurement, 19}(1), 39-49.
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # generate examinees abilities
#' set.seed(12)
#' theta <- rnorm(10)
#'
#' # simulate the item response data
#' data <- simdat(x, theta, D=1)
#'
#' \donttest{
#' # estimate the abilities using MLE
#' est_score(x, data, D=1, method="MLE", range=c(-4, 4), se=TRUE, ncore=2)
#'
#' # estimate the abilities using MLEF with default fences of item difficulty parameters
#' est_score(x, data, D=1, method="MLEF", fence.a=3.0, fence.b=NULL, se=TRUE, ncore=2)
#'
#' # estimate the abilities using MLEF with different fences of item difficulty parameters
#' est_score(x, data, D=1, method="MLEF", fence.a=3.0, fence.b=c(-5, 5), se=TRUE, ncore=2)
#'
#' # estimate the abilities using MAP
#' est_score(x, data, D=1, method="MAP", norm.prior=c(0, 1), nquad=30, se=TRUE, ncore=2)
#'
#' # estimate the abilities using EAP summed scoring
#' est_score(x, data, D=1, method="EAP.SUM", norm.prior=c(0, 1), nquad=30)
#'
#' # estimate the abilities using inverse TCC scoring
#' est_score(x, data, D=1, method="INV.TCC", constant=0.1, constraint=FALSE, range.tcc=c(-7, 7))
#'
#' }
#'
#' @export
est_score <- function(x, ...) UseMethod("est_score")

#' @describeIn est_score Default method to estimate examinees' latent ability parameters using a data.frame \code{x} containing the item meta data.
#'
#' @export
est_score.default <- function(x, data, D = 1, method = "MLE", range = c(-4, 4), norm.prior = c(0, 1),
                      nquad = 41, weights = NULL, fence.a = 3.0, fence.b = NULL, se = TRUE,
                      constant=0.1, constraint=FALSE, range.tcc=c(-7, 7), missing = NA, ncore=1, ...) {

  method <- toupper(method)

  # check if the data set is a vector of an examinee
  if(is.vector(data)) {
    data <- rbind(data)
  }

  # scoring of MLE, MAP, and EAP
  if(method %in% c("MLE", "MAP", "EAP", "MLEF")) {

    # check the number of examinees
    nstd <- nrow(data)

    # recode missing values
    if(!is.na(missing)) {
      data[data == missing] <- NA
    }

    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

    # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
    if(ncol(x[, -c(1, 2, 3)]) == 2) {
      x <- data.frame(x, par.3=NA)
    }

    # add two more items and data responses when "MLE" with Fences method is used
    if(method == "MLEF") {
      if(is.null(fence.b)) {
        # find the range of b-parameters in the item meta data
        range.b <- range(x[, 4])
        range.b[1] <- floor(range.b[1] - 0.001)
        range.b[2] <- ceiling(range.b[2] + 0.001)

        # adjust the range of b-parameters to be used as a fence
        fence.b[1] <- ifelse(range.b[1] >= -3.5, -3.5, range.b[1])
        fence.b[2] <- ifelse(range.b[2] >= 3.5, range.b[2], 3.5)
      }

      # add two more response columns for the two fence items
      data <- data.frame(data, f.lower=rep(1, nstd), f.upper=rep(0, nstd))

      # create a new item meta data for the two fence items
      x.fence <- shape_df(par.dc=list(a=rep(fence.a, 2), b=fence.b, g=0),
                          item.id=c("fence.lower", "fence.upper"), cats=rep(2, 2), model="3PLM")
      if(ncol(x) > ncol(x.fence)) {
        add.colnum <- ncol(x) - ncol(x.fence)
        x.fence <- data.frame(x.fence, matrix(NA, nrow=2, ncol=add.colnum))
        colnames(x.fence) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x.fence) - 3)))
      }

      # create the new item meta data by adding two fence items
      x <- rbind(x, x.fence)
    }

    # listrize the data.frame
    meta <- metalist2(x)

    # create equations for gradient vector and hessian matrix
    if(!is.null(meta$drm)) {
      # create equations for gradient vector and hessian matrix
      eq_grad_drm <- equation_drm(model="3PLM", use.pprior=FALSE, hessian=FALSE, type="ability")$params_fun
      eq_hess_drm <- equation_drm(model="3PLM", use.pprior=FALSE, hessian=TRUE, type="ability")$params_fun
    } else {
      eq_grad_drm <- NULL
      eq_hess_drm <- NULL
    }
    if(!is.null(meta$plm)) {
      eq_grad_plm <-
        purrr::map(.x=1:length(meta$plm$cats),
                   .f=function(i) equation_plm(cats=meta$plm$cats[i], pmodel=meta$plm$model[i], use.pprior=FALSE,
                                               hessian=FALSE, type="ability")$params_fun)
      eq_hess_plm <-
        purrr::map(.x=1:length(meta$plm$cats),
                   .f=function(i) equation_plm(cats=meta$plm$cats[i], pmodel=meta$plm$model[i], use.pprior=FALSE,
                                               hessian=TRUE, type="ability")$params_fun)
    } else {
      eq_grad_plm <- NULL
      eq_hess_plm <- NULL
    }
    if(method == "MAP") {
      eq_grad_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=FALSE, type="ability")$pprior_fun
      eq_hess_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=TRUE, type="ability")$pprior_fun
    } else {
      eq_grad_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=FALSE, type="ability")$pprior_fun
      eq_hess_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=TRUE, type="ability")$pprior_fun
    }

    # create lists of the equations
    FUN.grad=list(drm=eq_grad_drm, plm=eq_grad_plm, prior=eq_grad_prior)
    FUN.hess=list(drm=eq_hess_drm, plm=eq_hess_plm, prior=eq_hess_prior)

    # check the number of CPU cores
    if(ncore < 1) {
      stop("The number of logical CPU cores must not be less than 1.", call.=FALSE)
    }

    # estimation
    if(ncore == 1L) {

      # set a function for scoring
      f <- function(i) est_score_indiv(meta=meta, resp=data[i, ], D=D, method=method, range=range,
                                       norm.prior=norm.prior, nquad=nquad, weights=weights, se=se,
                                       FUN.grad=FUN.grad, FUN.hess=FUN.hess)

      # scoring
      est <- purrr::map(.x=1:nrow(data), .f=function(i) f(i))

      # assign estimated values
      est.theta <- purrr::map_dbl(est, .f=function(x) x$est.theta)
      if(se) {
        se.theta <- purrr::map_dbl(est, .f=function(x) x$se.theta)
      } else {
        se.theta <- NULL
      }

      rst <- list(est.theta=est.theta, se.theta=se.theta)

    } else {

      # specify the number of CPU cores
      numCores <- ncore

      # create a parallel processesing cluster
      cl = parallel::makeCluster(numCores, ...)

      # load some specific variable names into processing cluster
      parallel::clusterExport(cl, c("meta", "data", "D", "method", "range",
                                    "norm.prior", "nquad", "weights", "se",
                                    "FUN.grad", "FUN.hess",
                                    "est_score_indiv", "loglike_score", "grad_score", "hess_score",
                                    "ll_brute",
                                    "drm", "plm", "grm", "gpcm", "gen.weight"), envir = environment())
      parallel::clusterEvalQ(cl, library(dplyr))

      # set a function for scoring
      f <- function(i) est_score_indiv(meta=meta, resp=data[i, ], D=D, method=method, range=range,
                                       norm.prior=norm.prior, nquad=nquad, weights=weights, se=se,
                                       FUN.grad=FUN.grad, FUN.hess=FUN.hess)

      # parallel scoring
      est <- pbapply::pblapply(X=1:nrow(data), FUN=f, cl=cl) # to see the progress bar

      # finish
      parallel::stopCluster(cl)

      # assign estimated values
      est.theta <- purrr::map_dbl(est, .f=function(x) x$est.theta)
      if(se) {
        se.theta <- purrr::map_dbl(est, .f=function(x) x$se.theta)
      } else {
        se.theta <- NULL
      }

      rst <- list(est.theta=est.theta, se.theta=se.theta)

    }

  }

  if(method == "EAP.SUM") {
    if(is.null(weights)) {
      rst <- eap_sum(x, data, norm.prior=norm.prior, nquad=nquad, D=D)
    } else {
      rst <- eap_sum(x, data, weights=weights, D=D)
    }
  }

  if(method == "INV.TCC") {
    rst <- inv_tcc(x, data, D=D, constant=constant, constraint=constraint, range.tcc=range.tcc)
  }

  # return results
  rst

}

#' @describeIn est_score An object created by the function \code{\link{est_irt}}.
#'
#' @export
#'
est_score.est_irt <- function(x, method = "MLE", range = c(-4, 4), norm.prior = c(0, 1),
                              nquad = 41, weights = NULL, fence.a = 3.0, fence.b = NULL, se = TRUE,
                              constant=0.1, constraint=FALSE, range.tcc=c(-7, 7), missing = NA, ncore=1, ...) {

  method <- toupper(method)

  # extract information from an object
  data <- x$data
  D <- x$scale.D
  x <- x$par.est

  # check if the data set is a vector of an examinee
  if(is.vector(data)) {
    data <- rbind(data)
  }

  # scoring of MLE, MAP, and EAP
  if(method %in% c("MLE", "MAP", "EAP", "MLEF")) {

    # check the number of examinees
    nstd <- nrow(data)

    # recode missing values
    if(!is.na(missing)) {
      data[data == missing] <- NA
    }

    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

    # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
    if(ncol(x[, -c(1, 2, 3)]) == 2) {
      x <- data.frame(x, par.3=NA)
    }

    # add two more items and data responses when "MLE" with Fences method is used
    if(method == "MLEF") {
      if(is.null(fence.b)) {
        # find the range of b-parameters in the item meta data
        range.b <- range(x[, 4])
        range.b[1] <- floor(range.b[1] - 0.001)
        range.b[2] <- ceiling(range.b[2] + 0.001)

        # adjust the range of b-parameters to be used as a fence
        fence.b[1] <- ifelse(range.b[1] >= -3.5, -3.5, range.b[1])
        fence.b[2] <- ifelse(range.b[2] >= 3.5, range.b[2], 3.5)
      }

      # add two more response columns for the two fence items
      data <- data.frame(data, f.lower=rep(1, nstd), f.upper=rep(0, nstd))

      # create a new item meta data for the two fence items
      x.fence <- shape_df(par.dc=list(a=rep(fence.a, 2), b=fence.b, g=0),
                          item.id=c("fence.lower", "fence.upper"), cats=rep(2, 2), model="3PLM")
      if(ncol(x) > ncol(x.fence)) {
        add.colnum <- ncol(x) - ncol(x.fence)
        x.fence <- data.frame(x.fence, matrix(NA, nrow=2, ncol=add.colnum))
        colnames(x.fence) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x.fence) - 3)))
      }

      # create the new item meta data by adding two fence items
      x <- rbind(x, x.fence)
    }

    # listrize the data.frame
    meta <- metalist2(x)

    # create equations for gradient vector and hessian matrix
    if(!is.null(meta$drm)) {
      # create equations for gradient vector and hessian matrix
      eq_grad_drm <- equation_drm(model="3PLM", use.pprior=FALSE, hessian=FALSE, type="ability")$params_fun
      eq_hess_drm <- equation_drm(model="3PLM", use.pprior=FALSE, hessian=TRUE, type="ability")$params_fun
    } else {
      eq_grad_drm <- NULL
      eq_hess_drm <- NULL
    }
    if(!is.null(meta$plm)) {
      eq_grad_plm <-
        purrr::map(.x=1:length(meta$plm$cats),
                   .f=function(i) equation_plm(cats=meta$plm$cats[i], pmodel=meta$plm$model[i], use.pprior=FALSE,
                                               hessian=FALSE, type="ability")$params_fun)
      eq_hess_plm <-
        purrr::map(.x=1:length(meta$plm$cats),
                   .f=function(i) equation_plm(cats=meta$plm$cats[i], pmodel=meta$plm$model[i], use.pprior=FALSE,
                                               hessian=TRUE, type="ability")$params_fun)
    } else {
      eq_grad_plm <- NULL
      eq_hess_plm <- NULL
    }
    if(method == "MAP") {
      eq_grad_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=FALSE, type="ability")$pprior_fun
      eq_hess_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=TRUE, type="ability")$pprior_fun
    } else {
      eq_grad_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=FALSE, type="ability")$pprior_fun
      eq_hess_prior <- equation_drm(model="3PLM", pprior=list(dist="norm", params=norm.prior), use.pprior=TRUE,
                                    hessian=TRUE, type="ability")$pprior_fun
    }

    # create lists of the equations
    FUN.grad=list(drm=eq_grad_drm, plm=eq_grad_plm, prior=eq_grad_prior)
    FUN.hess=list(drm=eq_hess_drm, plm=eq_hess_plm, prior=eq_hess_prior)

    # check the number of CPU cores
    if(ncore < 1) {
      stop("The number of logical CPU cores must not be less than 1.", call.=FALSE)
    }

    # estimation
    if(ncore == 1L) {

      # set a function for scoring
      f <- function(i) est_score_indiv(meta=meta, resp=data[i, ], D=D, method=method, range=range,
                                       norm.prior=norm.prior, nquad=nquad, weights=weights, se=se,
                                       FUN.grad=FUN.grad, FUN.hess=FUN.hess)

      # scoring
      est <- purrr::map(.x=1:nrow(data), .f=function(i) f(i))

      # assign estimated values
      est.theta <- purrr::map_dbl(est, .f=function(x) x$est.theta)
      if(se) {
        se.theta <- purrr::map_dbl(est, .f=function(x) x$se.theta)
      } else {
        se.theta <- NULL
      }

      rst <- list(est.theta=est.theta, se.theta=se.theta)

    } else {

      # specify the number of CPU cores
      numCores <- ncore

      # create a parallel processesing cluster
      cl = parallel::makeCluster(numCores, ...)

      # load some specific variable names into processing cluster
      parallel::clusterExport(cl, c("meta", "data", "D", "method", "range",
                                    "norm.prior", "nquad", "weights", "se",
                                    "FUN.grad", "FUN.hess",
                                    "est_score_indiv", "loglike_score", "grad_score", "hess_score",
                                    "ll_brute",
                                    "drm", "plm", "grm", "gpcm", "gen.weight"), envir = environment())
      parallel::clusterEvalQ(cl, library(dplyr))

      # set a function for scoring
      f <- function(i) est_score_indiv(meta=meta, resp=data[i, ], D=D, method=method, range=range,
                                       norm.prior=norm.prior, nquad=nquad, weights=weights, se=se,
                                       FUN.grad=FUN.grad, FUN.hess=FUN.hess)

      # parallel scoring
      est <- pbapply::pblapply(X=1:nrow(data), FUN=f, cl=cl) # to see the progress bar

      # finish
      parallel::stopCluster(cl)

      # assign estimated values
      est.theta <- purrr::map_dbl(est, .f=function(x) x$est.theta)
      if(se) {
        se.theta <- purrr::map_dbl(est, .f=function(x) x$se.theta)
      } else {
        se.theta <- NULL
      }

      rst <- list(est.theta=est.theta, se.theta=se.theta)

    }

  }

  if(method == "EAP.SUM") {
    if(is.null(weights)) {
      rst <- eap_sum(x, data, norm.prior=norm.prior, nquad=nquad, D=D)
    } else {
      rst <- eap_sum(x, data, weights=weights, D=D)
    }
  }

  if(method == "INV.TCC") {
    rst <- inv_tcc(x, data, D=D, constant=constant, constraint=constraint, range.tcc=range.tcc)
  }

  # return results
  rst

}


# This function computes an abiltiy estimate for each examinee
est_score_indiv <- function(meta, resp, D = 1, method = "MLE", range = c(-4, 4), norm.prior = c(0, 1),
                            nquad = 41, weights=NULL, se = TRUE,
                            FUN.grad=list(drm=NULL, plm=NULL, prior=NULL),
                            FUN.hess=list(drm=NULL, plm=NULL, prior=NULL)) {

  # extract equations of a gradient and hessian
  eq_grad_drm <- FUN.grad$drm
  eq_grad_plm <- FUN.grad$plm
  eq_grad_prior <- FUN.grad$prior
  eq_hess_drm <- FUN.hess$drm
  eq_hess_plm <- FUN.hess$plm
  eq_hess_prior <- FUN.hess$prior

  # find missing data and
  # delete missing data from item data set and response data
  if(any(is.na(resp))) {
    loc.miss <- which(is.na(resp)) # check the locations of missing data

    # delete missing data from the item meta data
    # delete equations for the polytomous IRT models if missing data exist
    if(!is.null(meta$drm)) {
      meta$drm <- purrr::map(.x=meta$drm, .f=function(x) x[!meta$drm$loc %in% loc.miss])
    }
    if(!is.null(meta$plm)) {
      tmp <- purrr::map(.x=meta$plm, .f=function(x) x[!meta$plm$loc %in% loc.miss])
      eq_grad_plm <- eq_grad_plm[!meta$plm$loc %in% loc.miss]
      eq_hess_plm <- eq_hess_plm[!meta$plm$loc %in% loc.miss]
      meta$plm <- tmp
    }

    # create lists of the equations using the updated information for missing data
    FUN.grad <- list(drm=eq_grad_drm, plm=eq_grad_plm, prior=eq_grad_prior)
    FUN.hess <- list(drm=eq_hess_drm, plm=eq_hess_plm, prior=eq_hess_prior)

  }

  # factorize the response values
  max.cats <- max(c(meta$drm$cats, meta$plm$cats))
  resp.f <- factor(resp, levels=(seq_len(max.cats) - 1))

  # calculate the score categories
  tmp.id <- 1:length(resp.f)
  freq.cat <-
    stats::xtabs(~ tmp.id + resp.f, addNA=TRUE) %>%
    data.matrix()
  if(any(is.na(resp))) freq.cat <- freq.cat[, -ncol(freq.cat)]
  if(!is.null(meta$drm)) {
    if(length(meta$drm$id) == 1L) {
      freq.cat_drm <- rbind(freq.cat[meta$drm$loc, 1:2])
    } else {
      freq.cat_drm <- freq.cat[meta$drm$loc, 1:2]
    }
  } else {
    freq.cat_drm <- NULL
  }
  if(!is.null(meta$plm)) {
    if(length(meta$plm$id) == 1L) {
      freq.cat_plm <- rbind(freq.cat[meta$plm$loc, ])
    } else {
      freq.cat_plm <- freq.cat[meta$plm$loc, ]
    }
  } else {
    freq.cat_plm <- NULL
  }
  freq.cat <- list(freq.cat_drm=freq.cat_drm, freq.cat_plm=freq.cat_plm)

  ##----------------------------------------------------
  ## MLE and MAP scorings
  if(method %in% c("MLE", "MAP", "MLEF")) {

    # compute a perfect NC score
    total.nc <- sum(c(meta$drm$cats, meta$plm$cats) - 1)

    # compute a ininal value for scoring
    if(sum(resp, na.rm=TRUE) == 0) {
      startval <- log(1 / total.nc)
    }
    if(sum(resp, na.rm=TRUE) == total.nc) {
      startval <- log(total.nc / 1)
    }
    if(!sum(resp, na.rm=TRUE) %in% c(0, total.nc)) {
      startval <- log(sum(resp, na.rm=TRUE) / (total.nc - sum(resp, na.rm=TRUE)))
    }
    # if(method == "MLE") {
    #   startval_tmp <- c(seq(from=range[1], to=range[2], length.out=50), startval)
    #   ll_tmp <- ll_brute(theta=startval_tmp, meta=meta, freq.cat=freq.cat, method="MLE", D=D, norm.prior=norm.prior)
    #   startval <- startval_tmp[which.min(ll_tmp)]
    # }

    # estimate an abiltiy and SE
    if(method %in% c("MLE", "MLEF")) {

      # find a better starting value for MLE using a brute force method
      # prepare the discrete theta values
      startval_tmp <- c(seq(from=range[1], to=range[2], length.out=102))

      # compute the negative log-likelihood values for all the theta values
      ll_tmp <- ll_brute(theta=startval_tmp, meta=meta, freq.cat=freq.cat, method="MLE", D=D)

      # find the locations of thetas where the sign of slope changes
      # loc_change <- which(diff(sign(diff(ll_tmp))) != 0L) + 1
      loc_change <- which(diff(sign(diff(ll_tmp))) > 0L) + 1

      # select a theta value that has the minimum of negative log-likelihood value
      startval_tmp1 <- startval_tmp[loc_change][which.min(ll_tmp[loc_change])]

      # if there is no selected value from step 2, use the starting value from step 1
      startval <- ifelse(length(startval_tmp1) > 0L, startval_tmp1, startval)

      # estimation
      est.mle <-
        suppressWarnings(tryCatch({stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                                                 norm.prior=norm.prior, logL=TRUE,
                                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                                 gradient=grad_score, hessian=hess_score,
                                                 lower=range[1], upper=range[2])}, error = function(e) {NULL}))
      # when the estimation returns an error message
      if(is.null(est.mle)) {
        est.mle <- stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                                 norm.prior=norm.prior, logL=TRUE,
                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                 lower=range[1], upper=range[2])
      }
      if(est.mle$convergence == 0L) {
        est.theta <- est.mle$par
      } else {
        # this is to estimate an ability in a brute force way
        # when convergence is failed
        brute_theta <- seq(range[1], range[2], 0.001)
        brute_ll <- ll_brute(theta=brute_theta, meta=meta, freq.cat=freq.cat, method="MLE", D=D, norm.prior=norm.prior)
        est.theta <- brute_theta[which.min(brute_ll)]
      }

      if(se) {
        if(est.theta %in% range) {
          se.theta <- 99.9999
        } else {
          hess <- hess_score(theta=est.theta, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                             norm.prior=norm.prior, logL=TRUE,
                             FUN.grad=FUN.grad, FUN.hess=FUN.hess)
          # to prevent the case when the hessian has a negative value
          if(hess < 0L) {
            hess <- stats::optimHess(par=est.theta, fn=loglike_score, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                                     norm.prior=norm.prior, logL=TRUE,
                                     FUN.grad=FUN.grad, FUN.hess=FUN.hess)
          }
          se.theta <- as.numeric(sqrt(1 / hess))
        }
      } else {
        se.theta <- NULL
      }
    }
    if(method == "MAP") {
      est.map <-
        suppressWarnings(tryCatch({stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MAP", D=D,
                                                 norm.prior=norm.prior, logL=TRUE,
                                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                                 gradient=grad_score, hessian=hess_score,
                                                 lower=range[1], upper=range[2])}, error = function(e) {NULL}))
      # when the estimation returns an error message
      if(is.null(est.map)) {
        est.map <- stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MAP", D=D,
                                 norm.prior=norm.prior, logL=TRUE,
                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                 lower=range[1], upper=range[2])
      }
      if(est.map$convergence == 0L) {
        est.theta <- est.map$par
      } else {
        # this is to estimate an ability in a brute force way
        # when convergence is failed
        brute_theta <- seq(range[1], range[2], 0.001)
        brute_ll <- ll_brute(theta=brute_theta, meta=meta, freq.cat=freq.cat, method="MAP", D=D, norm.prior=norm.prior)
        est.theta <- brute_theta[which.min(brute_ll)]
      }

      if(se) {
        hess <- hess_score(theta=est.theta, meta=meta, freq.cat=freq.cat, method="MAP", D=D,
                           norm.prior=norm.prior, logL=TRUE,
                           FUN.grad=FUN.grad, FUN.hess=FUN.hess)
        # to prevent the case when the hessian has a negative value
        if(hess < 0L) {
          hess <- stats::optimHess(par=est.theta, fn=loglike_score, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                                   norm.prior=norm.prior, logL=TRUE,
                                   FUN.grad=FUN.grad, FUN.hess=FUN.hess)
        }
        se.theta <- as.numeric(sqrt(1 / hess))
      } else {
        se.theta <- NULL
      }
    }
  }

  ##----------------------------------------------------
  ## EAP scoring
  if(method == "EAP") {

    # generate quadrature points and weights
    if(is.null(weights)) {
      popdist <- gen.weight(n=nquad, dist="norm", mu=norm.prior[1], sigma=norm.prior[2])
    } else {
      popdist <- data.frame(weights)
    }

    # estimating posterior dist
    posterior <- purrr::map2_dbl(.x=popdist[, 1], .y=popdist[, 2],
                                 .f=function(k, y) loglike_score(theta=k, meta=meta, freq.cat=freq.cat, D=D, logL=FALSE) * y)

    # Expected A Posterior
    posterior <- posterior / sum(posterior)
    est.theta <- sum(popdist[, 1] * posterior)

    if(se) {
      # calculating standard error
      ex2 <- sum(popdist[, 1]^2 * posterior)
      var <- ex2 - (est.theta)^2
      se.theta <- sqrt(var)
    } else {
      se.theta <- NULL
    }

  }

  # return results
  rst <- list(est.theta=est.theta, se.theta=se.theta)

  rst

}
