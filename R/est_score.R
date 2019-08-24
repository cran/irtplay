#' Estimate examinees' ability (proficiency) parameters
#'
#' @description This function estimates examinees' latent ability parameters. Available scoring methods are maximum likelihood estimation (MLE),
#' maximum a posteriori estimation (MAP; Hambleton et al., 1991), expected a posteriori estimation (EAP; Bock & Mislevy, 1982),
#' and EAP summed scoring (Thissen et al., 1995; Thissen & Orlando, 2001).
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data.
#' This data.frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param method A character string indicating a scoring method. Available methods are "MLE" for
#' the maximum likelihood estimation, "MAP" for the maximum a posteriori estimation, "EAP" for the expected a posteriori estimation,
#' and "EAP.SUM" for the expected a posteriori summed scoring. Default method is "MLE".
#' @param range A numeric vector of two components to restrict the range of ability scale for the MLE. Default is c(-4, 4).
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "MLE".
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' Ignored if \code{method} is "MLE" or "MAP".
#' @param weights A two-column matrix or data.frame containing the theta values (in the first column) and the weights (in the second column)
#' for the prior distribution. The weights and theta values can be easily obtained using the function \code{\link{gen.weight}}.
#' If NULL and \code{method} is "EAP" or "EAP.SUM", default values are used (see the arguments of \code{norm.prior} and \code{nquad}). Ignored
#' if \code{method} is "MLE" or "MAP".
#' @param se A logical value. If TRUE, the standard errors of ability estimates are computed. However, if \code{method} is "EAP.SUM", the standard
#' errors are always returned.
#' @param missing A value indicating missing values in the response data set. Default is NA. See below for details.
#' @param ncore The number of logical CPU cores to use. Default is 1. See below for details.
#' @param ... additional arguments to pass to \code{parallel::makeCluster}.
#'
#' @details For MAP scoring method, only the normal prior distribution is available for the population distribution.
#'
#' When there are missing data in the response data set, the missing value must be specified in \code{missing}. The missing data are taken into account
#' when either of MLE, MAP, and EAP is used. However, there must be no missing data in the response data set when "EAP.SUM" is used.
#' One of possible ways to use "EAP.SUM" method when missing values exist is to remove rows with any missing values.
#'
#' To speed up the ability estimation for MLE, MAP, and EAP methods, this function applies a parallel process using multiple logical CPU cores.
#' You can set the number of logical CPU cores by specifying a positive integer value in the argument \code{ncore}. Default value is 1.
#'
#' Note that the standard errors of ability estimates are computed using observed information functions.
#'
#' @return A list including a vector of the ability estimates and a vector of the standard errors of ability estimates. When "EAP.SUM" is used,
#' raw sum scores of examinees and a table with the possible raw sum scores and corresponding ability estimates are returned as well.
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
#' Thissen, D. & Orlando, M. (2001). Item response theory for items scored in two categories. In D. Thissen & H. Wainer (Eds.),
#' \emph{Test scoring} (pp.73-140). Mahwah, NJ: Lawrence Erlbaum.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. (1995). Item Response Theory
#' for Scores on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
#' Measurement, 19}(1), 39-49.
#'
#' @export
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
#' # estimate the abilities using MAP
#' est_score(x, data, D=1, method="MAP", norm.prior=c(0, 1), nquad=30, se=TRUE, ncore=2)
#'
#' # estimate the abilities using EAP summed scoring
#' est_score(x, data, D=1, method="EAP.SUM", norm.prior=c(0, 1), nquad=30)
#' }
#'
est_score <- function(x, data, D = 1, method = "MLE", range = c(-4, 4), norm.prior = c(0, 1),
                      nquad = 41, weights = NULL, se = TRUE, missing = NA, ncore=1, ...) {

  method <- toupper(method)

  # check if the data set is a vector of an examinee
  if(is.vector(data)) {
    data <- rbind(data)
  }

  # scoring of MLE, MAP, and EAP
  if(method %in% c("MLE", "MAP", "EAP")) {

    # check the number of examinees
    nstd <- nrow(data)

    # recode missing values
    if(!is.na(missing)) {
      data[data == missing] <- NA
    }

    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

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
      eq_grad_prior <- NULL
      eq_hess_prior <- NULL
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
    if(length(meta$drm) == 1L) {
      freq.cat_drm <- rbind(freq.cat[meta$drm$loc, 1:2])
    } else {
      freq.cat_drm <- freq.cat[meta$drm$loc, 1:2]
    }
  } else {
    freq.cat_drm <- NULL
  }
  if(!is.null(meta$plm)) {
    if(length(meta$plm) == 1L) {
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
  if(method %in% c("MLE", "MAP")) {

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

    # estimate an abiltiy and SE
    if(method == "MLE") {
      est.theta <- stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                                 norm.prior=norm.prior, logL=TRUE,
                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                 gradient=grad_score, hessian=hess_score,
                                 lower=range[1], upper=range[2])$par
      if(se) {
        if(est.theta %in% range) {
          se.theta <- 99.9999
        } else {
          hess <- hess_score(theta=est.theta, meta=meta, freq.cat=freq.cat, method="MLE", D=D,
                             norm.prior=norm.prior, logL=TRUE,
                             FUN.grad=FUN.grad, FUN.hess=FUN.hess)
          se.theta <- as.numeric(sqrt(1 / hess))
        }
      } else {
        se.theta <- NULL
      }
    }
    if(method == "MAP") {
      est.theta <- stats::nlminb(startval, objective=loglike_score, meta=meta, freq.cat=freq.cat, method="MAP", D=D,
                                 norm.prior=norm.prior, logL=TRUE,
                                 FUN.grad=FUN.grad, FUN.hess=FUN.hess,
                                 gradient=grad_score, hessian=hess_score,
                                 lower=range[1], upper=range[2])$par
      if(se) {
        hess <- hess_score(theta=est.theta, meta=meta, freq.cat=freq.cat, method="MAP", D=D,
                           norm.prior=norm.prior, logL=TRUE,
                           FUN.grad=FUN.grad, FUN.hess=FUN.hess)
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
