#' Asymptotic variance-covariance matrices of item parameter estimates
#'
#' @description This function calculates the analytical asymptotic variance-covariance matrices (e.g., Li & Lissitz, 2004; Thissen & Wainer, 1982)
#' of item parameter estimates for dichotomous and polytomous IRT Models without examinee's responses to test items, 
#' given a set of item parameter estimates. The square roots of variance terms in the matrices can be used as the asymptotic 
#' standard errors of maximum likelihood item parameter estimates.
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
#' This data frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param nstd An integer value or a vector of integer values indicating a sample size. When a vector is specified, length of the vector must be
#' the same as the number of test items in the argument \code{x}. Default is 1,000. See below for details.
#' @param pcm.loc A vector of integer values indicating the locations of partial credit model (PCM) items. For the PCM items,
#' the variance-covariance matrices are computed only for the item category difficulty parameters. Default is NULL. See below for details.
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution.
#' Default is c(0,1).
#' @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
#' @param weights A two-column matrix or data frame containing the theta values (in the first column) and the weights (in the second column)
#' for the prior distribution. The weights and theta values can be easily obtained using the function \code{\link{gen.weight}}.
#' If NULL, default values are used for the prior distribution (see the arguments of \code{norm.prior} and \code{nquad}). Default is NULL.
#'
#' @details 
#' The standard errors obtained from the analytical approach are likely to represent lower bounds for the actual standard errors (Thissen & Wainer, 1982).
#' Therefore, they may be useful for assessing the degree of precision of a set of item parameter estimates when the corresponding standard errors of 
#' the estimates are presented in literature or research reports.
#' 
#' Sometimes item parameters need to be estimated using different sample size. If the item parameters in the argument \code{x} were
#' calibrated with different number of examinees, a vector of different sample sizes should be specified in the argument \code{nstd}. Suppose
#' that you want to compute the variance-covariance matrices of five IRT 3PLM items and the five items were calibrated with 500, 600, 1,000, 2,000,
#' and 700 examinees, respectively. Then, \code{nstd = c(500, 600, 1000, 2000, 700)} must be specified.
#'
#' Because you can specify only "GPCM" for both the partial credit model (PCM) or the generalized partial credit model (GPCM) in the item metadata,
#' you must indicate which items are the PCM items through the argument \code{pcm.loc}. This is because the item category difficulty parameters are estimated
#' from the PCM, meaning that the variance-covariance of item parameter estimates must be computed for the item category difficulty parameters. Suppose
#' that you want to compute the variance-covariance matrices of five polytomous items and the last two items were calibrated with the PCM. Then,
#' \code{pcm.loc = c(4, 5)} must be specified.
#'
#' @return A list of two internal objects. The first internal object contains a list of the variance-covariance matrices of item parameter estimates.
#' The second internal object contains a list of the standard errors of item parameter estimates.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Li, Y. & Lissitz, R. (2004). Applications of the analytically derived asymptotic standard errors of item response theory
#' item parameter estimates. \emph{Journal of educational measurement, 41}(2), 85-117.
#' 
#' Thissen, D. & Wainer, H. (1982). Weighted likelihood estimation of ability in item response theory. 
#' \emph{Psychometrika, 54}(3), 427-450.
#'
#' @seealso \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{gen.weight}}
#'
#' @export
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained sfrom a flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the first two dichotomous items and last polytomous item
#' x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # compute the var-covariance matrices with sample size of 2,000
#' covirt(x, D=1, nstd=2000, norm.prior=c(0, 1), nquad=40)
#'
covirt <- function(x, D=1, nstd=1000, pcm.loc=NULL, norm.prior=c(0, 1), nquad=41, weights=NULL) {


  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }

  # clear the item metadata set
  x <- back2df(metalist2(x))

  # consider DRM as 3PLM
  x[, 3] <- as.character(x[, 3])
  # consider DRM as 3PLM
  if("DRM" %in% x[, 3]) {
    x[, 3][x[, 3] == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items were considered as '3PLM' items in during the item parameter estimation."
    warning(memo, call.=TRUE)
  }

  # listlize the item metadata
  meta_list <- purrr::map(1:nrow(x), .f=function(i) metalist2(x[i, ]))

  # specify the locations of PCM items
  pcm.a.logic <- rep(FALSE, nrow(x))
  if(!is.null(pcm.loc)) {
    pcm.a.logic[pcm.loc] <- TRUE
  }

  # generate quadrature points and weights
  if(is.null(weights)) {
    wts <- gen.weight(n=nquad, dist="norm", mu=norm.prior[1], sigma=norm.prior[2])
  } else {
    wts <- data.frame(weights)
  }

  # sample size
  if(length(nstd) == 1L) nstd <- rep(nstd, nrow(x))
  if(length(nstd) != nrow(x)) {
    stop("Sample size must be an integer value or a vector with a length of total items.", call.=FALSE)
  }

  # compute the var-covariance matrix
  cov_par <- purrr::pmap(.l=list(x=meta_list, y=nstd, z=pcm.a.logic), .f=function(x, y, z) cov_mat(meta=x, D=D, nstd=y, wts=wts, pcm.a=z))
  se_par <- purrr::map(cov_par, .f=function(x) sqrt(diag(x)))
  names(cov_par) <- x[, 1]
  names(se_par) <- x[, 1]

  # return results
  rst <- list(cov=cov_par, se=se_par)
  rst

}



# This function computes the analytical var-covariance matrix of item parameter estimates for an item
cov_mat <- function(meta, D=1, nstd = 1000, wts, pcm.a=FALSE) {

  info_mat <- integrand(meta, theta=wts[, 1], dens=wts[, 2], D=D, pcm.a=pcm.a)

  rst <- tryCatch({solve(info_mat * nstd)}, error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
  if(is.null(rst)) {
    rst <- matrix(NA, nrow=dim(info_mat)[1], ncol=dim(info_mat)[2])
  }

  rst

}


# This function computes the kernel of interation when computing
# the analytical var-covariance matrix of item parameter estimates
integrand <- function(meta, theta, dens, D=1, pcm.a=FALSE) {

  # For a dichotomous item
  if(!is.null(meta$drm$a)) {

    # extract information from metadata
    a <- meta$drm$a
    b <- meta$drm$b
    g <- meta$drm$g
    item_par <- c(a, b, g)
    cats <- meta$drm$cats
    model <- meta$drm$model

    # compute a product of the probabilities p and q
    pq <- apply(trace3(meta, theta=theta, D=D, type="icc")$trace[[1]], 1, prod)

    # prevent that pq has zero probability
    pq <- ifelse(pq == 0L, 1e-20, pq)

    # create the gradient for score category probability
    funList <- equation_scocat(model=model, cats=NULL, hessian=FALSE, type="item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for(i in 1:3) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    # select a function for each score category
    fun.tmp <- funList[[1]]

    # implement the fuction for each score category
    tmp <- do.call("fun.tmp", args.list, envir=environment())

    # compute the kernel of integration
    p.d1 <- attributes(tmp)$gradient
    info.list  <- purrr::map(.x=1:length(theta),
                             .f=function(i) (1 / pq[i]) * outer(p.d1[i, ], p.d1[i, ]) * dens[i])
    info.mat <- Reduce(f='+', x=info.list)

  }

  # For a polytomous item
  if(!is.null(meta$plm$a)) {

    # extract information from metadata
    a <- meta$plm$a
    d <- meta$plm$d[[1]]
    cats <- meta$plm$cats
    model <- meta$plm$model
    item_par <- c(a, d)

    # compute a product of the probabilities p and q
    ps <-  trace3(meta, theta=theta, D=D, type="icc")$trace[[1]]
    pqs <- ps * (1 - ps)

    # prevent that pq has zero probability
    pqs[pqs == 0L] <- 1e-20

    # create the gradient for score category probability
    funList <- equation_scocat(model=model, cats=cats, fix.a.gpcm=pcm.a, hessian=FALSE, type="item")

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for(i in 1:cats) {
      args.pars[[i]] <- item_par[i]
    }
    args.list <- args.pars
    args.list$theta <- theta
    args.list$D <- D

    # compute the score category information
    info.mat <- vector('list', cats)
    for(i in 1:cats) {
      # select a function for each score category
      fun.tmp <- funList[[i]]

      # implement the fuction for each score category
      tmp <- do.call("fun.tmp", args.list, envir=environment())

      # compute the kernel of integration
      ps.d1 <- attributes(tmp)$gradient
      info.list  <- purrr::map(.x=1:length(theta),
                               .f=function(k) (1 / pqs[k, i]) * outer(ps.d1[k, ], ps.d1[k, ]) * dens[k])
      info.mat[[i]] <- Reduce(f='+', x=info.list)
    }

    info.mat <- Reduce(f='+', x=info.mat)

  }

  info.mat

}


