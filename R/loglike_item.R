#' Loglikelihood of Items
#'
#' @description This function computes the loglikelihood of individual items given the item parameters, ability values, and response data.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}} or \code{\link{simdat}} for more details about the item meta data.
#' This data.frame can be easily obtained using the function \code{\link{shape_df}}. If \code{prob = NULL}, this data.frame is
#' used in the recursion formula. See below for details.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param score A vector of examinees' ability estimates. Length of the vector must be the same as the number of rows in the
#' response data set.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param use.aprior A logical value. If TRUE, a prior distribution for the slope parameters is used when computing the loglikelihood values
#' across all items. Default is FALSE.
#' @param use.gprior A logical value. If TRUE, a prior distribution for the guessing parameters is used when computing the loglikelihood values
#' across all 3PLM items. Default is TRUE.
#' @param aprior A list containing the information of the prior distribution for item slope parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{\link[stats]{dbeta}}, \code{\link[stats]{dlnorm}},
#' and \code{\link[stats]{dnorm}} for more details.
#' @param gprior A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{\link[stats]{dbeta}}, \code{\link[stats]{dlnorm}},
#' and \code{\link[stats]{dnorm}} for more details.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#'
#' @return A vector of loglikelihood values
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the first two dichotomous items and last polytomous item
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(10)
#' score <- rnorm(10, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' # compute the loglikelihood values (no priors are used)
#' llike_item(x, data, score, D=1, use.aprior=FALSE, use.gprior=FALSE)
#'
#' @import purrr
#' @import dplyr
#'
#' @export
llike_item <- function(x, data, score, D=1, use.aprior=FALSE, use.gprior=FALSE,
                       aprior=list(dist="lnorm", params=c(0, 0.5)), gprior=list(dist="beta", params=c(5, 17)),
                       missing=NA) {

  ##-------------------------------------------------------------------------------------------------------
  ## 1. preperation of data
  # extract information about the number of score cetegories and models
  cats <- x[, 2]
  model <-
    as.character(x[, 3]) %>%
    toupper()

  # check wheter included data are correct
  if(nrow(x) != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)

  # consider DRM as 3PLM
  if("DRM" %in% model) {
    model[model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items were considered as '3PLM' items."
    warning(memo, call.=TRUE)
  }

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to data.frame
  if(nrow(data) == 1L) {
    data <- list(data)
  } else {
    data <- data.frame(data)
  }

  # transform scores to a vector form
  if(is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # copy scores
  score2 <- score

  # factorize the response values
  resp <- purrr::map2(.x=data, .y=cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))

  # calculate the score categories for each examinee
  freq.cat <- purrr::map(.x=resp, .f=function(k) stats::xtabs(~ score + k, addNA = FALSE))

  # delete 'resp' object
  rm(resp, envir=environment(), inherits = FALSE)

  # extract theta (score) values for each item
  score <- purrr::map(.x=freq.cat,
                      .f=function(k) attributes(k)$dimnames$score %>%
                        as.numeric())

  # transform the score category data.frame to a matrix
  freq.cat <- purrr::map(.x=freq.cat,
                         .f=function(k) data.matrix(k) %>%
                           unname())

  ##---------------------------------------------------------------
  ## 2. compute loglikelihood
  ##---------------------------------------------------------------
  # create empty vectors to contain results
  llike <- rep(NA, nrow(x))

  # listrize the item meta data to use the starting values
  meta <- metalist2(x)

  # compute the loglikelihood (or likelihood)
  for(i in 1:nrow(x)) {

    # prepare information to estimate item parameters
    mod <- model[i]
    theta <- score[[i]]
    score.cat <- cats[i]

    # in case of a dichotomous item
    if(score.cat == 2) {
      # response data
      f_i <-
        freq.cat[[i]] %>%
        rowSums()
      r_i <- freq.cat[[i]][, 2]

      # item parameters
      pos_item <- which(meta$drm$loc == i)
      a.val <- meta$drm$a[pos_item]
      b.val <- meta$drm$b[pos_item]
      g.val <- meta$drm$g[pos_item]
      item_par <- c(a.val, b.val, g.val)

      # negative loglikelihood
      llike[i] <- loglike_drm(item_par=item_par, f_i=f_i, r_i=r_i, theta=theta, model=mod, D=D,
                              fix.a=FALSE, fix.g=FALSE, aprior=aprior, gprior=gprior,
                              use.aprior=use.aprior,
                              use.gprior=use.gprior)
    }

    # in case of a polytomous item
    if(score.cat > 2) {
      # response data
      r_i <- freq.cat[[i]]

      # check the starting values
      pos_item <- which(meta$plm$loc == i)
      a.val <- meta$plm$a[pos_item]
      d.val <- meta$plm$d[[pos_item]]
      item_par <- c(a.val, d.val)

      # negative loglikelihood
      llike[i] <- loglike_plm(item_par=item_par, r_i=r_i, theta=theta, pmodel=mod, D=D, fix.a=FALSE,
                              aprior=aprior, use.aprior=use.aprior)
    }

  }

  ##---------------------------------------------------------------
  # loglikelihood values
  llike <- -llike
  names(llike) <- x[, 1]

  # return results
  llike

}


# Negetive Loglikelihood of GPCM and GRM items
# @description This function computes the negative loglikelihood of an item with the
# polytomous IRT model
# @param item_par A vector of item parameters. The first element is the item discrimination (or slope)
# parameter. From the second elements, all all parameters are item threshold (or step) parameters.
# @param r_i A matrix of the frequencies of score categories for thetas for an item.
# @param theta A vector of theta for an item.
# @param pmodel A vector of character strings specifying the polytomous model with which response data are simulated.
# For each polytomous model, "GRM" for the graded response model or "GPCM" for the (generalized) partial credit model can be
# specified.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal
# ogive function (if set to 1.7). Default is 1.
#
# @return A numeric value
loglike_plm <- function(item_par, r_i, theta, pmodel=c("GRM", "GPCM"), D=1, fix.a=FALSE, a.val=1,
                        aprior=list(dist="lnorm", params=c(1, 0.5)),
                        use.aprior=FALSE,
                        FUN.grad, FUN.hess) {

  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  ##-------------------------------------------------------------------------
  if(!fix.a) {
    # compute category probabilities for all thetas
    ps <- plm(theta, a=item_par[1], d=item_par[-1], D=D, pmodel=pmodel)

    # compute loglikelihood
    log_ps <- suppressWarnings(log(ps))

    # to prevent that log(p) and log(q) have -Inf values
    log_ps <- ifelse(is.nan(log_ps), log(1e-20), log_ps)
    log_ps <- ifelse(is.infinite(log_ps), log(1e-20), log_ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

    # when the slope parameter prior is used
    if(use.aprior) {
      aprior_call <- paste0("stats::d", aprior$dist, "(", item_par[1], ", ", aprior$params[1], ", ", aprior$params[2], ", log=TRUE)")
      ln.aprior <- eval(expr=parse(text=aprior_call), envir=environment())
      llike <- llike + ln.aprior
    }


  } else {
    # compute category probabilities for all thetas
    ps <- plm(theta, a=a.val, d=item_par, D=D, pmodel=pmodel)

    # compute loglikelihood
    log_ps <- log(ps)

    # to prevent that log(p) and log(q) have -Inf values
    log_ps <- ifelse(is.infinite(log_ps), log(1e-4), log_ps)

    # log-likelihood
    llike <- sum(r_i * log_ps)

  }

  # return negative loglikelihood
  - llike

}



# Negetive Loglikelihood of dichotomous item
# @description This function computes the negative loglikelihood of an item with the
# dichotomous IRT model
# @param item_par A vector of item parameters. The first element is the item discrimination (or slope)
# parameter, the second element is the item difficulty parameter, and the third element is the item guessing
# parameter. The third element is necessary only when the 3PLM is used.
# @param f_i A vector of the frequencies for thetas for an item.
# @param r_i A vector of the frequencies of score categories for thetas for an item.
# @param theta A vector of theta for an item.
# @param prior A list containing the information of piror distribution for item guessig parameters.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal
# ogive function (if set to 1.7). Default is 1.
# @param use.prior A logical value. If TRUE, a prior ditribution specified in the argument \code{prior} is used when
# estimating item parameters of the IRT 3PLM. Default is TRUE.
#
# @return A numeric value
loglike_drm <- function(item_par, f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                        fix.a=FALSE, fix.g=FALSE, a.val=1, g.val=.2, n.1PLM=NULL,
                        aprior=list(dist="lnorm", params=c(1, 0.5)),
                        gprior=list(dist="beta", params=c(5, 17)),
                        use.aprior=FALSE,
                        use.gprior=TRUE,
                        FUN.grad, FUN.hess) {

  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # compute loglikelihood
  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {

    # make a list of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # make a list of all arguments being used in drm function
    argus <- list(a=a, b=b, f_i=f_i, r_i=r_i, theta=theta)

    # compute the negative loglikelihood values for all 1PLM items
    llike <-
      purrr::pmap_dbl(argus, .f=llike_drm, g=0, D=D) %>%
      sum()

    # when the slope parameter prior is used
    if(use.aprior) {
      aprior_call <- paste0("stats::d", aprior$dist, "(", item_par[1], ", ", aprior$params[1], ", ", aprior$params[2], ", log=TRUE)")
      ln.aprior <- eval(expr=parse(text=aprior_call), envir=environment())
      llike <- llike + ln.aprior
    }

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {

    # sume of loglikelihood
    llike <- llike_drm(a=a.val, b=item_par, g=0, f_i=f_i, r_i=r_i, theta=theta, D=D)

  }

  # (3) 2PLM
  if(model == "2PLM") {

    # sum of loglikelihood
    llike <- llike_drm(a=item_par[1], b=item_par[2], g=0, f_i=f_i, r_i=r_i, theta=theta, D=D)

    # when the slope parameter prior is used
    if(use.aprior) {
      aprior_call <- paste0("stats::d", aprior$dist, "(", item_par[1], ", ", aprior$params[1], ", ", aprior$params[2], ", log=TRUE)")
      ln.aprior <- eval(expr=parse(text=aprior_call), envir=environment())
      llike <- llike + ln.aprior
    }

  }

  # (4) 3PLM
  if(!fix.g & model == "3PLM") {

    # sum of loglikelihood
    llike <- llike_drm(a=item_par[1], b=item_par[2], g=item_par[3], f_i=f_i, r_i=r_i, theta=theta, D=D)

    # when the slope parameter prior is used
    if(use.aprior) {
      aprior_call <- paste0("stats::d", aprior$dist, "(", item_par[1], ", ", aprior$params[1], ", ", aprior$params[2], ", log=TRUE)")
      ln.aprior <- eval(expr=parse(text=aprior_call), envir=environment())
      llike <- llike + ln.aprior
    }

    # when the guessing parameter prior is used
    if(use.gprior) {
      gprior_call <- paste0("stats::d", gprior$dist, "(", item_par[3], ", ", gprior$params[1], ", ", gprior$params[2], ", log=TRUE)")
      ln.gprior <- eval(expr=parse(text=gprior_call), envir=environment())
      llike <- llike + ln.gprior
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    # sum of loglikelihood
    llike <- llike_drm(a=item_par[1], b=item_par[2], g=g.val, f_i=f_i, r_i=r_i, theta=theta, D=D)

    # when the slope parameter prior is used
    if(use.aprior) {
      aprior_call <- paste0("stats::d", aprior$dist, "(", item_par[1], ", ", aprior$params[1], ", ", aprior$params[2], ", log=TRUE)")
      ln.aprior <- eval(expr=parse(text=aprior_call), envir=environment())
      llike <- llike + ln.aprior
    }

  }


  # return a negative loglikelihood value
  - llike

}


# compute a sum of the loglikelihood value for each dichotomous item
llike_drm <- function(a, b, g, f_i, r_i, theta, D=1) {

  # compute loglikelihood
  p <- drm(theta, a=a, b=b, g=g, D=D)
  q <- 1 - p
  log_p <- suppressWarnings(log(p))
  log_q <- suppressWarnings(log(q))

  # to prevent that log(p) and log(q) have -Inf values
  log_p <- ifelse(is.nan(log_p), log(1e-20), log_p)
  log_q <- ifelse(is.nan(log_q), log(1e-20), log_q)
  log_p <- ifelse(is.infinite(log_p), log(1e-20), log_p)
  log_q <- ifelse(is.infinite(log_q), log(1e-20), log_q)

  # log-likelihood
  L <- r_i * log_p + (f_i - r_i) * log_q

  # sum of loglikelihood
  llike <- sum(L)

  # return
  llike

}
