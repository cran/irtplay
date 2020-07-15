#' Loglikelihood of ability
#'
#' @description This function computes the loglikelihood of abilities for examinees given the item parameters and response data.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data.
#' This data.frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param data A matrix or vector containing examinees' response data for the items in the argument \code{x}. When a matrix is used, a row and column indicate
#' the examinees and items, respectively. When a vector is used, it should contains the item response data for an examinee.
#' @param theta A numeric vector of abilities of which loglikelihood values are computed.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param method A character string indicating a scoring method. Available methods are "MLE" for the maximum likelihood estimation,
#' "MLEF" for the maximum likelihood estimation with fences, "MAP" for the maximum a posteriori estimation. Default method is "MLE".
#' @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
#' These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution. Default is
#' c(0,1). Ignored if \code{method} is "MLE" or "MLEF".
#' @param fence.a A numeric value specifying the item slope parameter (i.e., \emph{a}-parameter) for the two imaginary items in MLEF. See below for details.
#' Default is 3.0.
#' @param fence.b A numeric vector of two components specifying the lower and upper fences of item difficulty parameters (i.e., \emph{b}-parameters)
#' for the two imaginary items, respectively, in MLEF. When \code{fence.b = NULL}, the lower and upper fences of item difficulty parameters were
#' automatically set. See below for details. Default is NULL.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#'
#' @details The loglikelihood function of ability for an examinee can be computed given the item parameters and the examinee's response data for the items.
#' For example, if you want to examine the loglikelihood functions of abilities for two examinees given the same test items specified in the argument \code{x},
#' then you should provide the item response data matrix with two rows in the argument \code{data} and a vector of ability points where the loglikelihood values
#' need to be computed in the argument \code{theta}. Or if you want to examine the loglikelihood function of ability for an examinee given the test items
#' specified in the argument \code{x}, then you should provide the item response data matrix with one row (or a vector of item response data) in the argument
#' \code{data} and a vector of ability points where the loglikelihood values need to be computed in the argument \code{theta}.
#'
#' @return A data.frame of loglikelihood values. A row indicates the ability value where the loglikelihood is computed and
#' a column represents a response pattern, respectively.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(10)
#' score <- rnorm(5, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' # set the ability values where the loglikelihood values are computed
#' theta <- seq(-3, 3, 0.5)
#'
#' # compute the loglikelihood values (When MLE method is used)
#' llike_score(x=x, data=data, theta=theta, D=1, method="MLE")
#'
#' @import purrr
#' @import dplyr
#' @export
#'
llike_score <- function(x, data, theta, D = 1, method = "MLE", norm.prior = c(0, 1), fence.a = 3.0, fence.b = NULL, missing = NA) {

  method <- toupper(method)

  # check if the data set is a vector of an examinee
  if(is.vector(data)) {
    data <- rbind(data)
  }

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

  # clear the item meta data set
  x <- back2df(metalist2(x))

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

  # create an empty list to contain the loglikelihood values
  llike_list <- vector('list', nrow(data))

  # compute the loglikelihood values for every examinee
  for(i in 1:nrow(data)) {

    # listrize the data.frame
    meta <- metalist2(x)

    # select a response vector for an examinee
    resp <- data[i, ]

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
        meta$plm <- purrr::map(.x=meta$plm, .f=function(x) x[!meta$plm$loc %in% loc.miss])
      }
    }

    # factorize the response values
    max.cats <- max(c(meta$drm$cats, meta$plm$cats))
    resp.f <- factor(resp, levels=(seq_len(max.cats) - 1))

    # calculate the score categories
    tmp.id <- 1:length(resp.f)
    if(length(resp.f) > 1L) {
      freq.cat <- data.matrix(stats::xtabs(~ tmp.id + resp.f, addNA=TRUE))
    } else {
      freq.cat <- rbind(stats::xtabs(~ tmp.id + resp.f, addNA=TRUE))
    }

    if(any(is.na(resp))) freq.cat <- freq.cat[, -ncol(freq.cat), drop=FALSE]
    if(!is.null(meta$drm)) {
      freq.cat_drm <- freq.cat[meta$drm$loc, 1:2, drop=FALSE]
    } else {
      freq.cat_drm <- NULL
    }
    if(!is.null(meta$plm)) {
      freq.cat_plm <- freq.cat[meta$plm$loc, , drop=FALSE]
    } else {
      freq.cat_plm <- NULL
    }
    freq.cat <- list(freq.cat_drm=freq.cat_drm, freq.cat_plm=freq.cat_plm)

    # computing the negative loglikelihood values
    llike_tmp <- ll_brute(theta=theta, meta=meta, freq.cat=freq.cat, method=method, D=D, norm.prior=norm.prior)

    llike_list[[i]] <- - llike_tmp

  }

  # transform the list to a data.frame
  llike <-
    dplyr::bind_cols(llike_list) %>%
    data.frame()

  # return results
  llike

}

# This function computes a negative log likelihood or likelihood value
# This function is used for scoring
loglike_score <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL),
                          method=c("MLE", "MAP", "MLEF"), D=1, norm.prior=c(0, 1), logL=TRUE) {


  method <- match.arg(method)

  # make the empty vector to contain probabilities
  prob <- c()

  # check the number of thetas
  nstd <- length(theta)

  if(nstd == 1L) {

    # when there are dichotomous items
    if(!is.null(meta$drm)) {

      # compute the probabilities
      prob.1 <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
      prob.0 <- 1 - prob.1
      freq.cat_drm <- freq.cat$freq.cat_drm
      prob.drm <- prob.0 * freq.cat_drm[, 1] + prob.1 * freq.cat_drm[, 2]
      prob <- c(prob, prob.drm)

    }

    # when there are polytomous items
    if(!is.null(meta$plm)) {

      # extract polytomous model info
      model <- meta$plm$model

      # compute the category probabilities of items
      prob.plm <- c()
      for(i in 1:length(model)) {
        r_i <- freq.cat$freq.cat_plm[i, 1:meta$plm$cats[i]]
        prob.tmp <- plm(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=model[i])
        prob.plm[i] <- sum(prob.tmp * r_i)
      }
      prob <- c(prob, prob.plm)

    }

    # negative log-likelihood
    if(logL) {

      # loglikelihood
      log_prob <- suppressWarnings(log(prob))

      # to prevent that log(prob) has -Inf values
      log_prob <- ifelse(is.infinite(log_prob), log(1e-20), log_prob)

      # sum of loglikelihood
      rst <- sum(log_prob)

      # return a negative loglikelihood value
      rst <- switch(method,
                    MLE = -rst,
                    MLEF = -rst,
                    MAP = -(rst + stats::dnorm(theta, mean=norm.prior[1], sd=norm.prior[2], log=TRUE)))

    } else {

      # likelihood
      rst <- prod(prob)

    }

  }

  if(nstd > 1L) {

    # when there are dichotomous items
    if(!is.null(meta$drm)) {

      # check the number of items
      nitem <- length(meta$drm$a)

      # compute the probabilities
      prob.1 <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
      if(nitem == 1L) {
        prob.1 <- rbind(prob.1)
      } else {
        prob.1 <- t(prob.1)
      }
      prob.0 <- 1 - prob.1
      freq.cat_drm <- freq.cat$freq.cat_drm
      prob.drm <- prob.0 * freq.cat_drm[, 1] + prob.1 * freq.cat_drm[, 2]

      prob <- rbind(prob, prob.drm)

    }

    # when there are polytomous items
    if(!is.null(meta$plm)) {

      # extract polytomous model info
      model <- meta$plm$model

      # compute the category probabilities of items
      prob.plm <- array(0, c(length(model), nstd))
      for(i in 1:length(model)) {
        r_i <- freq.cat$freq.cat_plm[i, 1:meta$plm$cats[i]]
        prob.tmp <- t(plm(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=model[i]))
        prob.plm[i, ] <- colSums(prob.tmp * r_i)
      }

      prob <- rbind(prob, prob.plm)

    }

    # likelihood
    rst <- apply(prob, 2, prod)

  }

  # return
  rst

}


# This function computes a negative log likelihood for scoring
# This function is used only to find an ability estimate in brute force way
# when MLE or MAP fails to find the ability estimates
ll_brute <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL),
                     method=c("MLE", "MAP", "MLEF"), D=1, norm.prior=c(0, 1)) {


  method <- match.arg(method)

  # the number of examinee
  nstd <- length(theta)

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # the number of dichotomous items
    nitem <- length(meta$drm$a)

    # compute the probabilities
    prob.1 <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    if(nstd == 1L) {
      prob.1 <- rbind(prob.1)
    }
    if(nitem == 1L & nstd > 1L) {
      prob.1 <- cbind(prob.1)
    }
    prob.0 <- 1 - prob.1
    freq.cat_drm <- freq.cat$freq.cat_drm
    prob.drm <-
      cbind(prob.0[, which(freq.cat_drm[, 1] == 1), drop=FALSE], prob.1[, which(freq.cat_drm[, 2] == 1), drop=FALSE])

  } else {
    prob.drm <- NULL
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # compute the category probabilities of items
    freq.cat_plm <- freq.cat$freq.cat_plm
    prob.plm <- matrix(0, nrow=nstd, ncol=length(model))
    for(i in 1:length(model)) {
      prob.tmp <- plm(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=model[i])
      if(nstd == 1L) {
        prob.tmp <- rbind(prob.tmp)
      }
      prob.plm[, i] <- prob.tmp[, which(freq.cat_plm[i, 1:(meta$plm$cats[i])] == 1)]
    }

  } else {
    prob.plm <- NULL
  }

  # cbind of probability matrices
  prob <- cbind(prob.drm, prob.plm)

  # loglikelihood
  log_prob <- suppressWarnings(log(prob))

  # to prevent that log(prob) has -Inf values
  log_prob <- ifelse(is.infinite(log_prob), log(1e-100), log_prob)

  # sum of loglikelihood
  logL <- rowSums(log_prob)

  # return a negative loglikelihood value
  switch(method,
         MLE = -logL,
         MLEF = -logL,
         MAP = -(logL + stats::dnorm(theta, mean=norm.prior[1], sd=norm.prior[2], log=TRUE))
  )

}
