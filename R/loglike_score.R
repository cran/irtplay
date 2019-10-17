# This function computes a negative log likelihood or likelihood value
# This function is used for scoring
#' @import dplyr
loglike_score <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL), method=c("MLE", "MAP", "MLEF"),
                          D=1, norm.prior=c(0, 1), logL=TRUE,
                          FUN.grad=list(drm=NULL, plm=NULL, prior=NULL),
                          FUN.hess=list(drm=NULL, plm=NULL, prior=NULL)) {


  method <- match.arg(method)

  # make the empty vector to contain probabilities
  prob <- c()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities
    prob.1 <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    prob.0 <- 1 - prob.1
    freq.cat_drm <- freq.cat$freq.cat_drm
    prob.drm <-
      data.frame(prob.0 * freq.cat_drm[, 1], prob.1 * freq.cat_drm[, 2]) %>%
      rowSums()
    prob <- c(prob, prob.drm)

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    freq.cat_plm <- freq.cat$freq.cat_plm
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    prob.plm <- purrr::map_dbl(.x=1:length(meta$plm$cats),
                               .f=function(i)
                                 sum(prob.plm[[i]] * freq.cat_plm[i, 1:(meta$plm$cats[i])]))
    prob <- c(prob, prob.plm)

  }

  # negative log-likelihood
  if(logL) {

    # loglikelihood
    log_prob <- suppressWarnings(log(prob))

    # to prevent that log(prob) has -Inf values
    log_prob <- ifelse(is.infinite(log_prob), log(1e-20), log_prob)

    # sume of loglikelihood
    logL <- sum(log_prob)

    # return a negative loglikelihood value
    switch(method,
           MLE = -logL,
           MLEF = -logL,
           MAP = -(logL + stats::dnorm(theta, mean=norm.prior[1], sd=norm.prior[2], log=TRUE))
    )

  } else {

    # likelihood
    L <- prod(prob)

    L

  }

}


# This function computes a negative log likelihood for scoring
# This function is used only to find an ability estimate in brute force way
# when MLE or MAP fails to find the ability estimates
#' @import dplyr
#' @import purrr
ll_brute <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL),
                     method=c("MLE", "MAP", "MLEF"), D=1, norm.prior=c(0, 1)) {


  method <- match.arg(method)

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities
    prob.1 <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    prob.0 <- 1 - prob.1
    freq.cat_drm <- freq.cat$freq.cat_drm
    prob.drm <-
      cbind(prob.0[, which(freq.cat_drm[, 1] == 1)], prob.1[, which(freq.cat_drm[, 2] == 1)])
  } else {
    prob.drm <- NULL
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    freq.cat_plm <- freq.cat$freq.cat_plm
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    prob.plm <-
      purrr::map_dfc(.x=1:length(meta$plm$cats),
                     .f=function(i)
                       prob.plm[[i]][, which(freq.cat_plm[i, 1:(meta$plm$cats[i])] == 1)]) %>%
      data.matrix()
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

