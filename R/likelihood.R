# This function computest a matrix of loglikelihoods for obtaining the observed item responses across examinees given theta values
# This function is used for MMLE-EM algorithm
likelihood <- function(meta, data1_drm=NULL, data2_drm=NULL, data_plm=NULL, theta, D=1) {

  # compute the log-likelihoods given the quadrature points
  # when there are dichotomous items
  if(!is.null(data1_drm)) {

    # compute the probabilities of answerting correctly on items
    ps <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    if(length(theta) == 1L) ps <- matrix(ps)

    # compute the probabilities of answerting incorrectly on items
    qs <- 1 - ps

    # to prevent that log(ps) and log(qs) have -Inf values
    log_ps <- suppressWarnings(log(ps))
    log_qs <- suppressWarnings(log(qs))
    log_ps <- ifelse(is.nan(log_ps), log(1e-20), log_ps)
    log_qs <- ifelse(is.nan(log_qs), log(1e-20), log_qs)
    log_ps <- ifelse(is.infinite(log_ps), log(1e-20), log_ps)
    log_qs <- ifelse(is.infinite(log_qs), log(1e-20), log_qs)

    # compte the likelihood values for all examinees at each quadrature point
    # a row indicate the examinee and column indicate the quad point
    llike_drm <- tcrossprod(x=data1_drm, y=log_ps) + tcrossprod(x=data2_drm, y=log_qs)
  } else {
    llike_drm <- 0L
  }

  # when there are polytomous items
  if(!is.null(data_plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <-
      purrr::pmap(.l=args, .f=plm, theta=theta, D=D) %>%
      do.call(what='cbind')

    # to prevent that log(prob.plm) have -Inf values
    log_prob.plm <- log(prob.plm)
    log_prob.plm <- ifelse(is.infinite(log_prob.plm), log(1e-20), log_prob.plm)

    # compte the likelihood values for all examinees at each quadrature point
    # a row indicate the examinee and column indicate the quad point
    llike_plm <- tcrossprod(x=data_plm, y=log_prob.plm)
  } else {
    llike_plm <- 0L
  }

  # sum of log-likelihood matrix
  LL <- llike_drm + llike_plm

  # transform to likelihood matrix
  L <- exp(LL)

  # return results
  rst <- list(L=L, LL=LL)
  rst


}

