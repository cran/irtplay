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

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {

    # sume of loglikelihood
    llike <- llike_drm(a=a.val, b=item_par, g=0, f_i=f_i, r_i=r_i, theta=theta, D=D)

  }

  # (3) 2PLM
  if(model == "2PLM") {

    # sume of loglikelihood
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

    # sume of loglikelihood
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

    # sume of loglikelihood
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


# compute a sume of the loglikelihood value for each dichotomous item
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
  L <- r_i * log_p +  (f_i - r_i) * log_q

  # sume of loglikelihood
  llike <- sum(L)

  # return
  llike

}
