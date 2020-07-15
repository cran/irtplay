# This function analytically computes a gradient vector for scoring
grad_score <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL),
                       method=c("MLE", "MAP", "MLEF"), D=1, norm.prior=c(0, 1), logL=TRUE) {


  method <- match.arg(method)

  # empty vectors
  grad_drm <- 0
  grad_plm <- 0

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # assign item parameters
    a <- meta$drm$a
    b <- meta$drm$b
    g <- meta$drm$g

    # assign item response for each score category
    r_i <- freq.cat$freq.cat_drm[, 2]

    # compute the probabilities for each score category
    p <- drm(theta=theta, a=a, b=b, g=g, D=D)
    p <- ifelse(p == 0L, 1e-20, p)

    # compute the gradient
    grad_drm[1] <- - D * sum((a * (p - g) * (r_i - p)) / ((1 - g) * p))

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    for(i in 1:length(meta$plm$cats)) {

      # check a polytomous model for an item
      model <- meta$plm$model[i]

      # assign item response for each score category
      r_i <- freq.cat$freq.cat_plm[i, 1:meta$plm$cats[i]]

      if(model == "GRM") {

        # assign item parameters
        a <- meta$plm$a[i]
        d <- meta$plm$d[[i]]

        # check the number of step parameters
        m <- length(d)

        # calculate all the probabilities greater than equal to each threshold
        allPst <- drm(theta=theta, a=rep(a, m), b=d, g=0, D=D)
        allQst <- 1 - allPst

        # calculate category probabilities
        P <- c(1, allPst) - c(allPst, 0)
        P <- ifelse(P == 0L, 1e-20, P)

        # compute the component values to get gradients
        Da <- D * a
        pq_st <- allPst * allQst
        deriv_Pstth <- Da * pq_st
        deriv_Pth <- c(0, deriv_Pstth) - c(deriv_Pstth, 0)
        frac_rp <- r_i / P

        # compute the gradient
        grad_plm[i] <- - sum(frac_rp * deriv_Pth)

      }

      if(model == "GPCM") {

        # assign item parameters
        a <- meta$plm$a[i]

        # include zero for the step parameter of the first category
        d <- c(0, meta$plm$d[[i]])

        # check the number of step parameters
        m <- length(d) - 1

        # calculate category probabilities
        Da <- D * a
        z <- Da * (theta - d)
        numer <- exp(cumsum(z)) # numerator
        denom <- sum(numer) # denominator
        P <- numer / denom
        P <- ifelse(P == 0L, 1e-20, P)

        # compute the component values to get a gradient vector
        frac_rp <- r_i / P
        denom2 <- denom^2
        d1th_z <- Da * (1:(m+1))
        d1th_denom <- sum(numer * d1th_z)
        deriv_Pth <- (numer / denom2) * (d1th_z * denom - d1th_denom)

        # compute the gradients of a and bs parameters
        grad_plm[i] <- - sum(frac_rp * deriv_Pth)

      }

    }

    # sum of all gradients for polytomous models
    grad_plm <- sum(grad_plm)

  }

  # sum of all gradients across all items
  grad <- sum(grad_drm, grad_plm)

  # extract the gradient vector when MAP method is used
  if(method == "MAP") {

    # compute a gradient of prior distribution
    rst.prior <-
      logprior_deriv(val=theta, is.aprior=FALSE, D=NULL, dist="norm",
                     par.1=norm.prior[1], par.2=norm.prior[2])

    # extract the gradient
    grad.prior <- attributes(rst.prior)$gradient

    # add the gradient
    grad <- sum(grad, grad.prior)

  }

  # return results
  grad

}



# This function computes the gradients vectors of the negative log likelihood for an item across all individuals
# This function is used to compute the standard errors of item parameter estimates using the cross-product method.
grad_llike <- function(item_par, f_i, r_i, quadpt, model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), D=1,
                       fix.a.1pl=TRUE, fix.a.gpcm=FALSE, fix.g=FALSE, a.val.1pl=1, a.val.gpcm=1, g.val=.2, n.1PLM=NULL) {

  # for dichotomous models
  if(model %in% c("1PLM", "2PLM", "3PLM")) {

    if(!fix.a.1pl & model == "1PLM") {

      # 1PLM: when the item slope parameters are not constrained to be equal across all items
      # compute the gradient vectors
      grad <- grad_item_drm_se(item_par=item_par, f_i=f_i, r_i=r_i, theta=quadpt, model=model, D=D,
                               fix.a=fix.a.1pl, n.1PLM=n.1PLM)

    } else {

      # for all other dichotomous models
      # compute the gradient vectors
      grad <- grad_item_drm_se(item_par=item_par, f_i=f_i, r_i=r_i, theta=quadpt, model=model, D=D,
                               fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL)

    }

  } else {

    # for polytomous models
    # compute the gradient vectors
    grad <- grad_item_plm_se(item_par=item_par, r_i=r_i, theta=quadpt, pmodel=model, D=D, fix.a=fix.a.gpcm,
                             a.val=a.val.gpcm)

  }

  # return results
  grad

}


# This function analytically computes a gradient vector of dichotomous item parameters
grad_item_drm <- function(item_par, f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                          fix.a=FALSE, fix.g=TRUE, a.val=1, g.val=.2, n.1PLM=NULL,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          bprior=list(dist="norm", params=c(0.0, 1.0)),
                          gprior=list(dist="beta", params=c(5, 17)),
                          use.aprior=FALSE,
                          use.bprior=FALSE,
                          use.gprior=TRUE) {


  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    r_p <- r_i  -  f_i * p
    bmat <- matrix(b, nrow=nstd, ncol=n.1PLM, byrow=TRUE)
    theta_b <- theta - bmat

    # compute the gradients of a and bs parameters
    gr_a <- - D * sum(theta_b * r_p)
    gr_b <- D * a[1] * colSums(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {

    # assign a and b parameters
    a <- a.val
    b <- item_par

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    r_p <- r_i  -  f_i * p

    # compute the gradients of a and bs parameters
    gr_b <- - as.numeric(-D * a * sum(r_p))

    # combine all gradients into a vector
    grad <- gr_b

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior <- attributes(rst.bprior)$gradient
    }

  }

  # (3) 2PLM
  if(model == "2PLM") {

    # assign a and b parameters
    a <- item_par[1]
    b <- item_par[2]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    Da <- D * a
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b

    # compute the gradients of a and bs parameters
    gr_a <- - D * sum(theta_b * r_p)
    gr_b <- Da * sum(r_p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }

  }

  # (4) 3PLM
  if(!fix.g & model == "3PLM") {

    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=g, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- - u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)
    gr_g <- - (1 / g_1) * sum(r_p / p)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b, gr_g)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }

    # extract the gradient vector when the guessing parameter prior is used
    if(use.gprior) {
      rst.gprior <-
        logprior_deriv(val=g, is.aprior=FALSE, D=NULL, dist=gprior$dist,
                       par.1=gprior$params[1], par.2=gprior$params[2])

      # extract the gradient vector
      grad.prior[3] <- attributes(rst.gprior)$gradient
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    # assign a, b, g parameters
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=g, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and bs parameters
    gr_a <- - u2 * sum(theta_b * u1)
    gr_b <- u3 * sum(u1)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior[2] <- attributes(rst.bprior)$gradient
    }

  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad

}


# This function analytically computes a gradient vector of polytomous item parameters
grad_item_plm <- function(item_par, r_i, theta, pmodel, D=1, fix.a=FALSE, a.val=1,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          bprior=list(dist="norm", params=c(0.0, 1.0)),
                          use.aprior=FALSE,
                          use.bprior=FALSE) {


  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ##-------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if(pmodel == "GRM") {

    # assign a, b parameters
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # check the numbers of examinees
    nstd <- length(theta)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta=theta, a=rep(a, m), b=d, g=0, D=D)
    if(nstd == 1L) {
      allPst <- matrix(allPst, nrow=1)
    }
    allQst <- 1 - allPst[, ,drop=FALSE]

    # calculate category probabilities
    P <- cbind(1, allPst) - cbind(allPst, 0)
    P <- ifelse(P == 0L, 1e-20, P)

    # compute the component values to get gradients
    # theta_b <- purrr::map_dfc(.x=d, .f=function(x) (theta - x))
    dmat <- matrix(d, nrow=nstd, ncol=m, byrow=TRUE)
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    # w1 <- D * theta_b
    w1 <- D * (theta - dmat)
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- frac_rp[, -(m + 1), drop=FALSE] - frac_rp[, -1, drop=FALSE]

    # compute the gradients of a and bs parameters
    gr_a <- - sum(frac_rp * w3)
    gr_b <- - Da * colSums(pq_st * w4)

    # combine all gradients into a vector
    grad <- c(gr_a, gr_b)

    # create a prior gradient vector
    grad.prior <- rep(0, n.par)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])

      # extract the gradient vector
      grad.prior[1] <- attributes(rst.aprior)$gradient
    }

    # extract the gradient vector when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=d, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])

      # extract the gradient vector
      grad.prior[-1] <- attributes(rst.bprior)$gradient
    }

  }

  # (2) GPCM and PCM
  if(pmodel == "GPCM") {

    if(!fix.a) {
      # For GPCM
      # assign a parameter
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # check the numbers of examinees and items
      nstd <- length(theta)

      # create a matrix for step parameters
      dmat <- matrix(d, nrow=nstd, ncol=(m+1), byrow=TRUE)

      # calculate category probabilities
      Da <- D * a
      z <- Da * (theta - dmat)
      numer <- exp(t(apply(z, 1, cumsum))) # numerator
      denom <- rowSums(numer) # denominator
      P <- numer / denom
      P <- ifelse(P == 0L, 1e-20, P)

      # compute the component values to get a gradient vector
      dsmat <- matrix(0, nrow=(m+1), ncol=(m+1))
      dsmat[-1, -1][upper.tri(x=dsmat[-1, -1], diag = TRUE)] <- 1
      tr_dsmat <- t(dsmat)
      frac_rp <- r_i / P
      w1 <- D * (theta - dmat)
      w1cum <- t(apply(w1, 1, cumsum))
      denom2 <- denom^2
      d1a_denom <- rowSums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- - (Da * numer) %*% tr_dsmat
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- - sum(frac_rp * deriv_Pa)
      gr_b <- c()
      for(k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (- numer * (DaDenom_vec %*% dsmat[k + 1, , drop=FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- c(gr_a, gr_b)

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the slope parameter prior is used
      if(use.aprior) {
        rst.aprior <-
          logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                         par.1=aprior$params[1], par.2=aprior$params[2])

        # extract the gradient vector
        grad.prior[1] <- attributes(rst.aprior)$gradient
      }

      # extract the gradient vector when the difficulty parameter prior is used
      if(use.bprior) {
        rst.bprior <-
          logprior_deriv(val=d[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                         par.1=bprior$params[1], par.2=bprior$params[2])

        # extract the gradient vector
        grad.prior[-1] <- attributes(rst.bprior)$gradient
      }

    } else {
      # for PCM
      # assign a parameter
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # check the numbers of examinees and items
      nstd <- length(theta)

      # create a matrix for step parameters
      dmat <- matrix(d, nrow=nstd, ncol=(m+1), byrow=TRUE)

      # calculate category probabilities
      Da <- D * a
      z <- Da * (theta - dmat)
      numer <- exp(t(apply(z, 1, cumsum))) # numerator
      denom <- rowSums(numer) # denominator
      P <- numer / denom
      P <- ifelse(P == 0L, 1e-20, P)

      # compute the component values to get a gradient vector
      dsmat <- matrix(0, nrow=(m+1), ncol=(m+1))
      dsmat[-1, -1][upper.tri(x=dsmat[-1, -1], diag = TRUE)] <- 1
      tr_dsmat <- t(dsmat)
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- - (Da * numer) %*% tr_dsmat
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- c()
      for(k in 1:m) {
        gr_b[k] <- -sum(frac_rp * (- numer * (DaDenom_vec %*% dsmat[k + 1, , drop=FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # combine all gradients into a vector
      grad <- gr_b

      # create a prior gradient vector
      grad.prior <- rep(0, n.par)

      # extract the gradient vector when the difficulty parameter prior is used
      if(use.bprior) {
        rst.bprior <-
          logprior_deriv(val=d[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                         par.1=bprior$params[1], par.2=bprior$params[2])

        # extract the gradient vector
        grad.prior <- attributes(rst.bprior)$gradient
      }

    }

  }

  # add the prior gradient vector
  grad <- grad + grad.prior

  # return results
  grad

}


# This function analytically computes a matrix of gradients of dichotomous item parameters across all examinees
# This function is used to compute the cross-product information matrix
grad_item_drm_se <- function(item_par, f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                             fix.a=FALSE, fix.g=TRUE, a.val=1, g.val=.2, n.1PLM=NULL) {


  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- rep(item_par[1], n.1PLM)
    b <- item_par[-1]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    r_p <- r_i  -  f_i * p
    bmat <- matrix(b, nrow=nstd, ncol=n.1PLM, byrow=TRUE)
    theta_b <- theta - bmat

    # compute the gradients of a and bs parameters
    gr_a <- - D * rowSums(theta_b * r_p)
    gr_b <- D * a[1] * r_p

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_a
    grad[, -1] <- gr_b

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {

    # assign a and b parameters
    a <- a.val
    b <- item_par

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    r_p <- r_i  -  f_i * p

    # compute the gradients of b parameter
    gr_b <- D * a * r_p

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_b

  }

  # (3) 2PLM
  if(model == "2PLM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=0, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    Da <- D * a
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b

    # compute the gradients of a and b parameters
    gr_a <- - D * theta_b * r_p
    gr_b <- Da * r_p

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_a
    grad[, 2] <- gr_b

  }

  # (4) 3PLM
  if(!fix.g & model == "3PLM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]
    g <- item_par[3]

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=g, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a, b, and g parameters
    gr_a <- - u2 * theta_b * u1
    gr_b <- u3 * u1
    gr_g <- - (1 / g_1) * r_p / p

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_a
    grad[, 2] <- gr_b
    grad[, 3] <- gr_g

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    b <- item_par[2]
    g <- g.val

    # check the numbers of examinees
    nstd <- length(theta)

    # compute the probabilities of correct and incorrect
    p <- drm(theta=theta, a=a, b=b, g=g, D=D)
    if(nstd == 1L) {
      p <- rbind(p)
    }
    p <- ifelse(p == 0L, 1e-20, p)
    q <- 1 - p

    # compute the component values
    g_1 <- 1 - g
    p_g <- p - g
    r_p <- r_i  -  f_i * p
    theta_b <- theta - b
    u1 <- p_g * r_p / p
    u2 <- D / g_1
    u3 <- a * u2

    # compute the gradients of a and b parameters
    gr_a <- - u2 * theta_b * u1
    gr_b <- u3 * u1

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_a
    grad[, 2] <- gr_b

  }

  # return results
  grad

}


# This function analytically computes a matrix of gradients of polytomous item parameters across all examinees
# This function is used to compute the cross-product information matrix
grad_item_plm_se <- function(item_par, r_i, theta, pmodel, D=1, fix.a=FALSE, a.val=1) {


  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # count the number of item parameters to be estimated
  n.par <- length(item_par)

  ##-------------------------------------------------------------------------
  # compute the gradients
  # (1) GRM
  if(pmodel == "GRM") {

    # make vectors of a and b parameters for all 1PLM items
    a <- item_par[1]
    d <- item_par[-1]

    # check the number of step parameters
    m <- length(d)

    # check the numbers of examinees
    nstd <- length(theta)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta=theta, a=rep(a, m), b=d, g=0, D=D)
    if(nstd == 1L) {
      allPst <- matrix(allPst, nrow=1)
    }
    allQst <- 1 - allPst[, ,drop=FALSE]

    # calculate category probabilities
    P <- cbind(1, allPst) - cbind(allPst, 0)
    P <- ifelse(P == 0L, 1e-20, P)

    # compute the component values to get gradients
    # theta_b <- purrr::map_dfc(.x=d, .f=function(x) (theta - x))
    dmat <- matrix(d, nrow=nstd, ncol=m, byrow=TRUE)
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    # w1 <- D * theta_b
    w1 <- D * (theta - dmat)
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- frac_rp[, -(m + 1), drop=FALSE] - frac_rp[, -1, drop=FALSE]

    # compute the gradients of a and bs parameters
    gr_a <- - rowSums(frac_rp * w3)
    gr_b <- - Da * pq_st * w4

    # create a matrix of gradients across all examinees
    grad <- matrix(0, nrow=nstd, ncol=n.par)
    grad[, 1] <- gr_a
    grad[, -1] <- gr_b

  }

  # (2) GPCM and PCM
  if(pmodel == "GPCM") {

    if(!fix.a) {
      # For GPCM
      # count the number of item parameters to be estimated
      n.par <- length(item_par)

      # make vectors of a and b parameters for all 1PLM items
      a <- item_par[1]

      # include zero for the step parameter of the first category
      d <- c(0, item_par[-1])

      # check the number of step parameters
      m <- length(d) - 1

      # check the numbers of examinees and items
      nstd <- length(theta)

      # create a matrix for step parameters
      dmat <- matrix(d, nrow=nstd, ncol=(m+1), byrow=TRUE)

      # calculate category probabilities
      Da <- D * a
      z <- Da * (theta - dmat)
      numer <- exp(t(apply(z, 1, cumsum))) # numerator
      denom <- rowSums(numer) # denominator
      P <- numer / denom
      P <- ifelse(P == 0L, 1e-20, P)

      # compute the component values to get a gradient vector
      dsmat <- matrix(0, nrow=(m+1), ncol=(m+1))
      dsmat[-1, -1][upper.tri(x=dsmat[-1, -1], diag = TRUE)] <- 1
      tr_dsmat <- t(dsmat)
      frac_rp <- r_i / P
      w1 <- D * (theta - dmat)
      w1cum <- t(apply(w1, 1, cumsum))
      denom2 <- denom^2
      d1a_denom <- rowSums(numer * w1cum)
      deriv_Pa <- (numer * (w1cum * denom - d1a_denom)) / denom2
      d1b_denom <- - (Da * numer) %*% tr_dsmat
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_a <- - rowSums(frac_rp * deriv_Pa)
      gr_b <- matrix(0, nrow=nstd, ncol=m)
      for(k in 1:m) {
        gr_b[, k] <- - rowSums(frac_rp * (- numer * (DaDenom_vec %*% dsmat[k + 1, , drop=FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- matrix(0, nrow=nstd, ncol=n.par)
      grad[, 1] <- gr_a
      grad[, -1] <- gr_b


    } else {
      # for PCM
      # make vectors of a and b parameters for all 1PLM items
      a <- a.val

      # include zero for the step parameter of the first category
      d <- c(0, item_par)

      # check the number of step parameters
      m <- length(d) - 1

      # check the numbers of examinees and items
      nstd <- length(theta)

      # create a matrix for step parameters
      dmat <- matrix(d, nrow=nstd, ncol=(m+1), byrow=TRUE)

      # calculate category probabilities
      Da <- D * a
      z <- Da * (theta - dmat)
      numer <- exp(t(apply(z, 1, cumsum))) # numerator
      denom <- rowSums(numer) # denominator
      P <- numer / denom
      P <- ifelse(P == 0L, 1e-20, P)

      # compute the component values to get a gradient vector
      dsmat <- matrix(0, nrow=(m+1), ncol=(m+1))
      dsmat[-1, -1][upper.tri(x=dsmat[-1, -1], diag = TRUE)] <- 1
      tr_dsmat <- t(dsmat)
      frac_rp <- r_i / P
      denom2 <- denom^2
      d1b_denom <- - (Da * numer) %*% tr_dsmat
      denom_vec <- cbind(denom)
      DaDenom_vec <- Da * denom_vec

      # compute the gradients of a and bs parameters
      gr_b <- matrix(0, nrow=nstd, ncol=m)
      for(k in 1:m) {
        gr_b[, k] <- - rowSums(frac_rp * (- numer * (DaDenom_vec %*% dsmat[k + 1, , drop=FALSE] + d1b_denom[, k + 1]) / denom2))
      }

      # create a matrix of gradients across all examinees
      grad <- gr_b

    }

  }

  # return results
  grad

}

