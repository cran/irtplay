hess_score <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL),
                       method=c("MLE", "MAP", "MLEF"), D=1, norm.prior=c(0, 1), logL=TRUE) {


  method <- match.arg(method)

  # empty vectors
  hess_drm <- 0
  hess_plm <- 0

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
    q <- 1 - p

    # compute the hessian
    hess_drm[1] <- - D^2 * sum(a^2 * (((p - g) * q) / ((1 - g)^2 * p)) * (((r_i * g) / p) - p))

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    for(i in 1:length(meta$plm$cats)) {

      # check a polytomous model for an item
      model <- meta$plm$model[i]

      # assign item response for each score category
      r_i <- freq.cat$freq.cat_plm[i, 1:meta$plm$cats[i]]

      if(model == "GRM") {

        # assign a, b parameters
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
        frac_rp2 <- frac_rp * (1 / P)
        w1 <- (1 - 2 * allPst) * deriv_Pstth
        deriv2_pth <- Da * (c(0, w1) - c(w1, 0))

        # compute the hessian matrix
        hess_thth1 <- - frac_rp2 * deriv_Pth^2
        hess_thth2 <- frac_rp * deriv2_pth
        hess_plm[i] <- - sum(hess_thth1 + hess_thth2)

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
        frac_rp2 <- frac_rp * (1 / P)
        denom2 <- denom^2
        denom4 <- denom2^2
        d1th_z <- Da * (1:(m+1))
        d1th_denom <- sum(numer * d1th_z)
        d2th_denom <- sum(numer * d1th_z^2)
        deriv_Pth <- (numer / denom2) * (d1th_z * denom - d1th_denom)
        d1th_z_den <- d1th_z * denom
        part1 <- d1th_z_den * denom * (d1th_z_den - d1th_denom)
        part2 <- denom2 * (d1th_z * d1th_denom + d2th_denom) - 2 * denom * d1th_denom^2
        deriv2_Pth <- (numer / denom4) * (part1 - part2)

        # compute the hessian matrix
        hess_thth1 <- - frac_rp2 * deriv_Pth^2
        hess_thth2 <- frac_rp * deriv2_Pth
        hess_plm[i] <- - sum(hess_thth1 + hess_thth2)

      }

    }

    # sum of all hessians for polytomous models
    hess_plm <- sum(hess_plm)

  }

  # sum of all hessians across all items
  hess <- sum(hess_drm, hess_plm)

  # extract the hessiam matrix when MAP method is used
  if(method == "MAP") {

    # compute a hessian of prior distribution
    rst.prior <-
      logprior_deriv(val=theta, is.aprior=FALSE, D=NULL, dist="norm",
                     par.1=norm.prior[1], par.2=norm.prior[2])

    # extract the hessian
    hess.prior <- attributes(rst.prior)$hessian

    # add the hessian
    hess <- sum(hess, hess.prior)

  }

  # return results
  hess <- data.matrix(hess)
  hess

}


# This function analytically computes a hessian matrix of dichotomous item parameters
# Also, adjust the hessian matrix if the matrix is singular by adding small random values
hess_item_drm <- function(item_par, f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                          fix.a=FALSE, fix.g=TRUE, a.val=1, g.val=.2, n.1PLM=NULL,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          bprior=list(dist="norm", params=c(0.0, 1.0)),
                          gprior=list(dist="beta", params=c(5, 17)),
                          use.aprior=FALSE, use.bprior=FALSE, use.gprior=TRUE,
                          adjust=TRUE) {


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
    bmat <- matrix(b, nrow=nstd, ncol=n.1PLM, byrow=TRUE)
    theta_b <- theta - bmat

    # compute the elements of hessian matrix of a and bs parameters
    D2 <- D^2
    pqf <- f_i * p * q
    hs_aa <- D2 * sum(theta_b^2 * pqf)
    hs_bb <- D2 * a[1]^2 * colSums(pqf)
    hs_ab <- D * colSums(r_i - f_i * p + (-D * a) * (theta_b) * pqf)

    # create a hessian matrix
    hess <- diag(c(hs_aa, hs_bb))
    hess[2:n.par, 1] <- hs_ab
    hess[1, 2:n.par] <- hs_ab

    # extract a prior hessian matrix
    hess.prior <- array(0, c(n.par, n.par))

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess.prior[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      diag(hess.prior)[-1] <- attributes(rst.bprior)$hessian
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

    # compute the elements of hessian matrix of a and bs parameters
    hs_bb <- D^2 * a^2 * sum(q * f_i * p)

    # create a hessian matrix
    hess <- matrix(hs_bb, nrow=1, ncol=1)

    # extract a prior hessian matrix
    hess.prior <- array(0, c(1, 1))

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess.prior[1, 1] <- attributes(rst.bprior)$hessian
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
    theta_b <- theta - b

    # compute the elements of hessian matrix of a and bs parameters
    D2 <- D^2
    pqf <- f_i * p * q
    hs_aa <- D2 * sum(theta_b^2 * pqf)
    hs_bb <- D2 * a^2 * sum(pqf)
    hs_ab <- D * sum(r_i - f_i * p + (-Da) * (theta_b) * pqf)

    # create a hessian matrix
    hess <- diag(c(hs_aa, hs_bb))
    hess[2:n.par, 1] <- hs_ab
    hess[1, 2:n.par] <- hs_ab

    # extract a prior hessian matrix
    hess.prior <- array(0, c(n.par, n.par))

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess.prior[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess.prior[2, 2] <- attributes(rst.bprior)$hessian
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
    theta_b <- theta - b
    u2 <- D / g_1
    u3 <- a * u2

    # compute the component values
    u4 <- r_i / p
    frac_qp <- q / p
    w1 <- u4 * g - f_i * p
    w2 <- u4 - f_i
    w3 <- frac_qp * u4
    w4 <- p_g * w3
    u5 <- p_g * frac_qp * w1

    # compute the elements of hessian matrix of a and bs parameters
    hs_aa <- - u2^2 * sum(theta_b^2 * u5)
    hs_bb <- - u3^2 * sum(u5)
    hs_gg <- - (1 / g_1^2) * sum(w2 - w3)
    hs_ab <- u2 * sum(p_g * (w2 + u3 * theta_b * frac_qp * w1))
    hs_ag <- (D / g_1^2) * sum(theta_b * w4)
    hs_bg <- - (D * a / g_1^2) * sum(w4)

    # create a hessian matrix
    hess <- diag(c(hs_aa, hs_bb, hs_gg))
    hess[2:n.par, 1] <- c(hs_ab, hs_ag)
    hess[1, 2:n.par] <- c(hs_ab, hs_ag)
    hess[3, 2] <- hs_bg
    hess[2, 3] <- hs_bg

    # extract a prior hessian matrix
    hess.prior <- array(0, c(n.par, n.par))

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess.prior[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess.prior[2, 2] <- attributes(rst.bprior)$hessian
    }

    # extract the hessian matrix when the guessing parameter prior is used
    if(use.gprior) {
      rst.gprior <-
        logprior_deriv(val=g, is.aprior=FALSE, D=NULL, dist=gprior$dist,
                       par.1=gprior$params[1], par.2=gprior$params[2])
      hess.prior[3, 3] <- attributes(rst.gprior)$hessian
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    # assign a and b parameters
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
    theta_b <- theta - b
    u2 <- D / g_1
    u3 <- a * u2

    # compute the component values
    u4 <- r_i / p
    frac_qp <- q / p
    w1 <- u4 * g - f_i * p
    w2 <- u4 - f_i
    u5 <- p_g * frac_qp * w1

    # compute the elements of hessian matrix of a and bs parameters
    hs_aa <- - u2^2 * sum(theta_b^2 * u5)
    hs_bb <- - u3^2 * sum(u5)
    hs_ab <- u2 * sum(p_g * (w2 + u3 * theta_b * frac_qp * w1))

    # create a hessian matrix
    hess <- diag(c(hs_aa, hs_bb))
    hess[1, 2] <- hs_ab
    hess[2, 1] <- hs_ab

    # extract a prior hessian matrix
    hess.prior <- array(0, c(n.par, n.par))

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess.prior[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=b, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess.prior[2, 2] <- attributes(rst.bprior)$hessian
    }

  }


  # add the prior hessian matrix
  hess <- hess + hess.prior

  # check if the hess is invertable
  if(adjust) {
    tmp <- suppressWarnings(tryCatch({solve(hess, tol=1e-200)}, error = function(e) {NULL}))
    if(is.null(tmp)) {
      item_par <- item_par + 0.05
      hess <- hess_item_drm(item_par=item_par, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                            fix.a=fix.a, fix.g=fix.g, a.val=a.val, g.val=g.val, n.1PLM=n.1PLM,
                            aprior=aprior, bprior=bprior, gprior=gprior, use.aprior=use.aprior,
                            use.bprior=use.bprior, use.gprior=use.gprior, adjust=TRUE)
    }
  }

  # return results
  hess

}



# This function analytically computes a hessian matrix of polytomous item parameters
# Also, adjust the hessian matrix if the matrix is singular by adding small random values
hess_item_plm <- function(item_par, r_i, theta, pmodel, D=1, fix.a=FALSE, a.val=1,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          bprior=list(dist="norm", params=c(0.0, 1.0)),
                          use.aprior=FALSE, use.bprior=FALSE, adjust=TRUE) {


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

    # make vectors of a and b parameters
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
    dmat <- matrix(d, nrow=nstd, ncol=m, byrow=TRUE)
    Da <- D * a
    pq_st <- allPst * allQst
    q_p_st <- allQst - allPst
    w1 <- D * (theta - dmat)
    w2 <- w1 * pq_st
    w3 <- cbind(0, w2) - cbind(w2, 0)
    frac_rp <- r_i / P
    w4 <- frac_rp[, -(m + 1), drop=FALSE] - frac_rp[, -1, drop=FALSE]

    # compute the more component values to get a hessian
    Da2 <- Da^2
    frac_rp2 <- frac_rp * (1 / P)
    z1 <- - frac_rp2 * w3^2
    w5 <- w1 * w2 * q_p_st
    z2 <- frac_rp * (cbind(0, w5) - cbind(w5, 0))
    z3 <- frac_rp[, -(m + 1), drop=FALSE] * (q_p_st +  pq_st / P[, -(m + 1), drop=FALSE])
    z4 <- frac_rp[, -1, drop=FALSE] * (q_p_st -  pq_st / P[, -1, drop=FALSE])
    z5 <- frac_rp2 * (P - a * w3)
    z6 <- z5[, -(m + 1), drop=FALSE] - z5[, -1, drop=FALSE]

    # compute the hessian matrix
    hess_aa <- - sum(z1 + z2)
    hess_bb <- colSums(Da2 * pq_st * (z3- z4))
    hess_ab <- - D * colSums(pq_st * (z6 + a * w1 * q_p_st * w4))
    hess <- diag(c(hess_aa, hess_bb))
    hess[1, 2:n.par] <- hess_ab
    hess[2:n.par, 1] <- hess_ab
    if(n.par > 2) {
      hess_b1b2 <-
        - Da2 * colSums(frac_rp2[, -c(1, (m + 1)), drop=FALSE] *
                          pq_st[, -m, drop=FALSE] * pq_st[, -1, drop=FALSE])
      diag(hess[2:m, 3:(m+1)]) <- hess_b1b2
      diag(hess[3:(m+1), 2:m]) <- hess_b1b2
    }

    # extract a prior hessian matrix
    hess.prior <- array(0, c(n.par, n.par))

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess.prior[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=d, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      diag(hess.prior)[-1] <- attributes(rst.bprior)$hessian
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

      # compute the more component values to get a hessian matrix
      frac_rp2 <- frac_rp * (1 / P)
      w1cum2 <- w1cum^2
      d2aa_denom <- rowSums(numer * w1cum2)
      denom4 <- denom2^2
      z1 <- (numer * denom / denom4)
      z2 <- w1cum2 * denom2
      deriv_aa1 <- - frac_rp2 * (deriv_Pa)^2
      deriv_Paa <- z1 * (z2 - 2 * w1cum * denom * d1a_denom - d2aa_denom * denom + 2 * d1a_denom^2)
      deriv_aa2 <- frac_rp * deriv_Paa
      d2bb_denom <- ((Da)^2 * numer) %*% tr_dsmat
      d2ab_denom <- - (D * numer * (1 + a * w1cum)) %*% tr_dsmat
      d1b_denom2 <- 2 * denom * d1b_denom
      Dvec <- array(D, c(nstd, 1))

      # compute the hessian matrix
      hess_aa <- - sum(deriv_aa1 + deriv_aa2)
      hess_bb <- c()
      DaDenom <- vector('list', m)
      deriv_Pb <- vector('list', m)
      for(k in 1:m) {
        DaDenom[[k]] <- DaDenom_vec %*% dsmat[k + 1, , drop=FALSE]
        deriv_Pb[[k]] <- - numer * (DaDenom[[k]] + d1b_denom[, k + 1]) / denom2
        deriv_bb1 <- - frac_rp2 * (deriv_Pb[[k]])^2
        part1 <- numer * denom2 * ((Da^2 * denom_vec) %*% dsmat[k + 1, , drop=FALSE] - d2bb_denom[, k + 1])
        part2 <- numer * (DaDenom[[k]] + d1b_denom[, k + 1]) * d1b_denom2[, k + 1]
        deriv_Pbb <- (part1 + part2) / denom4
        deriv_bb2 <- frac_rp * deriv_Pbb
        hess_bb[k] <- - sum(deriv_bb1 + deriv_bb2)
      }
      hess_ab <- c()
      d1b_w1cum <- vector('list', m)
      d1b_zcum <- vector('list', m)
      for(k in 1:m) {
        deriv_ab1 <- - frac_rp2 * deriv_Pb[[k]] * deriv_Pa
        d1b_w1cum[[k]] <- -Dvec %*% dsmat[k + 1, , drop=FALSE]
        d1b_zcum[[k]] <- -(Dvec * a) %*% dsmat[k + 1, , drop=FALSE]
        part1 <- denom2 * (d1b_w1cum[[k]] * denom + w1cum * (d1b_zcum[[k]] * denom - d1b_denom[, k+1]))
        part2 <- denom2 * (d1b_zcum[[k]] * d1a_denom + d2ab_denom[, k+1]) - 2 * denom * d1a_denom * d1b_denom[, k+1]
        deriv_Pab <- (numer / denom4) * (part1 - part2)
        deriv_ab2 <- frac_rp * deriv_Pab
        hess_ab[k] <- - sum(deriv_ab1 + deriv_ab2)
      }
      hess <- diag(c(hess_aa, hess_bb))
      hess[1, 2:n.par] <- hess_ab
      hess[2:n.par, 1] <- hess_ab
      if(n.par > 2) {
        hess_b1b2 <- c()
        for(k in 1:(m-1)) {
          for(j in (k+1):m) {
            deriv_b1b21 <- - frac_rp2 * deriv_Pb[[j]] * deriv_Pb[[k]]
            part1 <- denom2 * d1b_zcum[[k]] * (- DaDenom[[j]] - d1b_denom[, j + 1])
            part2 <-
              denom2 * (d1b_zcum[[j]] * d1b_denom[, k + 1] + d2bb_denom[, j + 1]) -
              2 * denom * d1b_denom[, k + 1] * d1b_denom[, j + 1]
            deriv_Pb1b2 <- (numer / denom4) * (part1 - part2)
            deriv_b1b22 <- frac_rp * deriv_Pb1b2
            tmp_val <- - sum(deriv_b1b21 + deriv_b1b22)
            hess_b1b2 <- c(hess_b1b2, tmp_val)
          }
        }
        hess[2:(m+1), 2:(m+1)][lower.tri(hess[2:(m+1), 2:(m+1)])] <- hess_b1b2
        hess_tr <- t(hess[2:(m+1), 2:(m+1)])
        hess[2:(m+1), 2:(m+1)][upper.tri(hess[2:(m+1), 2:(m+1)])] <- hess_tr[upper.tri(hess_tr)]
        # hess <-
        #   Matrix::forceSymmetric(x=hess, uplo="L") %>%
        #   as.matrix()
      }

      # extract a prior hessian matrix
      hess.prior <- array(0, c(n.par, n.par))

      # extract the hessian matrix when the slope parameter prior is used
      if(use.aprior) {
        rst.aprior <-
          logprior_deriv(val=a, is.aprior=TRUE, D=D, dist=aprior$dist,
                         par.1=aprior$params[1], par.2=aprior$params[2])
        hess.prior[1, 1] <- attributes(rst.aprior)$hessian
      }

      # extract the hessian matrix when the difficulty parameter prior is used
      if(use.bprior) {
        rst.bprior <-
          logprior_deriv(val=d[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                         par.1=bprior$params[1], par.2=bprior$params[2])
        diag(hess.prior)[-1] <- attributes(rst.bprior)$hessian
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

      # compute the more component values to get a hessian matrix
      frac_rp2 <- frac_rp * (1 / P)
      denom4 <- denom2^2
      d2bb_denom <- ((Da)^2 * numer) %*% tr_dsmat
      d1b_denom2 <- 2 * denom * d1b_denom
      Dvec <- array(D, c(nstd, 1))

      # compute the hessian matrix
      hess_bb <- c()
      DaDenom <- vector('list', m)
      deriv_Pb <- vector('list', m)
      for(k in 1:m) {
        DaDenom[[k]] <- DaDenom_vec %*% dsmat[k + 1, , drop=FALSE]
        deriv_Pb[[k]] <- - numer * (DaDenom[[k]] + d1b_denom[, k + 1]) / denom2
        deriv_bb1 <- - frac_rp2 * (deriv_Pb[[k]])^2
        part1 <- numer * denom2 * ((Da^2 * denom_vec) %*% dsmat[k + 1, , drop=FALSE] - d2bb_denom[, k + 1])
        part2 <- numer * (DaDenom[[k]] + d1b_denom[, k + 1]) * d1b_denom2[, k + 1]
        deriv_Pbb <- (part1 + part2) / denom4
        deriv_bb2 <- frac_rp * deriv_Pbb
        hess_bb[k] <- - sum(deriv_bb1 + deriv_bb2)
      }
      d1b_zcum <- vector('list', m)
      for(k in 1:m) {
        d1b_zcum[[k]] <- -(Dvec * a) %*% dsmat[k + 1, , drop=FALSE]
      }
      if(n.par > 1) {
        hess <- diag(hess_bb)
        hess_b1b2 <- c()
        for(k in 1:(m-1)) {
          for(j in (k+1):m) {
            deriv_b1b21 <- - frac_rp2 * deriv_Pb[[j]] * deriv_Pb[[k]]
            part1 <- denom2 * d1b_zcum[[k]] * (- DaDenom[[j]] - d1b_denom[, j + 1])
            part2 <-
              denom2 * (d1b_zcum[[j]] * d1b_denom[, k + 1] + d2bb_denom[, j + 1]) -
              2 * denom * d1b_denom[, k + 1] * d1b_denom[, j + 1]
            deriv_Pb1b2 <- (numer / denom4) * (part1 - part2)
            deriv_b1b22 <- frac_rp * deriv_Pb1b2
            tmp_val <- - sum(deriv_b1b21 + deriv_b1b22)
            hess_b1b2 <- c(hess_b1b2, tmp_val)
          }
        }
        hess[lower.tri(hess)] <- hess_b1b2
        hess_tr <- t(hess)
        hess[upper.tri(hess)] <- hess_tr[upper.tri(hess_tr)]
      } else {
        hess <- matrix(hess_bb, nrow=1, ncol=1)
      }

      # extract a prior hessian matrix
      hess.prior <- array(0, c(n.par, n.par))

      # extract the hessian matrix when the difficulty parameter prior is used
      if(use.bprior) {
        rst.bprior <-
          logprior_deriv(val=d[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                         par.1=bprior$params[1], par.2=bprior$params[2])
        diag(hess.prior) <- attributes(rst.bprior)$hessian
      }

    }

  }

  # add the prior hessian matrix
  hess <- hess + hess.prior

  # check if the hess is invertable
  if(adjust) {
    tmp <- suppressWarnings(tryCatch({solve(hess, tol=1e-200)}, error = function(e) {NULL}))
    if(is.null(tmp)) {
      item_par <- item_par + 0.05
      hess <- hess_item_plm(item_par=item_par, r_i=r_i, theta=theta, pmodel=pmodel, D=D, fix.a=fix.a, a.val=a.val,
                            aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior, adjust=TRUE)
    }
  }

  # return results
  hess

}



# This function computes the hessian matrix of the negative log likelihood of priors for an item
# This function is used to compute the cross-product information matrix
hess_prior <- function(item_par, model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), D=1,
                       fix.a.1pl=TRUE, fix.a.gpcm=FALSE, fix.g=FALSE,
                       aprior=list(dist="lnorm", params=c(1, 0.5)),
                       bprior=list(dist="norm", params=c(0.0, 1.0)),
                       gprior=list(dist="beta", params=c(5, 16)),
                       use.aprior=FALSE, use.bprior=FALSE, use.gprior=FALSE) {

  # for dichotomous models
  if(model %in% c("1PLM", "2PLM", "3PLM")) {

    # for all dichotomous models
    # compute the hessian matrix
    hess <- hess_prior_drm(item_par=item_par, model=model, D=D, fix.a=fix.a.1pl, fix.g=fix.g,
                           aprior=aprior, bprior=bprior, gprior=gprior, use.aprior=use.aprior,
                           use.bprior=use.bprior, use.gprior=use.gprior)

  } else {

    # for polytomous models
    # compute the hessian matrix
    hess <- hess_prior_plm(item_par=item_par, pmodel=model, D=D, fix.a=fix.a.gpcm,
                           aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior)

  }

  # return results
  hess

}


# This function computes the hessian matrix of the negative log likelihood of priors for an dichotomous item
# This function is used to compute the cross-product information matrix
hess_prior_drm <- function(item_par, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                           fix.a=FALSE, fix.g=TRUE,
                           aprior=list(dist="lnorm", params=c(1, 0.5)),
                           bprior=list(dist="norm", params=c(0.0, 1.0)),
                           gprior=list(dist="beta", params=c(5, 16)),
                           use.aprior=FALSE, use.bprior=FALSE, use.gprior=FALSE) {

  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)
  n.par <- length(item_par)

  # create an empty hessian matrix
  hess <- array(0, c(n.par, n.par))

  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=item_par[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      diag(hess)[-1] <- attributes(rst.bprior)$hessian
    }

  }

  # (2) 1PLM: the slope parameter is fixed
  if(fix.a & model == "1PLM") {

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess[1, 1] <- attributes(rst.bprior)$hessian
    }

  }

  # (3) 2PLM
  if(model == "2PLM") {

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=item_par[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par[2], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess[2, 2] <- attributes(rst.bprior)$hessian
    }

  }

  # (4) 3PLM: the guessing parameters are estimated
  if(!fix.g & model == "3PLM") {

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=item_par[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par[2], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess[2, 2] <- attributes(rst.bprior)$hessian
    }

    # extract the hessian matrix when the guessing parameter prior is used
    if(use.gprior) {
      rst.gprior <-
        logprior_deriv(val=item_par[3], is.aprior=FALSE, D=NULL, dist=gprior$dist,
                       par.1=gprior$params[1], par.2=gprior$params[2])
      hess[3, 3] <- attributes(rst.gprior)$hessian
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=item_par[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par[2], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      hess[2, 2] <- attributes(rst.bprior)$hessian
    }

  }

  # return results
  hess

}


# This function computes the hessian matrix of the negative log likelihood of priors for an polytomous item
# This function is used to compute the cross-product information matrix
hess_prior_plm <- function(item_par, pmodel, D=1, fix.a=FALSE,
                           aprior=list(dist="lnorm", params=c(1, 0.5)),
                           bprior=list(dist="norm", params=c(0.0, 1.0)),
                           use.aprior=FALSE, use.bprior=FALSE) {


  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)
  n.par <- length(item_par)

  # create an empty hessian matrix
  hess <- array(0, c(n.par, n.par))

  ##-------------------------------------------------------------------------
  # (1) GRM & GPCM
  if(!fix.a) {

    # extract the hessian matrix when the slope parameter prior is used
    if(use.aprior) {
      rst.aprior <-
        logprior_deriv(val=item_par[1], is.aprior=TRUE, D=D, dist=aprior$dist,
                       par.1=aprior$params[1], par.2=aprior$params[2])
      hess[1, 1] <- attributes(rst.aprior)$hessian
    }

    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par[-1], is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      diag(hess)[-1] <- attributes(rst.bprior)$hessian
    }

  } else {
    # (2) PCM
    # extract the hessian matrix when the difficulty parameter prior is used
    if(use.bprior) {
      rst.bprior <-
        logprior_deriv(val=item_par, is.aprior=FALSE, D=NULL, dist=bprior$dist,
                       par.1=bprior$params[1], par.2=bprior$params[2])
      diag(hess) <- attributes(rst.bprior)$hessian
    }

  }

  # return results
  hess

}


