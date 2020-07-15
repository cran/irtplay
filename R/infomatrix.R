# This function computes the information matrix and SEs of the complete data of an item obtained from M step of EM algorithm (M step information matrix and SEs)
info_mstep <- function(item_par, f_i, r_i, quadpt, model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), D=1,
                       fix.a.1pl=TRUE, fix.a.gpcm=FALSE, fix.g=FALSE, a.val.1pl=1, a.val.gpcm=1, g.val=.2, n.1PLM=NULL) {

  if(model %in% c("1PLM", "2PLM", "3PLM")) {

    # (1) compute the negative hessian matrix of log likelihood for a dichotomous model
    hess <- hess_item_drm(item_par=item_par, f_i=f_i, r_i=r_i, theta=quadpt, model=model, D=D,
                          fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=n.1PLM,
                          use.aprior=FALSE, use.bprior=FALSE, use.gprior=FALSE)

  } else {

    # (2) compute the negative hessian matrix of log likelihood for a polytomous model
    hess <- hess_item_plm(item_par=item_par, r_i=r_i, theta=quadpt, pmodel=model, D=D,
                          fix.a=fix.a.gpcm, a.val=a.val.gpcm, use.aprior=FALSE, use.bprior=FALSE)

  }

  # check if the negative hessian matrix can be inversed
  cov_mat <- suppressWarnings(tryCatch({solve(hess, tol=1e-200)}, error = function(e) {NULL}))

  # compute the standard errors of item parameter estimates
  if(is.null(cov_mat)) {
    se <- rep(99999, length(diag(hess)))
  } else {
    se <- sqrt(diag(cov_mat))
  }

  # prevent showing NaN values of standard errors
  if(any(is.nan(se))) {
    se[is.nan(se)] <- 99999
  }

  # set an upper bound of standard error
  se <- ifelse(se > 99999, 99999, se)


  # return results
  rst <- list(info.mat=hess, cov.mat=cov_mat, se=se)

  rst

}


# This function computes the information matrix of the item parameter estimates using the cross-product method
info_xpd <- function(meta, freq.cat, post_dist, cats, model, quadpt, D=1, loc_1p_const, loc_else, nstd,
                     fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val, reloc.par) {

  # a create empty matrix to contain the kernel of fishre identify equation
  # across all item response patterns
  kernel_fisher <- 0L

  for(r in 1:length(quadpt)) {

    # a create empty matrix to contain the gradient matrix of joint log likelihood functions
    # across all item response patterns
    grad_mat <- NULL

    # the dichotomous items: 1PLM with constrained slope values
    if("1PLM" %in% model & !fix.a.1pl) {

      # check the number of 1PLM items
      n.1PLM <- length(loc_1p_const)

      # prepare input files to estimate the 1PLM item parameters
      # theta_val <- purrr::map(.x=1:n.1PLM, .f=function(x) rep(quadpt[r], nstd))
      # f_i <- purrr::map(.x=freq.cat[loc_1p_const], .f=function(k) rowSums(k))
      # r_i <- purrr::map(.x=freq.cat[loc_1p_const], .f=function(k) k[, 2])
      theta_val <- rep(quadpt[r], nstd)
      f_i <- r_i <- array(0, c(nstd, n.1PLM))
      for(k in 1:n.1PLM) {
        f_i[, k] <- rowSums(freq.cat[loc_1p_const][[k]])
        r_i[, k] <- freq.cat[loc_1p_const][[k]][, 2]
      }

      # use the final item parameter estimates
      pos_1p_const <- which(meta$drm$loc %in% loc_1p_const)
      a.val <- meta$drm$a[pos_1p_const][1]
      b.val <- meta$drm$b[pos_1p_const]
      item_par <- c(a.val, b.val)

      # compute the gradient vectors
      grad <-
        unname(grad_llike(item_par=item_par, f_i=f_i, r_i=r_i, quadpt=theta_val, model="1PLM", D=D,
                          fix.a.1pl=fix.a.1pl, n.1PLM=n.1PLM))

      # cbind the gradient matrix
      grad_mat <- cbind(grad_mat, grad)

    }

    # all other items
    if(length(loc_else) >= 1) {
      for(i in 1:length(loc_else)) {

        # prepare information to estimate item parameters
        mod <- model[loc_else][i]
        score.cat <- cats[loc_else][i]

        # in case of a dichotomous item
        if(score.cat == 2) {
          f_i <- rowSums(freq.cat[loc_else][[i]])
          r_i <- freq.cat[loc_else][[i]][, 2]

          # use the final item parameter estimates
          pos_item <- which(meta$drm$loc == loc_else[i])
          a.val <- meta$drm$a[pos_item]
          b.val <- meta$drm$b[pos_item]
          g.val <- meta$drm$g[pos_item]
          if(mod == "1PLM") {
            item_par <- b.val
          }
          if(mod == "2PLM") {
            item_par <- c(a.val, b.val)
          }
          if(mod == "3PLM") {
            if(fix.g) {
              item_par <- c(a.val, b.val)
            } else {
              item_par <- c(a.val, b.val, g.val)
            }
          }

          # compute the gradient vectors
          theta_val <- rep(quadpt[r], nstd)
          grad <-
            unname(grad_llike(item_par=item_par, f_i=f_i, r_i=r_i, quadpt=theta_val, model=mod, D=D,
                              fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE), fix.g=fix.g, a.val.1pl=a.val.1pl,
                              g.val=.2, n.1PLM=NULL))

          # cbind the gradient matrix
          grad_mat <- cbind(grad_mat, grad)

        }

        # in case of a polytomous item
        if(score.cat > 2) {
          r_i <- freq.cat[loc_else][[i]]

          # use the final item parameter estimates
          pos_item <- which(meta$plm$loc == loc_else[i])
          a.val <- meta$plm$a[pos_item]
          d.val <- meta$plm$d[[pos_item]]
          if(mod == "GRM") {
            item_par <- c(a.val, d.val)
          }
          if(mod == "GPCM") {
            if(fix.a.gpcm) {
              item_par <- d.val
            } else {
              item_par <- c(a.val, d.val)
            }
          }

          # compute the gradient vectors
          theta_val <- rep(quadpt[r], nstd)
          grad <-
            unname(grad_llike(item_par=item_par, f_i=f_i, r_i=r_i, quadpt=theta_val, model=mod, D=1,
                              fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE), a.val.gpcm=a.val.gpcm,
                              n.1PLM=NULL))

          # cbind the gradient matrix
          grad_mat <- cbind(grad_mat, grad)

        }
      }
    }

    # add the kernel_fisher matrixes
    kernel_fisher <- kernel_fisher + grad_mat * post_dist[, r]

  }

  # relocate the columns of matrix to locate the standard errors on the original position of item parameters
  kernel_fisher <- kernel_fisher[, order(reloc.par)]

  # compute the observed information matrix using the cross-product methods
  info_mat <- 0L
  for(i in 1:nstd) {
    info_mat <- info_mat + crossprod(x=rbind(kernel_fisher[i, ]))
  }

  # return the results
  info_mat

}

# This function computes the information matrix of priors using the second derivatives of item parameter estimates
#' @importFrom Matrix bdiag
info_prior <- function(meta, cats, model, D=1, loc_1p_const, loc_else, nstd,
                       fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val,
                       aprior, bprior, gprior, use.aprior, use.bprior, use.gprior, reloc.par) {

  # a create a vector containing the gradient values of priors
  # a create empty matrix to contain the gradient
  hess_list <- NULL

  # the dichotomous items: 1PLM with constrained slope values
  if("1PLM" %in% model & !fix.a.1pl) {

    # use the final item parameter estimates
    pos_1p_const <- which(meta$drm$loc %in% loc_1p_const)
    a.val <- meta$drm$a[pos_1p_const][1]
    b.val <- meta$drm$b[pos_1p_const]
    item_par <- c(a.val, b.val)

    # compute the hessian matrix
    hess <- hess_prior(item_par=item_par, model="1PLM", D=D, fix.a.1pl=fix.a.1pl,
                       aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior)
    hess_list <- c(hess_list, list(hess))
  }

  # all other items
  if(length(loc_else) >= 1) {
    for(i in 1:length(loc_else)) {

      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if(score.cat == 2) {

        # use the final item parameter estimates
        pos_item <- which(meta$drm$loc == loc_else[i])
        a.val <- meta$drm$a[pos_item]
        b.val <- meta$drm$b[pos_item]
        g.val <- meta$drm$g[pos_item]
        if(mod == "1PLM") {
          item_par <- b.val
        }
        if(mod == "2PLM") {
          item_par <- c(a.val, b.val)
        }
        if(mod == "3PLM") {
          if(fix.g) {
            item_par <- c(a.val, b.val)
          } else {
            item_par <- c(a.val, b.val, g.val)
          }
        }

        # compute the gradient vectors
        hess <- hess_prior(item_par=item_par, model=mod, D=D, fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE),
                           fix.g=fix.g, aprior=aprior, bprior=bprior, gprior=gprior,
                           use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior)
        hess_list <- c(hess_list, list(hess))

      }

      # in case of a polytomous item
      if(score.cat > 2) {

        # use the final item parameter estimates
        pos_item <- which(meta$plm$loc == loc_else[i])
        a.val <- meta$plm$a[pos_item]
        d.val <- meta$plm$d[[pos_item]]
        if(mod == "GRM") {
          item_par <- c(a.val, d.val)
        }
        if(mod == "GPCM") {
          if(fix.a.gpcm) {
            item_par <- d.val
          } else {
            item_par <- c(a.val, d.val)
          }
        }

        # compute the gradient vectors
        hess <- hess_prior(item_par=item_par, model=mod, D=D, fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE),
                           aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior)
        hess_list <- c(hess_list, list(hess))
      }
    }
  }

  # make a block diagonal matrix
  info_mat <- as.matrix(Matrix::bdiag(hess_list))

  # relocate the diagonal parts of the matrix to locate the standard errors on the original position of item parameters
  diag(info_mat) <- diag(info_mat)[order(reloc.par)]

  # return the results
  info_mat

}


