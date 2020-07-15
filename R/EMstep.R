# E-step function
Estep <- function(x, data1_drm, data2_drm, data_plm, freq.cat, weights, D=1) {

  # listrize the data.frame
  meta <- metalist2(x)

  # compute the likelihood and log-likelihood matrix
  L_LL <- likelihood(meta, data1_drm=data1_drm, data2_drm=data2_drm, data_plm=data_plm, theta=weights[, 1], D=D)
  likehd <- L_LL$L
  loglikehd <- L_LL$LL

  # posterior distribution
  post_dist <- posterior(likehd, weights)

  # compute the expected frequency of scores categories across all items
  # this is the conditional expectation of item responses with respect to posterior likelihood distribution
  post_tp <- t(post_dist)
  freq.exp <- purrr::map(.x=freq.cat, .f=function(x) post_tp %*% x)

  # return results
  rst <- list(meta=meta, post_dist=post_dist, freq.exp=freq.exp, loglikehd=loglikehd, likehd=likehd)
  rst

}

# E-step function when FIPC method is used
Estep_fipc <- function(x1, x2,
                       data1_drm1, data2_drm1, data_plm1,
                       data1_drm2, data2_drm2, data_plm2,
                       freq.cat, weights, D=1) {

  # listrize the data.frame
  meta1 <- metalist2(x1) # meta information of the new items
  meta2 <- metalist2(x2) # meta information of the fixed items or all items

  # compute the likelihood and log-likelihood matrix of the new items to be estimates
  L_LL1 <- likelihood(meta1, data1_drm=data1_drm1, data2_drm=data2_drm1, data_plm=data_plm1, theta=weights[, 1], D=D)
  likehd1 <- L_LL1$L
  loglikehd1 <- L_LL1$LL

  # compute the likelihood and log-likelihood matrix of the fixed items (in the first iteration of EM) or all items (in the rest of the iteration of EM)
  # this (log) likelihood matrix is used only for computing the posterior density of ability
  L_LL2 <- likelihood(meta2, data1_drm=data1_drm2, data2_drm=data2_drm2, data_plm=data_plm2, theta=weights[, 1], D=D)
  likehd2 <- L_LL2$L

  # posterior distribution
  post_dist <- posterior(likehd2, weights)

  # compute the expected frequency of scores categories across all items
  # this is the conditional expectation of item responses with respect to posterior likelihood distribution
  post_tp <- t(post_dist)
  freq.exp <- purrr::map(.x=freq.cat, .f=function(x) post_tp %*% x)

  # return results
  rst <- list(meta=meta1, post_dist=post_dist, freq.exp=freq.exp, loglikehd=loglikehd1, likehd=likehd1)
  rst

}

# M-step function
#' @importFrom Matrix bdiag
Mstep <- function(estep, id, cats, model, quadpt, equation, D=1, loc_1p_const, loc_else, EmpHist, weights,
                  fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val,
                  use.aprior, use.bprior, use.gprior, aprior, bprior, gprior, group.mean, group.var, nstd,
                  Quadrature, use.startval, control, iter=NULL, fipc=FALSE, reloc.par, info.mstep=FALSE) {

  # extract the results of E-step
  meta <- estep$meta
  post_dist <- estep$post_dist
  freq.exp <- estep$freq.exp
  loglikehd <- estep$loglikehd
  likehd <- estep$likehd

  ##----------------------------------------------------------------------
  # (1) item parameter estimation
  ##----------------------------------------------------------------------
  # create empty vectors to contain results
  est_par <- NULL
  est_pure <- NULL
  convergence <- NULL
  noconv_items <- NULL
  info_list <- NULL
  se <- NULL

  # the dichotomous items: 1PLM with constrained slope values
  if("1PLM" %in% model & !fix.a.1pl) {

    # check the number of 1PLM items
    n.1PLM <- length(loc_1p_const)

    # prepare input files to estimate the 1PLM item parameters
    f_i <- r_i <- array(0, c(length(quadpt), n.1PLM))
    for(k in 1:n.1PLM) {
      f_i[, k] <- rowSums(freq.exp[loc_1p_const][[k]])
      r_i[, k] <- freq.exp[loc_1p_const][[k]][, 2]
    }

    # check the starting values
    if(use.startval) {
      pos_1p_const <- which(meta$drm$loc %in% loc_1p_const)
      a.val <- meta$drm$a[pos_1p_const][1]
      b.val <- meta$drm$b[pos_1p_const]
      startval <- c(a.val, b.val)
    } else {
      startval <- NULL
    }

    # item parameter estimation or compute the M step informaton matrix
    if(!info.mstep) {
      # parameter estimation
      est <- estimation2(f_i=f_i, r_i=r_i, quadpt=quadpt, model="1PLM", D=D, fix.a.1pl=FALSE, n.1PLM=n.1PLM,
                         aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior,
                         control=control, startval=startval, iter=iter)

      # extract the results
      # item parameter estimates
      a <- est$pars[1]
      b <- est$pars[-1]
      pars <- purrr::map(1:n.1PLM, .f=function(x) c(a, b[x], NA))
      est_par <- c(est_par, pars)
      est_pure <- c(est_pure, list(est$pars))

      # convergence indicator
      convergence <- c(convergence, est$convergence)
      if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_1p_const)

    } else {
      info_m <- info_mstep(item_par=startval, f_i=f_i, r_i=r_i, quadpt=quadpt, model="1PLM", D=D,
                           fix.a.1pl=FALSE, n.1PLM=n.1PLM)
      info_list <- c(info_list, list(info_m$info.mat))
      se <- c(se, info_m$se)

    }

  }

  # all other items
  if(length(loc_else) >= 1) {
    for(i in 1:length(loc_else)) {

      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if(score.cat == 2) {
        f_i <- rowSums(freq.exp[loc_else][[i]])
        r_i <- freq.exp[loc_else][[i]][, 2]

        # check the starting values
        if(use.startval) {
          pos_item <- which(meta$drm$loc == loc_else[i])
          a.val <- meta$drm$a[pos_item]
          b.val <- meta$drm$b[pos_item]
          g.val <- meta$drm$g[pos_item]
          if(mod == "1PLM") {
            startval <- b.val
          }
          if(mod == "2PLM") {
            startval <- c(a.val, b.val)
          }
          if(mod == "3PLM") {
            if(fix.g) {
              startval <- c(a.val, b.val)
            } else {
              startval <- c(a.val, b.val, g.val)
            }
          }
        } else {
          startval <- NULL
        }

        # item parameter estimation or compute the M step informaton matrix
        if(!info.mstep) {
          # parameter estimation
          est <- estimation2(f_i=f_i, r_i=r_i, quadpt=quadpt, model=mod, D=D,
                             fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE),
                             fix.g=fix.g, a.val.1pl=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                             aprior=aprior, bprior=bprior, gprior=gprior,
                             use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                             control=control, startval=startval, iter=iter)

          # extract the results
          # item parameter estimates
          a <- ifelse(mod == "1PLM", a.val.1pl, est$pars[1])
          b <- ifelse(mod == "1PLM", est$pars[1], est$pars[2])
          g <- ifelse(mod == "3PLM", ifelse(fix.g, g.val, est$pars[3]), NA)
          pars <- c(a, b, g)
          est_par <- c(est_par, list(pars))
          est_pure <- c(est_pure, list(est$pars))

          # convergence indicator
          convergence <- c(convergence, est$convergence)
          if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        } else {
          info_m <- info_mstep(item_par=startval, f_i=f_i, r_i=r_i, quadpt=quadpt, model=mod, D=D,
                               fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE), fix.g=fix.g, a.val.1pl=a.val.1pl,
                               g.val=g.val, n.1PLM=NULL)
          info_list <- c(info_list, list(info_m$info.mat))
          se <- c(se, info_m$se)
        }

      }

      # in case of a polytomous item
      if(score.cat > 2) {
        r_i <- freq.exp[loc_else][[i]]

        # check the starting values
        if(use.startval) {
          pos_item <- which(meta$plm$loc == loc_else[i])
          a.val <- meta$plm$a[pos_item]
          d.val <- meta$plm$d[[pos_item]]
          if(mod == "GRM") {
            startval <- c(a.val, d.val)
          }
          if(mod == "GPCM") {
            if(fix.a.gpcm) {
              startval <- d.val
            } else {
              startval <- c(a.val, d.val)
            }
          }
        } else {
          startval <- NULL
        }

        # item parameter estimation or compute the M step informaton matrix
        if(!info.mstep) {
          # parameter estimation
          est <- estimation2(r_i=r_i, quadpt=quadpt, model=mod, cats=score.cat, D=D,
                             fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE), a.val.gpcm=a.val.gpcm, n.1PLM=NULL,
                             aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior,
                             control=control, startval=startval, iter=iter)

          # extract the results
          # item parameter estimates
          a <- ifelse(mod == "GRM", est$pars[1], ifelse(fix.a.gpcm, a.val.gpcm, est$pars[1]))
          if(mod == "GRM") {
            bs <- est$pars[-1]
          } else {
            if(fix.a.gpcm) {
              bs <- est$pars
            } else{
              bs <- est$pars[-1]
            }
          }
          pars <- c(a, bs)
          est_par <- c(est_par, list(pars))
          est_pure <- c(est_pure, list(est$pars))

          # convergence indicator
          convergence <- c(convergence, est$convergence)
          if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        } else {
          info_m <- info_mstep(item_par=startval, r_i=r_i, quadpt=quadpt, model=mod, D=D,
                               fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE), a.val.gpcm=a.val.gpcm)
          info_list <- c(info_list, list(info_m$info.mat))
          se <- c(se, info_m$se)
        }

      }
    }
  }

  if(!info.mstep) {
    ##---------------------------------------------------------------
    # arrange the estimated item parameters and standard errors
    par_df <- data.frame(bind.fill(est_par, type="rbind"))
    par_df$loc <- c(loc_1p_const, loc_else)
    par_vec <- unlist(est_pure)[order(reloc.par)]
    # par_df <-
    #   par_df %>%
    #   dplyr::arrange(.data$loc) %>%
    #   dplyr::select(-.data$loc)
    par_df <- par_df[order(par_df$loc), ]
    par_df <- par_df[, -ncol(par_df)]

    # create a full data.frame for the item parameter estimates
    full_par_df <- data.frame(id=id, cats=cats, model=model, par_df, stringsAsFactors=FALSE)
    full_par_df$id <- as.character(full_par_df$id)
    colnames(full_par_df) <- c("id", "cats", "model", paste0("par.", 1:ncol(par_df)))

    # add par.3 column when there is no par.3 column (just in case that all items are 1PLMs and/or 2PLMs)
    if(ncol(full_par_df[, -c(1, 2, 3)]) == 2) {
      full_par_df <- data.frame(full_par_df, par.3=NA)
    }

    ##----------------------------------------------------------------------
    # (2) update the prior ability distribution
    ##----------------------------------------------------------------------
    if(EmpHist) {

      # update the prior frequencies
      prior_freq <- unname(colSums(post_dist))

      # prevent that the frequency has less than 1e-20
      prior_freq[prior_freq < 1e-20] <- 1e-20

      # normalize the updated prior freqency to obtain prior density function
      prior_dense <- unname(colSums(post_dist) / nstd)

      # update the prior densities
      if(fipc) {
        # when FIPC is used, no recaling is applied
        weights <- data.frame(theta=quadpt, weight=prior_dense)

      } else {
        # rescale the prior density distribution using the same quadrature point by applying Woods (2007) method
        prior_dense2 <- scale_prior(prior_freq=prior_freq, prior_dense=prior_dense, quadpt=quadpt,
                                    scale.par=c(group.mean, group.var), Quadrature=Quadrature)
        weights <- data.frame(theta=quadpt, weight=prior_dense2)

      }

    } else {

      if(fipc) {

        # update the prior frequencies
        prior_freq <- unname(colSums(post_dist))

        # prevent that the frequency has less than 1e-20
        prior_freq[prior_freq < 1e-20] <- 1e-20

        # normalize the updated prior freqency to obtain prior density function
        prior_dense <- unname(colSums(post_dist) / nstd)

        # compute the mean and sd of the updated prior distribution
        moments <- cal_moment(node=quadpt, weight=prior_dense)
        mu <- moments[1]
        sigma <- sqrt(moments[2])

        # obtain the updated prior densities from the normal distribution
        weights <- gen.weight(dist="norm", mu=mu, sigma=sigma, theta=quadpt)

      } else {
        weights <- weights
      }

    }

    ##---------------------------------------------------------------
    # compute the sum of loglikelihood values
    llike <- sum(log(likehd %*% matrix(weights[, 2])))

    # organize the the results
    rst <- list(par_df=full_par_df, par_vec=par_vec, convergence=convergence, noconv_items=noconv_items,
                weights=weights, loglike=llike)

  } else {

    # make a block diagonal information matrix
    info_mat <- as.matrix(Matrix::bdiag(info_list))

    # reorder the information matrix
    info_mat <- info_mat[, order(reloc.par)]
    info_mat <- info_mat[order(reloc.par), ]

    # reoerder the values of SEs
    se <- se[order(reloc.par)]

    # organize the the results
    rst <- list(info.mat=info_mat, se=se)

  }

  # return the results
  rst

}
