# This function is an internal function used in the 'est_item' function.
estimation <- function(f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), cats, D=1,
                       fix.a.1pl=TRUE, fix.a.gpcm=FALSE, fix.g=FALSE, a.val.1pl=1, a.val.gpcm=1, g.val=.2, n.1PLM=NULL,
                       aprior=list(dist="lnorm", params=c(1, 0.5)),
                       gprior=list(dist="beta", params=c(5, 17)),
                       use.aprior=FALSE, use.gprior=TRUE,
                       control, startval=NULL) {


  # create an initial item parameter vector
  # and set the bouds of the item parameters
  if(model == "1PLM") {
    if(fix.a.1pl) {
      if(is.null(startval)) {item_par <- 0} else {item_par <- startval}
      lower <- -Inf
      upper <- Inf
    } else {
      if(is.null(startval)) {item_par <- c(1, rep(0, n.1PLM))} else {item_par <- startval}
      lower <- c(0, rep(-Inf, n.1PLM))
      upper <- c(Inf, rep(Inf, n.1PLM))
    }
  }
  if(model == "2PLM") {
    if(is.null(startval)) {item_par <- c(1, 0)} else {item_par <- startval}
    lower <- c(0.001, -Inf)
    upper <- c(Inf, Inf)
  }
  if(model == "3PLM") {
    if(fix.g) {
      if(is.null(startval)) {item_par <- c(1, 0)} else {item_par <- startval}
      lower <- c(0.001, -Inf)
      upper <- c(Inf, Inf)
    } else {
      if(is.null(startval)) {item_par <- c(1, 0, 0.2)} else {item_par <- startval}
      lower <- c(0.001, -Inf, 0)
      upper <- c(Inf, Inf, 1)
    }
  }
  if(model == "GRM") {
    if(is.null(startval)) {
      item_par <- c(1, seq(-1.0, 1.0, length.out=(cats-1)))
    } else {
      item_par <- startval
    }
    lower <- c(0.001, rep(-Inf, (cats-1)))
    upper <- c(Inf, rep(Inf, (cats-1)))
  }
  if(model == "GPCM") {
    if(fix.a.gpcm) {
      if(is.null(startval)) {
        item_par <- seq(-1.0, 1.0, length.out=(cats-1))
      } else {
        item_par <- startval
      }
      lower <- rep(-Inf, (cats-1))
      upper <- rep(Inf, (cats-1))
    } else {
      if(is.null(startval)) {
        item_par <- c(1, seq(-1.0, 1.0, length.out=(cats-1)))
      } else {
        item_par <- startval
      }
      lower <- c(0.001, rep(-Inf, (cats-1)))
      upper <- c(Inf, rep(Inf, (cats-1)))
    }
  }


  # estimation of item parameter & standard error
  if(!fix.a.1pl & model == "1PLM") {

    # when the item slope parameters are constrained to be equal across all 1PLM items
    # create the gradient vector and hessian matrix
    FUN.gh <- equation_drm(model=model, fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=n.1PLM,
                           aprior=aprior, gprior=gprior,
                           use.aprior=use.aprior,
                           use.gprior=use.gprior,
                           hessian=TRUE,
                           type="item")

    # estimate the item parameters
    est <- stats::nlminb(item_par, objective=loglike_drm, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                         fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=n.1PLM,
                         aprior=aprior, gprior=gprior,
                         use.aprior=use.aprior, use.gprior=use.gprior,
                         FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                         gradient=grad_item_drm,
                         # hessian=hess_item_drm,
                         control=control, lower=lower, upper=upper)

    # estimate the standard error of estimates
    hess <- hess_item_drm2(item_par=est$par, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                           fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=n.1PLM,
                           aprior=aprior, gprior=gprior,
                           use.aprior=use.aprior, use.gprior=use.gprior,
                           FUN.grad=FUN.gh, FUN.hess=FUN.gh)

  } else {

    # when the item slope parameters are not constrained to be across all items
    if(model %in% c("1PLM", "2PLM", "3PLM")) {

      # create the gradient vector and hessian matrix
      FUN.gh <- equation_drm(model=model, fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                             aprior=aprior, gprior=gprior,
                             use.aprior=use.aprior,
                             use.gprior=use.gprior,
                             hessian=TRUE,
                             type="item")

      # initial estimation to find better starting values
      tmp_est <- vector('list', 3)
      tmp_est[[1]] <- stats::nlminb(item_par, objective=loglike_drm, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                                    fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                                    aprior=aprior, gprior=gprior,
                                    use.aprior=use.aprior, use.gprior=use.gprior,
                                    FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                    gradient=grad_item_drm,
                                    # hessian=hess_item_drm,
                                    control=list(eval.max=50, iter.max=30, trace=0, step.min=1), lower=lower, upper=upper)
      tmp_est[[2]] <- stats::nlminb(item_par, objective=loglike_drm, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                                    fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                                    aprior=aprior, gprior=gprior,
                                    use.aprior=use.aprior, use.gprior=use.gprior,
                                    FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                    gradient=grad_item_drm,
                                    # hessian=hess_item_drm,
                                    control=list(eval.max=50, iter.max=30, trace=0, step.min=2), lower=lower, upper=upper)
      tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective))
      item_par <- tmp_est[[tmp_num]]$par

      # estimate the item parameters
      est <- tryCatch({stats::nlminb(item_par, objective=loglike_drm, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                                     fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                                     aprior=aprior, gprior=gprior,
                                     use.aprior=use.aprior, use.gprior=use.gprior,
                                     FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                     gradient=grad_item_drm,
                                     hessian=hess_item_drm,
                                     control=control, lower=lower, upper=upper)}, error = function(e) {NULL})
      if(is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(item_par, objective=loglike_drm, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                             fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                             aprior=aprior, gprior=gprior,
                             use.aprior=use.aprior, use.gprior=use.gprior,
                             FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                             gradient=grad_item_drm,
                             # hessian=hess_item_drm,
                             control=control, lower=lower, upper=upper)
      }

      # estimate the standard error of estimates
      hess <- hess_item_drm2(item_par=est$par, f_i=f_i, r_i=r_i, theta=theta, model=model, D=D,
                             fix.a=fix.a.1pl, fix.g=fix.g, a.val=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                             aprior=aprior, gprior=gprior,
                             use.aprior=use.aprior, use.gprior=use.gprior,
                             FUN.grad=FUN.gh, FUN.hess=FUN.gh)


    } else {

      # create the gradient vector and hessian matrix
      FUN.gh <- equation_plm(cats, pmodel=model, fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                             aprior=aprior, use.aprior=use.aprior,
                             hessian=TRUE,
                             type="item")

      # initial estimation to find better starting values
      tmp_est <- vector('list', 3)
      tmp_est[[1]] <- stats::nlminb(item_par, objective=loglike_plm, r_i=r_i, theta=theta, pmodel=model, D=D,
                                    fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                                    aprior=aprior, use.aprior=use.aprior,
                                    FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                    gradient=grad_item_plm,
                                    # hessian=hess_item_plm,
                                    control=list(eval.max=50, iter.max=30, step.min=1, trace=0), lower=lower, upper=upper)
      tmp_est[[2]] <- stats::nlminb(item_par, objective=loglike_plm, r_i=r_i, theta=theta, pmodel=model, D=D,
                                    fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                                    aprior=aprior, use.aprior=use.aprior,
                                    FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                    gradient=grad_item_plm,
                                    # hessian=hess_item_plm,
                                    control=list(eval.max=50, iter.max=30, step.min=2, trace=0), lower=lower, upper=upper)
      tmp_num <- which.min(c(tmp_est[[1]]$objective, tmp_est[[2]]$objective))
      item_par <- tmp_est[[tmp_num]]$par

      # estimate the item parameters
      est <- tryCatch({stats::nlminb(item_par, objective=loglike_plm, r_i=r_i, theta=theta, pmodel=model, D=D,
                                     fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                                     aprior=aprior, use.aprior=use.aprior,
                                     FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                                     gradient=grad_item_plm,
                                     hessian=hess_item_plm,
                                     control=control, lower=lower, upper=upper)}, error = function(e) {NULL})
      if(is.null(est) || est$convergence > 0L) {
        # if error or non-convergence, only use the gradient
        est <- stats::nlminb(item_par, objective=loglike_plm, r_i=r_i, theta=theta, pmodel=model, D=D,
                             fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                             aprior=aprior, use.aprior=use.aprior,
                             FUN.grad=FUN.gh, FUN.hess=FUN.gh,
                             gradient=grad_item_plm,
                             # hessian=hess_item_plm,
                             control=control, lower=lower, upper=upper)
      }

      # estimate the standard error of estimates
      hess <- hess_item_plm2(item_par=est$par, r_i=r_i, theta=theta, pmodel=model, D=D,
                             fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                             aprior=aprior, use.aprior=use.aprior,
                             FUN.grad=FUN.gh, FUN.hess=FUN.gh)

    }

  }

  # check if the hessian matrix can be inversed
  se <- suppressWarnings(tryCatch({sqrt(diag(solve(hess, tol=1e-200)))}, error = function(e) {NULL}))
  if(is.null(se)) {
    se <- rep(99999, length(diag(hess)))
  }
  if(any(is.nan(se))) {
    se[is.nan(se)] <- 99999
  }
  se <- ifelse(se > 99999, 99999, se)

  # return results
  rst <- list(pars=est$par, se=se, convergence=est$convergence, objective=est$objective)

  rst

}

