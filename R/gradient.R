# This function analytically computes a gradient vector for scoring
grad_score <- function(theta, meta, freq.cat=list(freq.cat_drm=NULL, freq.cat_plm=NULL), method=c("MLE", "MAP", "MLEF"),
                       D=1, norm.prior=c(0, 1), logL=TRUE,
                       FUN.grad=list(drm=NULL, plm=NULL, prior=NULL),
                       FUN.hess=list(drm=NULL, plm=NULL, prior=NULL)) {


  method <- match.arg(method)

  # extract the equations
  eq_grad_drm <- FUN.grad$drm
  eq_grad_plm <- FUN.grad$plm

  # empty vectors
  grad <- c()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {
    par.1 <- meta$drm$a
    par.2 <- meta$drm$b
    par.3 <- meta$drm$g
    n.1 <- freq.cat$freq.cat_drm[, 2]
    n.2 <- freq.cat$freq.cat_drm[, 1]
    rst <- eq_grad_drm(par.1=par.1, par.2=par.2, par.3=par.3, n.1=n.1, n.2=n.2, theta=theta, D=D)
    grad <- sum(c(grad, attributes(rst)$gradient))
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    freq.cat_plm <- freq.cat$freq.cat_plm

    for(i in 1:length(meta$plm$cats)) {

      # create a list containing the arguments to be used in the equation function
      args.pars <- list()
      args.pars[[1]] <- meta$plm$a[i]
      args.pars <- c(args.pars, as.list(meta$plm$d[[i]])[1:(meta$plm$cats[i] - 1)])
      args.ns <- as.list(freq.cat_plm[i, ])[1:meta$plm$cats[i]]
      args.list <- c(args.pars, args.ns)
      args.list$theta <- theta
      args.list$D <- D
      args.list <- unname(args.list)

      # implement the equation function
      params_fun <- eq_grad_plm[[i]]
      rst <- do.call("params_fun", args.list, envir=environment())
      grad <- sum(c(grad, attributes(rst)$gradient[1]))

    }

  }


  # extract the gradient vector when MAP method is used
  if(method == "MAP") {

    # equation
    eq_grad_prior <- FUN.grad$prior

    # implement the equation
    rst.prior <- eq_grad_prior(theta)
    grad <- sum(c(grad, attributes(rst.prior)$gradient[1]))

  }

  # return results
  grad

}


# This function analytically computes a gradient vector of dichotomous item parameters
grad_item_drm <- function(item_par, f_i, r_i, theta, model=c("1PLM", "2PLM", "3PLM", "DRM"), D=1,
                          fix.a=FALSE, fix.g=TRUE, a.val=1, g.val=.2, n.1PLM=NULL,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          gprior=list(dist="beta", params=c(5, 17)),
                          use.aprior=FALSE,
                          use.gprior=TRUE,
                          FUN.grad, FUN.hess) {


  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {
    n.par <- length(item_par)

    args.pars <- vector('list', n.par)
    for(i in 1:n.par) {
      args.pars[[i]] <- item_par[i]
    }
    args.ns <- vector('list', 2 * n.1PLM)
    for(i in 1:n.1PLM) {
      args.ns[[(2*i - 1)]] <- r_i[[i]]
    }
    for(i in 1:n.1PLM) {
      args.ns[[(2*i)]] <- f_i[[i]] - r_i[[i]]
    }
    args.list <- c(args.pars, args.ns)
    args.theta <- vector('list', n.1PLM)
    for(i in 1:n.1PLM) {
      args.theta[[i]] <- theta[[i]]
    }
    args.list <- c(args.list, args.theta)
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      aprior_fun <- FUN.grad$aprior_fun
      rst.aprior <- do.call("aprior_fun", list(item_par[1]), envir=environment())

      # extract the gradient vector
      grad.aprior <- attributes(rst.aprior)$gradient
      grad.aprior[is.nan(grad.aprior)] <- 0
      grad <- grad + colSums(grad.aprior)
    }

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {
    n.par <- length(item_par)

    args.list <- list()
    for(i in 1:n.par) {
      args.list[[i]] <- item_par[i]
    }
    args.list$n.1 <- r_i
    args.list$n.2 <- f_i - r_i
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

  }

  # (3) 2PLM
  if(model == "2PLM") {
    n.par <- length(item_par)

    args.list <- list()
    for(i in 1:n.par) {
      args.list[[i]] <- item_par[i]
    }
    args.list$n.1 <- r_i
    args.list$n.2 <- f_i - r_i
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      aprior_fun <- FUN.grad$aprior_fun
      rst.aprior <- do.call("aprior_fun", list(item_par[1]), envir=environment())

      # extract the gradient vector
      grad.aprior <- attributes(rst.aprior)$gradient
      grad.aprior[is.nan(grad.aprior)] <- 0
      grad <- grad + colSums(grad.aprior)
    }

  }

  # (4) 3PLM
  if(!fix.g & model == "3PLM") {
    n.par <- length(item_par)

    args.list <- list()
    for(i in 1:n.par) {
      args.list[[i]] <- item_par[i]
    }
    args.list$n.1 <- r_i
    args.list$n.2 <- f_i - r_i
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      aprior_fun <- FUN.grad$aprior_fun
      rst.aprior <- do.call("aprior_fun", list(item_par[1]), envir=environment())

      # extract the gradient vector
      grad.aprior <- attributes(rst.aprior)$gradient
      grad.aprior[is.nan(grad.aprior)] <- 0
      grad <- grad + colSums(grad.aprior)
    }

    # extract the gradient vector when the guessing parameter prior is used
    if(use.gprior) {
      gprior_fun <- FUN.grad$gprior_fun
      rst.gprior <- do.call("gprior_fun", list(item_par[3]), envir=environment())

      # extract the gradient vector
      grad.gprior <- attributes(rst.gprior)$gradient
      grad.gprior[is.nan(grad.gprior)] <- 0
      grad <- grad + colSums(grad.gprior)
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {
    n.par <- length(item_par)

    args.list <- list()
    for(i in 1:n.par) {
      args.list[[i]] <- item_par[i]
    }
    args.list$n.1 <- r_i
    args.list$n.2 <- f_i - r_i
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      aprior_fun <- FUN.grad$aprior_fun
      rst.aprior <- do.call("aprior_fun", list(item_par[1]), envir=environment())

      # extract the gradient vector
      grad.aprior <- attributes(rst.aprior)$gradient
      grad.aprior[is.nan(grad.aprior)] <- 0
      grad <- grad + colSums(grad.aprior)
    }

  }

  # return results
  grad

}


# This function analytically computes a gradient vector of polytomous item parameters
grad_item_plm <- function(item_par, r_i, theta, pmodel, D=1, fix.a=FALSE, a.val=1,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          use.aprior=FALSE,
                          FUN.grad, FUN.hess) {


  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  # transform item parameters as numeric values
  item_par <- as.numeric(item_par)

  ##-------------------------------------------------------------------------
  # check the number of categories and parameters to be estimated
  if(!fix.a) {
    cats <- length(item_par)
    n.par <- length(item_par)

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for(i in 1:n.par) {
      args.pars[[i]] <- item_par[i]
    }
    args.ns <- list()
    for(i in 1:cats) {
      args.ns[[i]] <- r_i[, i]
    }
    args.list <- c(args.pars, args.ns)
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

    # extract the gradient vector when the slope parameter prior is used
    if(use.aprior) {
      aprior_fun <- FUN.grad$aprior_fun
      rst.aprior <- do.call("aprior_fun", list(item_par[1]), envir=environment())

      # extract the gradient vector
      grad.aprior <- attributes(rst.aprior)$gradient
      grad.aprior[is.nan(grad.aprior)] <- 0
      grad <- grad + colSums(grad.aprior)
    }

  } else {
    cats <- length(item_par) + 1
    n.par <- length(item_par)

    # create a list containing the arguments to be used in the equation function
    args.pars <- list()
    for(i in 1:n.par) {
      args.pars[[i]] <- item_par[i]
    }
    args.ns <- list()
    for(i in 1:cats) {
      args.ns[[i]] <- r_i[, i]
    }
    args.list <- c(args.pars, args.ns)
    args.list$theta <- theta
    args.list$D <- D

    ##-------------------------------------------------------------------------
    # implement the equation function
    params_fun <- FUN.grad$params_fun
    rst <- do.call("params_fun", args.list, envir=environment())

    # extract the gradient vector
    grad <- attributes(rst)$gradient
    grad[is.nan(grad)] <- 0
    grad <- colSums(grad)

  }

  # return results
  grad

}

