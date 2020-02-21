# This function returns an equation of negative log of a prior distribution
equation_prior <- function(x, dist=c("norm", "lnorm", "beta"), params, pars, hessian=FALSE) {

  dist <- match.arg(dist)

  # select a distrubtion
  equation <-
    switch(dist,
           norm = paste0("(1 / (", params[2], " * sqrt(2 * pi))) * exp(-(", x, " - ", params[1], ")^2 / (2 * ", params[2], "^2))"),
           lnorm = paste0("(1 / (", x, " * ", params[2], " * sqrt(2 * pi))) * exp(-(log(", x, ") - ", params[1], ")^2 / (2 * ", params[2], "^2))"),
           beta = paste0("(gamma(", params[1], " + ", params[2], ") / (gamma(", params[1], ") * gamma(", params[2], "))) * ",
                         x, "^(", params[1], " - 1) * (1 - ", x, ")^(", params[2], " - 1)")
    )

  # take a log
  equation <- paste0("log(", equation, ")")

  # negarive loglikelihood
  equation <- paste0("-(", equation, ")")

  # return results
  stats::deriv(expr=parse(text=equation), namevec=pars, hessian=hessian,
               function.arg=x)

}

# This function returs an equation of loglikelihood, gradient, and hessian for dichotomous models (GRM and GPCM)
#' @import purrr
equation_drm <- function(model=c("1PLM", "2PLM", "3PLM", "DRM"), fix.a=FALSE, fix.g=FALSE, a.val=1, g.val=.2, n.1PLM=NULL,
                         aprior=list(dist="lnorm", params=c(1, 0.5)),
                         gprior=list(dist="beta", params=c(5, 17)),
                         pprior=list(dist="norm", params=c(0, 1)),
                         use.aprior=FALSE,
                         use.gprior=TRUE,
                         use.pprior=FALSE,
                         hessian=FALSE,
                         type=c("item", "ability"),
                         negative=TRUE) {


  if((model %in% c("2PLM", "3PLM", "DRM")) & fix.a) {
    stop("The slope parameter can't be fixed for 2PLM, 3PLM, and DRM.", call.=FALSE)
  }

  if(model == "1PLM" & !fix.a & is.null(n.1PLM)) {
    stop("The number of 1PLM items in which slope parameters are constrained to be equal must be specified.", call.=FALSE)
  }

  # consider DRM as 3PLM
  if(model == "DRM") model <- "3PLM"

  # create empty log prior equation functions
  aprior_fun <- NULL
  gprior_fun <- NULL

  # (1) 1PLM: the slope parameters are contrained to be equal across the 1PLM items
  if(!fix.a & model == "1PLM") {

    pars <- paste0("par.", 1:(n.1PLM + 1))
    ns <- purrr::map(.x=1:n.1PLM, .f=function(i) paste0("n.", i, ".", 1:2))
    thetas <- paste0("theta.", 1:n.1PLM)

    log.sum <- c()
    for(i in 1:n.1PLM) {
      log.1 <- paste0(ns[[i]][1], " * log(1 / (1 + exp(-D * ", pars[1], " * (theta.", i, " - ", pars[i + 1], "))) + 1e-20)")
      log.2 <- paste0(ns[[i]][2], " * log(1 - (1 / (1 + exp(-D * ", pars[1], " * (theta.", i, " - ", pars[i + 1], ")))) + 1e-20)")
      log.sum <- c(log.sum, paste(log.1, log.2, sep=" + "))
    }
    equation <- paste0("(", paste(log.sum, collapse = " + "), ")")
    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }

  }

  # (2) 1PLM: the slope parameters are fixed to be a specified value
  if(fix.a & model == "1PLM") {

    pars <- "par.1"
    ns <- c("n.1", "n.2")

    log.1 <- paste0(ns[1], " * log(1 / (1 + exp(-D * ", a.val, " * (theta - ", pars, "))) + 1e-20)")
    log.2 <- paste0(ns[2], " * log(1 - (1 / (1 + exp(-D * ", a.val, " * (theta - ", pars, ")))) + 1e-20)")
    equation <- paste(log.1, log.2, sep=" + ")

  }

  # (3) 2PLM
  if(model == "2PLM") {

    pars <- c("par.1", "par.2")
    ns <- c("n.1", "n.2")

    log.1 <- paste0(ns[1], " * log(1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], "))) + 1e-20)")
    log.2 <- paste0(ns[2], " * log(1 - (1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], ")))) + 1e-20)")
    equation <- paste(log.1, log.2, sep=" + ")
    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }

  }

  # (4) 3PLM
  if(!fix.g & model == "3PLM") {

    pars <- c("par.1", "par.2", "par.3")
    ns <- c("n.1", "n.2")

    log.1 <- paste0(ns[1], " * log(", pars[3], " + (1 - ", pars[3], ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], "))) + 1e-20)")
    log.2 <- paste0(ns[2], " * log(1 - (", pars[3], " + (1 - ", pars[3], ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], ")))) + 1e-20)")
    equation <- paste(log.1, log.2, sep=" + ")
    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }
    if(use.gprior) {
      gprior_fun <- equation_prior(x=pars[3], dist=gprior$dist, params=gprior$params, pars=pars, hessian=hessian)
    }

  }

  # (5) 3PLM: the guessing parameters are fixed to be specified value
  if(fix.g & model == "3PLM") {

    pars <- c("par.1", "par.2")
    ns <- c("n.1", "n.2")

    log.1 <- paste0(ns[1], " * log(", g.val, " + (1 - ", g.val, ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], "))) + 1e-20)")
    log.2 <- paste0(ns[2], " * log(1 - (", g.val, " + (1 - ", g.val, ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], ")))) + 1e-20)")
    equation <- paste(log.1, log.2, sep=" + ")
    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }

  }

  # negative loglikelihood
  if(negative) {
    equation <- paste0("-(", equation, ")")
  } else {
    equation <- paste0("(", equation, ")")
  }


  ##----------------------------------------------------------------------------
  # create a function for loglikelihood, gradient, hessian
  type <- match.arg(type)
  if(!fix.a & model == "1PLM") {
    params_fun <-
      switch(type,
             item = stats::deriv(expr=parse(text=equation), namevec=pars, hessian=hessian,
                                 function.arg=c(get("pars", inherits=FALSE), unlist(get("ns", inherits=FALSE)), get("thetas", inherits=FALSE), "D"))
      )
  } else {
    params_fun <-
      switch(type,
             item = stats::deriv(expr=parse(text=equation), namevec=pars, hessian=hessian,
                                 function.arg=c(get("pars", inherits=FALSE), get("ns", inherits=FALSE), "theta", "D")),
             ability = stats::deriv(expr=parse(text=equation), namevec=c("theta"), hessian=hessian,
                                    function.arg=c(get("pars", inherits=FALSE), get("ns", inherits=FALSE), "theta", "D"))
      )
  }

  # create a log prior function for a populationn
  if(use.pprior) {
    pprior_fun <- equation_prior(x="theta", dist=pprior$dist, params=pprior$params, pars="theta", hessian=TRUE)
  } else {
    pprior_fun <- NULL
  }

  ##----------------------------------------------------------------------------
  # return a function
  rst <-
    switch(type,
           item = list(params_fun=params_fun, aprior_fun=aprior_fun, gprior_fun=gprior_fun),
           ability = list(params_fun=params_fun, pprior_fun=pprior_fun)
    )

  rst

}



# This function returs an equation of loglikelihood, gradient, and hessian for polytomous models (GRM and GPCM)
equation_plm <- function(cats, pmodel=c("GRM", "GPCM"), fix.a=FALSE, a.val=1,
                         aprior=list(dist="lnorm", params=c(1, 0.5)),
                         pprior=list(dist="norm", params=c(0, 1)),
                         use.aprior=FALSE,
                         use.pprior=FALSE,
                         hessian=FALSE,
                         type=c("item", "ability"),
                         negative=TRUE) {

  if(pmodel == "GRM" & fix.a) {
    stop("The slope parameter can't be fixed for GRM.", call.=FALSE)
  }

  if(!fix.a) {
    pars <- paste0("par.", 1:cats)
    ns <-  paste0("n.", 1:cats)
  } else {
    pars <- paste0("par.", 1:(cats-1))
    ns <-  paste0("n.", 1:(cats))
  }

  # create empty log prior equation functions
  aprior_fun <- NULL

  ##----------------------------------------------------------------------------
  # loglikelihood equation
  # (1) GRM
  if(pmodel == "GRM") {
    equation <- c()
    for(i in 1:cats) {
      if(i == 1) {
        equation <- paste0(ns[1], " * log(1 - (1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], ")))) + 1e-20) + ")
      }
      if(i > 1 & i < cats) {

        equation <- paste0(equation,
                           paste0(ns[i], " * log((1 / (1 + exp(-D * ", pars[1],  " * (theta - ", pars[i], ")))) -
                         (1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[i + 1], ")))) + 1e-20) + "))

      }
      if(i == cats) {

        equation <- paste0(equation, paste0(ns[i], " * log(1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[i], "))) + 1e-20)"))

      }
    }
    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }

  }

  ##----------------------------------------------------------------------------
  # loglikelihood equation
  # (2) GPCM
  if(pmodel == "GPCM" & !fix.a) {

    # denominator
    denom <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - ", pars[k+1], "))"))
        }

      }
      denom <- c(denom, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }
    denom <- paste0("(", paste(denom, collapse = " + "), ")")

    # numerator
    numer <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - ", pars[k+1], "))"))
        }

      }
      numer <- c(numer, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }

    # loglikelihood equation
    llike_wts <- paste0(ns, " * log((", paste(numer, denom, sep=" / "), ") + 1e-20)")
    equation <- paste(llike_wts, collapse = " + ")

    # if priors are specified is used
    if(use.aprior) {
      aprior_fun <- equation_prior(x=pars[1], dist=aprior$dist, params=aprior$params, pars=pars, hessian=hessian)
    }

  }

  ##----------------------------------------------------------------------------
  # loglikelihood equation
  # (3) GPCM with a fixed slope value
  if(pmodel == "GPCM" & fix.a) {

    # denominator
    denom <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", a.val, " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", a.val, " * (theta - ", pars[k], "))"))
        }

      }
      denom <- c(denom, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }
    denom <- paste0("(", paste(denom, collapse = " + "), ")")

    # numerator
    numer <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", a.val, " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", a.val, " * (theta - ", pars[k], "))"))
        }

      }
      numer <- c(numer, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }

    # loglikelihood equation
    llike_wts <- paste0(ns, " * log((", paste(numer, denom, sep=" / "), ") + 1e-20)")
    equation <- paste(llike_wts, collapse = " + ")

  }

  # negative loglikelihood
  if(negative) {
    equation <- paste0("-(", equation, ")")
  } else {
    equation <- paste0("(", equation, ")")
  }

  ##----------------------------------------------------------------------------
  # create a function for loglikelihood, gradient, hessian
  type <- match.arg(type)
  params_fun <-
    switch(type,
           item = stats::deriv(expr=parse(text=equation), namevec=pars, hessian=hessian,
                               function.arg=c(get("pars", inherits=FALSE), get("ns", inherits=FALSE), "theta", "D")),
           ability = stats::deriv(expr=parse(text=equation), namevec=c("theta"), hessian=hessian,
                                  function.arg=c(get("pars", inherits=FALSE), get("ns", inherits=FALSE), "theta", "D"))
    )

  # create a log prior function for a populationn
  if(use.pprior) {
    pprior_fun <- equation_prior(x="theta", dist=pprior$dist, params=pprior$params, pars="theta", hessian=TRUE)
  } else {
    pprior_fun <- NULL
  }

  ##----------------------------------------------------------------------------
  # return a function
  rst <-
    switch(type,
           item = list(params_fun=params_fun, aprior_fun=aprior_fun),
           ability = list(params_fun=params_fun, pprior_fun=pprior_fun)
    )

  rst

}



# This function returns the gradient vector and hessian matrix for the score category probability equation
equation_scocat <- function(model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), cats=NULL, fix.a.gpcm=FALSE, hessian=TRUE, type=c("item", "ability")) {


  ##-------------------------------
  # set the item parameters to be used in the equation
  if(model %in% c("1PLM", "2PLM", "3PLM")) {
    pars <- c("par.1", "par.2", "par.3")
  }
  if(model %in% c("GRM", "GPCM")) {
    pars <- paste0("par.", 1:cats)
  }

  # set the number of score categories
  if(model %in% c("1PLM", "2PLM", "3PLM")) cats <- 2

  ##----------------------------------------------------------------------------
  # score category probability equation
  equation <- c()

  # (1) DRM
  if(model %in% c("1PLM", "2PLM", "3PLM")) {
    equation <- c(equation, paste0(pars[3], " + (1 - ", pars[3], ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], ")))"))
    equation <- c(equation, paste0("1 - (", pars[3], " + (1 - ", pars[3], ") / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], "))))"))
  }

  # (2) GRM
  if(model == "GRM") {
    for(i in 1:cats) {
      if(i == 1) {
        equation <- c(equation, paste0("1 - (1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[2], "))))"))
      }
      if(i > 1 & i < cats) {
        equation <- c(equation, paste0("(1 / (1 + exp(-D * ", pars[1],  " * (theta - ", pars[i], ")))) -
                                     (1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[i + 1], "))))"))
      }
      if(i == cats) {
        equation <- c(equation, paste0("1 / (1 + exp(-D * ", pars[1], " * (theta - ", pars[i], ")))"))
      }
    }
  }

  # (2) GPCM
  if(model == "GPCM") {

    # denominator
    denom <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - ", pars[k+1], "))"))
        }

      }
      denom <- c(denom, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }
    denom <- paste0("(", paste(denom, collapse = " + "), ")")

    # numerator
    numer <- c()
    for(i in 0:(cats-1)) {

      tmp <- c()
      for(k in 0:i) {

        if(k == 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - 0))"))
        }
        if(k > 0) {
          tmp <- c(tmp, paste0("(D * ", pars[1], " * (theta - ", pars[k+1], "))"))
        }

      }
      numer <- c(numer, paste0("exp(", paste(tmp, collapse = " + "), ")"))
    }

    # loglikelihood equation
    equation <- paste(numer, denom, sep=" / ")

  }

  ##----------------------------------------------------------------------------
  # create a function for a gradient and hessian
  funList <- vector('list', cats)
  if(type == "item") {

    # set the evaluated item parameters
    if(model == "1PLM") namevec <- pars[-c(1, 3)]
    if(model == "2PLM") namevec <- pars[-3]
    if(model == "3PLM") namevec <- pars
    if(model == "GRM") namevec <- pars
    if(model == "GPCM") {
      if(fix.a.gpcm) {
        namevec <- pars[-1]
      } else {
        namevec <- pars
      }
    }

    for(i in 1:cats) {
      funList[[i]] <- stats::deriv(expr=parse(text=equation[i]), namevec=namevec, hessian=hessian,
                                   function.arg=c(get("pars", inherits=FALSE), "theta", "D"))
    }

  } else {

    for(i in 1:cats) {
      funList[[i]] <- stats::deriv(expr=parse(text=equation[i]), namevec="theta", hessian=hessian,
                                   function.arg=c(get("pars", inherits=FALSE), "theta", "D"))
    }

  }

  # return results
  funList

}




