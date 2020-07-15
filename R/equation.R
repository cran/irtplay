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

    # likelihood equation
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




