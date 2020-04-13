#  This function creates a list containing the equations of gredient vector and hessian matrix of negative log likelihood
# across all items in a test form
eqlist <- function(model, cats, loc_1p_const, loc_else,
                   fix.a.1pl, fix.a.gpcm, fix.g, a.val.1pl, a.val.gpcm, g.val,
                   use.aprior, use.gprior, aprior, gprior) {

  # empty list to contain the equations
  equation <- vector('list', length(model))

  ##----------------------------------------------------------------------
  # create the equations
  ##----------------------------------------------------------------------
  # 1) the dichotomous items: 1PLM with constrained slope values
  if("1PLM" %in% model & !fix.a.1pl) {

    # check the number of 1PLM items
    n.1PLM <- length(loc_1p_const)

    # create the equations
    eq_1p_const <- equation_item(model="1PLM", fix.a.1pl=FALSE, n.1PLM=n.1PLM, aprior=aprior, use.aprior=use.aprior)
    for(i in 1:n.1PLM) {
      equation[[loc_1p_const[i]]] <- eq_1p_const
    }

  }

  # 2) all other items
  if(length(loc_else) >= 1) {
    for(i in 1:length(loc_else)) {

      # prepare model and score category information
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if(score.cat == 2) {

        # create the equations
        equation[[loc_else[i]]] <- equation_item(model=mod, fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE), fix.g=fix.g, a.val.1pl=a.val.1pl,
                                                 g.val=g.val, n.1PLM=NULL, aprior=aprior, gprior=gprior, use.aprior=use.aprior, use.gprior=use.gprior)

      }

      # in case of a polytomous item
      if(score.cat > 2) {

        # create the equations
        equation[[loc_else[i]]] <- equation_item(model=mod, cats=score.cat,
                                                 fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE), a.val.gpcm=a.val.gpcm, n.1PLM=NULL,
                                                 aprior=aprior, use.aprior=use.aprior)

      }
    }
  }

  ##---------------------------------------------------------------
  # return the results
  equation

}


# This function returns the equations for gradient vector and hessian matrix for parameter estimation of an item.
equation_item <- function(model=c("1PLM", "2PLM", "3PLM", "GRM", "GPCM"), cats,
                          fix.a.1pl=TRUE, fix.a.gpcm=FALSE, fix.g=FALSE, a.val.1pl=1, a.val.gpcm=1, g.val=.2, n.1PLM=NULL,
                          aprior=list(dist="lnorm", params=c(1, 0.5)),
                          gprior=list(dist="beta", params=c(5, 17)),
                          use.aprior=FALSE, use.gprior=TRUE) {


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

    } else {

      # create the gradient vector and hessian matrix
      FUN.gh <- equation_plm(cats, pmodel=model, fix.a=fix.a.gpcm, a.val=a.val.gpcm,
                             aprior=aprior, use.aprior=use.aprior,
                             hessian=TRUE,
                             type="item")

    }

  }

  # return results
  FUN.gh

}
