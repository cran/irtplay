# This function compute the standard error of ability estimates based on the expected item information 
se_expinfo <- function(meta, theta, D=1, method="MLE", norm.prior=c(0, 1)) {
  
  # extract information
  any.dc <- !is.null(meta$drm)
  any.py <- !is.null(meta$plm)
  
  # if there are dichotomous items
  if(any.dc) {
    
    # item information matrix
    infomat_dc <- info.dich(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    
  } else {
    infomat_dc <- NULL
  }
  
  # if there are polytomous items
  if(any.py) {
    
    # check the number of polytomous items
    n.plm <- length(meta$plm$a)
    
    # item information matrix
    infomat_py <- array(0, c(n.plm, 1))
    for(i in 1:n.plm) {
      infomat_py[i, ] <- info.poly(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=meta$plm$model[i])
    }
    
  } else {
    infomat_py <- NULL
  }
  
  # create an item infomation matrix for all items
  infomat <- rbind(infomat_dc, infomat_py)
  
  # create a vector for test infomation
  t_info <- colSums(infomat)
  
  # compute the hessian matrix of the normal distribution
  if(method == "MAP") {
    
    # compute a hessian of prior distribution
    rst.prior <-
      logprior_deriv(val=theta, is.aprior=FALSE, D=NULL, dist="norm",
                     par.1=norm.prior[1], par.2=norm.prior[2])
    
    # extract the hessian
    hess.prior <- attributes(rst.prior)$hessian
    
    # add the hessian
    t_info <- sum(t_info, hess.prior)
    
  }
  
  # compute the standard error 
  se <- 1/sqrt(t_info)
  
  # return the results
  se
  
}
