#' Extract various elements from 'est_irt' or 'est_item' objects
#'
#' @description This function extracts various internal objects from an object of class \code{\link{est_irt}}
#' or \code{\link{est_item}}
#'.
#' @param x An object of class \code{\link{est_irt}} or \code{\link{est_item}}.
#' @param what A character string indicating what to extract.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details
#' Objects which can be extracted from the object of class \code{\link{est_irt}} include:
#'
#' \describe{
#' \item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
#' \item{par.est}{A data frame containing the item parameter estimates.}
#' \item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
#' observed information functions. The standard errors are estimated using the cross-production approximation method (Meilijson, 1989).}
#' \item{pos.par}{A data frame containing the position number of item parameters being estimated. The position information is useful
#' when interpreting the variance-covariance matrix of item parameter estimates.}
#' \item{covariance}{A matrix of variance-covariance matrix of item parameter estimates.}
#' \item{loglikelihood}{A sum of the log-likelihood values of the observed data set (marginal log-likelihood) across all estimated items.}
#' \item{aic}{A model fit statistic of Akaike information criterion based on the loglikelihood.}
#' \item{bic}{A model fit statistic of Bayesian information criterion based on the loglikelihood.}
#' \item{group.par}{A data frame containing the mean, variance, and standard deviation of latent variable prior distribution.}
#' \item{weights}{A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the (updated) latent variable prior distribution.}
#' \item{posterior.dist}{A matrix of normalized posterior densities for all the response patterns at each of the quadrature points.
#' The row and column indicate the response pattern and the quadrature point, respectively.}
#' \item{data}{A data frame of the examinees' response data set.}
#' \item{scale.D}{A scaling factor in IRT models.}
#' \item{ncase}{A total number of response patterns.}
#' \item{nitem}{A total number of items included in the response data.}
#' \item{Etol}{A convergence criteria for E steps of the EM algorithm.}
#' \item{MaxE}{The maximum number of E steps in the EM algorithm.}
#' \item{aprior}{A list containing the information of the prior distribution for item slope parameters.}
#' \item{bprior}{A list containing the information of the prior distribution for item difficulty (or threshold) parameters.}
#' \item{gprior}{A list containing the information of the prior distribution for item guessing parameters.}
#' \item{npar.est}{A total number of the estimated parameters.}
#' \item{niter}{The number of EM cycles completed.}
#' \item{maxpar.diff}{A maximum item parameter change when the EM cycles were completed.}
#' \item{EMtime}{Time (in seconds) spent for the EM cycles.}
#' \item{SEtime}{Time (in seconds) spent for computing the standard errors of the item parameter estimates.}
#' \item{TotalTime}{Time (in seconds) spent for total compuatation.}
#' \item{test.1}{Status of the first-order test to report if the gradients has vanished sufficiently for the solution to be stable.}
#' \item{test.2}{Status of the second-order test to report if the information matrix is positive definite, which is a prerequisite
#' for the solution to be a possible maximum.}
#' \item{var.note}{A note to report if the variance-covariance matrix of item parameter estimates is obtainable from the information matrix.}
#' \item{fipc}{A logical value to indicate if FIPC was used.}
#' \item{fipc.method}{A method used for the FIPC.}
#' \item{fix.loc}{A vector of integer values specifying the location of the fixed items when the FIPC was implemented.}
#' }
#'
#' Objects which can be extracted from the object of class \code{\link{est_item}} include:
#'
#' \describe{
#' \item{estimates}{A data frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
#' \item{par.est}{A data frame containing the item parameter estimates.}
#' \item{se.est}{A data frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
#' observed information functions.}
#' \item{loglikelihood}{A sum of the log-likelihood values of the complete data set across all estimated items.}
#' \item{data}{A data frame of the examinees' response data set.}
#' \item{score}{A vector of the examinees' ability values used as the fixed effects.}
#' \item{scale.D}{A scaling factor in IRT models.}
#' \item{convergence}{A string indicating the convergence status of the item parameter estimation.}
#' \item{nitem}{A total number of items included in the response data.}
#' \item{deleted.item}{The items which have no item response data. Those items are excluded from the item parameter estimation.}
#' \item{npar.est}{A total number of the estimated parameters.}
#' \item{n.response}{An integer vector indicating the number of item responses for each item used to estimate the item parameters.}
#' \item{TotalTime}{Time (in seconds) spent for total compuatation.}
#' }
#'
#' See \code{\link{est_irt}} and \code{\link{est_item}} for more details.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{est_irt}}, \code{\link{est_item}}
#'
#' @examples
#' \donttest{
#' # fit the 2PL model to LSAT6 data
#' mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2)
#'
#' # extract the item parameter estimates
#' (est.par <- getirt(mod.2pl, what="par.est"))
#'
#' # extract the standard error estimates
#' (est.se <- getirt(mod.2pl, what="se.est"))
#'
#' # extract the variance-covariance matrix of item parameter estimates
#' (cov.mat <- getirt(mod.2pl, what="covariance"))
#' }
#'
#' @export
getirt <- function(x, ...) UseMethod("getirt")

#' @describeIn getirt An object created by the function \code{\link{est_irt}}.
#' @export
getirt.est_irt <- function(x, what, ...){
  rst <- switch(what,
               estimates = x$estimates,
               par.est = x$par.est,
               se.est = x$se.est,
               pos.par = x$pos.par,
               covariance = x$covariance,
               loglikelihood = x$loglikelihood,
               aic = x$aic,
               bic = x$bic,
               group.par = x$group.par,
               weights = x$weights,
               data = x$data,
               ncase = x$ncase,
               nitem = x$nitem,
               Etol = x$Etol,
               MaxE = x$MaxE,
               aprior = x$aprior,
               bprior = x$bprior,
               gprior = x$gprior,
               npar.est = x$npar.est,
               niter = x$niter,
               maxpar.diff = x$maxpar.diff,
               EMtime = x$EMtime,
               SEtime = x$SEtime,
               TotalTime = x$TotalTime,
               test.1 = x$test.1,
               test.2 = x$test.2,
               var.note = x$var.note,
               fipc = x$fipc,
               fipc.method = x$fipc.method,
               fix.loc = x$fix.loc,
               stop(sprintf("Could not extract element \'%s\'", what), call.=FALSE)
  )

  rst

}

#' @describeIn getirt An object created by the function \code{\link{est_item}}.
#' @export
getirt.est_item <- function(x, what, ...){
  rst <- switch(what,
                estimates = x$estimates,
                par.est = x$par.est,
                se.est = x$se.est,
                loglikelihood = x$loglikelihood,
                group.par = x$group.par,
                data = x$data,
                score = x$score,
                scale.D = x$scale.D,
                convergence = x$convergence,
                nitem = x$nitem,
                deleted.item = x$deleted.item,
                npar.est = x$npar.est,
                n.response = x$n.response,
                TotalTime = x$TotalTime,
                stop(sprintf("Could not extract element \'%s\'", what), call.=FALSE)
  )

  rst

}

