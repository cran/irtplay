#' Dichotomous Response Model Probabilities
#'
#' @description This function computes the probability of correct answers for one or more items for a given set of theta values
#' using the IRT 1PL, 2PL, and 3PL models.
#'
#' @param theta A vector of ability values.
#' @param a A vector of item discrimination (or slope) parameters.
#' @param b A vector of item difficulty (or threshold) parameters.
#' @param g A vector of item guessing parameters.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#'          Default is 1.
#'
#' @details \code{g} does not need to be specified when the response probabilities of the 1PL and 2PL models are computed.
#'
#' @return This function returns a vector or matrix. When a matrix is returned, rows indicate theta values and columns represent items.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{plm}}
#'
#' @examples
#' ## when vectors are used for both theta values and item parameters (3PLM)
#' drm(c(-0.1, 0.0, 1.5), a=c(1, 2), b=c(0, 1), g=c(0.2, 0.1), D=1)
#'
#' ## when vectors are only used for item parameters (2PLM)
#' drm(0.0, a=c(1, 2), b=c(0, 1), D=1)
#'
#' ## when vectors are only used for theta values (3PLM)
#' drm(c(-0.1, 0.0, 1.5), a=1, b=1, g=0.2, D=1)
#'
#' @export
drm <- function(theta, a, b, g=NULL, D=1) {

  # check the numbers of examinees and items
  nstd <- length(theta)
  nitem <- length(a)

  # check the guessing parmaters
  if(is.null(g)) g <- rep(0, nitem)

  # when both the numbers of examiness and items are greater than 1
  if(nstd > 1 & nitem > 1) {
    a <- matrix(a, nrow=nstd, ncol=nitem, byrow=TRUE)
    b <- matrix(b, nrow=nstd, ncol=nitem, byrow=TRUE)
    g <- matrix(g, nrow=nstd, ncol=nitem, byrow=TRUE)
  }

  # calculate probability of correct answer
  Da <- D * a
  z <- Da * (theta - b)
  P <- g + (1 - g) / (1 + exp(-z))

  P

}

#' Polytomous Response Model Probabilities (GRM and GPCM)
#'
#' @description This function computes the probability of selecting a specific category for an item
#' for a given set of theta values using the graded response model, partial credit model, and generalized
#' partial credit model.
#'
#' @param theta A vector of ability values.
#' @param a A numeric value of item discrimination (or slope) parameter.
#' @param d A vector of item threshold (or step) parameters.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function  (if set to 1.7).
#' Default is 1.
#' @param pmodel A character string indicating the polytomous model being used. Available models are "GRM" for
#' the the graded response model and "GPCM" for the (generalized) partial credit model.
#'
#' @details When the category probabilities are computed for an item with the partial credit model, \code{a = 1} for that item.
#' When \code{pmodel = "GPCM"}, \code{d} should include step parameters. Item step parameters are the overall item
#' difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for each category. Thus, the number of step parameters
#' for an item with m categories is m-1 because a step parameter for the first category does not affect the category probabilities. For example,
#' if an item has five categories under the (generalized) partial credit model, four step parameters should be specified. For more details about
#' the parameterization of the (generalized) partial credit model, see \code{\link{irtfit}}.
#'
#' @return This function returns a vector or matrix. When a matrix is returned, rows indicate theta values and columns represent
#' categories of an item.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{drm}}, \code{\link{irtfit}}
#'
#' @examples
#' ## Category probabilities for an item with four categories
#' ## using a generalized partial credit model
#' plm(theta=c(-0.2, 0, 0.5), a=1.4, d=c(-0.2, 0, 0.5), D=1, pmodel='GPCM')
#'
#' ## Category probabilities for an item with five categories
#' ## using a graded response model
#' plm(theta=c(-0.2, 0, 0.5), a=1.2, d=c(-0.4, -0.2, 0.4, 1.5), D=1, pmodel='GRM')
#'
#' @export
plm <- function(theta, a, d, D=1, pmodel=c("GRM", "GPCM")) {

  pModel <- toupper(pmodel)
  model <- match.arg(pmodel)
  switch(model,
         GRM = grm(theta=theta, a=a, d=d, D=D),
         GPCM = gpcm(theta=theta, a=a, d=d, D=D)
  )

}


# IRT GPC model
gpcm <- function (theta, a, d, D = 1) {

  # include zero for the step parameter of the first category
  d <- c(0, d)

  # check the numbers of examinees and items
  nstd <- length(theta)

  # create a matrix for step parameters
  d <- matrix(d, nrow=nstd, ncol=length(d), byrow=TRUE)

  # calculate category probabilities
  z <- D * a * (theta - d)
  numer <- exp(t(apply(z, 1, cumsum))) # numerator
  denom <- rowSums(numer) # denominator
  P <- numer / denom
  if(nstd == 1) P <- as.numeric(P)

  # return
  P

}

# IRT GRM model
#' @import purrr
grm <- function(theta, a, d, D = 1) {

  # check the number of step parameters
  m <- length(d)

  # check the numbers of examinees
  nstd <- length(theta)

  # calculate all the probabilities greater than equal to each threshold
  allP <- drm(theta=theta, a=rep(a, m), b=d, g=0, D=D)
  if(nstd == 1L) {
    allP <- matrix(allP, nrow=1)
  }

  # calculate category probabilities
  P <- cbind(1, allP) - cbind(allP, 0)
  if(nstd == 1) P <- as.numeric(P)

  # return
  P

}


# This function computes the probability of correct answers for a dichotomous item
# This is used only when computing the analytical var-covariance matrix of item parameter estimates
drm2 <- function(par, theta, D=1, model=c("1PLM", "2PLM", "3PLM"), a.1pl=1) {

  # assign item parameters
  if(model == "1PLM") {
    a <- a.1pl; b <- par; g <- 0
  }
  if(model == "2PLM") {
    a <- par[1]; b <- par[2]; g <- 0
  }
  if(model == "3PLM") {
    a <- par[1]; b <- par[2]; g <- par[3]
  }

  # calculate probability of correct answer
  Da <- D * a
  z <- Da * (theta - b)
  P <- g + (1 - g) / (1 + exp(-z))

  P

}

# This function computes the probability of answering each score category for a dichotomous item
# This is used only when computing the analytical var-covariance matrix of item parameter estimates
plm2 <- function(par, theta, score, D=1, pmodel=c("GRM", "GPCM"), pcm.a=FALSE) {

  if(!pcm.a) {
    a <- par[1]
    d <- par[2:length(par)]
  } else {
    a <- 1
    d <- par[1:length(par)]
  }

  pModel <- toupper(pmodel)
  model <- match.arg(pmodel)
  ps <-
    switch(model,
           GRM = grm(theta=theta, a=a, d=d, D=D),
           GPCM = gpcm(theta=theta, a=a, d=d, D=D)
    )

  ps[score + 1]

}




