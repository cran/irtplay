#' Generate Weights
#'
#' @description This function generates a set of weights based on a set of theta values to be used in the functions \code{\link{est_score}}
#' and \code{\link{sx2_fit}}.
#'
#' @param n An integer identifying the number of theta (or node) values for which weights are generated. Default is 41.
#' @param dist A character string specifying a probability distribution from which the weights are generated. Available distributions are
#' "norm" for a normal distribution, "unif" for a uniform distribution, and "emp" for an empirical distribution.
#' When \code{dist = "norm"}, either \code{n} or \code{theta} can be specified, when \code{dist = "unif"},
#' only \code{n} can be used, and when \code{dist = "emp"}, only \code{theta} can be used.
#' @param mu,sigma A mean and standard deviation of a normal distribution.
#' @param l,u Lower and upper limits of a uniform distribution.
#' @param theta A vector of empirical theta (or node) values for which weights are generated.
#'
#' @details When the argument \code{theta} is missing, \emph{n} weights can be generated from either the normal distribution or the uniform distribution.
#' Note that if \code{dist = "norm"}, gaussian quadrature points and weights from the normal distribution are generated. See
#' \code{gauss.quad.prob()} in the \pkg{statmod} package for more details.
#'
#' When the argument \code{theta} is not missing, the weights corresponding to the provided theta values are generated. Specifically, if
#' \code{dist = "norm"}, normalized weights from the normal distribution are returned. If \code{dist = "emp"}, every specified theta value has the equal
#' values of normalized weights.
#'
#' @return This function returns a data frame with two columns, where the first column has theta values (nodes) and the second column provides weights.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{est_score}}, \code{\link{sx2_fit}}
#'
#' @examples
#' ## example 1
#' ## generate 41 gaussian quadrature points and weights of normal distribution
#' gen.weight(n=41, dist = "norm", mu = 0, sigma = 1)
#'
#' ## example 2
#' ## generate 41 theta values and weights from the uniform normal distribution,
#' ## given the mininum value of -4 and the maximum value of 4
#' gen.weight(n=41, dist = "unif", l = -4, u = 4)
#'
#' ## example 3
#' ## generate the normalized weights from the standardized normal distribution,
#' ## given a set of theta values
#' theta <- seq(-4, 4, by=0.1)
#' gen.weight(dist = "norm", mu = 0, sigma = 1, theta = theta)
#'
#' ## example 4
#' ## generate the same values of normalized weights for the theta values that are
#' ## randomly sampled from the standardized normal distribution
#' theta <- rnorm(100)
#' gen.weight(dist = "emp", theta = theta)
#'
#' @export
gen.weight <- function(n=41, dist = "norm", mu = 0, sigma = 1, l = -4, u = 4, theta) {

  dist <- tolower(dist)

  if(missing(theta)) {

    if(dist == "emp") {
      stop("To use actual option in a distribution argument, theta values are needed.", call.=FALSE)
    }
    if(dist == "norm") {
      wts.nd <- statmod::gauss.quad.prob(n, dist=dist, mu=mu, sigma=sigma)
      nodes <- wts.nd$nodes
      whts <- wts.nd$weights
    }
    if(dist == "unif") {
      nodes <- seq(l, u, length.out=n)
      dens <- stats::dunif(nodes, min=l, max=u)
      whts <- dens / sum(dens)
    }

  }

  if(!missing(theta)) {

    if(dist == "emp") {
      nodes <- theta
      whts <- rep(1/length(theta), length(theta))
    }
    if(dist == "norm") {
      nodes <- theta
      dens <- stats::dnorm(theta, mean=mu, sd=sigma)
      whts <- dens / sum(dens)
    }

  }

  return(data.frame(theta=nodes, weight=whts))

}
