# Compute EAP Summed Score
#
# @description This function computes the expected a posterior (EAP) summed score (Thissen et al., 1995; Thissen & Orlando, 2001)
# for each examinee. The EAP summed score is the mean of the posterior density for the summed score (or observed score) given
# the item parameter estimates.
#
# @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
# See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data.
# This data.frame can be easily obtained using the function \code{\link{shape_df}}.
# @param data A matrix containing examinees' response data for the test items. A row and column indicate examinees and items, respectively.
# @param norm.prior A numeric vector of two components specifying a mean and standard deviation of the normal prior distribution.
# These two parameters are used to obtain the gaussian quadrature points and the corresponding weights from the normal distribution.
# Default is c(0,1).
# @param nquad An integer value specifying the number of gaussian quadrature points from the normal prior distribution. Default is 41.
# @param weights A two-column matrix or data.frame containing the theta values (in the first column) and the weights (in the second column)
# for the prior distribution. If missing, default values are used (see \code{norm.prior} and \code{nquad}).
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
# Default is 1.
#
# @details
#
# @return A vector of the EAP summed scores
#
# @author Hwanggyu Lim \email{hglim83@@gmail.com}
# @export
# @examples
# ## the use of a "-prm.txt" file obtained from a flexMIRT
# flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
# x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#
# # simulate the item responses of examinees
# set.seed(15)
# theta <- rnorm(20)
# data <- simdat(x, theta, D=1)
#
# # estimate the abilities
# eap_sum(x, data, norm.prior=c(0, 1), nquad=41, D=1)
#
eap_sum <- function(x, data, norm.prior = c(0, 1), nquad = 41, weights, D=1) {

  ##------------------------------------------------------------------------------------------------
  # check missing data
  if(any(is.na(data))) stop("There should be no missing data in the data set.", call.=FALSE)

  # select catogory variable
  cats <- x[, 2]

  # generate quadrature points and weights
  if(missing(weights)) {
    weights <- gen.weight(n=nquad, dist="norm", mu=norm.prior[1], sigma=norm.prior[2])
  } else {
    weights <- data.frame(weights)
  }

  # estimate likelihoods using lord-wingersky argorithem
  lkhd <- lwrc(x=x, theta=weights[, 1], prob=NULL, cats=cats, D=D)

  # estimate EAP for sum scores
  eap.est <- rep(NA, nrow(lkhd))
  se.est <- rep(NA, nrow(lkhd))
  for(i in 1:nrow(lkhd)) {
    ss.prob <- sum(lkhd[i,] * weights[ ,2])
    post <- lkhd[i, ] * weights[ ,2] / ss.prob
    ss.eap <- sum(post * weights[ ,1])
    ss.eap2 <- sum(post * weights[ ,1]^2)
    eap.est[i] <- ss.eap
    se.est[i] <- sqrt(ss.eap2 - ss.eap^2)
  }

  # assign the EAP summed scores to each examinee
  names(eap.est) <- 0:sum(cats-1)
  names(se.est) <- 0:sum(cats-1)
  sumScore <- apply(data, 1, sum)
  est_score <- eap.est[as.character(sumScore)]
  est_se <- se.est[as.character(sumScore)]
  score_table <- data.frame(sum.score=names(eap.est), est.theta=eap.est, se.theta=se.est)

  # return results
  rst <- list(est.theta=est_score, se.theta=est_se, sum.score=sumScore, score.table=score_table)
  rst

}
