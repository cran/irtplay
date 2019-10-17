# Evaluation of Inverse test characteristic curve (TCC) scoring
#
# @description This function analytically computes the conditional standard error of ability estimate (CSEE)
# for the inverse test characteristic curve (TCC) scoring method (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).
#
# @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
# See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data.
# This data.frame can be easily obtained using the function \code{\link{shape_df}}.
# @param eval.score A numeric vector containing theta values where the CSEEs need to be computed or
# a two column matrix containing theta values in the first column and the corresponding scale scores in the second column.
# @param scale.score A numeric vector indicating the corresponding scale scores for the raw scores.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
# Default is 1.
# @param constant A numeric value used to adjust zero and perfect raw sum scores, or the raw sum score equal to the sum of item guessing parameters,
# if necessary, to find estimable solutions for those raw sum scores. The zero raw score is forced to become the score of "zero raw score + constant"
# and the perfect raw score is forced to become the score of "perfect raw score - constant". If the 3PLM items are included in the item meta data,
# the raw sum score equal to the sum of item guessing parameters is forced to become the score of "the raw score + constant". Default is .1.
# @param constraint A logical value indicating whether the ability estimates will be restricted within a specific ability range
# specified in the argument \code{range.tcc}. If \code{constraint = TRUE}, all ability estimates less than the first value in the vector specified in
# the argument \code{range.tcc} are transformed to the first value and all ability estimates greather than the second value in the vector specified in
# the argument \code{range.tcc} are transformed to the second value. Also, when \code{constraint = TRUE} and the 3PLM items are contained
# in the item meta data, linear interpolation method is used to find the ability estimates for the raw sum scores less than the sum of item guessing
# parameters. When \code{constraint = FALSE} and the 3PLM items are contained in the item meta data, linear extrapolation method is used to find
# the ability estimates for the raw sum scores less than the sum of item guessing parameters. See below for details. Default is FALSE.
# @param range.tcc A numeric vector of two components to be used as the lower and upper bounds of ability estimates when \code{constraint = TRUE}.
# Default is c(-7, 7).
#
# @details When employing the IRT 3-parameter logistic model (3PLM) in a test, ability estimates for the raw sum scores less than the sum of item
# guessing parameters are not attainable. In this case, either of linear interpolation and linear extrapolation can be applied. Note that
# if \code{constraint = TRUE}, linear interpolation method is used. Otherwise, linear extrapolation method is used. Now, let \eqn{\theta_{min}} and
# \eqn{\theta_{max}} be the minimum and maximum ability estimates and \eqn{\theta_{X}} be the smallest raw score, X, greater than or equal to the sum of item
# guessing parameters. When linear interpolation method is used, a linear line is constructed between two points of (x=\eqn{\theta_{min}}, y=0) and
# (x=\eqn{\theta_{X}}, y=X). Because \code{constraint = TRUE}, \eqn{\theta_{min}} is the first value in the argument \code{range.tcc}.
# When linear extrapolation method is used, a linear line is constructed using two points of (x=\eqn{\theta_{X}}, y=X) and
# (x=\eqn{\theta_{max}}, y=maximum raw score). Then, ability estimates for the raw sum scores between zero and the smallest raw score greater than or equal
# to the sum of item guessing parameters are found using the constructed linear line.
#
# @return A data.frame
#
# @author Hwanggyu Lim \email{hglim83@@gmail.com}
#
# @references
# Kolen, M. J. & Brennan, R. L. (2004). \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
# Springer
#
# Kolen, M. J. & Tong, Y. (2010). Psychometric properties of IRT proficiency estimates.
# \emph{Educational Measurement: Issues and Practice, 29}(3), 8-14.
#
# Stocking, M. L. (1996). An alternative method for scoring adaptive tests.
# \emph{Journal of Educational and Behavioral Statistics, 21}(4), 365-389.
#
#' @import purrr
#' @import dplyr
#
se_tcc <- function(x, eval.score, scale.score=NULL, D=1, constant=0.1, constraint=TRUE, range.tcc=c(-5, 5)) {

  ##########################################
  ## collect all information for the test form
  cats <- x[, 2]
  drm <- any(cats == 2)

  # Find a maximum raw score
  max.score <- sum(cats - 1)

  # Find a minimum raw score
  if(drm) {
    gs <- x[which(cats == 2), 6]
    min.g <- sum(gs, na.rm=TRUE)
    min.score <- ifelse(min.g != 0L, ceiling(min.g), 0L)
  } else {
    min.g <- 0L
    min.score <- 0L
  }

  ##########################################
  # set the range of raw scores whose theta values are estimable
  scores <- seq(min.score, max.score, 1)

  # find the theta values correspoding to the raw scores
  theta <- raw2theta(raw.score=scores, x=x, D=D, constant=constant)

  # if the constraint is applied,
  # linear interpolation method is used to find ability estimates for the raw scores
  # less than the sum of item guessing parameters using the lower and upper bounds of abilities
  if(constraint) {

    # create a score conversion table
    if(min.g != 0L) {

      if(range.tcc[1] > theta[1]) {
        stop(paste0("A lower bound of theta must be less than ", round(theta[1], 3),
                    ", which is the theta estimate of the smallest raw score greater than or equal to the sum of item guessing parameters. Use the different lower bound of theta."),
             call.=FALSE)
      }

      # calculate slope and intercept for linear interpolation
      lessg.score <- 0:(min.score-1)
      slope <- scores[1] / (theta[1] - range.tcc[1])
      intercept <- - slope * range.tcc[1]
      lessg.theta <- (lessg.score - intercept) / slope

      theta <- c(lessg.theta, theta)
      theta <- ifelse(is.nan(theta), range.tcc[1], theta)
    }

    # restrict the range of thetas
    theta <- ifelse(theta < range.tcc[1], range.tcc[1], theta)
    theta <- ifelse(theta > range.tcc[2], range.tcc[2], theta)

    # score conversion table
    score.table <- data.frame(sum.score=0:max.score, est.theta=theta)

  }

  # if the constraint is not applied,
  # linear extrapolation method is used to find ability estimates for the raw scores
  # less than the sum of item guessing parameters using the maximum raw score and corresponding ability estimate
  # and the smallest raw score greater than or equal to the sum of item guessing parameters and corresponding ability estimate
  if(!constraint) {

    # create a score conversion table
    if(min.g != 0L) {
      # calculate slope and intercept for linear interpolation
      lessg.score <- 0:(min.score-1)
      slope <- (max.score - scores[1]) / (theta[length(theta)] - theta[1])
      intercept <- max.score - slope * theta[length(theta)]
      lessg.theta <- (lessg.score - intercept) / slope
      theta <- c(lessg.theta, theta)
    }
    score.table <- data.frame(sum.score=0:max.score, est.theta=theta)

  }

  # estimate observed score function using lord-wingersky argorithem
  if(is.null(scale.score)) {
    lkhd <- lwrc(x=x, theta=eval.score, D=D)
  } else {
    lkhd <- lwrc(x=x, theta=eval.score[, 1], D=D)
  }

  # calculate the standard error of ability estimates
  se <-
    purrr::map_dfc(.x=1:ncol(lkhd),
                   .f=function(i) sqrt(cal_moment(node=theta, weight=lkhd[, i])[2])) %>%
    as.numeric()
  if(!is.null(scale.score)) {
    se.scaled <-
      purrr::map_dfc(.x=1:ncol(lkhd),
                     .f=function(i) sqrt(cal_moment(node=scale.score, weight=lkhd[, i])[2])) %>%
      as.numeric()
  }

  # return results
  if(is.null(scale.score)) {
    rst <- data.frame(score=eval.score, scaled=NA, se.score=se, se.scaled=NA)
  } else {
    rst <- data.frame(score=eval.score[, 1], scaled=eval.score[, 2], se.score=se, se.scaled=se.scaled)
  }

  rst

}
