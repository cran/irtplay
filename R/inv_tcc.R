# Inverse test characteristic curve (TCC) scoring
#
# @description This function estimates examinees' abilities corresponding to thire raw sum scores using
# the inverse test characteristic curve (TCC) scoring method (e.g., Kolen & Brennan, 2004; Kolen & Tong, 2010; Stocking, 1996).
#
# @param x A data.frame containing the item metadata (e.g., item parameters, number of categories, models ...).
# See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item metadata.
# This data.frame can be easily obtained using the function \code{\link{shape_df}}.
# @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
# the examinees and items, respectively.
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
# Default is 1.
# @param constant A numeric value used to adjust zero and perfect raw sum scores, or the raw sum score equal to the sum of item guessing parameters,
# if necessary, to find estimable solutions for those raw sum scores. The zero raw score is forced to become the score of "zero raw score + constant"
# and the perfect raw score is forced to become the score of "perfect raw score - constant". If the 3PLM items are included in the item metadata,
# the raw sum score equal to the sum of item guessing parameters is forced to become the score of "the raw score + constant". Default is .1.
# @param constraint A logical value indicating whether the ability estimates will be restricted within a specific ability range
# specified in the argument \code{range.tcc}. If \code{constraint = TRUE}, all ability estimates less than the first value in the vector specified in
# the argument \code{range.tcc} are transformed to the first value and all ability estimates greather than the second value in the vector specified in
# the argument \code{range.tcc} are transformed to the second value. Also, when \code{constraint = TRUE} and the 3PLM items are contained
# in the item metadata, linear interpolation method is used to find the ability estimates for the raw sum scores less than the sum of item guessing
# parameters. When \code{constraint = FALSE} and the 3PLM items are contained in the item metadata, linear extrapolation method is used to find
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
# @return A list including a vector of the ability estimates and a vector of the standard errors of ability estimates. Also,
# raw sum scores of examinees and a table with the possible raw sum scores and corresponding ability estimates are returned as well.
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
inv_tcc <- function(x, data, D=1, constant=0.1, constraint=FALSE, range.tcc=c(-7, 7)) {
  
  # check missing data
  if(any(is.na(data))) stop("There should be no missing data in the data set.", call.=FALSE)
  
  ##########################################
  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))
  
  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }
  
  # clear the item metadata set
  x <- back2df(metalist2(x))
  
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
  lkhd <- lwrc(x=x, theta=theta, D=D)
  
  # calculate the standard error of ability estimates
  se.est <-
    purrr::map_dbl(.x=1:length(theta),
                   .f=function(i) sqrt(cal_moment(node=theta, weight=lkhd[, i])[2]))
  
  # assign the inverse TCC scores and the corresponding SE valsues to each examinee
  names(theta) <- 0:max.score
  names(se.est) <- 0:max.score
  sumScore <- apply(data, 1, sum)
  est_score <- theta[as.character(sumScore)]
  est_se <- se.est[as.character(sumScore)]
  
  # add SEs to the score conversion table
  score.table <- data.frame(score.table, se.theta=se.est)
  
  # return results
  rst <- list(est.theta=est_score, se.theta=est_se, sum.score=sumScore, score.table=score.table)
  rst
  
}

# "raw2theta" function
# This is a function to find theta corresponding to raw score using TCC.
#' @import purrr
#' @importFrom rlang .data
raw2theta <- function(raw.score, x, D = 1, constant = 0.1){

  drm <- any(x[, 2] == 2) # if there are dichotomous items

  # Find a maximum raw score
  max.score <- sum(x[, 2] - 1)

  # Adjust perfect and zero raw scores by adding or subtracting a constant
  raw.score[raw.score == max.score] <- max.score - constant
  raw.score[raw.score == 0] <- constant

  # Listrize the item meta data.frame
  meta <- metalist2(x)

  # Set low raw scores greater than sum of guessing parameters
  # to something that can be estimated
  if(drm) {
    gs <- x[which(x[, 2] == 2), 6]
    min.g <- sum(gs, na.rm=TRUE)
  } else {
    min.g <- 0
  }

  # Adjust raw scores less than or equal to the sum of item guessing parameters
  # if(min.g != 0L) {
  #   if(ceiling(min.g) %% min.g != 0L) {
  #     raw.score[raw.score <= min.g] <- ceiling(min.g)
  #   } else {
  #     raw.score[raw.score <= min.g] <- ceiling(min.g) + constant
  #   }
  # }
  if(min.g != 0L) {
    raw.score[raw.score <= min.g] <- min.g + constant
  }

  # set a loss function to estimate theta
  loss <- function(theta, score) trace2(meta, theta, D, "tcc")$trace - score

  # find the thetas corresponding to the raw scores
  root <- purrr::map_dbl(.x=raw.score, .f=function(x) stats::uniroot(loss, interval=c(-20, 20),
                                                                     score=x, maxiter = 2000, extendInt="yes")$root)
  # return results
  root

}

