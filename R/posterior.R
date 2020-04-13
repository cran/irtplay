# This function computes the posterior distribution of abilities for examinees
# given their item response data, item parameters, and population distribution
posterior <- function(likehd, weights) {

  # joint likelihood matrix of likelihoods and population distribution
  joint_like <- purrr::map2_dfc(.x=weights[, 2], .y=1:nrow(weights), .f=function(x, y) likehd[, y] * x)

  # denominator of the posterior distribution
  denom <- rowSums(joint_like)

  # compute the posterior distribution of examinees across the node values
  # a row and column indicate the examinees and nodes, respectively
  posterior <- joint_like / denom

  # return results
  posterior

}


