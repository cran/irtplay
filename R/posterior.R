# This function computes the posterior distribution of abilities for examinees
# given their item response data, item parameters, and population distribution
posterior <- function(likehd, weights) {
  
  # joint likelihood matrix of likelihoods and population distribution
  ntheta <- nrow(weights)
  joint_like <- array(0, c(nrow(likehd), ntheta))
  for(i in 1:ntheta) {
    joint_like[, i] <- likehd[, i] * weights[i, 2]
  }
  
  # denominator of the posterior distribution
  denom <- rowSums(joint_like)
  
  # compute the posterior distribution of examinees across the node values
  # a row and column indicate the examinees and nodes, respectively
  posterior <- joint_like / denom
  
  # return results
  posterior
  
}



