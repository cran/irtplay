#' Lord-Wingersky Recursion Formula
#'
#' @description This function computes the conditional distributions of number-correct (or observed) scores
#' given probabilities of category responses to items or given a set of theta values using the Lord and
#' Wingersky recursion formula (1984).
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}} or \code{\link{simdat}} for more details about the item meta data.
#' This data.frame can be easily obtained using the function \code{\link{shape_df}}. If \code{prob = NULL}, this data.frame is
#' used in the recursion formula. See below for details.
#' @param theta A vector of theta values where the conditional distribution of observed scores are computed.
#' The theta values are only required when a data.frame is specified in the argument \code{x}.
#' @param prob A matrix containing the probability of answering each category of an item. Each row indicates an item and
#' each column represents each category of the item. When the number of categories differs between items, the empty cells
#' should be filled with zeros or NA values. If \code{x = NULL}, this probability matrix is used in the recursion Formula.
#' @param cats A numeric vector specifying the number of categories for each item. For example, a dichotomous
#' item has two categories. This information is only required when a probability matrix is specified in the argument
#' \code{prob}.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function
#' (if set to 1.7). Default is 1.
#'
#' @details The Lord and Wingersky recursive algorithm is an efficient way of calculating the compound probabilities
#' of any number-correct scores on a test based on IRT models. This algorithm is particularly useful when computing
#' the IRT model-based observed score distribution for a test.
#'
#' To compute the conditional distributions of observed scores, either the item meta data set specified in \code{x} or
#' the probability matrix specified in \code{prob} can be used.
#'
#' @return This function returns a vector of the probabilities of obtaining every observed score on a test.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Kolen, M. J. & Brennan, R. L. (2004) \emph{Test Equating, Scaling, and Linking} (2nd ed.). New York:
#' Springer.
#'
#' Lord, F. & Wingersky, M. (1984). Comparison of IRT true score and equipercentile observed score equatings.
#' \emph{Applied Psychological Measurement, 8}(4), 453-461.
#'
#' @examples
#' ## example 1: when a matrix of probabilities is used as a data set
#' ## this is an example from Kolen and Brennan (2004, p. 183)
#' # create a matrix of probabilities of getting correct and incorrect answers for three items
#' probs <- matrix(c(.74, .73, .82, .26, .27, .18), nrow=3, ncol=2, byrow = FALSE)
#'
#' # creat a vector of score cagtegories for the three items
#' cats <- c(2,2,2)
#'
#' # compute the conditional distributions of observed scores
#' lwrc(prob=probs, cats=cats)
#'
#' ## example 2: when a matrix of probabilities is used as a data set
#' ## with a mixed-format test
#' # category probabilies for a dichotomous item
#' p1 <- c(0.2, 0.8, 0, 0, 0)
#' # category probabilies for a dichotomous item
#' p2 <- c(0.4, 0.6, NA, NA, NA)
#' # category probabilies for a polytomous item with five categories
#' p3 <- c(0.1, 0.2, 0.2, 0.4, 0.1)
#' # category probabilies for a polytomous item with three categories
#' p4 <- c(0.5, 0.3, 0.2, NA, NA)
#'
#' # rbind the probability vectors
#' p <- rbind(p1, p2, p3, p4)
#'
#' # creat a vector of score cagtegories for the four items
#' cats <- c(2, 2, 5, 3)
#'
#' # compute the conditional distributions of observed scores
#' lwrc(prob=p, cats=cats)
#'
#' ## example 3: when a data.frame for the item meta data is used instead of a probabiliy matrix
#' ## with a mixed-format test
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' x <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # compute the conditional distributions of observed scores
#' lwrc(x=x, theta=seq(-1, 1, 0.1), D=1)
#'
#' @export
lwrc <- function(x=NULL, theta, prob=NULL, cats, D=1) {

  if(is.null(x)) {

    logic <- sapply(1:ncol(prob), function(i) is.numeric(prob[[i]]))
    if(!all(logic)) stop("All numbers in 'prob' should be numeric.", call.=FALSE)

    if(nrow(prob) != length(cats)) {
      stop(paste0("There are ", nrow(prob), " items in the probability matrix (or data.frame), ",
                  "whereas there are ", length(cats), " items in the category vector."), call.=FALSE)
    }

    if(min(cats) < 2) {
      stop("Minimum number of categories for each item is 2", call.=FALSE)
    }

    # Probabilities for each category at 1st item
    p <- as.double(prob[1, 1:cats[1]])

    # Create a temporary vector to contain probabilities
    tmp <- c()

    # Possible observed score range for a first item
    obs.range <- 0:(cats[1]-1)

    # Caculate probabilities to earn an observed score by accumulating over remaining items
    for(j in 2:nrow(prob)) {  # should start from a 2nd item

      # Probability to earn zero score. This is a special case
      tmp[1] <- p[1]*prob[j, 1]

      # The range of category scores for an added item
      cat.range <- 0:(cats[j]-1)

      # Possible minimum and maximum observed scores when the item is added
      # but, except zero and perfect score
      min.s <- min(obs.range) + min(cat.range) + 1
      max.s <- max(obs.range) + max(cat.range) - 1

      # Calculate probability to obtain possible observed score after the item is added
      # but, still except zero and perfect score
      for(score in min.s:max.s) {

        # Difference between the possible observed score and each category score if the added item
        poss.diff <- score - min(cat.range):max(cat.range)

        # The difference above should be greater than or equal to the minimum observed score
        # where the item is added, and less than or equal to the maximum observed score where
        # the item is added
        cols <- poss.diff >= min(obs.range) & poss.diff <= max(obs.range)
        poss.diff <- poss.diff[cols]  # Now, these are remained difference

        # Find columns including the probability of getting the category score of added item
        # in which magnitude is the same as differences obtained from above
        cols.t <- rep(FALSE, ncol(prob))
        cols.t[1:length(cols)] <- cols

        # Final probability to earn the observed score
        # "p[rev(poss.diff + 1)]" is the probabilities for possible observed score
        # before the item is added
        prob.score <- sum(p[rev(poss.diff + 1)] * rev(prob[j, cols.t]))
        tmp[1+score] <- prob.score
      }

      # Probability to earn perfect score. This is a special case.
      tmp[score + 2] <- p[length(p)]*prob[j, cats[j]]

      # Update probabilities
      p <- tmp

      # Update the range of possible observed scores
      obs.range <- (min.s-1):(max.s+1)

      # Reset the temporary vector
      tmp <- c()

    }

    p

  } else {

    prepdat <- prep4lw2(x, theta, D)
    p <- lwRecurive(prepdat$prob.list, prepdat$score.cat)

    p

  }

}



###--------------------------------------------------------------------------------------------------------------------
# "lwRecursive" function
# Arg:
# prob.list: (list) each element of a list includes a probability matrix (or data.frame) for an item.
#					In each probability matrix, rows indicate theta values (or weights) and columns indicate
#					score categories. Thus, each elements in the matrix represents the probability that a person
#					who has a certain ability earns a certain score.
# score.cats: (vector) containes score categories for all items
lwRecurive <- function(prob.list, score.cats) {

  if(length(unique(sapply(prob.list, nrow))) != 1L) {
    stop("At least, one probability matrix (or data.frame) has the difference number of rows across all marices in a list", call.=FALSE)
  }

  if(length(prob.list) != length(score.cats)) {
    stop(paste0("There are ", nrow(prob.list), " items in the probability matrix (or data.frame) at each element of a list, ",
                "whereas there are ", length(score.cats), " items in the category vector."), call.=FALSE)
  }

  if(min(score.cats) < 2) {
    stop("Minimum number of categories for each item is 2", call.=FALSE)
  }

  # Probabilities for each category at 1st item
  p <- prob.list[[1]]

  # Possible observed score range for a first item
  obs.range <- 0:(score.cats[1]-1)

  # Create a temporary matrix to contain all probabilities
  tScore.range <- 0:sum(score.cats-1)
  tmp <- matrix(0, nrow=nrow(prob.list[[1]]), ncol=length(tScore.range))

  # Caculate probabilities to earn an observed score by accumulating over remaining items
  for(j in 2:length(prob.list)) {

    # Probability to earn zero score. This is a special case.
    tmp[,1] <- p[,1]*prob.list[[j]][,1]

    # The range of category scores for an added item
    cat.range <- 0:(score.cats[j]-1)

    # Possible minimum and maximum observed scores when the item is added
    # but, except zero and perfect score
    min.s <- min(obs.range) + min(cat.range) + 1
    max.s <- max(obs.range) + max(cat.range) - 1

    # Calculate probability to obtain possible observed score after the item is added
    # but, still except zero and perfect score
    for(score in min.s:max.s) {

      # Difference between the possible observed score and each category score if the added item
      poss.diff <- score - min(cat.range):max(cat.range)

      # The difference above should be greater than or equal to the minimum observed score
      # where the item is added, and less than or equal to the maximum observed score where
      # the item is added
      cols <- poss.diff >= min(obs.range) & poss.diff <= max(obs.range)
      poss.diff <- poss.diff[cols]  # Now, these are remained difference

      # Final probability to earn the observed score
      # "p[ ,(poss.diff + 1)]" is the probabilities for possible observed score
      prob.score <- apply(p[ ,(poss.diff + 1)] * prob.list[[j]][, cols], 1, sum)
      tmp[,1+score] <- prob.score

    }


    # Probability to earn perfect score. This is a special case.
    tmp[ ,score+2] <- p[ ,ncol(p)]*prob.list[[j]][ ,score.cats[j]]

    # Update the range of possible observed scores
    obs.range <- (min.s-1):(max.s+1)

    # Update probabilities
    p <- tmp[, 1:length(obs.range)]

  }

  colnames(p) <- paste0("score.", 0:sum(score.cats-1))
  rownames(p) <- paste0("theta.", 1:nrow(prob.list[[1]]))

  t(p)

}


# "prep4lw2" function
# This function is used only for "lwrc" function
prep4lw2 <- function(dframe, theta, D) {

  meta <- metalist2(dframe)

  # compute category probabilities for all items
  prob_cats <- trace2(meta, theta, D, type="icc")$trace
  prob_cats <- lapply(prob_cats, function(x) data.frame(x))

  # if there are dichotomous items
  if(!is.null(meta$drm)) {
    loc_drm <- meta$drm$loc
    score.1 <- prob_cats[loc_drm]
    prob_cats[loc_drm] <- lapply(score.1, function(x) data.frame(s0=1 - x, s1=x))
  }

  # extract score categories for all items
  score.cats <- c(meta$drm$cats, meta$plm$cats)
  pos <- c(meta$drm$loc, meta$plm$loc)
  score.cats <- score.cats[order(pos)] # reorder

  # Return results
  list(prob.list=prob_cats, score.cats=score.cats)

}

