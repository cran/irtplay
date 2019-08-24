# Count Frequency and Proportion of score categories
#
# @description This function computes frequency/proportion of each category scores for each theta
#
# @param item A integer value indicating an item to be analyzed.
# @param cats A vector of score categories for items.
# @param score A list containing abiltiy estimates obtained from the IRT summed scoring method across all evaluated items.
# @param res_dat A matrix (or data.frame) of size (# of items * # of examinees) containing response data
# @param theta A list containing unique theta values for each item
# @param score_freq A list containing frequencies of examinees who have the each theta values for each item
# @param freq Logical whether frequencies or proportions of score categories are computed
#
# @return A list
#
#' @importFrom rlang .data
#'
countFreq <- function(item, cats, score, res_dat, theta, score_freq, freq=TRUE) {


  count_mat <-
    array(NA, c(length(theta[[item]]), cats[item])) %>%
    data.frame()

  for(i in 0:(cats[item]-1)) {

    tmp_mat <-
      data.frame(score=score[[item]], response=res_dat[, item]) %>%
      dplyr::filter(.data$response == i)
    sel_score <- factor(tmp_mat[, 1], levels=theta[[item]])

    if(freq) {
      res <- as.numeric(table(sel_score))
    } else {
      res <- as.numeric(table(sel_score)) / score_freq[[item]]
    }

    count_mat[, (i+ 1)] <- res

  }

  colnames(count_mat) <- paste0("score.", 0:(cats[item]-1))
  count_mat

}


# Count Frequency and Proportion of score categories under CAT
#
# @description This function computes frequency/proportion of each category scores for each theta under CAT
#
# @param item A integer value indicating an item to be analyzed.
# @param cats A vector of score categories for items.
# @param score A maxtrix of size (# of items * # of examinees) containing abiltiy estimates
# obtained from the IRT summed scoring method.
# @param fitItem_dat A list containing response data for each studied item
# @param theta A list containing unique theta values for each item
# @param score_freq A list containing frequencies of examinees who have the each theta values for each item
# @param freq Logical whether frequencies or proportions of score categories are computed
#
# @return A list
#
#' @importFrom rlang .data
#'
countFreq2 <- function(item, cats, score, fitItem_dat, theta, score_freq, freq=TRUE) {

  count_mat <-
    array(NA, c(length(theta[[item]]), cats[item])) %>%
    data.frame()

  for(i in 0:(cats[item]-1)) {

    tmp_mat <-
      data.frame(score=score[[item]], response=fitItem_dat[[item]]) %>%
      dplyr::filter(.data$response == i)
    sel_score <- factor(tmp_mat[, 1], levels=theta[[item]])

    if(freq) {
      res <- as.numeric(table(sel_score))
    } else {
      res <- as.numeric(table(sel_score)) / score_freq[[item]]
    }

    count_mat[, (i+ 1)] <- res

  }

  colnames(count_mat) <- paste0("score.", 0:(cats[item]-1))
  count_mat

}

