#' Bind Fill
#'
#' @description This function creates a cbind matrix or rbind matrix using a list containing different length
#' of numeric vectors.
#' @param List A list containing different length of numeric vectors
#' @param type A character string specifying whether rbind is used or cbind is used.
#'
#' @return A matrix.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @examples
#' # sample list
#' score_list <- list(item1=c(0:3), item2=c(0:2), item3=c(0:5), item3=c(0:4))
#'
#' # examples
#' # 1) create a rbind with the sample score list
#' bind.fill(score_list, type="rbind")
#'
#' # 2) create a cbind with the sample score list
#' bind.fill(score_list, type="cbind")
#'
#' @import dplyr
#' @export

bind.fill <- function(List, type=c("rbind", "cbind")){
  type <- tolower(type)
  type <- match.arg(type)
  nm <- List
  nm <- purrr::map(nm, as.matrix)
  names(nm) <- 1:length(nm)
  n <- max(purrr::map_dbl(nm, nrow))
  df <-
    purrr::map_dfc(nm, function(x) {rbind(x, matrix(NA, n-nrow(x), ncol(x)))}) %>%
    as.matrix()
  switch(type,
         cbind = unname(df),
         rbind = unname(t(df))
  )
  
}

# a function to calculate a mean and variance at each theta point
cal_moment <- function(node, weight) {
  mu <- sum(node * weight)
  sigma2 <- sum(node^2 * weight) - mu^2
  rst <- c(mu=mu, sigma2=sigma2)
  rst
}


# This function divides the item response data sets into the two dichotomous (correct and incorrect)
# and one polytomous item parts.
divide_data <- function(data, cats, freq.cat) {

  # divide the data set for the mixed-item format
  if(any(cats == 2L)) {
    data1_drm <- data[, cats == 2L]
    data2_drm <- 1 - data1_drm
    data1_drm[is.na(data1_drm)] <- 0
    data2_drm[is.na(data2_drm)] <- 0
  } else {
    data1_drm <- NULL
    data2_drm <- NULL
  }
  if(any(cats > 2)) {
    data_plm <-
      freq.cat[cats > 2L] %>%
      do.call(what='cbind')
  } else {
    data_plm <- NULL
  }

  # return the results
  list(data1_drm=data1_drm, data2_drm=data2_drm, data_plm=data_plm)

}

