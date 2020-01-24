#' Compute Item/Test Characteristic Functions
#'
#' @description This function computes the item category probabilities, item characteristic function, and
#' test characteristic function given a set of theta values. The returned object of this function can be used
#' to draw the item or test characteristic curve using the function \code{\link{plot.traceline}}.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...) or an object
#' of class \code{\link{est_item}} obtained from the function \code{\link{est_item}}. See \code{\link{irtfit}}, \code{\link{test.info}},
#' or \code{\link{simdat}} for more details about the item meta data. The data.frame of item meta data can be easily obtained using
#' the function \code{\link{shape_df}}.
#' @param theta A vector of theta values.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#'          Default is 1.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return This function returns an object of class \code{\link{traceline}}. This object contains a list containing
#' the item category probabilities, item characteristic function, and test characteristic function.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{plot.traceline}}, \code{\link{est_item}}
#'
#' @examples
#' ## example
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-3, 3, 0.5)
#'
#' # compute the item category probabilities and item/test
#' # characteristic functions given the theta values
#' traceline(x=test_flex, theta, D=1)
#'
#' @export
traceline <- function(x, ...) UseMethod("traceline")

#' @describeIn traceline Default method to compute the item category probabilities, item characteristic function, and
#' test characteristic function for a data.frame \code{x} containing the item meta data.
#' @import purrr
#' @import dplyr
#' @export
traceline.default <- function(x, theta, D = 1, ...) {

  meta <- metalist2(x)

  # make the empty list and data.frame to contain probabilities
  icc <- NULL
  prob.cat <- list()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    prob.drm <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    f <- function(theta, a, b, g, D) {
      p <- drm(theta, a, b, g, D)
      data.frame(score.0=(1 - p), score.1=p)
    }
    args <- list(a=meta$drm$a, b=meta$drm$b, g=meta$drm$g)
    prob <- purrr::pmap(.l=args, .f=f, theta=theta, D=D)

    if(length(theta) == 1L) {
      prob.drm <- as.data.frame(t(prob.drm))
    } else{
      prob.drm <- as.data.frame(prob.drm)
    }
    colnames(prob.drm) <- meta$drm$id
    names(prob) <- meta$drm$id
    # fill the empty list and data.frame
    prob.cat <- c(prob.cat, prob)
    icc <- dplyr::bind_cols(icc, prob.drm)

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    if(length(theta) == 1L) {
      prob.plm <- purrr::map(prob.plm, function(x) matrix(x, nrow=1)) %>%
        purrr::map(as.data.frame) %>%
        purrr::map(function(x) stats::setNames(object=x, nm=paste0("score.", c(0:(ncol(x)-1)))))
    } else {
      prob.plm <- purrr::map(prob.plm, data.frame) %>%
        purrr::map(function(x) stats::setNames(object=x, nm=paste0("score.", c(0:(ncol(x)-1)))))
    }
    names(prob.plm) <- meta$plm$id

    # compute the TCCs for items
    prob.icc <- purrr::map_dfc(prob.plm, function(x) as.matrix(x) %*% c(0:(ncol(x)-1)))
    names(prob.icc) <- meta$plm$id

    # fill the empty list and data.frame
    prob.cat <- c(prob.cat, prob.plm)
    icc <- dplyr::bind_cols(icc, prob.icc)

  }

  # relocate columns as the original order of the data
  loc <- order(c(meta$drm$loc, meta$plm$loc))
  icc <- as.matrix(icc[, loc])
  if(ncol(icc) == 1L & !is.null(meta$drm)) colnames(icc) <- meta$drm$id
  prob.cat <- prob.cat[loc]

  # compute TCCs for theta values
  tcc <- rowSums(icc)

  # return results
  rst <- list(prob.cat=prob.cat, icc=icc, tcc=tcc, theta=theta)
  class(rst) <- "traceline"
  rst

}


#' @describeIn traceline An object created by the function \code{\link{est_item}}.
#' @import purrr
#' @import dplyr
#' @export
traceline.est_item <- function(x, theta, ...) {

  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  # listrize the meta data
  meta <- metalist2(x)

  # make the empty list and data.frame to contain probabilities
  icc <- NULL
  prob.cat <- list()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    prob.drm <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    f <- function(theta, a, b, g, D) {
      p <- drm(theta, a, b, g, D)
      data.frame(score.0=(1 - p), score.1=p)
    }
    args <- list(a=meta$drm$a, b=meta$drm$b, g=meta$drm$g)
    prob <- purrr::pmap(.l=args, .f=f, theta=theta, D=D)

    if(length(theta) == 1L) {
      prob.drm <- as.data.frame(t(prob.drm))
    } else{
      prob.drm <- as.data.frame(prob.drm)
    }
    colnames(prob.drm) <- meta$drm$id
    names(prob) <- meta$drm$id
    # fill the empty list and data.frame
    prob.cat <- c(prob.cat, prob)
    icc <- dplyr::bind_cols(icc, prob.drm)

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    if(length(theta) == 1L) {
      prob.plm <- purrr::map(prob.plm, function(x) matrix(x, nrow=1)) %>%
        purrr::map(as.data.frame) %>%
        purrr::map(function(x) stats::setNames(object=x, nm=paste0("score.", c(0:(ncol(x)-1)))))
    } else {
      prob.plm <- purrr::map(prob.plm, data.frame) %>%
        purrr::map(function(x) stats::setNames(object=x, nm=paste0("score.", c(0:(ncol(x)-1)))))
    }
    names(prob.plm) <- meta$plm$id

    # compute the TCCs for items
    prob.icc <- purrr::map_dfc(prob.plm, function(x) as.matrix(x) %*% c(0:(ncol(x)-1)))
    names(prob.icc) <- meta$plm$id

    # fill the empty list and data.frame
    prob.cat <- c(prob.cat, prob.plm)
    icc <- dplyr::bind_cols(icc, prob.icc)

  }

  # relocate columns as the original order of the data
  loc <- order(c(meta$drm$loc, meta$plm$loc))
  icc <- as.matrix(icc[, loc])
  if(ncol(icc) == 1L & !is.null(meta$drm)) colnames(icc) <- meta$drm$id
  prob.cat <- prob.cat[loc]

  # compute TCCs for theta values
  tcc <- rowSums(icc)

  # return results
  rst <- list(prob.cat=prob.cat, icc=icc, tcc=tcc, theta=theta)
  class(rst) <- "traceline"
  rst

}




# Item and Test Characteristic Functions
# This function is used for several functions in the packages
#' @import purrr
#' @import dplyr
trace <- function(meta, theta, D = 1) {

  # make the empty list and data.frame to contain probabilities
  tcc_df <- NULL
  icc_list <- list()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    prob.drm <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    if(length(theta) == 1L) {
      prob.drm <- as.data.frame(t(prob.drm))
    } else{
      prob.drm <- as.data.frame(prob.drm)
    }
    colnames(prob.drm) <- meta$drm$id

    # fill the empty list and data.frame
    icc_list <- c(icc_list, as.list(prob.drm))
    tcc_df <- dplyr::bind_cols(tcc_df, prob.drm)

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    names(prob.plm) <- meta$plm$id

    # compute the TCCs for items
    if(length(theta) == 1L) {
      prob.tcc <- purrr::map_dfc(prob.plm, function(x) x %*% c(0:(length(x)-1)))
    } else {
      prob.tcc <- purrr::map_dfc(prob.plm, function(x) as.matrix(x) %*% c(0:(ncol(x)-1)))
    }
    names(prob.tcc) <- meta$plm$id

    # fill the empty list and data.frame
    icc_list <- c(icc_list, prob.plm)
    tcc_df <- dplyr::bind_cols(tcc_df, prob.tcc)

  }

  # relocate columns as the original order of the data
  loc <- order(c(meta$drm$loc, meta$plm$loc))
  tcc_df <- as.matrix(tcc_df[, loc])
  icc_list <- icc_list[loc]

  # compute TCCs for theta values
  tcc <- rowSums(tcc_df)

  # return results
  rr <- list(trace_list=icc_list, trace_df=tcc_df, tcc=tcc, theta=theta)
  rr

}



# Item and Test Characteristic Functions
# This function is used only for the estimation of the IRT linking coefficients
#' @import purrr
#' @import dplyr
trace2 <- function(meta, theta, D = 1, type=c("icc", "tcc")) {

  # make the empty list and data.frame to contain probabilities
  tcc_df <- NULL
  icc_list <- list()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    prob.drm <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    if(length(theta) == 1L) {
      prob.drm <- as.data.frame(t(prob.drm))
    } else{
      prob.drm <- as.data.frame(prob.drm)
    }
    colnames(prob.drm) <- meta$drm$id

    # fill the empty list and data.frame
    if(type == "icc") {
      icc_list <- c(icc_list, as.list(prob.drm))
    }
    if(type == "tcc") {
      tcc_df <- dplyr::bind_cols(tcc_df, prob.drm)
    }
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    names(prob.plm) <- meta$plm$id

    # fill the empty list and data.frame
    if(type == "icc") {
      icc_list <- c(icc_list, prob.plm)
    }
    if(type == "tcc") {
      # compute the TCCs for items
      if(length(theta) == 1L) {
        prob.tcc <- purrr::map_dfc(prob.plm, function(x) x %*% c(0:(length(x)-1)))
      } else {
        prob.tcc <- purrr::map_dfc(prob.plm, function(x) as.matrix(x) %*% c(0:(ncol(x)-1)))
      }
      names(prob.tcc) <- meta$plm$id
      tcc_df <- dplyr::bind_cols(tcc_df, prob.tcc)
    }
  }

  # relocate columns as the original order of the data
  loc <- order(c(meta$drm$loc, meta$plm$loc))
  if(type == "icc") {
    icc_list <- icc_list[loc]
    tr <- icc_list
  }
  if(type == "tcc") {
    tcc_df <- as.matrix(tcc_df[, loc])

    # compute TCCs for theta values
    tcc <- rowSums(tcc_df)
    tr <- tcc
  }

  # return results
  rr <- list(trace=tr, theta=theta)
  rr

}


# Item and Test Characteristic Functions
# This function is used only for Haebara method in the IRT linking
#' @import purrr
#' @import dplyr
trace3 <- function(meta, theta, D = 1, type=c("icc", "tcc")) {

  # make the empty list and data.frame to contain probabilities
  tcc_df <- NULL
  icc_list <- list()

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    prob.drm <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    f <- function(theta, a, b, g, D) {
      p <- drm(theta, a, b, g, D)
      data.frame(score.0=(1 - p), score.1=p)
    }
    args <- list(a=meta$drm$a, b=meta$drm$b, g=meta$drm$g)
    prob <- purrr::pmap(.l=args, .f=f, theta=theta, D=D)

    if(length(theta) == 1L) {
      prob.drm <- as.data.frame(t(prob.drm))
    } else{
      prob.drm <- as.data.frame(prob.drm)
    }
    colnames(prob.drm) <- meta$drm$id
    names(prob) <- meta$drm$id

    # fill the empty list and data.frame
    if(type == "icc") {
      icc_list <- c(icc_list, prob)
    }
    if(type == "tcc") {
      tcc_df <- dplyr::bind_cols(tcc_df, prob.drm)
    }
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    names(prob.plm) <- meta$plm$id

    # fill the empty list and data.frame
    if(type == "icc") {
      icc_list <- c(icc_list, prob.plm)
    }
    if(type == "tcc") {
      # compute the TCCs for items
      if(length(theta) == 1L) {
        prob.tcc <- purrr::map_dfc(prob.plm, function(x) x %*% c(0:(length(x)-1)))
      } else {
        prob.tcc <- purrr::map_dfc(prob.plm, function(x) as.matrix(x) %*% c(0:(ncol(x)-1)))
      }
      names(prob.tcc) <- meta$plm$id
      tcc_df <- dplyr::bind_cols(tcc_df, prob.tcc)
    }
  }

  # relocate columns as the original order of the data
  loc <- order(c(meta$drm$loc, meta$plm$loc))
  if(type == "icc") {
    icc_list <- icc_list[loc]
    tr <- icc_list
  }
  if(type == "tcc") {
    tcc_df <- as.matrix(tcc_df[, loc])

    # compute TCCs for theta values
    tcc <- rowSums(tcc_df)
    tr <- tcc
  }

  # return results
  rr <- list(trace=tr, theta=theta)
  rr

}

# This function compuest item characteristic functions
# This function is used in item parameter estimation
#' @import purrr
trace4 <- function(meta, theta, D = 1) {

  # make the empty list and data.frame to contain probabilities
  icc_list <- vector('list', sum(length(meta$drm$a), length(meta$plm$a)))

  # when there are dichotomous items
  if(!is.null(meta$drm)) {

    # compute the probabilities of answerting correclty on items
    ps <- drm(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)
    qs <- 1 - ps
    if(length(theta) == 1L) {
      prob <- purrr::map(.x=1:length(meta$drm$a), .f=function(i) c(qs[i], ps[i]))
    } else {
      if(length(meta$drm$a) == 1L) {
        prob <- purrr::map(.x=1:length(meta$drm$a), .f=function(i) unname(cbind(qs, ps)))
      } else {
        prob <- purrr::map(.x=1:length(meta$drm$a), .f=function(i) cbind(qs[, i], ps[, i]))
      }
    }

    # fill the empty list
    icc_list[meta$drm$loc] <- prob
  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # extract polytomous model info
    model <- meta$plm$model

    # make a list of arguments
    args <- list(meta$plm$a, meta$plm$d, model)

    # compute the category probabilities of items
    prob.plm <- purrr::pmap(.l=args, .f=plm, theta=theta, D=D)
    # if(length(theta) == 1L) {
    #   prob.plm <- purrr::map(.x=prob.plm, .f=rbind)
    # }

    # fill the empty list and data.frame
    icc_list[meta$plm$loc] <- prob.plm
  }

  # return results
  rr <- list(trace=icc_list, theta=theta)
  rr

}
