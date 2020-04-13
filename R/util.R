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
  n <- max(purrr::map_dbl(nm, nrow))
  df <-
    purrr::map_dfc(nm, function(x) rbind(x, matrix(NA, n-nrow(x), ncol(x)))) %>%
    as.matrix()
  switch(type,
         cbind = unname(df),
         rbind = unname(t(df))
  )

}



#' @export
print.irtfit <- function(x, ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Significance level for chi-square fit statistic:", x$ancillary$alpha, "\n\n")
  cat("Item fit statistics: \n")
  print(x$fit_stat, print.gap=2, quote=FALSE)
  cat("\n")
  cat("Caution is needed in interpreting infit and outfit statistics for non-Rasch models. \n")
  invisible(x)

}


#' @export
print.est_item <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("1. Convergence:\n")
  cat(x$convergence, "\n\n")
  cat("2. -2loglikelihood:", (-2 * x$loglikelihood), "\n\n")
  cat("3. Estimation Results:\n")
  cat("(1) Item Parameters \n")
  item.par <- purrr::modify_if(.x=x$estimates, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat("\n")
  cat("(2) Group Parameters \n")
  group.par <- round(x$group.par, digits=digits)
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)

}

#' @export
print.est_irt <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Summary of the Data \n")
  cat(" Number of Items: ", x$nitem, "\n", sep="")
  cat(" Number of Cases: ", x$ncase, "\n\n", sep="")

  cat("Summary of Estimation Process \n")
  cat(" Maximum number of EM cycles: ", x$MaxE, "\n", sep="")
  cat(" Convergence criterion of E-step: ", x$Etol, "\n", sep="")
  cat(" Number of rectangular quadrature points: ", nrow(x$weights), "\n", sep="")
  cat(" Minimum & Maximum quadrature points: ", x$weights[1, 1], ", ", -x$weights[1, 1], "\n", sep="")
  cat(" Number of free parameters: ", x$npar.est, "\n", sep="")
  cat(" Number of fixed items: ", length(x$fix.loc), "\n", sep="")
  cat(" Number of E-step cycles completed: ", x$niter, "\n", sep="")
  cat(" Maximum parameter change: ", x$maxpar.diff, "\n\n", sep="")

  cat("Processing time (in seconds) \n")
  cat(" EM algorithm: ", x$EMtime, "\n", sep="")
  cat(" Standard error computation: ", x$SEtime, "\n\n", sep="")

  cat("Convergence and Stability of Solution \n")
  cat(" First-order test: ", x$test.1, "\n", sep="")
  cat(" Second-order test: ", x$test.2, "\n", sep="")
  cat(" Computation of variance-covariance matrix: ", x$var.note, "\n\n", sep="")

  cat("Summary of Estimation Results \n")
  cat(" -2loglikelihood: ", (-2 * x$loglikelihood), "\n", sep="")
  cat(" Item Parameters: \n")
  item.par <- purrr::modify_if(.x=x$estimates, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat(" Group Parameters: \n")
  group.par <- round(x$group.par, digits=digits)
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)

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

