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


