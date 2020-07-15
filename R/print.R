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
print.est_irt <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Item parameter estimation using MMLE-EM. \n")
  cat(x$niter, " E-step cycles were completed using ", nrow(x$weights), " quadrature points.", "\n", sep="")
  cat("First-order test: ", x$test.1, "\n", sep="")
  cat("Second-order test: ", x$test.2, "\n", sep="")
  cat("Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")
  cat("Log-likelihood: ", (x$loglikelihood), "\n", sep="")

  cat("\n")
  invisible(x)

}


#' @export
print.summary.est_irt <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  cat("\nCall:\n", paste(x$call.expr, sep = "\n", collapse = "\n"),
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
  cat(" Standard error computation: ", x$SEtime, "\n", sep="")
  cat(" Total computation: ", x$TotalTime, "\n\n", sep="")

  cat("Convergence and Stability of Solution \n")
  cat(" First-order test: ", x$test.1, "\n", sep="")
  cat(" Second-order test: ", x$test.2, "\n", sep="")
  cat(" Computation of variance-covariance matrix: \n", "  ", x$var.note, "\n\n", sep="")

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


#' @export
print.est_item <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  call.expr <- deparse(x$call)
  cat("\nCall:\n", paste(call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Item calibration using Method-A. \n")
  cat(x$convergence, "\n\n")
  cat("Log-likelihood: ", (x$loglikelihood), "\n", sep="")

  cat("\n")
  invisible(x)

}


#' @export
print.summary.est_item <- function(x, digits = max(2L, getOption("digits") - 5L), ...) {

  cat("\nCall:\n", paste(x$call.expr, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Summary of the Data \n")
  cat(" Number of Items in Response Data: ", x$nitem, "\n", sep="")
  cat(" Number of Excluded Items: ", x$nitem.del, "\n", sep="")
  cat(" Number of free parameters: ", x$npar.est, "\n", sep="")
  cat(" Number of Responses for Each Item: \n")
  print(data.frame(id=x$estimates$id, n=x$n.response), print.gap=2, quote=FALSE)
  cat("\n")

  cat("Processing time (in seconds) \n")
  cat(" Total computation: ", x$TotalTime, "\n\n", sep="")

  cat("Convergence of Solution \n")
  cat(" ", x$convergence, "\n\n", sep="")

  cat("Summary of Estimation Results \n")
  cat(" -2loglikelihood: ", (-2 * x$loglikelihood), "\n", sep="")
  cat(" Item Parameters: \n")
  item.par <- purrr::modify_if(.x=x$estimates, .p=is.numeric, .f=round, digits=digits)
  print(item.par, print.gap=2, quote=FALSE)
  cat("\n")
  cat(" Group Parameters: \n")
  group.par <- round(x$group.par, digits=digits)
  print(group.par, print.gap=2, quote=FALSE)
  cat("\n")
  invisible(x)

}
