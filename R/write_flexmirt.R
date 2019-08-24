#' Write a "-prm.txt" file for flexMIRT
#'
#' @description This function writes an output file of "-prm.txt" for flexMIRT. The current version of this function can be used
#' only for the unidimensional IRT models.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...).
#' See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}} for more details about the item meta data.
#' This data.frame can be easily obtained using the function \code{\link{shape_df}}.
#' @param file The destination file name.
#' @param norm.pop A numeric vector of two components specifying a mean and standard deviation of the normal
#' population distribution. Default is c(0,1).
#' @param rePrm A logical value indicating whether the item parameters in the item meta data
#' are the reparameterized item parameters. If TRUE, the item intercepts and logits of item guessing parameters
#' should be included in the item meta data. If FALSE, the item difficulty and item guessing parameters
#' should be included in the item meta data.
#'
#' @return A "-prm.txt" file.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Cai, L. (2015). flexMIRT 3.0 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
#' Chapel Hill, NC: Vector Psychometric Group.
#'
#' @export
#' @examples
#' ## use the simulated CAT data
#' # extract the item meta data
#' x <- simCAT_MX$item.prm
#'
#' # set a name of "-prm.txt" file
#' temp_prm <- file.path(tempdir(), "temp-prm.txt")
#'
#' # write out the "-prm.txt" file
#' write.flexmirt(x, file=temp_prm, norm.pop=c(0, 1), rePrm=FALSE)
#'
write.flexmirt <- function(x, file=NULL, norm.pop=c(0, 1), rePrm=TRUE) {

  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # warning message
  if (is.null(file)) {
    stop("You must specify the destination of the file.")
  }

  # open a new file
  prm_file <- file(file, open = "w")

  # select parameter columns
  param <- dplyr::select(x, dplyr::starts_with("par."))

  # write prm file one item by one item
  for(i in 1:nrow(x)) {

    # check meta data for each item
    id <- as.character(x$id[i])
    cats <- x$cats[i]
    model <- as.character(x$model[i])
    group <- 1
    n.factor <- 1

    # (1) Dichotomous item
    if(model %in% c("3PLM", "DRM")) {
      a <- param[i, 1]
      if(!rePrm) {
        c <- - a * param[i, 2]
        logitg <- log(param[i, 3] / (1 - param[i, 3]))
      } else {
        c <- param[i, 2]
        logitg <- param[i, 3]
      }
      a <- format(a, nsmall = 7)
      c <- format(c, nsmall = 7)
      logitg <- format(logitg, nsmall = 7)
      cat(c(1, id, group, n.factor, 1, cats, logitg, c, a), sep="\t",
          file=prm_file, fill=TRUE)
    }

    if(model %in% c("1PLM", "2PLM")) {
      a <- param[i, 1]
      if(!rePrm) {
        c <- - a * param[i, 2]
      } else {
        c <- param[i, 2]
      }
      a <- format(a, nsmall = 7)
      c <- format(c, nsmall = 7)
      cat(c(1, id, group, n.factor, 2, cats, c, a), sep="\t",
          file=prm_file, fill=TRUE)
    }

    # (2) Polytomous item: GRM
    if(model %in% c("GRM")) {
      a <- param[i, 1]
      if(!rePrm) {
        cs <- - a * as.numeric(param[i, 2:cats])
      } else {
        cs <- as.numeric(param[i, 2:cats])
      }
      a <- format(a, nsmall = 7)
      cs <- format(cs, nsmall = 7)
      cat(c(1, id, group, n.factor, 2, cats, cs, a), sep="\t",
          file=prm_file, fill=TRUE)
    }

    # (3) Polytomous item: PCM and GPCM
    if(model %in% c("PCM", "GPCM")) {
      a <- param[i, 1]
      if(!rePrm) {

        alpha <- c(1, rep(0, cats - 2))
        ds.new <- param[i, 2:cats]
        b <- sum(ds.new) / 4
        ds <- b - ds.new

        Tmat <- matrix(0, nrow=cats, ncol=(cats-1))
        Tmat[, 1] <- 0:(cats-1)
        for(k in 2:(cats-1)) {
          for(j in 2:(cats-1)) {
            Tmat[k, j] <- sin( pi*(j-1)*(k-1) / (cats-1) )
          }
        }

        c.vec <- rep(0, cats)
        c.vec[cats] <- - a * (cats - 1) * b
        for(k in length(ds):2) {
          c.vec[k] <- c.vec[k+1] - a * (as.numeric(ds[k]) - b)
        }

        gam <- as.numeric(solve(Tmat[-1, ]) %*% c.vec[-1])

      } else {
        alpha <- as.numeric(param[i, 2:(cats-1)])
        gam <- as.numeric(param[i, (1+cats-2+1+1):(1+cats-2+1+1+cats-2)])

      }
      a <- format(a, nsmall = 7)
      alpha <- format(alpha, nsmall = 7)
      gam <- format(gam, nsmall = 7)
      cat(c(1, id, group, n.factor, 3, cats, 0, alpha, a, 0, gam), sep="\t",
          file=prm_file, fill=TRUE)

    }

  }

  # group information
  cat(c(0, "Group1", 1, n.factor, 0, format(norm.pop[1], nsmall = 7), format(norm.pop[2], nsmall = 7)),
      sep="\t", file=prm_file, fill=TRUE)

  # close the file
  close(prm_file)

}
