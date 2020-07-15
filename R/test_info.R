#' Item and Test Information Function
#'
#' @description This function computes both item and test information functions (Hambleton et al., 1991) given a set of theta values.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...), an object of class \code{\link{est_item}}
#' obtained from the function \code{\link{est_item}}, or an object of class \code{\link{est_irt}} obtained from the function \code{\link{est_irt}}.
#' The data.frame of item meta data can be easily obtained using the function \code{\link{shape_df}}. See below for details.
#' @param theta A vector of theta values where item and test information values are computed.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details A specific form of a data.frame should be used for the argument \code{x}. The first column should have item IDs,
#' the second column should contain the number of score categories of the items, and the third column should include IRT models.
#' The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items.
#' Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. From the fourth column, item parameters should be included.
#' For dichotomous items, the fourth, fifth, and sixth columns represent the item discrimination (or slope), item difficulty, and
#' item guessing parameters, respectively. When "1PLM" or "2PLM" is specified for any items in the third column, NAs should be inserted
#' for the item guessing parameters. For polytomous items, the item discrimination (or slope) parameters should be contained in the
#' fourth column and the item threshold (or step) parameters should be included from the fifth to the last columns.
#' When the number of categories differs between items, the empty cells of item parameters should be filled with NAs.
#' In this package, item step parameters should be used for the (generalized) partial credit model. The item step parameter is the
#' overall item difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for each category.
#' Thus, the number of step parameters for an item with m categories is m-1 because a step parameter for the first category does not
#' affect the category probabilities. For example, if an item has five categories under the (generalized) partial credit model,
#' four step parameters should be specified. An example of a data.frame with a single-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 1.736 \tab  1.501 \tab  0.203 \cr
#'   ITEM4  \tab 2 \tab 3PLM \tab 0.835 \tab -1.049 \tab  0.182 \cr
#'   ITEM5  \tab 2 \tab DRM \tab 0.926 \tab  0.394 \tab  0.099
#' }
#' And an example of a data.frame for a mixed-format test is as follows:
#' \tabular{lrlrrrrr}{
#'   ITEM1  \tab 2 \tab 1PLM \tab 1.000 \tab  1.461 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM2  \tab 2 \tab 2PLM \tab 1.921 \tab -1.049 \tab         NA \tab         NA \tab         NA\cr
#'   ITEM3  \tab 2 \tab 3PLM \tab 0.926 \tab  0.394 \tab  0.099 \tab         NA \tab         NA\cr
#'   ITEM4  \tab 2 \tab DRM \tab 1.052 \tab -0.407 \tab  0.201 \tab         NA \tab         NA\cr
#'   ITEM5  \tab 4 \tab GRM  \tab 1.913 \tab -1.869 \tab -1.238 \tab -0.714 \tab         NA \cr
#'   ITEM6  \tab 5 \tab GRM  \tab 1.278 \tab -0.724 \tab -0.068 \tab  0.568 \tab  1.072\cr
#'   ITEM7  \tab 4 \tab GPCM  \tab 1.137 \tab -0.374 \tab  0.215 \tab  0.848 \tab         NA \cr
#'   ITEM8  \tab 5 \tab GPCM  \tab 1.233 \tab -2.078 \tab -1.347 \tab -0.705 \tab -0.116
#' }
#' For more details about the parameterization of the (generalized) partial credit model, see \code{IRT Models} section in
#' the page of \code{\link{irtplay-package}} for more details about the IRT models. An easier way to create a data.frame for
#' the argument \code{x} is by using the function \code{\link{shape_df}}.
#'
#' @return This function returns an object of class \code{\link{test.info}}. This object contains item and test information values
#' given the specified theta values.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Hambleton, R. K., & Swaminathan, H., & Rogers, H. J. (1991) \emph{Fundamentals of item response theory}.
#' Newbury Park, CA: Sage.
#'
#' @seealso \code{\link{plot.test.info}}, \code{\link{shape_df}}, \code{\link{est_item}}
#'
#' @examples
#' ## example 1.
#' ## using the function "shape_df" to create a data.frame of test meta data
#' # create a list containing the dichotomous item parameters
#' par.dc <- list(a=c(1.1, 1.2, 0.9, 1.8, 1.4),
#'                b=c(0.1, -1.6, -0.2, 1.0, 1.2),
#'                g=rep(0.2, 5))
#'
#' # create a list containing the polytomous item parameters
#' par.py <- list(a=c(1.4, 0.6),
#'                d=list(c(0.0, -1.9, 1.2), c(0.4, -1.1, 1.5, 0.2)))
#'
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 2, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("DRM", "GRM", "DRM", "DRM", "GPCM", "DRM", "DRM")
#'
#' # create an item meta data set
#' test <- shape_df(par.dc=par.dc, par.py=par.py, cats=cats, model=model) # create a data.frame
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' test.info(x=test, theta=theta, D=1)
#'
#'
#' ## example 2.
#' ## using a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-2, 2, 0.1)
#'
#' # compute item and test information values given the theta values
#' test.info(x=test_flex, theta=theta, D=1)
#'
#' @export
test.info <- function(x, ...) UseMethod("test.info")

#' @describeIn test.info Default method to compute item and test information functions for a data.frame \code{x} containing the item meta data.
#' @export
test.info.default <- function(x, theta, D=1, ...) {

  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }

  # clear the item meta data set
  x <- back2df(metalist2(x))

  # extract information
  any.dc <- any(x[, 2] == 2)
  any.py <- any(x[, 2] > 2)
  id <- x[, 1]
  meta <- metalist2(x)

  # if there are dichotomous items
  if(any.dc) {

    # item information matrix
    infomat_dc <- info.dich(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)

  } else {
    infomat_dc <- NULL
  }

  # if there are polytomous items
  if(any.py) {

    # check the number of polytomous items
    n.plm <- length(meta$plm$a)

    # item information matrix
    infomat_py <- array(0, c(n.plm, length(theta)))
    for(i in 1:n.plm) {
      infomat_py[i, ] <- info.poly(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=meta$plm$model[i])
    }

  } else {
    infomat_py <- NULL
  }

  # creat a item infomation matrix for all items
  infomat <- rbind(infomat_dc, infomat_py)

  # re-order the item information maxtirx along with the original order of items
  pos <- c(meta$drm$loc, meta$plm$loc)
  if(length(pos) > 1) {
    infomat <- cbind(infomat[order(pos), ])
  } else {
    infomat <- infomat
  }
  rownames(infomat) <- id
  colnames(infomat) <- paste0("theta.", 1:length(theta))

  # create a vector for test infomation
  testInfo <- colSums(infomat)

  rr <- list(itemInfo=infomat, testInfo=testInfo, theta=theta)
  class(rr) <- c("test.info")

  rr

}


#' @describeIn test.info An object created by the function \code{\link{est_item}}.
#' @export
test.info.est_item <- function(x, theta, ...) {

  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  any.dc <- any(x[, 2] == 2)
  any.py <- any(x[, 2] > 2)
  id <- x[, 1]
  meta <- metalist2(x)

  # if there are dichotomous items
  if(any.dc) {

    # item information matrix
    infomat_dc <- info.dich(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)

  } else {
    infomat_dc <- NULL
  }

  # if there are polytomous items
  if(any.py) {

    # check the number of polytomous items
    n.plm <- length(meta$plm$a)

    # item information matrix
    infomat_py <- array(0, c(n.plm, length(theta)))
    for(i in 1:n.plm) {
      infomat_py[i, ] <- info.poly(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=meta$plm$model[i])
    }

  } else {
    infomat_py <- NULL
  }

  # creat a item infomation matrix for all items
  infomat <- rbind(infomat_dc, infomat_py)

  # re-order the item information maxtirx along with the original order of items
  pos <- c(meta$drm$loc, meta$plm$loc)
  if(length(pos) > 1) {
    infomat <- cbind(infomat[order(pos), ])
  } else {
    infomat <- infomat
  }
  rownames(infomat) <- id
  colnames(infomat) <- paste0("theta.", 1:length(theta))

  # create a vector for test infomation
  testInfo <- colSums(infomat)

  rr <- list(itemInfo=infomat, testInfo=testInfo, theta=theta)
  class(rr) <- c("test.info")

  rr

}

#' @describeIn test.info An object created by the function \code{\link{est_irt}}.
#' @export
test.info.est_irt <- function(x, theta, ...) {

  # extract information from an object
  D <- x$scale.D
  x <- x$par.est

  any.dc <- any(x[, 2] == 2)
  any.py <- any(x[, 2] > 2)
  id <- x[, 1]
  meta <- metalist2(x)

  # if there are dichotomous items
  if(any.dc) {

    # item information matrix
    infomat_dc <- info.dich(theta=theta, a=meta$drm$a, b=meta$drm$b, g=meta$drm$g, D=D)

  } else {
    infomat_dc <- NULL
  }

  # if there are polytomous items
  if(any.py) {

    # check the number of polytomous items
    n.plm <- length(meta$plm$a)

    # item information matrix
    infomat_py <- array(0, c(n.plm, length(theta)))
    for(i in 1:n.plm) {
      infomat_py[i, ] <- info.poly(theta=theta, a=meta$plm$a[i], d=meta$plm$d[[i]], D=D, pmodel=meta$plm$model[i])
    }

  } else {
    infomat_py <- NULL
  }

  # creat a item infomation matrix for all items
  infomat <- rbind(infomat_dc, infomat_py)

  # re-order the item information maxtirx along with the original order of items
  pos <- c(meta$drm$loc, meta$plm$loc)
  if(length(pos) > 1) {
    infomat <- cbind(infomat[order(pos), ])
  } else {
    infomat <- infomat
  }
  rownames(infomat) <- id
  colnames(infomat) <- paste0("theta.", 1:length(theta))

  # create a vector for test infomation
  testInfo <- colSums(infomat)

  rr <- list(itemInfo=infomat, testInfo=testInfo, theta=theta)
  class(rr) <- c("test.info")

  rr

}


# item information function for dichotomous data
info.dich <- function(theta, a, b, g, D=1) {


  # check the numbers of examinees and items
  nstd <- length(theta)
  nitem <- length(a)

  # when both the numbers of examiness and items are greater than 1
  a <- matrix(a, nrow=nstd, ncol=nitem, byrow=TRUE)
  b <- matrix(b, nrow=nstd, ncol=nitem, byrow=TRUE)
  g <- matrix(g, nrow=nstd, ncol=nitem, byrow=TRUE)

  # calculate item information
  z <- D * a * (theta - b)
  numer <- D^2 * a^2 * (1 - g)
  denom <- (g + exp(z)) * (1 + exp(-z))^2
  info <- numer / denom

  # return
  t(info)

}


# item information function for polytomous data
info.poly <- function(theta, a, d, D=1, pmodel) {

  pmodel <- toupper(pmodel)
  if(!pmodel %in% c("GRM", "GPCM")) stop("'pmodel' must be either 'GRM' or 'GPCM'.", call.=FALSE)

  if(pmodel == "GRM") {

    # check the number of step parameters
    m <- length(d)

    # check the numbers of examinees
    nstd <- length(theta)

    # calculate all the probabilities greater than equal to each threshold
    allPst <- drm(theta=theta, a=rep(a, m), b=d, g=0, D=D)
    if(nstd == 1L) {
      allPst <- matrix(allPst, nrow=1)
    }
    allQst <- 1 - allPst[, ,drop=FALSE]

    # calculate category probabilities
    P <- cbind(1, allPst) - cbind(allPst, 0)
    P <- ifelse(P == 0L, 1e-20, P)

    # compute the component values to obtain information
    Da <- D * a
    pq_st <- allPst * allQst
    deriv_Pstth <- Da * pq_st
    deriv_Pth <- cbind(0, deriv_Pstth) - cbind(deriv_Pstth, 0)
    w1 <- (1 - 2 * allPst) * deriv_Pstth
    deriv2_Pth <- Da * (cbind(0, w1) - cbind(w1, 0))

    # compute the informaton for all score categories
    ip <- (deriv_Pth^2 / P) - deriv2_Pth

  }

  if(pmodel == "GPCM") {

    # include zero for the step parameter of the first category
    d <- c(0, d)

    # check the number of step parameters
    m <- length(d) - 1

    # check the numbers of examinees and items
    nstd <- length(theta)

    # create a matrix for step parameters
    dmat <- matrix(d, nrow=nstd, ncol=(m+1), byrow=TRUE)

    # calculate category probabilities
    Da <- D * a
    z <- Da * (theta - dmat)
    numer <- exp(t(apply(z, 1, cumsum))) # numerator
    denom <- rowSums(numer) # denominator
    P <- numer / denom
    P <- ifelse(P == 0L, 1e-20, P)

    # compute the component values to obtain information
    denom2 <- denom^2
    d1th_z <- matrix(Da * (1:(m+1)), nrow=nstd, ncol=(m+1), byrow=TRUE)
    d1th_denom <- rowSums(numer * d1th_z)
    deriv_Pth <- (numer / denom2) * (d1th_z * denom - d1th_denom)
    denom4 <- denom2^2
    d2th_denom <- rowSums(numer * d1th_z^2)
    d1th_z_den <- d1th_z * denom
    part1 <- d1th_z_den * denom * (d1th_z_den - d1th_denom)
    part2 <- denom2 * (d1th_z * d1th_denom + d2th_denom) - 2 * denom * d1th_denom^2
    deriv2_Pth <- (numer / denom4) * (part1 - part2)

    # compute the informaton for all score categories
    ip <- (deriv_Pth^2 / P) - deriv2_Pth

  }

  # sum of all score category information
  info <- rowSums(ip)

  info

}
