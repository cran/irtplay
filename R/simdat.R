#' Simulated Response Data
#'
#' @description This function generates a simulated response data for a single- or a mixed-format test forms. For dichotomous
#' item response data, the IRT 1PL, 2PL, and 3PL models are available. For polytomous item response data, the graded response model,
#' the partial credit model, and the generalized partial credit model are available.
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...). This data.frame
#' can be easily obtained using the function \code{\link{shape_df}}. See below for details.
#' @param theta A vector of theta values.
#' @param a.dc A vector of item discrimination (or slope) parameters for dichotomous IRT models.
#' @param b.dc A vector of item difficulty (or threshold) parameters for dichotomous IRT models.
#' @param g.dc A vector of item guessing parameters for dichotomous IRT models.
#' @param a.py A vector of item discrimination (or slope) parameters for polytomous IRT models.
#' @param d.py A list containing vectors of item threshold (or step) parameters for polytomous IRT models.
#' @param cats A vector containing the number of score categories for items.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param pmodel A vector of character strings specifying the polytomous model with which response data are simulated.
#' For each polytomous model, "GRM" for the graded response model or "GPCM" for the (generalized) partial credit model can be
#' specified.
#'
#' @details There are two ways of generating the simulated response data.
#' The first way is by using the argument \code{x} to read in a data.frame
#' of item meta data. In the data.frame, the first column should have item IDs, the second column should contain the number of score
#' categories of the items, and the third column should include IRT models. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM"
#' for dichotomous items, and "GRM" and "GPCM" for polytomous items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM",
#' "2PLM", and "3PLM") and "GRM" and "GPCM" represent the graded response model and (generalized) partial credit model, respectively.
#' From the fourth column, item parameters should be included. For dichotomous items, the fourth, fifth, and sixth columns represent
#' the item discrimination (or slope), item difficulty, and item guessing parameters, respectively. When "1PLM" or "2PLM" is specified
#' for any items in the third column, NAs should be inserted for the item guessing parameters. For polytomous items, the item
#' discrimination (or slope) parameters should be contained in the fourth column and the item threshold (or step) parameters
#' should be included from the fifth to the last columns. When the number of categories differs between items, the empty cells of item
#' parameters should be filled with NAs. In this package, item step parameters should be used for the (generalized) partial credit model.
#' The item step parameters are the overall item difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for
#' each category. Thus, the number of step parameters for an item with m categories is m-1 because a step parameter for the first category does
#' not affect the category probabilities. For example, if an item has five categories under the (generalized) partial credit model, four step
#' parameters should be specified. An example of a data.frame with a single-format test is as follows:
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
#' The second way is by directly specifying item parameters for each item for which response data should be simulated
#' (i.e., without using a data.frame, as shown in the examples that follow). In addition to item parameters,
#' \code{theta}, \code{cats}, \code{pmodel}, and  \code{D} should be specified as well. \code{g.dc} does not need to be specified when only
#' the 1PL and 2PL models are used for dichotomous item response data. For dichotomous items, 2s should be specified in \code{cats}.
#' For polytomous items, the number of score categories should be specified in \code{cats}. When a response data set is generated with
#' a mixed-format test, it is important to clearly specify \code{cats} according to the order of items in the test form. Suppose that the response
#' data of ten examinees are simulated with five items, including three dichotomous items and two polytomous items with three categories.
#' Also, suppose that the second and the forth items are the polytomous items. Then, \code{cats = c(2, 3, 2, 3, 2)} should be used.
#' Additionally, among those two polytomous items, if the first and second item response data are simulated from the graded response model
#' and generalized partial credit model, respectively, then \code{pmodel = c('GRM', 'GPCM')}.
#'
#' @return This function returns a vector or a matrix. When a matrix is returned, rows indicate theta values and columns represent items.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{drm}}, \code{\link{plm}}
#'
#' @examples
#' ## example 1.
#' ## simulates response data with a mixed-format test.
#' ## for the first two polytomous items, the generalized partial credit model is used
#' ## for the last polytomous item, the graded response model is used
#' # 100 examinees are sampled
#' theta <- rnorm(100)
#'
#' # set item parameters for three dichotomous items with the 3PL model
#' a.dc <- c(1, 1.2, 1.3); b.dc <- c(-1, 0, 1); g.dc <- rep(0.2, 3)
#'
#' # set item parameters for three polytomous item parameters
#' # note that 4, 4, and 5 categories are used for polytomous items
#' a.py <- c(1.3, 1.2, 1.7)
#' d.py <- list(c(-1.2, -0.3, 0.4), c(-0.2, 0.5, 1.6), c(-1.7, 0.2, 1.1, 2.0))
#'
#' # create a numeric vector of score categoires for both dichotomous and polytomous item data
#' # this score category vector is used to specify the location of the polytomous items
#' cats <- c(2, 2, 4, 4, 5, 2)
#'
#' # create a character vector of the IRT model for the polytomous items
#' pmodel <- c('GPCM', 'GPCM', 'GRM')
#'
#' # simulate the response data
#' simdat(theta=theta, a.dc=a.dc, b.dc=b.dc, g.dc=NULL,
#'        a.py=a.py, d.py=d.py, cats=cats, pmodel=pmodel, D=1)
#'
#'
#' ## example 2.
#' ## simulates response data with a sigle-format test with the 2PL model.
#' # create a numeric vector of score categoires for the three 2PL model items
#' cats <- rep(2, 3)
#'
#' # simulate the response data
#' simdat(theta=theta, a.dc=a.dc, b.dc=b.dc, cats=cats, D=1)
#'
#' ## example 3.
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # simulate the response data
#' simdat(x=test_flex, theta=theta, D=1) # use a data.farame of item meta information
#'
#' @export
simdat <- function(x=NULL, theta, a.dc, b.dc, g.dc=NULL, a.py, d.py, cats, pmodel, D=1) {

  if(!is.null(x)) { # if a data.frame is inserted for the item parameter infomation

    prm <- metalist2(x)

    # Set conditions
    nstd <- length(theta)
    nitem <- nrow(x)

    # Create an empty matrix
    res <- matrix(NA, nrow=nstd, ncol=nitem)

    # Set initial numbers
    idx.dc <- prm$drm$loc
    idx.py <- prm$plm$loc

    # Simulate data
    # (1) for dichotomous items
    if(!is.null(idx.dc)) {
      res.dc <- simdat_dc(theta, a=prm$drm$a, b=prm$drm$b, g=prm$drm$g, D=D)
      res[, idx.dc] <- res.dc
    }

    # (2) for polytomous items
    if(!is.null(idx.py)) {
      for(i in 1:length(idx.py)) {
        res[, idx.py[i]] <- simdat_py(theta, a=prm$plm$a[i], d=prm$plm$d[[i]], D=D,
                                      pmodel=prm$plm$model[i])
      }
    }

    if(nstd == 1 | nitem == 1) res <- as.numeric(res)


  } else {

    # check whether argument is correctly specified
    if(missing(cats)) stop("Category of each item is missing", call.=FALSE)

    # Set conditions
    nstd <- length(theta)
    nitem <- length(cats)

    # Create an empty matrix
    res <- matrix(NA, nrow=nstd, ncol=nitem)

    # Set initial numbers
    idx.dc <- which(cats == 2)
    idx.py <- which(cats > 2)

    # Simulate data
    # (1) for dichotomous items
    if(any(cats == 2)) {
      res.dc <- simdat_dc(theta, a=a.dc, b=b.dc, g=g.dc, D=D)
      res[, idx.dc] <- res.dc
    }

    # (2) for polytomous items
    if(any(cats > 2)) {
      for(i in 1:length(idx.py)) {
        res[, idx.py[i]] <- simdat_py(theta, a=a.py[i], d=d.py[[i]], D=D, pmodel=pmodel[i])
      }
    }

    if(nstd == 1 | nitem == 1) res <- as.numeric(res)

  }

  res

}


# A function for generating binary data for one item
simdat_dc <- function(theta, a, b, g, D) {

  # Number of examinees
  nstd <- length(theta)

  # Number of items
  nitem <- length(a)

  # Calculate true probability for each category
  fit <- drm(theta, a, b, g, D)

  # Sample random variables from uniform dist
  tmp <- stats::runif(nstd * nitem, 0, 1)
  rv_unif <- matrix(tmp, nrow=nstd, ncol=nitem)

  # Simulated Response data for one item
  res <- ifelse(fit >= rv_unif, 1, 0)
  if(nstd == 1 | nitem == 1) res <- as.numeric(res)

  res

}


# A function for generating categorical data for one item
simdat_py <- function(theta, a, d, D, pmodel) {

  # Number of examinees
  nstd <- length(theta)

  # Calculate true probability for each category
  fit <- plm(theta, a, d, D, pmodel)
  if(length(nstd) == 1L) {
    fit <- rbind(fit)
  }

  # Sample random variables from uniform dist
  rv_unif <- stats::runif(nstd, 0, 1)

  # Simulated Response data for one item
  cumprob <- t(apply(fit, 1, cumsum))
  matTF <- rv_unif >= cumprob
  res <- apply(matTF, 1, sum)

  res

}


