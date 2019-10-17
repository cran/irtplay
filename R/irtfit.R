#' Traditional IRT item fit statistics
#'
#' @description This function computes traditional IRT item fit statistics (i.e., \eqn{\chi^{2}} fit statistic (e.g., Bock, 1960; Yen, 1981),
#' loglikelihood ratio \eqn{\chi^{2}} fit statistic (\eqn{G^{2}}; McKinley & Mills, 1985), and infit and outfit statistics (Ames et al., 2015)) and returns
#' contingency tables to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics. Note that caution is needed in interpreting the infit and
#' outfit statistics for non-Rasch models. The saved object of this function, especially the object of contingency tables,
#' is used in the function of \code{\link{plot.irtfit}} to draw a raw and standardized residual plots (Hambleton et al., 1991).
#'
#' @param x A data.frame containing the item meta data (e.g., item parameters, number of categories, models ...) or an object of class \code{\link{est_item}}
#' obtained from the function \code{\link{est_item}}. The data.frame of item meta data can be easily obtained using the function \code{\link{shape_df}}.
#' See below for details.
#' @param score A vector of examinees' ability estimates.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param group.method A character string indicating how to group examinees along the ability scale for computing the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
#' Available methods are "equal.width" for grouping examinees by dividing the ability scale into intervals of equal width and "equal.freq"
#' for grouping examinees by dividing the ability scale into intervals with equal frequencies of examinees. However, "equal.freq" does not
#' always guarantee exactly the same frequency of examinees for all groups. Default is "equal.width". To divide the ability scale, the range
#' of ability scale and the number of divided groups must be specified in the arguments of \code{range.score} and \code{n.width}, respectively.
#' See below for details.
#' @param n.width An integer value to specify the number of divided groups along the ability scale. Default is 10. See below for details.
#' @param loc.theta A character string to indicate the location of ability point at each group (or interval) where the expected probabilities
#' of score categories are calculated using the IRT models. Available locations are "average" for computing the expected probability
#' at the average point of examinees' ability estimates in each group and "middle" for computing the expected probability at the midpoint of each group.
#' Default is "average".
#' @param range.score A vector of two numeric values to restrict the range of ability scale. All ability estimates less than
#' the first value are transformed to the first value. All ability estimates greater than the second value are transformed to the second value.
#' If NULL, the minimum and maximum values of ability estimates in the argument \code{score} is used as the range of ability scale. Note that
#' selection of grouping method in the argument \code{group.method} has nothing to do with the range of ability scale. Default is NULL.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param alpha A numeric value to specify significance \eqn{\alpha}-level of the hypothesis test for the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics.
#' Default is .05.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param overSR A numeric value to specify a criterion to find ability groups (or intervals) which have standardized residuals
#' greater than the specified value. Default is 2.
#' @param min.collapse An integer value to indicate the minimum frequency of cells to be collapsed when computing the \eqn{\chi^{2}} and \eqn{G^{2}}
#' fit statistics. Neighboring interval groups will be collapsed to avoid expected interval frequencies less than the specified minimum cell frequency.
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
#' To calculate the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics, two methods are used in the argument \code{group.method} to divide the ability scale
#' into several groups. If \code{group.method = "equal.width"}, the examinees are grouped based on equal length of intervals.
#' If \code{group.method = "equal.freq"}, the examinees are grouped so that all groups have equal frequencies. However, the grouping method
#' of "equal.freq" does guarantee that every group has the exactly same frequency of examinees. This is because the examinees are divided by
#' the same size of quantile.
#'
#' When dividing the ability scale into intervals to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics, the intervals should be wide enough not to include
#' too small number of examinees. On the other hand, the interval should be narrow enough to include homogeneous examinees in terms of ability
#' (Hambleton et al, 1991). Thus, if you want to divide the ability scale into other than ten groups, you need to specify the number of groups
#' in the argument \code{n.width}. Yen (1981) fixed the number of groups to 10, whereas Bock (1960) allowed for any number of groups.
#'
#' Regarding degrees of freedom (\emph{df}), the \eqn{\chi^{2}} is assumed to be distributed approximately as a chi-square with \emph{df} equal to
#' the number of groups less the number of the IRT model parameters (Ames et al., 2015) whereas the \eqn{G^{2}} is assumed to be distributed approximately
#' as a chi-square with \emph{df} equal to the number of groups (Ames et al., 2015; Muraki & Bock, 2003)
#'
#' Note that if "DRM" is specified for an item in the item meta data set, the item is considered as "3PLM" to compute degrees of freedom of
#' the \eqn{\chi^{2}} fit statistic.
#'
#' @return This function returns an object of class \code{\link{irtfit}}. Within this object, several internal objects are contained such as:
#' \item{fit_stat}{A data.frame containing the results of three IRT fit statistics (i.e., \eqn{\chi^{2}} and \eqn{G^{2}}, infit, outfit statistics) across
#' all evaluated items. In the data.frame, the columns indicate item's ID, \eqn{\chi^{2}} fit statistic, \eqn{G^{2}} fit statistic, degrees of freedom for the \eqn{\chi^{2}},
#' degrees of freedom for the \eqn{G^{2}}, critical value for the \eqn{\chi^{2}}, critical value for the \eqn{G^{2}}, p-value for the \eqn{\chi^{2}},
#' p-value for the \eqn{G^{2}}, outfit statistic, infit statistic, the number of examinees used to compute the five fit statistics, and the proportion of
#' ability groups (or intervals), before collapsing the cells, that have standardized residuals greater than the specified criterion in the argument \code{overSR},
#' respectively.}
#' \item{contingency.fitstat}{A list of contingency tables used to compute the \eqn{\chi^{2}} and \eqn{G^{2}} fit statistics for all items.
#' Note that the collapsing cell strategy is implemented to these contingency tables.}
#' \item{contingency.plot}{A list of contingency tables used to draw a raw and standardized residual plots (Hambleton et al., 1991) in the function of
#' \code{\link{plot.irtfit}}. Note that the collapsing cell strategy is \emph{not} implemented to these contingency tables.}
#' \item{individual.info}{A list of data.frames including individual residual and variance values. Those information are used to compute
#' infit and outfit statistics.}
#' \item{item_df}{The item meta data specified in the argument \code{x}.}
#' \item{ancillary}{A list of ancillary information used in the item fit analysis.}
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{plot.irtfit}}, \code{\link{shape_df}}, \code{\link{est_item}}
#'
#' @references
#' Ames, A. J., & Penfield, R. D. (2015). An NCME Instructional Module on Item-Fit Statistics for Item Response Theory Models.
#' \emph{Educational Measurement: Issues and Practice, 34}(3), 39-48.
#'
#' Bock, R.D. (1960), \emph{Methods and applications of optimal scaling}. Chapel Hill, NC: L.L. Thurstone Psychometric Laboratory.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.
#'
#' McKinley, R., & Mills, C. (1985). A comparison of several goodness-of-fit statistics.
#' \emph{Applied Psychological Measurement, 9}, 49-57.
#'
#' Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
#' scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com
#'
#' Wells, C. S., & Bolt, D. M. (2008). Investigation of a nonparametric procedure for assessing goodness-of-fit in
#' item response theory. \emph{Applied Measurement in Education, 21}(1), 22-40.
#'
#' Yen, W. M. (1981). Using simulation results to choose a latent trait model. \emph{Applied Psychological Measurement, 5}, 245-262.
#'
#' @examples
#' ## example 1
#' ## use the simulated CAT data
#' # find the location of items that have more than 10,000 responses
#' over10000 <- which(colSums(simCAT_MX$res.dat, na.rm=TRUE) > 10000)
#'
#' # select the items that have more than 10,000 responses
#' x <- simCAT_MX$item.prm[over10000, ]
#'
#' # select the response data for the items
#' data <- simCAT_MX$res.dat[, over10000]
#'
#' # select the examinees' abilities
#' score <- simCAT_MX$score
#'
#' \donttest{
#' # compute fit statistics
#' fit1 <- irtfit(x=x, score=score, data=data, group.method="equal.width",
#'                n.width=10, loc.theta="average", range.score=NULL, D=1, alpha=0.05,
#'                missing=NA, overSR=2)
#'
#' # fit statistics
#' fit1$fit_stat
#'
#' # contingency tables
#' fit1$contingency.fitstat
#' }
#'
#' ## example 2
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the first two dichotomous items and last polytomous item
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(10)
#' score <- rnorm(1000, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' # compute fit statistics
#' fit2 <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
#'                n.width=11, loc.theta="average", range.score=c(-4, 4), D=1, alpha=0.05)
#'
#' # fit statistics
#' fit2$fit_stat
#'
#' # contingency tables
#' fit2$contingency.fitstat
#'
#' # residual plots for the first item (dichotomous item)
#' plot(x=fit2, item.loc=1, type = "both", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)
#'
#' # residual plots for the third item (polytomous item)
#' plot(x=fit2, item.loc=3, type = "both", ci.method = "wald", show.table=FALSE, ylim.sr.adjust=TRUE)
#'
#' @export
irtfit <- function(x, ...) UseMethod("irtfit")

#' @describeIn irtfit Default method to compute the traditional IRT item fit statistics for a data.frame \code{x} containing the item meta data.
#'
#' @export
#'
irtfit.default <- function(x, score, data, group.method=c("equal.width", "equal.freq"),
                           n.width=10, loc.theta="average", range.score=NULL, D=1, alpha=0.05,
                           missing=NA, overSR=2, min.collapse=1, ...) {


  # match.call
  cl <- match.call()

  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # consider DRM as 3PLM
  x$model <- as.character(x$model)
  # consider DRM as 3PLM
  if("DRM" %in% x$model) {
    x$model[x$model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items were considered as '3PLM' items in during the item parameter estimation."
    warning(memo, call.=TRUE)
  }

  # transform scores to a vector form
  if(is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check if there are items which have zero or one response frequency
  n.score <-  colSums(!is.na(data))
  if(all(n.score %in% c(0L, 1L))) {
    stop("Every item has frequency of zero or one for the item response data. Each item must have more than two item responses.", call.=TRUE)
  }

  if(any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # delete the items which have no frequency of scores from the data set
    x <- x[-del_item, ]
    data <- data[, -del_item]

    # warning message
    memo <- paste0(paste0("item ", del_item, collapse = ", "),
                   " is/are deleted in the analysis. Because the item(s) has/have frequency of zero or one for the item response data.")
    warning(memo, call.=TRUE)
  }

  # restrict the range of scores if required
  if(!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm=TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    purrr::map(1:nrow(x), .f=function(i)
      itemfit(item_meta=x[i, ], score=score, resp=data[, i], group.method=group.method,
              n.width=n.width, loc.theta=loc.theta, D=D, alpha=alpha, overSR=overSR, min.collapse=min.collapse))

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f=function(i) i$fit.stats) %>%
    do.call(what="rbind")
  fit_stat <- data.frame(id=x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f=function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f=function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f=function(i) i$individual.info)
  names(individual.info) <- x$id

  # return results
  rst <- list(fit_stat=fit_stat, contingency.fitstat=contingency.fitstat, contingency.plot=contingency.plot,
              item_df=x, individual.info=individual.info, ancillary=list(range.score=range.score, alpha=alpha, overSR=overSR, scale.D=D))
  class(rst) <- "irtfit"
  rst$call <- cl

  rst

}



#' @describeIn irtfit An object created by the function \code{\link{est_item}}.
#'
#' @export
#'
irtfit.est_item <- function(x, group.method=c("equal.width", "equal.freq"),
                            n.width=10, loc.theta="average", range.score=NULL, alpha=0.05,
                            missing=NA, overSR=2, min.collapse=1, ...) {


  # match.call
  cl <- match.call()

  # extract information from an object
  data <- x$data
  score <- x$score
  D <- x$scale.D
  x <- x$par.est

  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # consider DRM as 3PLM
  x$model <- as.character(x$model)
  # consider DRM as 3PLM
  if("DRM" %in% x$model) {
    x$model[x$model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items were considered as '3PLM' items in during the item parameter estimation."
    warning(memo, call.=TRUE)
  }

  # transform scores to a vector form
  if(is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # transform the response data to a matrix form
  data <- data.matrix(data)

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # check if there are items which have zero or one response frequency
  n.score <-  colSums(!is.na(data))
  if(all(n.score %in% c(0L, 1L))) {
    stop("Every item has frequency of zero or one for the item response data. Each item must have more than two item responses.", call.=TRUE)
  }

  if(any(n.score %in% c(0L, 1L))) {
    del_item <- which(n.score %in% c(0L, 1L))

    # delete the items which have no frequency of scores from the data set
    x <- x[-del_item, ]
    data <- data[, -del_item]

    # warning message
    memo <- paste0(paste0("item ", del_item, collapse = ", "),
                   " is/are deleted in the analysis. Because the item(s) has/have frequency of zero or one for the item response data.")
    warning(memo, call.=TRUE)
  }

  # restrict the range of scores if required
  if(!is.null(range.score)) {
    score <- ifelse(score < range.score[1], range.score[1], score)
    score <- ifelse(score > range.score[2], range.score[2], score)
  } else {
    tmp.val <- max(ceiling(abs(range(score, na.rm=TRUE))))
    range.score <- c(-tmp.val, tmp.val)
  }

  # compute item fit statistics and obtain contingency tables across all items
  fits <-
    purrr::map(1:nrow(x), .f=function(i)
      itemfit(item_meta=x[i, ], score=score, resp=data[, i], group.method=group.method,
              n.width=n.width, loc.theta=loc.theta, D=D, alpha=alpha, overSR=overSR, min.collapse=min.collapse))

  # extract fit statistics
  fit_stat <-
    purrr::map(fits, .f=function(i) i$fit.stats) %>%
    do.call(what="rbind")
  fit_stat <- data.frame(id=x$id, fit_stat)

  # extract the contingency tables used to compute the fit statistics
  contingency.fitstat <-
    purrr::map(fits, .f=function(i) i$contingency.fitstat)
  names(contingency.fitstat) <- x$id

  # extract the contingency tables to be used to draw residual plots
  contingency.plot <-
    purrr::map(fits, .f=function(i) i$contingency.plot)
  names(contingency.plot) <- x$id

  # extract the individual residuals and variances
  individual.info <-
    purrr::map(fits, .f=function(i) i$individual.info)
  names(individual.info) <- x$id

  # return results
  rst <- list(fit_stat=fit_stat, contingency.fitstat=contingency.fitstat, contingency.plot=contingency.plot,
              item_df=x, individual.info=individual.info, ancillary=list(range.score=range.score, alpha=alpha, overSR=overSR, scale.D=D))
  class(rst) <- "irtfit"
  rst$call <- cl

  rst

}



