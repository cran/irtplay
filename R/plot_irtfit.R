#' Draw raw and standardized residual plots
#'
#' @description This function provides graphical displays to look at residuals between the observed data
#' and model-based predictions (Hambleton, Swaminathan, & Rogers, 1991). This function gives two residual plots for
#' each score category of an item: (a) the raw residual plot and (b) the standardized residual plot. Note that
#' for dichotomous items the residual plots are drawn only for the score category of 1.
#'
#' @param x An object of class \code{\link{irtfit}}.
#' @param item.loc An integer value indicating that the \emph{n}th item (or the location of the item) is plotted. See below for
#' details.
#' @param type A character string indicating what type of residual plot is returned. Available options
#' are "icc" for the raw residual plot, "sr" for the standardized residual plot, and "both" for both of them.
#' Default is "both".
#' @param ci.method A character string indicating what method is used to estimate the confidence interval for the raw residual plot.
#' Available options are "wald" for Wald method, "cp" for Clopper-Pearson interval, "wilson" for Wilson score interval, and
#' "wilson.cr" for Wilson score interval with continuity correction. Default is "wald". See below for details.
#' @param show.table A logical value. If TRUE, a contingency table containing the information used to draw the residual
#' plots for the studied item is returned. Default is TRUE.
#' @param layout.col An integer value indicating the number of columns in the panel when a polytomous item is used.
#' Default is 2.
#' @param xlab.text A title for the x axis. If missing, the default string is used.
#' @param ylab.text A title for the y axis. If \code{type = "both"}, two character strings can be
#' specified for the raw residual and standardized residual plots, respectively. If missing,
#' the default strings are used.
#' @param main.text An overall title for the plot. If \code{type = "both"}, two character strings
#' can be specified for the raw residual and standardized residual plots, respectively. If missing,
#' the default strings are used.
#' @param lab.size The size of xlab and ylab. Default is 15.
#' @param main.size The size of \code{main.text}. Default is 15.
#' @param axis.size The size of labels along the x and y axes. Default is 15.
#' @param line.size The size of lines. Default is 1.
#' @param point.size The size of points. Default is 2.5.
#' @param strip.size The size of facet labels. Default is 12.
#' @param ylim.icc A vector of two numeric values specifying the range of y axis for the raw residual plot. Default is c(0, 1).
#' @param ylim.sr.adjust A logical value. If TRUE, the range of y axis for the standardized residual plot is adjusted for each item.
#' If FALSE, the range of y axis for the standardized residual plot is fixed to the values specified in the argument \code{ylim.sr}.
#' @param ylim.sr A vector of two numeric values specifying the range of y axis for the standardized residual plot.
#' Default is c(-4, 4).
#' @param ... Further arguments passed from the function \code{\link[ggplot2]{ggplot}}.
#'
#' @details All of the plots are drawn using the ggplot2 package.
#'
#' Once the results of the IRT model fit analysis are obtained from the function \code{\link{irtfit}},
#' an object of class \code{\link{irtfit}} can be used to draw the IRT raw residual and standardized residual plots.
#'
#' Because the residual plots are drawn for an item at a time, you have to indicate which item will be evaluated. For this,
#' you should specify an integer value, which is the location of the studied item, in the argument \code{item.loc}.
#' For example, if you want to draw the residual plots for the third item, then \code{item.loc = 3}.
#'
#' In terms of the raw residual plot, the argument \code{ci.method} is used to select a method to estimate the confidence intervals
#' among four methods. Those methods are "wald" for the Wald interval, which is based on the normal approximation (Laplace, 1812),
#' "cp" for Clopper-Pearson interval (Clopper & Pearson, 1934), "wilson" for Wilson score interval (Wilson, 1927), and
#' "wilson.cr" for Wilson score interval with continuity correction (Newcombe, 1998).
#' See \url{https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval} for more details about
#' the binomial proportion confidence intervals. Note that the width of confidence interval is determined by the \eqn{\alpha}-level
#' specified in the argument \code{alpha} of the function \code{\link{irtfit}}.
#'
#' Regarding the standardized residual plot, any standardized residuals greater than the specified criterion value
#' in the argument {\code{overSR}} of the function \code{\link{irtfit}} are displayed as triangles. Otherwise,
#' they are displayed as circles.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}
#'
#' @references
#' Clopper, C. J., & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial.
#' \emph{Biometrika, 26}(4), 404-413.
#'
#' Hambleton, R. K., Swaminathan, H., & Rogers, H. J. (1991).\emph{Fundamentals of item response theory}. Newbury Park, CA: Sage.
#'
#' Laplace, P. S. (1820).\emph{Theorie analytique des probabilites} (in French). Courcier.
#'
#' Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of seven methods.
#' \emph{Statistics in medicine, 17}(8), 857-872.
#'
#' Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference.
#' \emph{Journal of the American Statistical Association, 22}(158), 209-212.
#'
#' @examples
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the first two dichotomous items and last polytomous item
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df[c(1:2, 55), ]
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(1000, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' \donttest{
#' # compute fit statistics
#' fit <- irtfit(x=x, score=score, data=data, group.method="equal.freq",
#'                n.width=11, loc.theta="average", range.score=c(-4, 4), D=1, alpha=0.05, overSR=1.5)
#'
#' # residual plots for the first item (dichotomous item)
#' plot(x=fit, item.loc=1, type = "both", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)
#'
#' # residual plots for the third item (polytomous item)
#' plot(x=fit, item.loc=3, type = "both", ci.method = "wald", show.table=FALSE, ylim.sr.adjust=TRUE)
#'
#' # raw residual plot for the third item (polytomous item)
#' plot(x=fit, item.loc=3, type = "icc", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)
#'
#' # standardized residual plot for the third item (polytomous item)
#' plot(x=fit, item.loc=3, type = "sr", ci.method = "wald", show.table=TRUE, ylim.sr.adjust=TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @describeIn plot An object created by the function \code{\link{irtfit}}.
#' @export
#'
plot.irtfit <- function(x, item.loc=NULL, type = "both",
                        ci.method = c("wald", "cp", "wilson", "wilson.cr"), show.table = TRUE,
                        layout.col=2, xlab.text, ylab.text, main.text, lab.size=15, main.size=15,
                        axis.size=15, line.size=1.0, point.size=2.5, strip.size=12,
                        ylim.icc=c(0, 1), ylim.sr.adjust=FALSE, ylim.sr=c(-4, 4), ...) {

  ##-------------------------------------------------------------------------
  if(is.null(item.loc)) {
    stop("The location of item being plotted should be specified in 'item.loc'.", call.=FALSE)
  }

  if(!is.null(item.loc)) {
    if(length(item.loc) > 1) {
      stop("A length of 'item.loc' must be 1.", call.=FALSE)
    }
  }

  # prepare data
  # extract an meta information of an item in which residual plot is drawn
  item_meta <- x$item_df[item.loc, ]

  # obtain the number of score categories of the item
  score.cats <- item_meta$cats

  # extract a score range for the plot
  range.score <- x$ancillary$range.score

  # extract a scaling factor
  D <- x$ancillary$scale.D

  # extract an alpha level
  alpha <- x$ancillary$alpha

  # a criterion of standardized residuals
  overSR <- x$ancillary$overSR

  # extract contingency table for the item
  contingency <- x$contingency[[item.loc]]

  # extract total frequencies for each score group
  N <- contingency$N

  # subset of the frequency table for score categories
  freq <- dplyr::select(contingency, dplyr::starts_with("freq"))

  # subset of the expected probability table for score categories
  exp.prob <- dplyr::select(contingency, dplyr::starts_with("exp.prob"))

  # subset of the observed probability table for score categories
  obs.prop <-
    dplyr::select(contingency, .data$theta, dplyr::starts_with("obs.prop")) %>%
    dplyr::rename_all(.funs=function(k) gsub(pattern="obs.prop.", replacement="", x=k))

  ##-------------------------------------------------------------------------
  # standardized Residual
  # find a z-score corresponding to significance level
  zscore <- stats::qnorm(1-alpha)

  # compute raw residuals
  raw_resid <- obs.prop[, -1] - exp.prob

  # compute the standard errors
  se <- sqrt((exp.prob * (1 - exp.prob)) / N)

  # standardize the raw residuals
  std_resid <- raw_resid/se

  # a data.frame including the standardized residuals and and information
  # if the standardized residuals are greater than two SE.
  resid_df <-
    std_resid %>%
    dplyr::mutate(theta=contingency$theta) %>%
    reshape2::melt(id.vars="theta", variable.name="score", value.name = "std_resid") %>%
    dplyr::mutate(lessSR= factor(ifelse(abs(std_resid) > overSR, 1, 0), levels=c(0, 1)))
  resid_df$score <- gsub(pattern="^*", replacement = "Score: ", x=resid_df$score)

  # when the item is dichotomous
  if(score.cats == 2) {
    resid_df <-
      resid_df %>%
      dplyr::filter(.data$score == "Score: 1")
  }

  ##-------------------------------------------------------------------------
  # observed data information for ICC plot
  # reorganize the "obs.prop" data.frame
  df1 <- reshape2::melt(data=obs.prop, variable.name="score", id.vars="theta", value.name="obs.prop")

  # compute the confidence interval
  ci.method <- tolower(ci.method)
  ci.method <- match.arg(ci.method)
  ci <- suppressWarnings(
    switch(ci.method,
           wald =
             purrr::map2(.x=obs.prop[, -1], .y=zscore * se, .f=function(k, j) data.frame(lower=k - j, upper=k + j)) %>%
             do.call(what="rbind"),

           cp =
             purrr::map(freq, .f=function(k)
               purrr::map2_dfc(.x=k, .y=N, .f=function(i, j)
                 stats::binom.test(i, j, alternative="two.sided", conf.level=1-alpha)$conf.int) %>%
                 t() %>%
                 data.frame() %>%
                 dplyr::rename("lower"=.data$X1, "upper"=.data$X2)) %>%
             do.call(what="rbind"),

           wilson =
             purrr::map(freq, .f=function(k)
               purrr::map2_dfc(.x=k, .y=N, .f=function(i, j)
                 stats::prop.test(i, j, alternative="two.sided", conf.level=1-alpha, correct=FALSE)$conf.int) %>%
                 t() %>%
                 data.frame() %>%
                 dplyr::rename("lower"=.data$X1, "upper"=.data$X2)) %>%
             do.call(what="rbind"),

           wilson.cr =
             purrr::map(freq, .f=function(k)
               purrr::map2_dfc(.x=k, .y=N, .f=function(i, j)
                 stats::prop.test(i, j, alternative="two.sided", conf.level=1-alpha, correct=TRUE)$conf.int) %>%
                 t() %>%
                 data.frame() %>%
                 dplyr::rename("lower"=.data$X1, "upper"=.data$X2)) %>%
             do.call(what="rbind")
    )
  )

  # restrict the confidence intervesl from 0 to 1 when Wald statistic is used
  if(ci.method == "wald") {
    ci[ci < 0] <- 0
    ci[ci > 1] <- 1
  }

  # a data.frame including the observed proportion and the confidence interval
  obs_df <- cbind(df1, ci)
  obs_df$score <- gsub(pattern="^*", replacement = "Score: ", x=obs_df$score)

  # define theta values for x-axis
  theta <- seq(range.score[1], range.score[2], 0.1)

  # a data.frame including the ICC acros all score categories
  icc_df <-
    data.frame(theta, traceline(item_meta, theta, D=D)$prob.cat[[1]]) %>%
    dplyr::rename_all(.funs=function(k) gsub(pattern="score.", replacement="", x=k)) %>%
    reshape2::melt(id.vars="theta", variable.name="score", value.name = "icc")
  icc_df$score <- gsub(pattern="^*", replacement = "Score: ", x=icc_df$score)

  # when the item is dichotomous
  if(score.cats == 2) {
    obs_df <-
      obs_df %>%
      dplyr::filter(.data$score == "Score: 1")
    icc_df <-
      icc_df %>%
      dplyr::filter(.data$score == "Score: 1")
  }

  ##-------------------------------------------------------------------------
  # draw plots
  type <- tolower(type)
  # (1) draw ICC plots
  if(type == "icc") {

    if(missing(xlab.text)) xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- 'Probability'
    if(missing(main.text)) main.text <- paste0('Raw Residuals: ', item_meta$id)

    p <-
      ggplot2::ggplot(data=icc_df, mapping=aes_string(x="theta", y="icc"), ...) +
      ggplot2::geom_line(color="blue", size=line.size) +
      ggplot2::geom_point(data=obs_df, mapping=aes_string(x="theta", y="obs.prop"),
                          shape=1, size=point.size, color="red") +
      ggplot2::geom_segment(data=obs_df,
                            mapping=aes_string(x="theta", y="upper", xend="theta", yend="lower"),
                            lineend = "square", size=0.7) +
      ggplot2::labs(title=main.text, x=xlab.text, y=ylab.text) +
      ggplot2::ylim(ylim.icc[1], ylim.icc[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol=layout.col) +
      ggplot2::theme(plot.title = element_text(size=main.size),
                     axis.title = element_text(size=lab.size),
                     axis.text = element_text(size=axis.size)) +
      ggplot2::theme(strip.text.x = element_text(size = strip.size, face='bold'))

    print(p)

  }

  ##-------------------------------------------------------------------------
  if(type == "sr") {

    # (2) draw standardized residual plots
    if(missing(xlab.text)) xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- 'Standardized Residual'
    if(missing(main.text)) main.text <- paste0('Standardized Residuals: ', item_meta$id)
    point.color <- c("blue", "red")
    point.shape <- c(1, 2)

    # find a maximum value of the absolute standardized residuals
    # use the maximum value as the ylim value
    if(ylim.sr.adjust)  {
      y.max <- round(max(abs(resid_df$std_resid)) + 1, 0)
      ylim.sr <- c(-y.max, y.max)
    }

    p <-
      ggplot2::ggplot(data=resid_df, mapping=aes_string(x="theta", y="std_resid"), ...) +
      ggplot2::geom_point(mapping=aes_string(color="lessSR", shape="lessSR"), size=point.size) +
      ggplot2::labs(title=main.text, x=xlab.text, y=ylab.text) +
      ggplot2::ylim(ylim.sr[1], ylim.sr[2]) +
      ggplot2::xlim(range.score[1], range.score[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol=layout.col) +
      ggplot2::theme(plot.title = element_text(size=main.size),
                     axis.title = element_text(size=lab.size),
                     axis.text = element_text(size=axis.size)) +
      ggplot2::theme(strip.text.x = element_text(size = strip.size, face='bold')) +
      ggplot2::scale_shape_manual(values=point.shape) +
      ggplot2::scale_colour_manual(values=point.color) +
      ggplot2::theme(legend.position="none") +
      ggplot2::geom_hline(yintercept = c(-overSR, overSR), linetype="dashed")

    print(p)

  }

  ##-------------------------------------------------------------------------
  if(type == "both") {

    xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- c('Probability', 'Standardized Residual')
    if(missing(main.text)) {
      main.text <- c(paste0('Raw Residuals: ', item_meta$id),
                     paste0('Standardized Residuals: ', item_meta$id))
    }

    # (1) draw ICC plots
    p1 <-
      ggplot2::ggplot(data=icc_df, mapping=aes_string(x="theta", y="icc"), ...) +
      ggplot2::geom_line(color="blue", size=line.size) +
      ggplot2::geom_point(data=obs_df, mapping=aes_string(x="theta", y="obs.prop"),
                          shape=1, size=point.size, color="red") +
      ggplot2::geom_segment(data=obs_df,
                            mapping=aes_string(x="theta", y="upper", xend="theta", yend="lower"),
                            lineend = "square", size=0.7) +
      ggplot2::labs(title=main.text[1], x=xlab.text, y=ylab.text[1]) +
      ggplot2::ylim(ylim.icc[1], ylim.icc[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol=1) +
      ggplot2::theme(plot.title = element_text(size=main.size),
                     axis.title = element_text(size=lab.size),
                     axis.text = element_text(size=axis.size)) +
      ggplot2::theme(strip.text.x = element_text(size = strip.size, face='bold'))

    # (2) draw standardized residual plots
    point.color <- c("blue", "red")
    point.shape <- c(1, 2)

    # find a maximum value of the absolute standardized residuals
    # use the maximum value as the ylim value
    if(ylim.sr.adjust)  {
      y.max <- round(max(abs(resid_df$std_resid)) + 1, 0)
      ylim.sr <- c(-y.max, y.max)
    }

    p2 <-
      ggplot2::ggplot(data=resid_df, mapping=aes_string(x="theta", y="std_resid"), ...) +
      ggplot2::geom_point(mapping=aes_string(color="lessSR", shape="lessSR"), size=point.size) +
      ggplot2::labs(title=main.text[2], x=xlab.text, y=ylab.text[2]) +
      ggplot2::ylim(ylim.sr[1], ylim.sr[2]) +
      ggplot2::xlim(range.score[1], range.score[2]) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~score, ncol=1) +
      ggplot2::theme(plot.title = element_text(size=main.size),
                     axis.title = element_text(size=lab.size),
                     axis.text = element_text(size=axis.size)) +
      ggplot2::theme(strip.text.x = element_text(size = strip.size, face='bold')) +
      ggplot2::scale_shape_manual(values=point.shape) +
      ggplot2::scale_colour_manual(values=point.color) +
      ggplot2::theme(legend.position="none") +
      ggplot2::geom_hline(yintercept = c(-overSR, overSR), linetype="dashed")

    # combine the two plots in a window
    gridExtra::grid.arrange(p1, p2, ncol=2) # For multitple plots

  }

  # add more information to the contingency table
  raw_resid <- dplyr::rename_all(raw_resid, .funs=function(k) paste0("raw.resid.", k))
  se  <- dplyr::rename_all(se, .funs=function(k) gsub(pattern="exp.prob.", replacement="se.", x=k))
  std_resid <- dplyr::rename_all(std_resid, .funs=function(k) paste0("std.resid.", k))
  contingency <- data.frame(contingency, raw_resid, se, std_resid)

  # return the contingency table
  if(show.table) {
    return(contingency)
  }

}
