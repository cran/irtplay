#' Plot ICC and TCC
#'
#' @description This function plots item or test characteristic curve using the ggplot2 package. The item characteristic
#' (or category) curve (ICC) or item score curve is drawn for an individual item. The test characteristic curve (TCC) is drawn
#' based on a total test form.
#'
#' @param x An object of class \code{\link{traceline}}.
#' @param item.loc A numeric value indicating that the \emph{n}th item (or the location of item) is plotted.
#' If NULL, the TCC based on a total test form is drawn. Default is NULL.
#' @param score.curve Logical value. If TRUE, item score curve (i.e., a weighted sum of item category probabilities over the item scores) is plotted
#' in a panel. Otherwise, ICCs for all score categories are plotted in separate panels. For a dichotomous item, the item score curve is the same as
#' the ICC of score category 1. Ignored when \code{item.loc = NULL}. Default is FALSE.
#' @param layout.col An integer value indicating the number of columns in the panel when displaying ICCs for an item.
#' @param xlab.text,ylab.text A title for the x and y axes.
#' @param main.text An overall title for the plot.
#' @param lab.size The size of xlab and ylab. Default is 15.
#' @param main.size The size of \code{main.text}. Default is 15.
#' @param axis.size The size of labels along the x and y axes. Default is 15.
#' @param line.color A character string specifying the color for a line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/} for more details
#' about colors used in ggplot2.
#' @param line.size The size of lines. Default is 1.
#' @param strip.size The size of facet labels when ICCs for an item are plotted.
#' @param ... Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.
#'
#' @details All of the plots are drawn using the ggplot2 package.
#' If \code{item.loc = NULL}, the TCC based on the total test form is plotted. In the argument \code{item.loc},
#' a numeric value should be specified to indicate the \emph{n}th item among the total test form. For example,
#' if there are ten items in the test form and the ICCs of the third item should be plotted, then \code{item.loc = 3}.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{traceline}}
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
#' theta <- seq(-3, 3, 0.1)
#'
#' # compute the item category probabilities and item/test
#' # characteristic functions given the theta values
#' x <- traceline(x=test_flex, theta, D=1)
#'
#' # plot TCC based on the toal test form
#' plot(x, item.loc=NULL)
#'
#' # plot ICCs for the first item (dichotomous item)
#' plot(x, item.loc=1, score.curve=FALSE, layout.col=2)
#'
#' # plot item score curve for the first item (dichotomous item)
#' plot(x, item.loc=1, score.curve=TRUE)
#'
#' # plot ICCs for the last item (plolytomous item)
#' plot(x, item.loc=55, score.curve=FALSE, layout.col=2)
#'
#' # plot item score curve for the last item (plolytomous item)
#' plot(x, item.loc=55, score.curve=TRUE)
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
#'
plot.traceline <- function(x, item.loc=NULL,
                           score.curve=FALSE, layout.col=2,
                           xlab.text, ylab.text, main.text, lab.size=15, main.size=15, axis.size=15,
                           line.color, line.size=1, strip.size=12, ...) {

  ##----------------------------------------------------------
  if(!is.null(item.loc)) {
    if(length(item.loc) > 1) {
      stop("A length of 'item.loc' must be 1.", call.=FALSE)
    }
  }

  # extract theta values for x-axis
  theta <- x$theta

  # 1. plot TCCs
  if(is.null(item.loc)) {

    # data manipulation for plotting
    tcc.trace <- x$tcc
    df_tcc <- data.frame(tcc=tcc.trace, theta=theta)

    # plot
    # Set plot conditions
    if(missing(xlab.text)) xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- 'Expected Score'
    if(missing(main.text)) main.text <- 'Test Characteristic Curve'
    if(missing(line.color)) line.color <- "#F8766D" else line.color <- line.color
    max.score <-
      purrr::map_dbl(x$prob.cat, .f=function(k) ncol(k) - 1) %>%
      sum()

    # draw a plot
    p <-
      df_tcc %>%
      ggplot2::ggplot(mapping=aes_string(x="theta", y="tcc")) +
      ggplot2::geom_line(size=line.size, color=line.color, ...) +
      ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
      ggplot2::ylim(0, max.score) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = element_text(size=main.size),
                     axis.title = element_text(size=lab.size),
                     axis.text = element_text(size=axis.size))

  }

  # 2. plot ICCs
  if(!is.null(item.loc)) {

    # check the number of score categories
    cats <- ncol(x$prob.cat[[item.loc]])

    if(!score.curve) {

      # a data.frame including the ICC across all score categories
      icc_df <-
        data.frame(theta=theta, x$prob.cat[[item.loc]]) %>%
        dplyr::rename_all(.funs=function(k) gsub(pattern="score.", replacement="", x=k)) %>%
        reshape2::melt(id.vars="theta", variable.name="score", value.name = "icc")
      icc_df$score <- gsub(pattern="^*", replacement = "Score: ", x=icc_df$score)

      ##-------------------------------------------------------------------------
      # draw ICC plots
      if(missing(xlab.text)) xlab.text <- expression(theta)
      if(missing(ylab.text)) ylab.text <- 'Probability'
      if(missing(main.text)) main.text <- paste0('Item Characteristic Curve: ', names(x$prob.cat[item.loc]))
      if(missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

      p <-
        ggplot2::ggplot(data=icc_df, mapping=aes_string(x="theta", y="icc")) +
        ggplot2::geom_line(color=line.color, size=line.size, ...) +
        ggplot2::labs(title=main.text, x=xlab.text, y=ylab.text) +
        ggplot2::ylim(0, 1) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~score, ncol=layout.col) +
        ggplot2::theme(plot.title = element_text(size=main.size),
                       axis.title = element_text(size=lab.size),
                       axis.text = element_text(size=axis.size)) +
        ggplot2::theme(strip.text.x = element_text(size = strip.size, face='bold'))

    }

    if(score.curve) {

      # data manipulation for plotting
      tcc.trace <- x$icc[, item.loc]
      df_tcc <- data.frame(tcc=tcc.trace, theta=theta)

      # plot
      # Set plot conditions
      if(missing(xlab.text)) xlab.text <- expression(theta)
      if(missing(ylab.text)) ylab.text <- 'Expected Score'
      if(missing(main.text)) main.text <- paste0('Item Score Curve: ', names(x$prob.cat[item.loc]))
      if(missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

      # draw a plot
      p <-
        df_tcc %>%
        ggplot2::ggplot(mapping=aes_string(x="theta", y="tcc")) +
        ggplot2::geom_line(size=line.size, color=line.color, ...) +
        ggplot2::labs(title = main.text, x = xlab.text, y = ylab.text) +
        ggplot2::ylim(0, (cats - 1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = element_text(size=main.size),
                       axis.title = element_text(size=lab.size),
                       axis.text = element_text(size=axis.size))

    }

  }

  print(p)

}





