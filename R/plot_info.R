#' Plot Item and Test Information Functions
#'
#' @description This function plots item or test information function given a specified theta values.
#'
#' @param x x An object of class \code{\link{test.info}}.
#' @param item.loc A vector of numeric values indicating that the item information functions of the \emph{n}th items
#' (or the location of items in a test form) are plotted. If NULL, the test information function for the total test form is drawn.
#' Default is NULL.
#' @param xlab.text,ylab.text A title for the x and y axes.
#' @param main.text An overall title for the plot.
#' @param lab.size The size of xlab and ylab. Default is 15.
#' @param main.size The size of \code{main.text}. Default is 15.
#' @param axis.size The size of labels along the x and y axes. Default is 15.
#' @param line.color A character string specifying a color for the line. See \url{http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/}
#' for more details about colors used in ggplot2.
#' @param line.size The size of lines. Default is 1.
#' @param layout.col An integer value indicating the number of columns in the panel when displaying the item information functions of
#' the multiple items. Default is 4.
#' @param strip.size The size of facet labels when the item information functions of the multiple items are drawn.
#' @param ... Further arguments passed from the function \code{geom_line()} in the \pkg{ggplot2} package.
#'
#' @details All of the plots are drawn using the ggplot2 package.
#' The object of class \code{\link{test.info}} can be obtained from the function \code{\link{test.info}}.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{test.info}}
#'
#' @examples
#' ## the use of a "-prm.txt" file obtained from a flexMIRT
#' # import the "-prm.txt" output file from flexMIRT
#' flex_prm <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item metadata
#' test_flex <- bring.flexmirt(file=flex_prm, "par")$Group1$full_df
#'
#' # set theta values
#' theta <- seq(-4, 4, 0.1)
#'
#' # compute item and test information values given the theta values
#' x <- test.info(x=test_flex, theta=theta, D=1)
#'
#' # draw a plot of the test information function
#' plot(x)
#'
#' # draw a plot of the item information function for the second item
#' plot(x, item.loc=2)
#'
#' # draw a plot of the item information function for the mutiple items
#' plot(x, item.loc=1:8)
#'
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
plot.test.info <- function(x, item.loc=NULL, xlab.text, ylab.text, main.text, lab.size=15, main.size=15,
                           axis.size=15, line.color, line.size=1, layout.col=4, strip.size=12, ...) {


  # 1. plot test infomation
  if(is.null(item.loc)) {

    # data manipulation for plotting
    df_info <- data.frame(theta=x$theta, info = x$testInfo)

    # plot
    # Set plot conditions
    if(missing(xlab.text)) xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- 'Information'
    if(missing(main.text)) main.text <- 'Test Information'
    if(missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

    # draw a plot
    p <- 
      df_info %>% 
      ggplot(mapping=aes_string(x="theta", y="info")) +
      geom_line(size=line.size, color=line.color, ...) +
      ggplot2::theme_bw() +
      labs(title = main.text, x = xlab.text, y = ylab.text) +
      theme(plot.title = element_text(size=main.size),
            axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size))

  }

  # 2. plot item information
  if(!is.null(item.loc)) {

    # data manipulation for plotting
    df_info <- 
      data.frame(t(x$itemInfo), theta=x$theta) %>%
      stats::setNames(nm=c(1:nrow(x$itemInfo), "theta")) %>%
      reshape2::melt(variable.name="item", id.vars="theta", value.name="info")
    df_info$item <- as.numeric(df_info$item)
    df_info <- dplyr::filter(df_info, .data$item %in% item.loc)

    # plot
    # Set plot conditions
    if(missing(xlab.text)) xlab.text <- expression(theta)
    if(missing(ylab.text)) ylab.text <- 'Information'
    if(missing(main.text)) main.text <- 'Item Information'
    if(missing(line.color)) line.color <- "#F8766D" else line.color <- line.color

    p <- 
      df_info %>% 
      ggplot2::ggplot(mapping=aes_string(x="theta", y="info")) +
      ggplot2::geom_line(size=line.size, color=line.color, ...) +
      ggplot2::theme_bw() +
      labs(title = main.text, x = xlab.text, y = ylab.text) +
      theme(plot.title = element_text(size=main.size),
            axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size)) +
      facet_wrap(~item, ncol=layout.col) +
      theme(strip.text.x = element_text(size = strip.size, face = 'bold'))

  }

  p

}


