#' Run flexMIRT through R
#'
#' @description This function implements flexMIRT (Cai, 2017) to run a model specified in the syntax file of
#' flexMIRT (i.e., *.flexmirt) through R. To run this function, flexMIRT software must be installed in advance.
#' This function will be useful especially when conducting a simulation study using flexMIRT.
#'
#' @param file.syntax A single string or vector containing the file path(s) of a flexmirt syntax file(s) to be run.
#' An example is “C:/Users/Data/irtmodel.flexmirt".
#' @param dir.flex A path of directory where flexMIRT is installed. The path may include a directory name of "flexMIRT3".
#' If NULL, a default path of "C:/Program Files/flexMIRT3" will be used.
#' @param show.output.on.console A logical value to indicate whether to capture the output of the command and show it on the R console.
#' Default is FALSE. See \code{\link[base]{system}}.
#' @param ... Further arguments passed from the function \code{\link[base]{system}}.
#'
#' @details When a path of directory where flexMIRT is installed is provided in the argument \code{dir.flex}, the directory must
#' include following six files:
#' \itemize{
#'   \item WinFlexMIRT.exe
#'   \item FlexMIRT_x64.exe
#'   \item FlexMIRT_x86.exe
#'   \item vpg.dll
#'   \item vpg.licensing.client.dll
#'   \item vpg.licensing.dll
#' }
#'
#' @return output files of flexMIRT
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @references
#' Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
#' Chapel Hill, NC: Vector Psychometric Group.
#'
#' @examples
#'
#' # Emxaples below will run when flexMIRT software is installed
#' # in a default path of "C:/Program Files/flexMIRT3".
#' # Otherwise provide a path where flexMIRT software is installed
#' # in the argument 'dir.flex'.
#'
#' \dontrun{
#' # (1) run a single syntax file
#' # import an example of flexMIRT syntax file to run the item parameter estimation of IRT 3PL model
#' file.syntax <- system.file("extdata", "2PLM_example.flexmirt", package = "irtplay")
#'
#' # run flexMIRT to estimate the item parameters of IRT 3PL model
#' run_flexmirt(file.syntax=file.syntax, dir.flex=NULL, show.output=TRUE)
#'
#' # check the output file
#' out.file <- system.file("extdata", "2PLM_example-prm.txt", package = "irtplay")
#' bring.flexmirt(out.file, type="par")
#'
#' # (2) run multiple syntax files
#' # import two examples of flexMIRT syntax files
#' file.syntax1 <- system.file("extdata", "2PLM_example.flexmirt", package = "irtplay")
#' file.syntax2 <- system.file("extdata", "3PLM_example.flexmirt", package = "irtplay")
#'
#' # run flexMIRT to estimate the item parameters
#' run_flexmirt(file.syntax=c(file.syntax1, file.syntax2), dir.flex=NULL, show.output=FALSE)
#'
#' # check the output file
#' out.file1 <- system.file("extdata", "2PLM_example-prm.txt", package = "irtplay")
#' out.file2 <- system.file("extdata", "3PLM_example-prm.txt", package = "irtplay")
#' bring.flexmirt(out.file1, type="par")
#' bring.flexmirt(out.file2, type="par")
#' }
#'
#'
#' @export
#' @import utils
run_flexmirt <- function(file.syntax, dir.flex=NULL, show.output.on.console=FALSE, ...) {

  if(is.null(dir.flex)) {

    # check if a default directory where flexMIRT is installed exists
    flex_exist <- dir.exists(paths="C:/Program Files/flexMIRT3")

    # provide a warning message when a default directory where flexMIRT is installed does not exist
    if(!flex_exist) {
      stop("A default directory of 'C:/Program Files/flexMIRT3' does not exist.
           Please provide a full path of directory where flexMIRT is installed.", call.=TRUE)
    }

    # a path of flexMIRT execution file
    file.flex <- "C:/Program Files/flexMIRT3/WinFlexMIRT.exe"

  } else {

    # check if the provided directory where flexMIRT is installed exists
    flex_exist <- dir.exists(paths=dir.flex)

    # provide a warning message when the provided directory where flexMIRT is installed does not exist
    if(!flex_exist) {
      stop(paste0("The provided directory of '", dir.flex, "' does not exist.
           Please provide a full path of directory where flexMIRT is installed."), call.=TRUE)
    }

    # a path of flexMIRT execution file
    file.flex <- file.path(dir.flex, "WinFlexMIRT.exe")

  }

  # run flexMIRT
  # when a single syntax files are provided
  if(length(file.syntax) == 1) {
    system(command=paste("\"", file.flex, "\"", " -r ", "\"", file.syntax, "\"", sep=""),
           show.output.on.console=show.output.on.console, ...)
  }

  # when multiple syntax files are provided
  if(length(file.syntax) > 1) {
    # start the progress bar
    pb <- utils::txtProgressBar(min=0, max=100, initial=0, width=50, style=3)
    for(r in 1:length(file.syntax)) {
      system(command=paste("\"", file.flex, "\"", " -r ", "\"", file.syntax[r], "\"", sep=""),
             show.output.on.console = show.output.on.console, ...)

      # modify the progress bar
      info <- sprintf("Parameter Estimation Progress: %d%%", round((r/length(file.syntax))*100))
      utils::setTxtProgressBar(pb, value=r/length(file.syntax)*100, title=info, label=paste0("Estimation ", r, " is done"))
    }

    # closing the progress Bar
    close(pb)
  }

}

