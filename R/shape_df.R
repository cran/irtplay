#' Create a data.frame of Item Meta Data
#'
#' @description This function creates a data.frame which includes item meta data (e.g., item parameter, categories, models ...) to be
#' used for the IRT model-data fit analysis as well as other analyses.
#'
#' @param par.dc A list containing three vectors of dichotomous item parameters. Namely, the item discrimination (a), item difficulty (b),
#' and item guessing parameters.
#' @param par.py A list containing a vector of polytomous item discrimination (or slope) parameters and a list of polytomous item threshold
#' (or step) parameters. In the list, the argument \code{a} should have a vector of slope parameters and the argument \code{d} should include
#' a list of threshold (or step) parameters. See below for more details.
#' @param item.id A character vector of item IDs. If NULL, an ID is automatically given to each item.
#' @param cats A vector containing the number of score categories for items.
#' @param model A character vector of IRT models corresponding to items. The available IRT models are "1PLM", "2PLM", "3PLM", and "DRM" for
#' dichotomous items, and "GRM" and "GPCM" for polytomous items. Note that "DRM" covers all dichotomous IRT models (i.e, "1PLM", "2PLM", and
#' "3PLM") and "GRM" and "GPCM" represent the graded response model and (generalized) partial credit model, respectively.
#' @param empty.par A logical value to create an empty item meta data. If TRUE, the number of score categories and corresponding IRT models should be specified
#' in the arguments of \code{cats} and \code{model}, respectively. In the empty item meta data, the item slope parameter has a fixed value of 1, the item difficulty
#' (or threshold) parameter has a fixed value of 0, and the item guessing parameter has a fixed value of .2. Default is FALSE.
#'
#' @details For any item where "1PLM" or "2PLM" is specified in \code{model}, the item guessing parameter will be NA. If \code{model} is
#' a vector of \eqn{length = 1}, the specified model is replicated across all items. As in the function \code{\link{simdat}}, it is important
#' to clearly specify \code{cats} according to the order of items in the test form when a data.frame for a mixed-format test needs to be created.
#' See \code{\link{simdat}} for more details about how to specify \code{cats}.
#'
#' When specifying item parameters in \code{par.dc} and \code{par.dc}, keep the order of item parameter types. For example,
#' in the list of \code{par.dc}, the order of items parameters should be the slope, the difficulty, and the guessing parameters.
#'
#' When specifying item parameters in \code{par.dc}, note that in the list of the threshold (or step) parameters, each vector should contain
#' the threshold (or step) parameters for each item. When an item follows the (generalized) partial credit model, the item step parameters
#' are the overall item difficulty (or location) parameter subtracted by the difficulty (or threshold) parameter for each category. Thus, the number
#' of step parameters for item with m categories is m-1 because a step parameter for the first category does not affect the category probabilities.
#'
#' @return This function returns a data.frame.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{test.info}}
#'
#' @examples
#' ## a mixed-item format test form
#' ## with five dichotomous and two polytomous items
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
#' shape_df(par.dc=par.dc, par.py=par.py, cats=cats, model=model)
#'
#' ## an empty item meta data with five dichotomous and two polytomous items
#' # create a numeric vector of score categories for the items
#' cats <- c(2, 4, 3, 2, 5, 2, 2)
#'
#' # create a character vector of IRT models for the items
#' model <- c("1PLM", "GRM", "GRM", "2PLM", "GPCM", "DRM", "3PLM")
#'
#' # create an empty item meta data set
#' shape_df(cats=cats, model=model, empty.par=TRUE)
#'
#' ## an item meta data for a single-item format test form with five dichotomous
#' shape_df(par.dc=par.dc, cats=rep(2, 5), model="DRM")
#'
#'
#' @export
shape_df <- function(par.dc=list(a=NULL, b=NULL, g=NULL), par.py=list(a=NULL, d=NULL), item.id=NULL, cats, model, empty.par=FALSE) {

  model <- toupper(model)

  # only to create an empty item meta data
  if(empty.par) {

    if(missing(cats) | missing(model)) {
      stop("The number of score categories and IRT models must be specified.", call.=FALSE)
    }

    any.dc <- any(cats == 2)
    any.py <- any(cats > 2)

    which.drm <- which(cats == 2)
    which.plm <- which(cats > 2)

    if(any.dc) {
      par.dc <- list(a=rep(1, length(which.drm)), b=rep(0, length(which.drm)), g=rep(NA, length(which.drm)))
      par.dc$g[model[which.drm] == "3PLM" | model[which.drm] == "DRM"] <- 0.2
    } else {
      par.dc <- list(a=NULL, b=NULL, g=NULL)
    }

    if(any.py) {
      par.py <- list(a=rep(1, length(which.plm)), d=vector('list', length(which.plm)))
      for(i in 1:length(which.plm)) {
        par.py$d[[i]] <- rep(0, cats[which.plm[i]] - 1)
      }
    } else {
      par.py <- list(a=NULL, d=NULL)
    }

  }

  nitem <- length(par.dc$b) + length(par.py$d)
  max.cat <- max(cats)
  if(is.null(item.id)) item.id <- paste0("V", 1:nitem)
  if(length(cats) == 1) cats <- rep(cats, nitem)
  if(length(model) == 1) model <- rep(model, nitem)

  if(!all(model %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))) {
    stop("At least one model is not one of '1PLM', '2PLM', '3PLM', 'DRM', 'GRM', and 'GPCM'.", call.=FALSE)
  }

  any.dc <- any(cats == 2)
  any.py <- any(cats > 2)

  if(!any.py) {
    prm_mat <- array(NA, c(nitem, 3))
  } else{
    prm_mat <- array(NA, c(nitem, max.cat))
  }

  # if there are dichotomous items
  if(any.dc) {
    if(is.null(par.dc[[3]])) par.dc[[3]] <- 0
    prm_mat[cats == 2, 1] <- par.dc[[1]]
    prm_mat[cats == 2, 2] <- par.dc[[2]]
    prm_mat[cats == 2, 3] <- par.dc[[3]]
  }

  # if there are polytomous items
  if(any.py) {
    row.py <- which(cats > 2)
    prm_mat[row.py, 1] <- par.py[[1]]
    for(i in 1:length(row.py)) {
      prm_mat[row.py[i], 2:(length(par.py[[2]][[i]])+1)] <- par.py[[2]][[i]]
    }
  }

  # change item guessing parameters into NA when 1PLMs or 2PLMs are specified in 'model' argument
  prm_mat[, 3][model %in% c('1PLM', '2PLM')] <- NA

  if(!any.py) {
    colnames(prm_mat) <- paste0("par.", 1:3)
  } else {
    colnames(prm_mat) <- paste0("par.", 1:max.cat)
  }

  full_df <- data.frame(id=item.id, cats=cats, model=model, prm_mat, stringsAsFactors = FALSE)

  # last check
  if(any(full_df[full_df$cats == 2, 3] %in% c("GRM", "GPCM"))) {
    stop("Dichotomous items must have models among '1PLM', '2PLM', '3PLM', and 'DRM'.", call.=FALSE)
  }
  if(any(full_df[full_df$cats > 2, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
    stop("Polytomous items must have models among 'GRM' and 'GPCM'.", call.=FALSE)
  }


  full_df

}
