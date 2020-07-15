# Transform a data.frame without anchor information to a list
#' @import purrr
metalist2 <- function(x) {

  # chagne all factor variables into character variables
  x <- purrr::modify_if(x, is.factor, as.character)
  x[, 3] <- toupper(x[, 3])

  modelGood <- all(x[,3] %in% c('1PLM', '2PLM', '3PLM', 'DRM', 'GRM', 'GPCM'))
  catsGood <- all(x[,2] >= 1)
  if(!modelGood) stop("At least, one of model is mis-specified. Available models are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM")
  if(!catsGood) stop("At least, one of score category is less than 2. Score category should be greater than 1")

  # check the maximum of category numbers
  max.cat <- max(x[, 2])

  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }

  # give location numbers to each item
  x$loc <- 1:nrow(x)

  # give column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 4)), "loc")

  # create an empty list to containl all meta information
  meta <- list()

  # when dichotomous items are included in a test
  if(2 %in% x[, 2]) {

    # drm.df <-
    #   x %>%
    #   dplyr::filter(.data$cats <= 2) %>%
    #   dplyr::select(.data$id, .data$cats, .data$model, .data$par.1, .data$par.2, .data$par.3, .data$loc) %>%
    #   dplyr::rename(a="par.1", b="par.2", g="par.3")
    drm.df <- x[x$cats <= 2, ]
    drm.df <- drm.df[, c(1:6, ncol(drm.df))]
    colnames(drm.df) <- c("id", "cats", "model", "a", "b", "g", "loc")

    # give zero values to NA in the guessing parameter column
    drm.df$g[is.na(drm.df$g)] <- 0

    # give zero values to the guessing parameter column for 1PLM and 2PLM items
    drm.df$g[drm.df$model %in% c("1PLM", "2PLM")] <- 0

    # listrize a dafa.frame
    drm.list <- as.list(drm.df)

    # contain the list in the meta list
    meta$drm <- drm.list

  }

  # when polytomous items are included in a test
  if(max.cat > 2) {

    # plm.df <-
    #   x %>%
    #   dplyr::filter(.data$cats > 2) %>%
    #   dplyr::select(.data$id, .data$cats, .data$model, .data$par.1, .data$loc) %>%
    #   dplyr::rename(a="par.1")
    plm.df <- x[x$cats > 2, ]
    plm.df.a <- plm.df[, c(1:4, ncol(plm.df))]
    colnames(plm.df.a) <- c("id", "cats", "model", "a", "loc")

    # manipulate step parameters to make them as a list
    # d <-
    #   x %>%
    #   dplyr::filter(.data$cats > 2) %>%
    #   dplyr::select(paste0("par.", 2:max.cat)) %>%
    #   t() %>%
    #   data.frame() %>%
    #   purrr::map(data.frame) %>%
    #   purrr::map(tidyr::drop_na) %>%
    #   purrr::map(unlist) %>%
    #   purrr::map(unname)
    d <- t(plm.df[, paste0("par.", 2:max.cat)])
    d <- data.frame(d)
    d <- purrr::map(.x=d, .f=function(x) x[!is.na(x)])
    names(d) <- paste0("d", 1:nrow(plm.df))

    # compare the number of score categories and the number of item parameters
    num1 <- purrr::map_dbl(.x=d, .f=function(k) length(k) + 1)
    if(!all(plm.df$cats == num1)) {
      stop("The number of item parameters and the numeber of score categories do not match at least for one polytomous item.", call.=FALSE)
    }

    # listrize a dafa.frame
    plm.list <- as.list(plm.df.a)
    plm.list <- c(plm.list, d=list(d))
    plm.list <- plm.list[c("id","cats", "model", "a", "d", "loc")]

    meta$plm <- plm.list

  }

  meta

}

