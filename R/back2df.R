# a function "back2df"
# This function is used to transfrom a list back to a data.frame of item meta data
#' @import purrr
#' @import dplyr
#' @importFrom rlang .data
back2df <- function(meta) {

  # when there are dichotomous items
  if(!is.null(meta$drm)){

    drm_df <- data.frame(meta$drm)
    colnames(drm_df) <- c("id", "cats", "model", paste0("par.", 1:3), "loc")

    # chagne all factor variables into character variables
    drm_df <- purrr::modify_if(drm_df, is.factor, as.character)

  } else {

    drm_df <- NULL

  }

  # when there are polytomous items
  if(!is.null(meta$plm)) {

    # check the maximum of category numbers
    max.cat <- max(meta$plm$cats)

    # change a list of step parameters to a matrix
    d.mat <- bind.fill(meta$plm$d, type="rbind")

    pre_df <- data.frame(meta$plm[-5])
    plm_df <- data.frame(pre_df[, 1:4], d.mat, pre_df[, -c(1:4)])
    colnames(plm_df) <- c("id", "cats", "model", paste0("par.", 1:max.cat), "loc")

    # chagne all factor variables into character variables
    plm_df <- purrr::modify_if(plm_df, is.factor, as.character)

  } else {

    plm_df <- NULL

  }

  # rbind data.frames
  meta_df <- dplyr::bind_rows(drm_df, plm_df)

  # relocation rows and columns
  meta_df <-
    meta_df %>%
    dplyr::arrange(.data$loc) %>%
    dplyr::select(.data$id, .data$cats, .data$model, dplyr::starts_with("par."), dplyr::everything(), -.data$loc)

  ## return a data.frame
  meta_df

}




