#' Simulated single-item format CAT Data
#'
#' This data set contains an item pool information, response data, and examinee's ability estimates.
#'
#' @usage simCAT_DC
#'
#' @format This data includes a list of length three. The first internal object is a data.frame of
#' the item pool consisting of 100 dichotomous items. The item parameters of the first 90 items were generated with
#' the IRT 2PL model and calibrated with the same model. However, the item parameters of the last 10 items were
#' generated with the IRT 3PL model but calibrated with the IRT 2PL model. The second internal object is the response
#' data set including a sparse response data set of 10,000 examinees for the items in the item pool.
#' The third internal object is the examinee's ability estimates for 10,000 examinees.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
"simCAT_DC"
