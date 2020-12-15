#' Updated prior (a.k.a. posterior) latent ability distribution 
#'
#' @description This function computes updated prior (a.k.a. posterior) densities of the latent ability distribution given 
#' a prior ability distribution, item parameters, and item response data. 
#'
#' @param x A data frame containing the item metadata (e.g., item parameters, number of categories, models ...). 
#' See \code{\link{irtfit}}, \code{\link{test.info}}, or \code{\link{simdat}}  for more details about the item metadata.
#' This data frame can be easily obtained using the function \code{\link{shape_df}}. 
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param Quadrature A numeric vector of two components specifying the number of quadrature points (in the first component) and
#' the symmetric minimum and maximum values of these points (in the second component). For example, a vector of c(49, 6) indicates 49 rectangular
#' quadrature points over -6 and 6. Default is c(49, 6).
#' @param weights A two-column matrix or data frame containing the quadrature points (in the first column) and the corresponding weights
#' (in the second column) of the latent variable prior distribution. The weights and quadrature points can be easily obtained
#' using the function \code{\link{gen.weight}}. If NULL, a normal prior density is used based on the information provided in the arguments
#' of \code{Quadrature}, \code{group.mean}, and \code{group.var}). Default is NULL.
#' @param group.mean A numeric value to set the mean of latent variable prior distribution. Default is 0. 
#' @param group.var A positive numeric value to set the variance of latent variable prior distribution. Default is 1. 
#' @param missing A value indicating missing values in the response data set. Default is NA.
#'
#' @return This function returns a list containing two internal objects. The first internal object is a data frame with two columns, 
#' where the first column has theta values (nodes) and the second column provides the weights of the posterior latent ability distribution. 
#' The second internal object is a data frame containing the mean, variance, and standard deviation of the distribution.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{shape_df}}, \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}
#'
#' @examples
#' 
#' \donttest{
#' # fit the 2PL model to LSAT6 data
#' (mod.2pl <- est_irt(data=LSAT6, D=1, model="2PLM", cats=2))
#' 
#' # extract the item parameter estimates
#' (x <- getirt(x=mod.2pl, what="par.est"))
#' 
#' # update the standard normal prior deisnty of the ability distribution
#' # using the estimated item parameters
#' (upd_prior <- post_den(x=x, data=LSAT6, D=1, group.mean=0, group.var=1))
#' }
#' 
#' @export
post_den <- function(x, data, D=1, Quadrature=c(49, 6.0), weights=NULL, 
                     group.mean=0, group.var=1, missing=NA) {
  
  
  # give column names
  x <- data.frame(x)
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))
  
  # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
  if(ncol(x[, -c(1, 2, 3)]) == 2) {
    x <- data.frame(x, par.3=NA)
  }
  
  # clear the item metadata set
  x <- back2df(metalist2(x))
  cats <- x[, 2]
  model <-
    as.character(x[, 3]) %>%
    toupper()
  
  # check whether included data are correct
  if(nrow(x) != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)
  
  # consider DRM as 3PLM
  if("DRM" %in% model) {
    model[model == "DRM"] <- "3PLM"
    x$model <- model
    memo <- "All 'DRM' items were considered as '3PLM' items during the item parameter estimation. \n"
    warning(memo, call.=FALSE)
  }
  
  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }
  
  # transform a data set to matrix
  data <- data.matrix(data)
  
  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))
  
  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)
  if(length(loc_allmiss) > 0L) {
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "), " has/have no item response data. \n")
    stop(memo2, call.=TRUE)
  }
  
  # check the total number of examinees
  nstd <- nrow(data)
  
  # create initial weights of prior ability distribution when it is not specified
  if(is.null(weights)) {
    # create quadrature points
    quadpt <- seq(-Quadrature[2], Quadrature[2], length.out=Quadrature[1])
    
    # create the data.frame containing the quadurature points and weights
    weights <- gen.weight(dist="norm", mu=group.mean, sigma=sqrt(group.var), theta=quadpt)
    
  } else {
    quadpt <- weights[, 1]
  }
  
  # factorize the response values and create a list of item responses across items
  resp <- purrr::map2(.x=data.frame(data), .y=cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))
  
  # calculate the score categories for each examinee
  std.id <- 1:nstd
  freq.cat <- purrr::map(.x=resp, .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))
  
  # transform the score category data.frame to a matrix
  freq.cat <- purrr::map(.x=freq.cat,
                         .f=function(k) data.matrix(k) %>%
                           unname())
  
  # delete 'resp' object
  rm(resp, envir=environment(), inherits = FALSE)
  
  # divide the data set for the mixed-item format
  datlist <- divide_data(data=data, cats=cats, freq.cat=freq.cat)
  data1_drm <- datlist$data1_drm
  data2_drm <- datlist$data2_drm
  data_plm <- datlist$data_plm
  
  # delete 'datlist' object
  rm(datlist, envir=environment(), inherits = FALSE)
  
  # listrize the data.frame
  meta <- metalist2(x)
  
  # compute the likelihood and log-likelihood matrix
  L_LL <- likelihood(meta, data1_drm=data1_drm, data2_drm=data2_drm, data_plm=data_plm, theta=weights[, 1], D=D)
  likehd <- L_LL$L
  
  # posterior distribution
  post_dist <- posterior(likehd, weights)
  
  # compute the expected frequency of scores categories across all items
  # this is the conditional expectation of item responses with respect to posterior likelihood distribution
  post_tp <- t(post_dist)
  freq.exp <- purrr::map(.x=freq.cat, .f=function(x) post_tp %*% x)
  
  # update the frequency of the distribution
  prior_freq <- unname(colSums(post_dist))
  
  # normalize the updated freqency distribution to obtain the density distribution
  prior_dense <- prior_freq /nstd
  
  # updated density distribution
  pop_dense <- data.frame(theta=quadpt, weight=prior_dense)
  
  # mean and sd of the updated distribution
  moments <- cal_moment(node=pop_dense$theta, weight=pop_dense$weight)
  sigma <- sqrt(moments[2])
  
  # data.frame for the population density parameter estimates
  group.par <- data.frame(rbind(c(moments, sigma)))
  colnames(group.par) <- c("mu", "sigma2", "sigma")
  
  # return
  rst <- list(post.den=pop_dense, group.par=group.par)
  rst
  
  
}