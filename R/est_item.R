#' Fixed effect item parameter calibration
#'
#' @description This function performs the fixed effects (or fixed ability) item parameter calibration, called
#' Method A (Stocking, 1988), which is the maximum likelihood estimation of item parameters given the ability
#' estimates (Baker & Kim, 2004; Ban, Hanson, Wang, Yi, & Harris, 2001; Stocking, 1988). Also, this could be considered as a special
#' type of the joint maximum likelihood estimation where only one cycle of parameter estimation is implemented
#' given the ability estimates (Birnbaum, 1968). Method A is one of potentially useful online item calibration methods
#' for computerized adaptive testing (CAT) to put the parameter estimates of pretest items on the same scale of
#' operational item parameter estimates and recalibrate the operational items to evaluate the parameter drifts
#' of the operational items (Chen & Wang, 2016; Stocking, 1988).
#'
#' @param x A data.frame containing the item meta data. This meta data is necessary to obtain the information of
#' each item (i.e., number of score categories and IRT model) to be calibrated. You can easily create an empty
#' item meta data using the function \code{\link{shape_df}}. When \code{use.startval = TRUE}, the item parameters
#' specified in the item meta data are used as the starting values for the item parameter estimation.
#' If \code{x = NULL}, the arguments of \code{model} and \code{cats} must be specified. See \code{\link{irtfit}},
#' \code{\link{test.info}} or \code{\link{simdat}} for more details about the item meta data. See below for details.
#' @param data A matrix containing examinees' response data for the items in the argument \code{x}. A row and column indicate
#' the examinees and items, respectively.
#' @param score A vector of examinees' ability estimates. Length of the vector must be the same as the number of rows in the
#' response data set.
#' @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
#' Default is 1.
#' @param model A vector of character strings indicating what IRT model is used to calibrate each item. Available IRT models are
#' "1PLM", "2PLM", "3PLM", and "DRM" for dichotomous items, and "GRM" and "GPCM" for polytomous items. "GRM" and "GPCM" represent the graded
#' response model and (generalized) partial credit model, respectively. Note that "DRM" is considered as "3PLM" in this function.
#' If a single character of the IRT model is specified, that model will be recycled across all items. This information is only required
#' when \code{x = NULL}.
#' @param cats A numeric vector specifying the number of score categories for each item. For example, a dichotomous
#' item has two score categories. If a single numeric value is specified, that value will be recycled across all items. If NULL and all items
#' are binary items (i.e., dichotomous items), it assumes that all items have two score categories. This information is only required
#' when \code{x = NULL}.
#' @param fix.a.1pl A logical value. If TRUE, the slope parameters of the 1PLM items are fixed to a specific value specified in the argument
#' \code{a.val.1pl}. Otherwise, the slope parameters of all 1PLM items are constrained to be equal and estimated. Default is FALSE.
#' @param fix.a.gpcm A logical value. If TRUE, the GPCM items are calibrated with the partial credit model and the slope parameters of
#' the GPCM items are fixed to a specific value specified in the argument \code{a.val.gpcm}. Otherwise, the slope parameter of each GPCM item
#' is estimated. Default is FALSE.
#' @param fix.g A logical value. If TRUE, the guessing parameters of the 3PLM items are fixed to a specific value specified in the argument
#' \code{g.val}. Otherwise, the guessing parameter of each 3PLM item is estimated. Default is FALSE.
#' @param a.val.1pl A numeric value. This value is used to fixed the slope parameters of the 1PLM items.
#' @param a.val.gpcm A numeric value. This value is used to fixed the slope parameters of the GPCM items.
#' @param g.val A numeric value. This value is used to fixed the guessing parameters of the 3PLM items.
#' @param use.aprior A logical value. If TRUE, a prior distribution for the slope parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.bprior A logical value. If TRUE, a prior distribution for the difficulty (or threshold) parameters is used for the parameter calibration
#' across all items. Default is FALSE.
#' @param use.gprior A logical value. If TRUE, a prior distribution for the guessing parameters is used for the parameter calibration
#' across all 3PLM items. Default is TRUE.
#' @param aprior A list containing the information of the prior distribution for item slope parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param bprior A list containing the information of the prior distribution for item difficulty (or threshold) parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param gprior A list containing the information of the prior distribution for item guessing parameters. Three probability distributions
#' of Beta, Log-normal, and Normal distributions are available. In the list, a character string of the distribution name must be specified
#' in the first internal argument and a vector of two numeric values for the two parameters of the distribution must be specified in the
#' second internal argument. Specifically, when Beta distribution is used, "beta" should be specified in the first argument. When Log-normal
#' distribution is used, "lnorm" should be specified in the first argument. When Normal distribution is used, "norm" should be specified
#' in the first argument. In terms of the two parameters of the three distributions, see \code{dbeta()}, \code{dlnorm()},
#' and \code{dnorm()} in the \pkg{stats} package for more details.
#' @param missing A value indicating missing values in the response data set. Default is NA.
#' @param use.startval A logical value. If TRUE, the item parameters provided in the item meta data (i.e., the argument \code{x}) are used as
#' the starting values for the item parameter estimation. Otherwise, internal starting values of this function are used. Default is FALSE.
#' @param control A list of control parameters to be passed to the optimization function of \code{nlminb()} in the \pkg{stats} package. The control parameters
#' set the conditions of the item parameter estimation process such as the maximum number of iterations. See \code{nlminb()} in the \pkg{stats} package for details.
#'
#' @details In most cases, the function \code{\link{est_item}} will return successfully converged item parameter estimates using
#' the default internal starting values. However, if there is a convergence problem in the calibration, one possible solution is using
#' different starting values. When the item parameter values are specified in the item meta data (i.e., the argument \code{x}), those values
#' can be used as the starting values for the item parameter calibration by setting \code{use.startval = TRUE}.
#'
#' @return This function returns an object of class \code{\link{est_item}}. Within this object, several internal objects are contained such as:
#' \item{estimates}{A data.frame containing both the item parameter estimates and the corresponding standard errors of estimates.}
#' \item{par.est}{A data.frame containing the item parameter estimates.}
#' \item{se.est}{A data.frame containing the standard errors of the item parameter estimates. Note that the standard errors are estimated using
#' observed information functions.}
#' \item{loglikelihood}{A sum of the log-likelihood values of the complete data set across all estimated items.}
#' \item{data}{A data.frame of the examinees' response data set.}
#' \item{score}{A vector of the examinees' ability values used as the fixed effects.}
#' \item{scale.D}{A scaling factor in IRT models.}
#' \item{convergence}{A string indicating the convergence status of the item parameter estimation.}
#' \item{nitem}{A total number of items included in the response data.}
#' \item{deleted.item}{The items which have no item response data. Those items are excluded from the item parameter estimation.}
#' \item{npar.est}{A total number of the estimated parameters.}
#' \item{n.response}{An integer vector indicating the number of item responses for each item used to estimate the item parameters.}
#' \item{TotalTime}{Time (in seconds) spent for total compuatation.}
#'
#' The internal objects can be easily extracted using the function \code{\link{getirt}}.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}, \code{\link{test.info}}, \code{\link{simdat}}, \code{\link{shape_df}}, \code{\link{sx2_fit}},
#' \code{\link{traceline.est_item}}, \code{\link{getirt}}
#'
#' @references
#' Baker, F. B., & Kim, S. H. (2004). \emph{Item response theory: Parameter estimation techniques.} CRC Press.
#'
#' Ban, J. C., Hanson, B. A., Wang, T., Yi, Q., & Harris, D., J. (2001) A comparative study of on-line pretest item calibration/scaling methods
#' in computerized adaptive testing. \emph{Journal of Educational Measurement, 38}(3), 191-212.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F. M. Lord & M. R. Novick (Eds.),
#' \emph{Statistical theories of mental test scores} (pp. 397-479). Reading, MA: Addison-Wesley.
#'
#' Chen, P., & Wang, C. (2016). A new online calibration method for multidimensional computerized adaptive testing.
#' \emph{Psychometrika, 81}(3), 674-701.
#'
#' Stocking, M. L. (1988). \emph{Scale drift in on-line calibration} (Research Rep. 88-28). Princeton, NJ: ETS.
#'
#' @examples
#' ## import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # select the item meta data
#' x <- bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' # modify the item meta data so that some items follow 1PLM, 2PLM and GPCM
#' x[c(1:3, 5), 3] <- "1PLM"
#' x[c(1:3, 5), 4] <- 1
#' x[c(1:3, 5), 6] <- 0
#' x[c(4, 8:12), 3] <- "2PLM"
#' x[c(4, 8:12), 6] <- 0
#' x[54:55, 3] <- "GPCM"
#'
#' # generate examinees' abilities from N(0, 1)
#' set.seed(23)
#' score <- rnorm(500, mean=0, sd=1)
#'
#' # simulate the response data
#' data <- simdat(x=x, theta=score, D=1)
#'
#' \donttest{
#' # 1) item parameter estimation: constrain the slope parameters of the 1PLM to be equal
#' (mod1 <- est_item(x, data, score, D=1, fix.a.1pl=FALSE, use.gprior=TRUE,
#'                   gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
#' summary(mod1)
#'
#' # extract the item parameter estimates
#' getirt(mod1, what="par.est")
#'
#' # 2) item parameter estimation: fix the slope parameters of the 1PLM to 1
#' (mod2 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, a.val.1pl=1, use.gprior=TRUE,
#'                   gprior=list(dist="beta", params=c(5, 17)), use.startval=FALSE))
#' summary(mod2)
#'
#' # extract the standard error estimates
#' getirt(mod2, what="se.est")
#'
#' # 3) item parameter estimation: fix the guessing parameters of the 3PLM to 0.2
#' (mod3 <- est_item(x, data, score, D=1, fix.a.1pl=TRUE, fix.g=TRUE, a.val.1pl=1, g.val=.2,
#'                   use.startval=FALSE))
#' summary(mod3)
#'
#' # extract both item parameter and standard error estimates
#' getirt(mod2, what="estimates")
#'
#' }
#'
#' @import purrr
#' @import dplyr
#'
#' @export
#'
est_item <- function(x=NULL, data, score, D=1, model=NULL, cats=NULL, fix.a.1pl=FALSE, fix.a.gpcm=FALSE, fix.g=FALSE,
                     a.val.1pl=1, a.val.gpcm=1, g.val=.2, use.aprior=FALSE, use.bprior=FALSE, use.gprior=TRUE,
                     aprior=list(dist="lnorm", params=c(0, 0.5)), bprior=list(dist="norm", params=c(0.0, 1.0)),
                     gprior=list(dist="beta", params=c(5, 17)), missing=NA, use.startval=FALSE,
                     control=list(eval.max=500, iter.max=500)) {

  # check start time
  start.time <- Sys.time()

  # match.call
  cl <- match.call()

  ##-------------------------------------------------------------------------------------------------------
  ## 1. preperation of item parameter estimation
  cat("Starting...", '\n')

  # check if the starting values are available
  if(use.startval & is.null(x)) {
    stop("To use starting values for item parameter estimation, the item meta data must be specified in the argument 'x'.", call.=FALSE)
  }

  # extract information about the number of score cetegories and models
  cat("Parsing input...", '\n')
  if(!is.null(x)) {
    # give column names
    x <- data.frame(x)
    colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))
    cats <- x[, 2]

    # add par.3 column when there is no par.3 column (just in case that all items are 2PLMs)
    if(ncol(x[, -c(1, 2, 3)]) == 2) {
      x <- data.frame(x, par.3=NA)
    }

    # clear the item meta data set
    x <- back2df(metalist2(x))
    model <-
      as.character(x[, 3]) %>%
      toupper()
  } else {
    model <- toupper(model)
    if(length(model) == 1) {
      model <- rep(model, ncol(data))
    }
    if(is.null(cats)) {
      if(all(model %in% c("1PLM", "2PLM", "3PLM", "DRM"))) {
        cats <- rep(2, ncol(data))
      } else {
        stop("The number of score categories for the items should be specified in the argument 'cats'.", call.=FALSE)
      }
    }
    if(length(cats) == 1) {
      cats <- rep(cats, ncol(data))
    }
    x <- shape_df(cats=cats, model=model, empty.par=TRUE)
  }


  # check whether included data are correct
  if(nrow(x) != ncol(data)) stop("The number of items included in 'x' and 'data' must be the same.", call.=FALSE)

  # consider DRM as 3PLM
  if("DRM" %in% model) {
    model[model == "DRM"] <- "3PLM"
    x$model <- model
    memo <- "All 'DRM' items were considered as '3PLM' items during the item parameter estimation."
    warning(memo, call.=FALSE)
  }

  # recode missing values
  if(!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform a data set to data.frame
  data <- data.frame(data)

  # check the total number of item in the response data set
  nitem <- ncol(data)

  # check the number of item responses across all items
  n.resp <- colSums(!is.na(data))

  # check the items which have all missing responses
  loc_allmiss <- which(n.resp == 0L)

  # delete items which have all missing responses from the item meta data set
  if(length(loc_allmiss) > 0L) {
    loc_nomiss <- which(n.resp > 0L)
    x_pre <- x
    x <- x[-loc_allmiss, ]
    model <- x[, 3]
    cats <- x[, 2]
    data <- data[, loc_nomiss]
    memo2 <- paste0(paste0("item ", loc_allmiss, collapse = ", "),
                    " is/are excluded in the item parameter estimation because the item(s) has/have no item response data.")
    warning(memo2, call.=FALSE)
  } else {
    loc_allmiss <- NULL
  }

  # group parameter estimates
  group.par <- c(mu=mean(score), sigma=stats::sd(score))

  # transform scores to a vector form
  if(is.matrix(score) | is.data.frame(score)) {
    score <- as.numeric(data.matrix(score))
  }

  # copy scores to theta values
  theta <- score

  # find the location of 1PLM items in which slope parameters should be constrained to be equal
  # also, find the location of other items
  if("1PLM" %in% model & !fix.a.1pl) {
    loc_1p_const <- which(model == "1PLM")
    loc_else <- which(model != "1PLM")
  } else {
    loc_1p_const <- NULL
    loc_else <- 1:nrow(x)
  }

  # record the original location of item parameters to be estimated, and
  # the relocated position of item parameters when computing
  # the variance-covariance matrix of item parameter estimates
  param_loc <- parloc(x=x, loc_1p_const=loc_1p_const, loc_else=loc_else,
                      fix.a.1pl=fix.a.1pl, fix.a.gpcm=fix.a.gpcm, fix.g=fix.g)

  # factorize the response values
  resp <- purrr::map2(.x=data, .y=cats, .f=function(k, m) factor(k, levels=(seq_len(m) - 1)))

  # check the total number of examinees
  nstd <- nrow(data)

  # calculate the score categories for each examinee
  std.id <- 1:nstd
  freq.cat <- purrr::map(.x=resp,
                         .f=function(k) stats::xtabs(~ std.id + k, na.action=stats::na.pass, addNA = FALSE))

  # transform the score category data.frame to a matrix
  freq.cat <- purrr::map(.x=freq.cat,
                         .f=function(k) unname(data.matrix(k)))

  # delete 'resp' object
  rm(resp, envir=environment(), inherits = FALSE)

  # extract theta (score) values for each item
  # score <- purrr::map(.x=1:length(freq.cat), ~as.numeric(score))

  ##-------------------------------------------------------------------------------------------------------
  ## 2. item parameter estimation
  ##---------------------------------------------------------------
  # Item parameter estimation
  # create empty vectors to contain results
  est_par <- NULL
  est_se <- NULL
  loc_items <- NULL
  convergence <- NULL
  noconv_items <- NULL
  objective <- NULL

  # estimation
  cat("Estimating item parameters...", '\n')

  # listrize the item meta data to use the starting values
  meta <- metalist2(x)

  # the dichotomous items: 1PLM with constrained slope values
  if("1PLM" %in% model & !fix.a.1pl) {

    # check the number of 1PLM items
    n.1PLM <- length(loc_1p_const)

    # prepare input files to estimate the 1PLM item parameters
    f_i <- r_i <- array(0, c(nstd, n.1PLM))
    for(k in 1:n.1PLM) {
      f_i[, k] <- rowSums(freq.cat[loc_1p_const][[k]])
      r_i[, k] <- freq.cat[loc_1p_const][[k]][, 2]
    }

    # check the starting values
    if(use.startval) {
      pos_1p_const <- which(meta$drm$loc %in% loc_1p_const)
      a.val <- meta$drm$a[pos_1p_const][1]
      b.val <- meta$drm$b[pos_1p_const]
      startval <- c(a.val, b.val)
    } else {
      startval <- NULL
    }

    # parameter estimation
    est <- estimation(f_i=f_i, r_i=r_i, theta=theta, model="1PLM", D=D, fix.a.1pl=FALSE, n.1PLM=n.1PLM,
                      aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior,
                      control=control, startval=startval)

    # extract the results
    # item parameter estimates
    a <- est$pars[1]
    b <- est$par[-1]
    pars <- purrr::map(1:n.1PLM, .f=function(x) c(a, b[x], NA))
    est_par <- c(est_par, pars)

    # standrad errors
    a.se <- est$se[1]
    b.se <- est$se[-1]
    pars.se <- purrr::map(1:n.1PLM, .f=function(x) c(NA, b.se[x], NA))
    pars.se[[1]][1] <- a.se
    est_se <- c(est_se, pars.se)

    # convergence indicator
    convergence <- c(convergence, est$convergence)
    if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_1p_const)

    # negative loglikelihood value
    objective <- c(objective, est$objective)

  }

  # all other items
  if(length(loc_else) >= 1) {
    for(i in 1:length(loc_else)) {

      # prepare information to estimate item parameters
      mod <- model[loc_else][i]
      score.cat <- cats[loc_else][i]

      # in case of a dichotomous item
      if(score.cat == 2) {
        f_i <- rowSums(freq.cat[loc_else][[i]])
        r_i <- freq.cat[loc_else][[i]][, 2]

        # check the starting values
        if(use.startval) {
          pos_item <- which(meta$drm$loc == loc_else[i])
          a.val <- meta$drm$a[pos_item]
          b.val <- meta$drm$b[pos_item]
          g.val <- meta$drm$g[pos_item]
          if(mod == "1PLM") {
            startval <- b.val
          }
          if(mod == "2PLM") {
            startval <- c(a.val, b.val)
          }
          if(mod == "3PLM") {
            if(fix.g) {
              startval <- c(a.val, b.val)
            } else {
              startval <- c(a.val, b.val, g.val)
            }
          }
        } else {
          startval <- NULL
        }

        # parameter estimation
        est <- estimation(f_i=f_i, r_i=r_i, theta=theta, model=mod, D=D,
                          fix.a.1pl=ifelse(mod == "1PLM", TRUE, FALSE),
                          fix.g=fix.g, a.val.1pl=a.val.1pl, g.val=g.val, n.1PLM=NULL,
                          aprior=aprior, bprior=bprior, gprior=gprior,
                          use.aprior=use.aprior, use.bprior=use.bprior, use.gprior=use.gprior,
                          control=control, startval=startval)

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "1PLM", a.val.1pl, est$pars[1])
        b <- ifelse(mod == "1PLM", est$pars[1], est$pars[2])
        g <- ifelse(mod == "3PLM", ifelse(fix.g, g.val, est$pars[3]), NA)
        pars <- c(a, b, g)
        est_par <- c(est_par, list(pars))

        # standrad errors
        a.se <- ifelse(mod == "1PLM", NA, est$se[1])
        b.se <- ifelse(mod == "1PLM", est$se[1], est$se[2])
        g.se <- ifelse(mod == "3PLM", ifelse(fix.g, NA, est$se[3]), NA)
        pars.se <- c(a.se, b.se, g.se)
        est_se <- c(est_se, list(pars.se))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative loglikelihood value
        objective <- c(objective, est$objective)

      }

      # in case of a polytomous item
      if(score.cat > 2) {
        r_i <- freq.cat[loc_else][[i]]

        # check the starting values
        if(use.startval) {
          pos_item <- which(meta$plm$loc == loc_else[i])
          a.val <- meta$plm$a[pos_item]
          d.val <- meta$plm$d[[pos_item]]
          if(mod == "GRM") {
            startval <- c(a.val, d.val)
          }
          if(mod == "GPCM") {
            if(fix.a.gpcm) {
              startval <- d.val
            } else {
              startval <- c(a.val, d.val)
            }
          }
        } else {
          startval <- NULL
        }

        # parameter estimation
        est <- estimation(r_i=r_i, theta=theta, model=mod, cats=score.cat, D=D,
                          fix.a.gpcm=ifelse(mod == "GPCM", fix.a.gpcm, FALSE), a.val.gpcm=a.val.gpcm, n.1PLM=NULL,
                          aprior=aprior, bprior=bprior, use.aprior=use.aprior, use.bprior=use.bprior,
                          control=control, startval=startval)

        # extract the results
        # item parameter estimates
        a <- ifelse(mod == "GRM", est$pars[1], ifelse(fix.a.gpcm, a.val.gpcm, est$pars[1]))
        if(mod == "GRM") {
          bs <- est$pars[-1]
        } else {
          if(fix.a.gpcm) {
            bs <- est$pars
          } else{
            bs <- est$pars[-1]
          }
        }
        pars <- c(a, bs)
        est_par <- c(est_par, list(pars))

        # standrad errors
        a.se <- ifelse(mod == "GRM", est$se[1], ifelse(fix.a.gpcm, NA, est$se[1]))
        if(mod == "GRM") {
          bs.se <- est$se[-1]
        } else {
          if(fix.a.gpcm) {
            bs.se <- est$se
          } else{
            bs.se <- est$se[-1]
          }
        }
        pars.se <- c(a.se, bs.se)
        est_se <- c(est_se, list(pars.se))

        # convergence indicator
        convergence <- c(convergence, est$convergence)
        if(est$convergence > 0L) noconv_items <- c(noconv_items, loc_else[i])

        # negative loglikelihood value
        objective <- c(objective, est$objective)

      }
    }
  }

  ##---------------------------------------------------------------
  # check the convergence of parameter estimation
  if(sum(convergence) == 0L) {
    note <- "All item parameters were successfully converged."
  } else {
    note <- paste0(paste0("Item ", sort(noconv_items), collapse=", "), " was(were) not successfully converged.")
  }

  ##---------------------------------------------------------------
  # compute the sum of loglikelihood values
  llike <- -sum(objective)

  ##---------------------------------------------------------------
  # arrange the estimated item parameters and standard errors
  par_df <- data.frame(bind.fill(est_par, type="rbind"))
  par_df$loc <- c(loc_1p_const, loc_else)
  par_df <-
    par_df %>%
    dplyr::arrange(.data$loc) %>%
    dplyr::select(-.data$loc)
  se_df <- data.frame(bind.fill(est_se, type="rbind"))
  se_df$loc <- c(loc_1p_const, loc_else)
  se_df <-
    se_df %>%
    dplyr::arrange(.data$loc) %>%
    dplyr::select(-.data$loc)

  # arrange the estimated item parameters and standard errors
  # when there is any item which have all missing response data
  if(length(loc_allmiss) > 0L) {
    tmp <-
      x_pre[, -c(1:3)] %>%
      dplyr::mutate_all(.funs=function(x) NA) %>%
      purrr::map(.x=1:nrow(x_pre), .f=function(i, .data) as.numeric(.data[i, ]))
    par_tmp <- se_tmp <- tmp
    par_tmp[loc_nomiss] <-
      par_df %>%
      purrr::map(.x=1:nrow(par_df), .f=function(i, .data) as.numeric(.data[i, ]))
    se_tmp[loc_nomiss] <-
      se_df %>%
      purrr::map(.x=1:nrow(se_df), .f=function(i, .data) as.numeric(.data[i, ]))
    par_df <- data.frame(bind.fill(par_tmp, type="rbind"))
    se_df <- data.frame(bind.fill(se_tmp, type="rbind"))
  }

  # create a full data.frame for the item parameter estimates
  if(length(loc_allmiss) > 0L) x <- x_pre

  full_par_df <- data.frame(x[, 1:3], par_df)
  full_par_df$id <- as.character(full_par_df$id)
  colnames(full_par_df) <- c("id", "cats", "model", paste0("par.", 1:ncol(par_df)))

  # create a full data.frame for the standard error estimates
  full_se_df <- data.frame(x[, 1:3], se_df)
  full_se_df$id <- as.character(full_se_df$id)
  colnames(full_se_df) <- c("id", "cats", "model", paste0("se.", 1:ncol(se_df)))

  # create a full data.frame including both the item parameter estimates and standard error estimates
  all_df <- data.frame(matrix(NA, nrow=nrow(x), ncol=2*ncol(par_df)))
  all_df[, seq(1, 2*ncol(par_df), 2)] <- par_df
  all_df[, seq(2, 2*ncol(par_df), 2)] <- se_df
  col.names <- rep(NA, 2*ncol(par_df))
  col.names[seq(1, 2*ncol(par_df), 2)] <- paste0("par.", 1:ncol(par_df))
  col.names[seq(2, 2*ncol(par_df), 2)] <- paste0("se.", 1:ncol(se_df))
  colnames(all_df) <- col.names
  full_all_df <- data.frame(x[, 1:3], all_df)

  ##---------------------------------------------------------------
  # check end time
  end.time <- Sys.time()

  # record total computation time
  est_time <- round(as.numeric(difftime(end.time, start.time, units = "secs")), 2)

  # return results
  rst <- structure(list(estimates=full_all_df, par.est=full_par_df, se.est=full_se_df, loglikelihood=llike,
                        group.par=group.par, data=data, score=score, scale.D=D, convergence=note, nitem=nitem,
                        deleted.item=as.numeric(loc_allmiss), npar.est=length(param_loc$reloc.par),
                        n.response=as.numeric(n.resp), TotalTime=est_time),
                   class="est_item")
  rst$call <- cl

  cat('Estimation is finished.', '\n')
  return(rst)


}


