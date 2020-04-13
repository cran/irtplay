#' Import Item and Ability Parameters from IRT Software
#'
#' @description These functions import item and/or ability parameters from BILOG-MG 3, PARSCALE 4, flexMIRT, and
#' mirt (R package).
#'
#' @param file A file name (including a directory) containing the item or ability parameters.
#' @param type A character string indicating a type of output file. Available types are "par" for a file
#' containing item parameter estimates and "sco" for a file containing ability parameter estimates.
#' @param rePrm A logical value. If TRUE, the item intercept and logit of item guessing parameters are
#' reparameterized into the item difficulty and item guessing parameters, respectively. To estimate the IRT
#' linking coefficients in this package, the item difficulty and item guessing parameters must be used.
#' Default is TRUE.
#' @param n.factor A numeric value indicating the number of estimated factors. This argument should be specified
#' when \code{type = "sco"}. Default is 1.
#' @param x An output object obtained from the function \code{\link[mirt]{mirt}}.
#'
#' @details The \code{\link{bring.flexmirt}} was written by modifying the function \code{read.flexmirt}
#' (Pritikin, 2018). The functions \code{\link{bring.bilog}} and \code{\link{bring.parscale}}
#' were written by modifying the functions \code{read.bilog} and \code{read.parscale}
#' (Weeks, 2017), respectively.
#'
#' The file extensions for item parameter and ability files, respectively, are: ".par" and ".sco"
#' for BILOG-MG and PARSCALE, and "-prm.txt" and "-sco.txt" for flexMIRT. For mirt, the name of the output
#' object is specified by the user.
#'
#' Although \code{\link{bring.flexmirt}} is able to extract multidimensional item and ability parameter estimates,
#' this package only deals with unidimensional IRT methods.
#'
#' For polytomous item parameters, \code{\link{bring.flexmirt}} and \code{\link{bring.mirt}} are able to import
#' the item parameters of the graded response model and the (generalized) partial credit model.
#'
#' @return These functions return a list including several objects. Only for the output of flexMIRT, the results of
#' multiple group analysis can be returned. In that case, each element of the list contains the estimation results for
#' each group.
#'
#' @note Regarding the item parameter files for any IRT software, only the internal object "full_df" in the returned list is
#' necessary for the IRT linking. The object "full_df" is a data.frame containing the item meta data
#' in a test form (e.g., item parameters, number of categories, models). See \code{\link{test.info}}
#' or \code{\link{simdat}} for more details about the item meta data.
#'
#' Also, when item parameters are estimated using the partial credit or the generalized partial credit model,
#' item step parameters are returned in the object "full_df". Item step parameters are the overall item difficulty (or location)
#' parameter subtracted by the difficulty (or threshold) parameter for each category. See \code{\link{irtfit}} for more details
#' about the parameterization of the (generalized) partial credit model.
#'
#' @section Sample Output Files of IRT software:
#'
#' To illustrate how to import the item parameter estimate files of PARSCALE 4 and flexMIRT
#' using \code{\link{bring.parscale}} and \code{\link{bring.flexmirt}}, two item parameter
#' estimate output files are included in this package.
#'
#' Among the two output files, one of them is from PARSCALE 4 with a file extension of ".PAR"
#' (i.e., "parscale_sample.PAR") and another one is from flexMIRT
#' with a file extension of "-prm.txt" (i.e., "flexmirt_sample-prm.txt").
#'
#' For the two item parameter estimate output files, both are mixed-format tests with 55 items
#' consisting of fifty dichotomous items following the IRT 3PL model and five polytomous items with five
#' categories following the graded response model. The examples below show how to import those output files.
#'
#' @author Hwanggyu Lim \email{hglim83@@gmail.com}
#'
#' @seealso \code{\link{irtfit}}
#'
#' @references
#' Cai, L. (2017). flexMIRT 3.5 Flexible multilevel multidimensional item analysis and test scoring [Computer software].
#' Chapel Hill, NC: Vector Psychometric Group.
#'
#' Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment.
#' \emph{Journal of Statistical Software, 48}(6), 1-29.
#'
#' Weeks, J. P. (2010). plink: An R Package for Linking Mixed-Format Tests Using IRT-Based Methods.
#' \emph{Journal of Statistical Software, 35}(12), 1-33. URL http://www.jstatsoft.org/v35/i12/.
#'
#' Pritikin, J. (2018). \emph{rpf: Response Probability Functions}. R package version 0.59.
#' https://CRAN.R-project.org/package=rpf
#'
#' Muraki, E. & Bock, R. D. (2003). PARSCALE 4: IRT item analysis and test scoring for rating
#' scale data [Computer Program]. Chicago, IL: Scientific Software International. URL http://www.ssicentral.com
#'
#' Zimowski, M. F., Muraki, E., Mislevy, R. J., & Bock, R. D. (2003). BILOG-MG 3: Multiple-group
#' IRT analysis and test maintenance for binary items [Computer Program]. Chicago, IL: Scientific
#' Software International. URL http://www.ssicentral.com
#'
#' @examples
#' ## example 1
#' # import the "-prm.txt" output file from flexMIRT
#' flex_sam <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' bring.flexmirt(file=flex_sam, "par")$Group1$full_df
#'
#' ## example 2
#' ## import the ".par" output file from PARSCALE
#' pscale_sam <- system.file("extdata", "parscale_sample.PAR", package = "irtplay")
#'
#' # read item parameters and transform them to item meta data
#' bring.parscale(file=pscale_sam, "par")$full_df
#'
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.flexmirt <- function(file, type=c("par", "sco"), rePrm = TRUE, n.factor=1) {

  if(length(type)==0L) stop("Specify one of file types: 'par' or 'sco'", call.=FALSE)

  type <- toupper(type)

  # When "par" file is inserted
  if(type=="PAR") {
    results <- bring.flexPrm(file, rePrm)
  }

  # When "sco" file is inserted
  if(type=="SCO") {
    results <- bring.flexSco(file, n.factor)
  }
  results
}

# "bring.flexPrm" is a modified version of "read.flexmirt" in R package of "rpf" (Pritikin, Weeks, Cai, Houts, and Chalmers, 2016)
# Read a flexMIRT PRM file
# Load the item parameters from a flexMIRT PRM file.
#' @importFrom utils count.fields read.delim read.fwf read.table
bring.flexPrm <- function(file, rePrm = TRUE) {

  groups <- list()

  ncol <- max(utils::count.fields(file, sep = "\t"))
  prm <- utils::read.delim(file,header=FALSE, col.names=1:ncol)

  # find out how many rows
  nlines <- nrow(prm)

  # determine the number of groups
  ng <- length( unique(prm[,3]) )

  for (g in 1:ng) {
    # select all items and this group
    index <- (prm[,3] == g)
    thisGroup <- prm[index,]

    g.name <- "?"
    g.dist <- list()
    g.label <- list()
    g.category <- list()
    g.dim <- list()
    g.param <- list()
    g.model <- list()

    for (i in 1:nrow(thisGroup)) {

      param <- c()
      dims <- thisGroup[i,4]

      if (thisGroup[i,1] == 1) { # item
        catg <- thisGroup[i,6]

        if (thisGroup[i,5] == 1) { # 3PL
          # grab item parameters
          logitg <- thisGroup[i,7]
          c <- thisGroup[i,8]
          a <- thisGroup[i,9:(9+dims-1)]
          if(rePrm) {
            if(!dims == 1) stop("Re-parameterization works only for one-dimensional model", call.=FALSE)
            param <- c(unlist(a), -c/unlist(a), exp(logitg)/(1+exp(logitg)))
          } else {
            param <- c(unlist(a), c, logitg)
          }
          model <- "3PLM"
        }

        if (thisGroup[i,5] == 2) { # graded
          # grab item parameters
          c <- rep(0,catg-1)
          for (k in 1:(catg-1)) {
            c[k] <- thisGroup[i,6+k]
          }
          a <- thisGroup[i,(6+catg):(6+catg+dims-1)]
          if(rePrm) {
            if(!dims == 1) stop("Re-parameterization works only for one-dimensional model", call.=FALSE)
            param <- c(unlist(a), -c/unlist(a))
          } else {
            param <- c(unlist(a), c)
          }
          model <- ifelse(catg==2, "2PLM", "GRM")
        }

        if (thisGroup[i,5] == 3) { # gpc
          offset <- 7
          # grab item parameters
          Tmat <- thisGroup[i,offset]
          offset <- offset+1

          # scoring fn
          alf <- rep(0,catg-1)
          for (k in 1:(catg-1)) {
            alf[k] <- thisGroup[i,offset]
            offset <- offset+1
          }
          alfGpc <- c(1, rep(0, catg-2))
          isGpc <- all(alf == alfGpc)

          # slope
          a <- thisGroup[i,offset:(offset+dims-1)]
          offset <- offset+dims

          # intercept
          Lmat <- thisGroup[i,offset]
          offset <- offset+1
          gam <- rep(0,catg-1)
          for (k in 1:(catg-1)) {
            gam[k] <- thisGroup[i,offset]
            offset <- offset+1
          }

          # For GPC parameterization (only when dimension is one)
          if(dims == 1 & isGpc) {
            Tmat <- matrix(0, nrow=catg, ncol=(catg-1))
            Tmat[, 1] <- 0:(catg-1)
            if(catg > 2) {
              for(k in 2:(catg-1)) {
                for(j in 2:(catg-1)) {
                  Tmat[k, j] <- sin( pi*(j-1)*(k-1) / (catg-1) )
                }
              }
            }
            c.Vec <- Tmat%*%gam
            b <- -gam[1] / a
            d <- rep(0, length(c.Vec))
            for(k in 2:length(d)) {
              d[k] <- ((c.Vec[k] - c.Vec[k-1]) / a) + b
            }
            d.new <- b - d[-1]
          }

          if(rePrm) {
            if(!dims == 1 | !isGpc) stop("Re-parameterization works only for one-dimensional GPC model", call.=FALSE)
            if(dims == 1 & isGpc)  param <- c(a, d.new)
          } else {
            param <- c(unlist(a), alf, gam)
          }
          if(catg==2) model <- "2PLM"
          if(catg > 2 & isGpc) model <- "GPCM"
          if(catg > 2 & !isGpc) model <- "NOMM"
        }

        g.label[[i]] <- as.character(thisGroup[i,2])
        g.category[[i]] <- catg
        g.dim[[i]] <- dims
        g.param[[i]] <- param
        g.model[[i]] <- model

      } else { # group
        if(thisGroup[i, 5] == 0) { # EmpHis = No (default)
          Cov <- matrix(NA, nrow=dims, ncol=dims)
          col <- 6+dims
          for (rx in 1:dims) {
            for (cx in 1:rx) {
              Cov[rx,cx] <- thisGroup[i,col]
              if (rx != cx) Cov[cx,rx] <- thisGroup[i,col]
              col <- col+1
            }
          }
          Mean <- unlist(thisGroup[i,6:(6+dims-1)])
          names(Mean) <- paste0("Mu.", 1:length(Mean))
          colnames(Cov) <- paste0("Theta.", 1:ncol(Cov))
          rownames(Cov) <- paste0("Theta.", 1:nrow(Cov))
          g.dist <- list(Mean=Mean, Cov=Cov)
        }
        if(thisGroup[i, 5] == 1) { # EmpHis = Yes
          Cov <- matrix(NA, nrow=dims, ncol=dims)
          col <- 8+dims
          for (rx in 1:dims) {
            for (cx in 1:rx) {
              Cov[rx,cx] <- thisGroup[i,col]
              if (rx != cx) Cov[cx,rx] <- thisGroup[i,col]
              col <- col+1
            }
          }
          Mean <- unlist(thisGroup[i,8:(8+dims-1)])
          names(Mean) <- paste0("Mu.", 1:length(Mean))
          colnames(Cov) <- paste0("Theta.", 1:ncol(Cov))
          rownames(Cov) <- paste0("Theta.", 1:nrow(Cov))
          g.dist <- list(Mean=Mean, Cov=Cov)
        }
      }
    } # for every item

    pmat <- matrix(NA, nrow=length(g.param), ncol=max(vapply(g.param, length, 0)))
    for (i in 1:length(g.param)) {
      v <- g.param[[i]]
      pmat[i, 1:length(v)] <- v
    }
    # DF <- data.frame(ID=unlist(g.label), FACTOR=unlist(g.dim), CATEGORY=unlist(g.category),
    # 				MODEL=unlist(g.model), PARAM=pmat)
    DF <- data.frame(id=unlist(g.label), cats=unlist(g.category),
                     model=unlist(g.model), par=pmat, stringsAsFactors=FALSE)
    rownames(pmat) <- g.label
    groups[[g]] <- list(id=unlist(g.label), factor=unlist(g.dim), cats=unlist(g.category),
                        model=unlist(g.model), par=pmat, full_df=DF, mean=g.dist$Mean, cov=g.dist$Cov)
    names(groups)[[g]] <- paste0("Group", g)
  } # for every group

  groups
}


# "bring.flexSco" function
# This function reads a scoring output of flexMIRT, "-sco.txt",
# and returns a list
#' @importFrom utils count.fields read.delim read.fwf read.table
bring.flexSco <- function(file, n.factor=1) {

  groups <- list()

  ncol <- max(utils::count.fields(file, sep = ""))

  # read "-sco.txt" file
  sco <- utils::read.delim(file,header=FALSE, sep="", col.names=1:ncol)

  # (1) unidimensional model
  if(n.factor == 1) {
    if(ncol == 4) {
      if(typeof(sco[1,3]) == "integer") {
        # determine the number of groups
        ng <- length( unique(sco[,2]) )
        for(g in 1:ng) {
          index <- (sco[,2] == g)
          thisGroup <- sco[index,]

          g.imputation <- sco[index, 1]
          g.ability <- sco[index, 4]
          g.se <- NULL

          groups[[g]] <- list(theta = g.ability, se = g.se, imputation=g.imputation)
          names(groups)[[g]] <- paste0("Group", g)
        }
      } else {
        # determine the number of groups
        ng <- length(unique(sco[,1]))
        for(g in 1:ng) {
          index <- (sco[,1] == g)
          thisGroup <- sco[index,]

          g.ability <- sco[index, 3]
          g.se <- sco[index, 4]

          groups[[g]] <- list(theta = g.ability, se = g.se)
          names(groups)[[g]] <- paste0("Group", g)
        }
      }
    }

    if(ncol == 5) {
      # determine the number of groups
      ng <- length( unique(sco[,1]) )
      for(g in 1:ng) {
        index <- (sco[,1] == g)
        thisGroup <- sco[index,]

        g.iter <- sco[index, 3]
        g.ability <- sco[index, 4]
        g.se <- sco[index, 5]

        groups[[g]] <- list(theta = g.ability, se = g.se, iter = g.iter)
        names(groups)[[g]] <- paste0("Group", g)
      }
    }
  }

  # (2) multidimensional model
  if(n.factor > 1) {
    if(typeof(sco[1,3]) == "integer") {
      # determine the number of groups
      ng <- length( unique(sco[,2]) )
      for(g in 1:ng) {
        index <- (sco[,2] == g)
        thisGroup <- sco[index,]

        g.imputation <- sco[index, 1]
        g.ability <- sco[index, 4:(4+(n.factor-1))]
        g.se <- NULL

        groups[[g]] <- list(theta = g.ability, se = g.se, imputation=g.imputation)
        names(groups)[[g]] <- paste0("Group", g)
      }
    } else {
      # determine the number of groups
      ng <- length(unique(sco[,1]))
      for(g in 1:ng) {
        index <- (sco[,1] == g)
        thisGroup <- sco[index,]

        g.ability <- sco[index, 3:(3+(n.factor-1))]
        g.se <- NULL # will be updated

        groups[[g]] <- list(theta = g.ability, se = g.se)
        names(groups)[[g]] <- paste0("Group", g)
      }
    }
  }

  groups
}

######################################################################################################################################
# "bring.bilog" function
# This function is to read "par" or "sco" output files from Bilog-MG

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.bilog <- function(file, type=c("par", "sco")) {

  if(length(type)==0L) stop("Specify one of file types: 'par' or 'sco'.", call.=FALSE)

  type <- toupper(type)

  # When "par" file is inserted
  if(type=="PAR") {
    df <- utils::read.fwf(file, widths=c(8, 8, rep(10, 10)), skip=4)
    df[, 1] <- trimws(as.character(df[, 1]), )
    prm <- df[c(1, 5, 7, 11)]
    prm.df <- data.frame(id=prm[,1], cats=rep(2,nrow(prm)), model="DRM",
                         par.1=prm[, 2], par.2=prm[, 3], par.3=prm[, 4], stringsAsFactors=FALSE)
    results <- list(par=prm[-1], full_df=prm.df)
  }

  # When "sco" file is inserted
  if(type=="SCO") {
    # df <- utils::read.table(file, skip=2, fill=TRUE)
    df1 <- utils::read.fwf(file, widths=c(7, 8, 5, 5, 10, 12, 12), skip=2)
    df2 <- utils::read.fwf(file, widths=c(3, 2, 30), skip=2)
    num1 <- seq(2, nrow(df1), 2)
    num2 <- seq(1, (nrow(df2)-1), 2)

    # Delete rows including examinees' id
    sco.df1 <- df1[num1, 3:7]
    sco.df2 <- df2[num2, c(1, 3)]
    sco.df <- cbind(sco.df2, sco.df1)
    colnames(sco.df) <- c("group","id","nitem","ncorrect","p","theta","se")

    # Splite dataframe based on a group factor
    sco.df$group <- as.factor(sco.df$group)

    # Make a list including all info per a group
    scoList <- list()
    ng <- length(levels(sco.df$group))
    dfList <- split(sco.df, sco.df$group)

    for(g in 1:ng) {
      g.name <- levels(sco.df$group)[g]
      g.ability <- dfList[[g]]$theta
      g.se <- dfList[[g]]$se
      g.df <- dfList[[g]]

      scoList[[g]] <- list(theta = g.ability, se = g.se, full_df = g.df)
      names(scoList)[[g]] <- as.character(g.name)
    }

    results <- scoList
  }
  results
}



######################################################################################################################################
# "bring.parscale" function
# This function is to read "par" or "sco" output files from Parscale
# "bring.parscale" is a modified version of "read.parscale" in R package of "plink" (Jonathan P. Weeks, 2017)

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.parscale <- function(file, type=c("par", "sco")) {

  if(length(type)==0L) stop("Specify one of file types: 'par' or 'sco.", call.=FALSE)

  type <- toupper(type)

  ##   Import item parameters
  if (type=="PAR") {

    ##   Read the third line of the .PAR file which identifies the number
    ##   of parameter blocks, the total number of items, and the model code
    nums <- as.numeric(utils::read.fwf(file, c(8,rep(5,5)), skip=2, n=1)[-1])

    ##   Read the fourth line in the .PAR file which identifies
    ##   the number of items in each block
    line_quot <- nums[1] %/% 30 # quotient
    numLines <- 2 * line_quot - 1
    if(nums[1] %% 30 != 0) {
      numLines <- numLines + 2
    }
    block_str1 <- readLines(file, n = (3+numLines))[-c(1, 2, 3)] # strings containing block information
    block_str2 <- strsplit(block_str1, '\\s+')
    block_str3 <- as.numeric(do.call('c', block_str2))
    block <- block_str3[!is.na(block_str3)]

    ##   Initialize an object to store the item parameters
    pars <- NULL

    ##   Initialize a list to store the deviation threshold/step parameters for each block
    pars1 <- vector("list",length(block))

    ##   Initialize an object to store the number of response
    ##   categories for each item
    cat <- NULL

    ##   Identify the starting row for reading in the
    ##   parameters for a given block
    skip <- 4 + numLines

    ##   Identify the starting row for reading in the
    ##   category parameters for a given block
    skip1 <- skip+block[1]

    ##   Loop through each block of parameters
    for (i in 1:length(block)) {
      if (i>1) {
        ##   Update the starting rows for the given block
        skip <- (4 + numLines) + sum(block[1:(i-1)]) + (i-1)*2
        skip1 <- (4 + numLines) + sum(block[1:i]) + (i-1)*2
      }

      ##   Read the item parameters for the given block
      tmp <- utils::read.fwf(file, c(8,5,4,rep(10,6)), skip=skip, n=block[i], colClasses=c("character","numeric","character",rep("numeric",6)))

      ##   Add this block of parameters to the full set of parameters
      pars <- rbind(pars, tmp)

      ##   If there is more than one item in the block
      if (nrow(tmp) > 1) {
        cat <- c(cat, as.numeric(tmp[,2]))

        ##   If there is only one item in the block
      } else {
        cat <- c(cat, as.numeric(tmp[2]))
      }

      ##   Read in the deviation threshold/step parameters parameters for the given block
      tmp1 <- utils::read.fwf(file, rep(10,cat[length(cat)]), skip=skip1, n=2)

      ##   For the graded response models
      ##   Eliminate the last category parameter
      if (nums[3]<5) {
        tmp1 <- tmp1[,-ncol(tmp1)]

        ##   For the  partial credit models
        ##   Eliminate the first category parameter
      } else {
        tmp1 <- tmp1[,-1]
      }
      pars1[[i]] <- tmp1
    }
    colnames(pars) <- c("block.name","cat","item","slope","slope.se","location","loc.se","asymptote","asymp.se")

    ##   For polytomous items
    if (max(cat)>2) {
      ##   Create a matrix for the threshold/step parameters
      step <- matrix(NA, nrow(pars), 2*(max(cat)-1))

      ##   Because multiple items can all have the same
      ##   deviation threshold/step parameters, expand {pars1}
      ##   so that each item has a corresponding list element
      ##   containing threshold/step parameters
      pars1 <- rep(pars1,block)

      ##   Loop through all of the list elements for {pars1}
      ##   and insert the deviation threshold/step parameters
      ##   into the matrix {step}
      for (i in 1:length(pars1)) {
        ##   We only need deviation threshold/step parameters
        ##   for polytomous items. In the list, these will be
        ##   formatted as a data.frame whereas the parameters
        ##   for dichotomous items will be formatted as a
        ##   vector (of zeros)
        if (is.data.frame(pars1[[i]])) {
          step[i,] <- unlist(pars1[[i]])
        }
      }

      if (nums[3]<5) prefix <- "thresh" else prefix <- "step"

      ##   Create column names for the matrix of deviation threshold/step parameters
      step.names <- character()
      for (j in 1:(max(cat)-1)) {
        step.names <- c(step.names, paste(prefix,j,sep=""), paste(prefix,j,".se",sep=""))
      }
      colnames(step) <- step.names

      ##   Combine the matrix of item parameters with the matrix
      ##   of deviation threshold/step parameters
      pars <- cbind(pars, step)
    }
    pars_loc <- pars <- data.frame(pars)

    ##   For polytomous items, reformulate the parameters to
    ##   exclude the location parameter (if applicable)
    if (max(cat)>2) {
      if (nums[3] %in% c(2,4,6,8)) {
        pars[cat>2,6] <- mean(pars[cat>2,seq(10,ncol(pars),2)], na.rm=TRUE)
        pars[cat>2,7] <- NA
        pars[cat>2,seq(10,ncol(pars),2)] <- pars[cat>2,6] - pars[cat>2,seq(10,ncol(pars),2)]
      }
      if (nums[3] %in% c(1,3,5,7)) {
        pars[cat>2,seq(10,ncol(pars),2)] <- pars[cat>2,6] - pars[cat>2,seq(10,ncol(pars),2)]
        pars[cat>2,6:7] <- 0
      }
    }

    ##   Get rid of all the columns of containing standard errors
    ##   the item names, and the number of categories
    pars_nose <- pars[,seq(4,ncol(pars),2)]

    ##   reverse the order of threshold/step parameters for polytomous items
    # if(max(cat)>2) {
    #   pars_nose <- cbind(pars_nose[, 1:3], rev(pars_nose[, 4:ncol(pars_nose)]))
    # }

    ##   Make sure that the matrix of parameters is formatted properly
    ##   The threshold/step parameters in the matrix need to be shifted
    ##   one or two columns to the left
    if (max(cat)>2) {
      pars_nose[cat>2,2:(ncol(pars_nose)-2)] <- pars_nose[cat>2,4:ncol(pars_nose)]
      pars_nose <- pars_nose[,1:max(cat)]
    }
    colnames(pars_nose) <- paste0("par.", 1:ncol(pars_nose))

    ##   Number of items
    n <- nrow(pars)

    ##   Model of items
    if (nums[3]<5) {
      mod <- ifelse(cat == 2, "DRM", "GRM")
    } else {
      mod <- ifelse(cat == 2, "DRM", "GPCM")
    }

    ##   create a full data frame
    # pars_full <- data.frame(ID=pars$item, FACTOR=rep(1, n), CATEGORY=cat, MODEL=mod, pars_nose)
    pars_full <- data.frame(id=pars$item, cats=cat, model=mod, pars_nose, stringsAsFactors=FALSE)

    ##   a list to be returned
    rr <- list(id=pars$item, factor=rep(1, n), cats=cat, model=mod,
               par_parscale=pars_loc, par=pars_nose, full_df=pars_full)

  }

  ##   Import ability estimates
  if (type=="SCO") {
    pars <- read.fwf(file, c(21,32,12,12))
    id <- as.character(pars[seq(1,nrow(pars),2),1])
    ability <- pars[seq(2,nrow(pars),2),3:4]
    ability[ability==999] <- NA
    pars <- cbind(id,ability)
    colnames(pars) <- c("id","theta","se")
    rr <- pars
  }

  return(rr)
}


######################################################################################################################################
# "bring.parscale.ph2" function
# This function is to read the information of posterior distribution and the mean and sd of population dist'n.

bring.parscale.ph2 <- function(file) {

  phase2 <- readLines(file)

  # find the lines containing posterior information
  pts_nums <- grep("POINT", phase2, value=FALSE)
  wts_nums <- pts_nums + 1

  # select points of posterior distribution
  pts_str <- do.call('c', strsplit(phase2[pts_nums], '\\s+'))
  pts_str <- suppressWarnings(as.numeric(pts_str))

  # select weights of posterior distribution
  wts_str <- do.call('c', strsplit(phase2[wts_nums], '\\s+'))

  # create a data frame for posterior information
  posterior <- data.frame(POINT=pts_str, WEIGHT=wts_str, stringsAsFactors=FALSE)
  posterior <- posterior[!is.na(posterior$POINT), ]
  num.wts <- suppressWarnings(as.numeric(posterior$WEIGHT))
  if(any(is.na(num.wts))) num.wts[is.na(num.wts)] <- 0
  posterior$WEIGHT <- num.wts

  # find a mean and a sd of population distribution
  num <- grep("TOTAL WEIGHT", phase2, value=FALSE)
  mu_str <- phase2[num + 1]
  sigma_str <- phase2[num + 2]
  mu <- as.numeric(gsub(" MEAN        :", "", mu_str))
  sigma <- as.numeric(gsub(" S.D.        :", "", sigma_str))
  stat.pop <- c(MEAN=mu, SD=sigma)

  rr <- list(POSTERIOR=posterior, STAT.POPULATION=stat.pop)
  rr

}


######################################################################################################################################
# "bring.mirt" function
# This function is to read parameter estimates of an object obtained from MIRT

#' @rdname bring.flexmirt
#' @importFrom utils count.fields read.delim read.fwf read.table
#' @export
bring.mirt <- function(x) {

  # read paramter estimates from an object of mirt
  prm_all <- mirt::coef(x, simplify=TRUE, IRTpars=TRUE)
  cats <- x@Data$K
  any_py <- any(cats > 2)

  # extract other information from an object
  wts <- as.numeric(x@Internals$Prior[[1]])
  nods <- as.numeric(x@Model$Theta)
  EMPHIST <- cbind(nods=nods, weights=wts)
  mods <- x@Model$itemtype
  MODEL <- mods
  MODEL[MODEL %in% c("gpcm", "gpcmIRT")] <- "GPCM"
  MODEL[MODEL == "Rasch" & cats > 2] <- "GPCM"
  MODEL[MODEL == "graded"] <- "GRM"
  MODEL[MODEL == "2PL"| MODEL == "3PL"] <- "DRM"
  MODEL[MODEL == "Rasch" & cats == 2L ] <- "DRM"
  CATEGORY <- cats
  FACTOR <- mirt::extract.mirt(x, what="nfact")
  ID <- mirt::extract.mirt(x, what="itemnames")
  nitem <- mirt::extract.mirt(x, what="nitems")
  MEAN <- as.numeric(prm_all$mean)
  COV <- as.numeric(prm_all$cov)

  # check dimensionality
  if(FACTOR > 1) stop("This function only works for an unidimensional model", call.=FALSE)

  if(!any_py) {
    # check models for dichotomous data
    if(!all(mods %in% c("Rasch", "2PL", "3PL"))) stop("For dichotomous data, only Rasch, 2PL, 3PL models are available", call.=FALSE)
    prm_item <- prm_all$items
    prm_item <- prm_item[, which(colnames(prm_item) != "u")]
  }

  if(any_py) {
    mods_py <- mods[cats > 2]

    # check models for polytomous data
    if(!all(mods_py %in% c("Rasch", "graded", "gpcm", "gpcmIRT"))) stop("For polytomous data, only Rasch, graded, gpcm, gpcmIRT models are available", call.=FALSE)

    # when a single polytimous model is used
    if(length(unique(mods_py)) == 1) {
      prm_item <- prm_all$items

      # when a graded or gpcm model is used
      if(unique(mods_py) != "gpcmIRT") {
        prm_item <- prm_item[, which(colnames(prm_item) != "u")]
        if(any(cats == 2L)) {
          prm_item[cats > 2,2:(ncol(prm_item)-2)] <- prm_item[cats > 2,4:ncol(prm_item)]
        }
        prm_item <- prm_item[,1:max(cats)]

        # when a gpcmIRT model is used
      } else {
        prm_item <- prm_item[, which(!colnames(prm_item) %in% c("u", "c"))]
        if(any(cats == 2L)) {
          prm_item[cats > 2,1:(ncol(prm_item)-3)] <- prm_item[cats > 2,4:ncol(prm_item)]
        }
        prm_item <- prm_item[,1:max(cats)]
      }
    }

    # when two different polytimous models are used
    if(length(unique(mods_py)) > 1) {
      prm_item <- prm_all$items

      # when one of the models is gpcmIRT model
      if(any(mods_py %in% "gpcmIRT")) {
        prm_item <- prm_item[, which(!colnames(prm_item) %in% c("u", "c"))]

        # when there are dichotomous models
        if(any(cats == 2L)) {

          if(all(mods_py %in% "gpcmIRT")) {
            # when all polytomous models are gpcmIRT models
            prm_item[cats > 2 & mods == "gpcmIRT", 1:length(4:ncol(prm_item))] <- prm_item[cats > 2 & mods == "gpcmIRT", 4:ncol(prm_item)]

          } else {
            # when gpcmIRT models are used with other polytomous models
            prm_item[cats > 2 & mods == "gpcmIRT", 1] <- prm_item[cats > 2 & mods == "gpcmIRT", ncol(prm_item)]
            prm_item <- prm_item[, -ncol(prm_item)]
            prm_item[cats > 2,2:(ncol(prm_item)-2)] <- prm_item[cats > 2,4:ncol(prm_item)]
          }
        } else {

          # when there is no dichotomous model
          prm_item[cats > 2 & mods == "gpcmIRT", 1] <- prm_item[cats > 2 & mods == "gpcmIRT", ncol(prm_item)]
        }

        prm_item <- prm_item[,1:max(cats)]
      }

      # when there is no gpcmIRT model
      if(all(mods_py %in% c("graded", "gpcm", "Rasch"))) {
        prm_item <- prm_item[, which(colnames(prm_item) != "u")]
        if(any(cats == 2L)) {
          prm_item[cats > 2,2:(ncol(prm_item)-2)] <- prm_item[cats > 2,4:ncol(prm_item)]
        }
        prm_item <- prm_item[,1:max(cats)]
      }

    }

  }

  # give colnames for a prm matrix
  colnames(prm_item) <- paste0("par.", 1:ncol(prm_item))

  # original item parameters matrix obtained from mirt
  prm_mirt <- mirt::coef(x, simplify=TRUE, IRTpars=FALSE)$items

  # full data.frame of item parameters
  # prm_full <- data.frame(ID=ID, FACTOR=rep(FACTOR, nitem), CATEGORY=cats, MODEL=MODEL, prm_item)
  prm_full <- data.frame(id=ID, cats=cats, model=MODEL, prm_item, stringsAsFactors=FALSE)

  rr <- list(id=ID, factor=rep(FACTOR, nitem), cats=cats, model=list(old=mods, new=MODEL),
             emphist = EMPHIST, par_mirt=prm_mirt, par=prm_item, full_df=prm_full,
             ability=c(MEAN=MEAN, COV=COV))

  rr

}

