# Item fit
#
# @description This function returns IRT item fit statistics (i.e., chi-square fit, infit, outfit statistics) and
# contingency table to be used in the residual plot for each item
#
# @param item_meta A row of data.frame containing each item's metadata (e.g., item parameters, number of categories, models ...).
# @param score A vector of examinees' ability estimates.
# @param resp A vector of examinees' responses for the item. Missing data should be recorded as NA.
# @param group.method A character string indicating how to group examinees. Available methods are
# "equal.width" and "equal.freq".
# @param n.width An integer value to specify the number of groups.
# @param loc.theta A character string to indicate the location of abiltiy point at each group where the expected probability
# is calculated. Available locations are "average" and "middle".
# @param D A scaling factor in IRT models to make the logistic function as close as possible to the normal ogive function (if set to 1.7).
# Default is 1.
# @param alpha A numeric value to specify significance level of the hypothesis test.
# @param overSR A numeric value to specify a criterion to find ability groups (or intervals) which have standardized residuals
# greater than the specified value. Default is 2.
# @param min.collapse An integer value to indicate the minimum frequency of cells to be collapsed when computing the \eqn{\chi^{2}} and \eqn{G^{2}}
# fit statistics. Neighboring interval groups will be collapsed to avoid expected interval frequencies less than the specified minimum cell frequency.
# Default is 1.
# @return A list
#
#' @importFrom rlang .data
itemfit <- function(item_meta, score, resp, group.method=c("equal.width", "equal.freq"),
                    n.width=10, loc.theta="average", D=1, alpha=0.05, overSR=2, min.collapse=1) {

  # number of score categories
  cats <- item_meta$cats
  
  # create a data.frame after removing NA rows
  tmp_df <-
    data.frame(resp, score) %>%
    tidyr::drop_na(.data$resp)
  
  # assign score levels to the response variable
  tmp_df$resp <- factor(tmp_df$resp, levels=(seq_len(item_meta$cats) - 1))
  
  # compute cut scores to divide score groups
  group.method <- match.arg(group.method)
  cutscore <- switch(group.method,
                     equal.width = seq(from=min(tmp_df$score), to=max(tmp_df$score), length.out=n.width + 1) ,
                     equal.freq = stats::quantile(tmp_df$score, probs=seq(0, 1, length.out=n.width + 1), type=9)
  )
  
  # when there are same cutscores,
  # delete the overlapped cutscores
  if(length(cutscore) != length(unique(cutscore))) {
    
    len_cuts <- length(cutscore)
    record <- rep(NA, len_cuts)
    for(i in 1:(len_cuts-1)) {
      if(cutscore[i] == cutscore[i + 1]){
        record[i + 1] <- 1
      }
    }
    cutscore <- cutscore[is.na(record)]
    
  }
  
  # assign score group variable to each score
  tmp_df$group.val <- cut(tmp_df$score, breaks=cutscore, include.lowest=TRUE, dig.lab = 7)
  
  # create a contingency table for the category frequencies
  obs.freq <-
    table(tmp_df$group.val, tmp_df$resp)
  
  # remove a row where total sum of category frequency is zero
  obs.freq2 <- obs.freq[rowSums(obs.freq) > 0, ]
  
  # create a contingency table for the category proportions
  obs.prop <- prop.table(obs.freq2, margin=1)
  
  # transfrom the contingency tables to the data.frame
  obs.freq2 <- as.data.frame.matrix(obs.freq2)
  obs.prop <- as.data.frame.matrix(obs.prop)
  
  # find a theta point for each score group
  loc.theta <- tolower(loc.theta)
  if(loc.theta == "middle") {
    theta <- purrr::map_dbl(.x=1:length(cutscore), .f=function(i) mean(c(cutscore[i], cutscore[i + 1])))
    theta <- theta[!is.na(theta)]
    theta <- theta[rowSums(obs.freq) > 0]
  }
  if(loc.theta == "average") {
    theta <-
      dplyr::group_by(tmp_df, .data$group.val) %>%
      dplyr::summarize(ave=mean(.data$score)) %>%
      dplyr::select(.data$ave) %>%
      unlist() %>%
      unname()
  }
  
  # compute expected probabilities of endorseing each score category
  exp.prob <-
    traceline(x=item_meta, theta, D=D)$prob.cat[[1]] %>%
    dplyr::rename_all(.funs=gsub, pattern="score.", replacement="")
  
  ##-------------------------------------------------------------------------
  # The contingency table made from here is used only for drawing IRT residual plots
  # standardized Residual
  # find a z-score corresponding to significance level
  zscore <- stats::qnorm(1-alpha)
  
  # total frequencies for each score group
  N <- rowSums(obs.freq2)
  
  # compute raw residuals
  raw_resid <- obs.prop - exp.prob
  
  # compute the standard errors
  se <- sqrt((exp.prob * (1 - exp.prob)) / N)
  
  # standardize the raw residuals
  std_resid <- raw_resid/se
  
  # compute a proportion of groups (or intervals) that have the standardized residuals greater than a specified criterion
  over_sr <- sum(abs(std_resid) > overSR)
  over_sr_prop <- round(over_sr / (cats * nrow(std_resid)), 3)
  
  # create a full contingency table to draw IRT residual plots
  contingency <- data.frame(theta=theta, N=N, freq=obs.freq2, obs.prop=obs.prop, exp.prob=exp.prob,
                            raw_resid=raw_resid, se=se, std_resid=std_resid)
  
  ##------------------------------------------------------------------------------
  # collapsing the contingency talbes to compute the chi-square fit statistics
  # check the number of expected frequency for all cells
  exp.freq <- exp.prob * N
  
  # collapse the expected and observed frequency tables
  ftable_info <- data.frame(exp.freq=exp.freq, obs.freq2=obs.freq2)
  for(i in 1:cats) {
    ftable_info <- collapse_ftable(x=ftable_info, col=i, min.collapse=min.collapse)
  }
  
  # new contingency tables after collapsing
  exp.freq.cp <- dplyr::select(ftable_info, dplyr::contains("exp.freq"))
  obs.freq2.cp <- dplyr::select(ftable_info, dplyr::contains("obs.freq2"))
  N.cp <- rowSums(obs.freq2.cp)
  exp.prob.cp <- exp.freq.cp / N.cp
  obs.prop.cp <- obs.freq2.cp / N.cp
  
  ##-------------------------------------------------------------------------
  # The contingency table made from here is used only for computing the chi-square fit statistic
  # compute raw residuals
  raw_resid.cp <- obs.prop.cp - exp.prob.cp
  
  # create a full contingency table for chi-square fit statistic
  contingency.cp <- data.frame(N=N.cp, freq=obs.freq2.cp, obs.prop=obs.prop.cp,
                               exp.prob=exp.prob.cp, raw_resid=raw_resid.cp)
  rownames(contingency.cp) <- 1:nrow(contingency.cp)
  col.names <- grep(pattern="N|freq|obs|exp|raw", x=colnames(contingency), value=TRUE)
  colnames(contingency.cp) <- col.names
  
  # compute the chi-square statistic (X2)
  x2 <- sum(N.cp * (raw_resid.cp^2 / exp.prob.cp), na.rm=TRUE)
  
  # compute the likelihood ratio chi-square fit statistic (G2)
  g2 <- 2 * sum(obs.freq2.cp * log(obs.prop.cp / exp.prob.cp), na.rm = TRUE)
  
  # find the number of parameters for each item
  model <- as.character(item_meta$model)
  count_prm <- NA
  count_prm[model %in% "1PLM"] <- 1
  count_prm[model %in% "2PLM"] <- 2
  count_prm[model %in% c("3PLM", "DRM")] <- 3
  count_prm[model %in% "PCM"] <- item_meta[model %in% "PCM", 2] - 1
  count_prm[model %in% "GPCM"] <- item_meta[model %in% "GPCM", 2]
  count_prm[model %in% "GRM"] <- item_meta[model %in% "GRM", 2]
  
  # find a critical value and compute the p values
  df.x2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1) - count_prm
  df.g2 <- nrow(exp.freq.cp) * (ncol(exp.freq.cp) - 1)
  crtval.x2 <- stats::qchisq(1-alpha, df=df.x2, lower.tail=TRUE)
  crtval.g2 <- stats::qchisq(1-alpha, df=df.g2, lower.tail=TRUE)
  pval.x2 <- 1 - stats::pchisq(x2, df=df.x2, lower.tail = TRUE)
  pval.g2 <- 1 - stats::pchisq(g2, df=df.g2, lower.tail = TRUE)
  
  ##------------------------------------------------------------------------------
  # infit & outfit
  # individual expected probabilities for each score category
  indiv_exp.prob <- traceline(item_meta, tmp_df$score, D=D)$prob.cat[[1]]
  
  # individual observed proportion for each score category
  indiv_obs.prob <-
    table(1:nrow(tmp_df), tmp_df$resp) %>%
    as.data.frame.matrix()
  
  # a matrix of the expected scores for each category
  Emat <- matrix(0:(cats-1), nrow(indiv_exp.prob), ncol(indiv_exp.prob), byrow = TRUE)
  
  # residuals
  resid <- rowSums(indiv_obs.prob * Emat) - rowSums(Emat * indiv_exp.prob)
  
  # variance
  Var <- rowSums((Emat - rowSums(Emat * indiv_exp.prob))^2 * indiv_exp.prob)
  
  # compute outfit & infit
  outfit <- sum(resid^2/Var) / length(tmp_df$score)
  infit <- sum(resid^2) / sum(Var)
  
  ##------------------------------------------------------------------------------
  # summary of fit statistics
  fitstats <- data.frame(X2=round(x2, 3), G2=round(g2, 3), df.X2=df.x2, df.G2=df.g2,
                         crit.value.X2=round(crtval.x2, 2), crit.value.G2=round(crtval.g2, 2),
                         p.value.X2=round(pval.x2, 3), p.value.G2=round(pval.g2, 3),
                         outfit=round(outfit, 3), infit=round(infit, 3),
                         N=nrow(tmp_df), overSR.prop=over_sr_prop)

  ##------------------------------------------------------------------------------
  # return results
  list(fit.stats=fitstats, contingency.fitstat=contingency.cp, contingency.plot=contingency,
       individual.info=data.frame(resid=resid, Var=Var))

}
