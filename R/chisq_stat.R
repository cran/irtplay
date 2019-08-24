# This function calculates the central and non-central chi-square fit statistics 
chisq_stat <- function(exp.freq, obs.freq, count.prm, crt.delta, alpha) {
  
  # transform the two frequency tables to the matrix forms 
  exp.freq2 <- data.matrix(exp.freq)
  obs.freq2 <- data.matrix(obs.freq)
  
  # replace NA with 0 for the two frequency tables 
  exp.freq2[is.na(exp.freq2)] <- 0
  obs.freq2[is.na(obs.freq2)] <- 0
  
  # create the two proportion tables
  exp.prop <- prop.table(exp.freq2, margin=1)
  obs.prop <- prop.table(obs.freq2, margin=1)
  
  # compute the chi-square statistic
  chisq_fit <- sum(rowSums(obs.freq2) * ((obs.prop - exp.prop)^2 / exp.prop), na.rm=TRUE)
  
  # copy the expected proportion table and replace 0 with NA 
  exp.prop2 <- exp.prop
  obs.prop2 <- obs.prop
  exp.prop2[is.na(exp.freq)] <- NA
  obs.prop2[is.na(obs.freq)] <- NA
  
  # check if any cell exists whose value is greater than crt.delta
  # and and replace the value with the crt.delta
  diff_up <- 1 - exp.prop2
  diff_low <- exp.prop2 - 0
  diff_up[diff_up > crt.delta] <- crt.delta
  diff_low[diff_low > crt.delta] <- crt.delta
  
  # compute the two non-centrality parameters
  ncp_up <- sum(rowSums(obs.freq) * (diff_up^2 / exp.prop2), na.rm=TRUE)
  ncp_low <- sum(rowSums(obs.freq) * (diff_low^2 / exp.prop2), na.rm=TRUE)
  ncp <- max(ncp_up, ncp_low)
  
  # compute degrees of freedom
  # the number of collapsed cells
  counted_NA <- sum(is.na(exp.freq))
  
  # degrees of freedom
  df <- nrow(exp.freq) * (ncol(exp.freq) - 1) - count.prm - counted_NA
  
  # compute the central critical value
  crtval_cen <- stats::qchisq(1-alpha, df=df, lower.tail=TRUE)
  
  # compute the non-central critical value
  crtval_non <- stats::qchisq(1-alpha, df=df, ncp=ncp, lower.tail=TRUE)
  
  # return results    
  list(chisq_fit=chisq_fit, df=df, crtval_cen=crtval_cen, crtval_non=crtval_non, ncp=ncp, 
       exp.prop=exp.prop2, obs.prop=obs.prop2)
  
} 