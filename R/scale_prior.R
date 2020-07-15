# rescaling process by applying Woods's (2007) empirical histogram method
scale_prior <- function(prior_freq, prior_dense, quadpt, scale.par=c(0, 1), Quadrature) {

  # mean and sd of the updated prior distribution
  moments <- cal_moment(node=quadpt, weight=prior_dense)
  mu <- moments[1]
  sigma <- sqrt(moments[2])

  # standardize the updated prior distribution by adjusting the original quadrature points
  quadpt_star <- (quadpt - mu) / sigma

  # scale transforamtion using the specified mean and sd
  quadpt_star <- quadpt_star * sqrt(scale.par[2]) + scale.par[1]

  # distance betwenn two new quadrature points
  delta <- quadpt_star[2] - quadpt_star[1]

  # first point of the new quad points
  qstar_0 <- quadpt_star[1]

  # end point of the new quad points
  qstar_Q <- quadpt_star[Quadrature[1]]

  # translate back the standardized prior distribution to the original quadrature points
  # (a) extrapolate the frequencies of the original quad points less than or equal to the first new quad point
  quad_tmp1 <- quadpt[quadpt <= qstar_0]
  if(length(quad_tmp1) > 0) {
    freq_1 <- ((prior_freq[1] / prior_freq[2])^((qstar_0 - quad_tmp1) / delta)) * prior_freq[1]
  } else {
    freq_1 <- NULL
  }

  # (b) extrapolate the frequencies of the original quad points greater than or equal to the end new quad point
  quad_tmp3 <- quadpt[quadpt >= qstar_Q]
  if(length(quad_tmp3) > 0) {
    freq_3 <- ((prior_freq[Quadrature[1]] / prior_freq[Quadrature[1] - 1])^((quad_tmp3 - qstar_Q) / delta)) * prior_freq[Quadrature[1]]
  } else {
    freq_3 <- NULL
  }

  # (c) interpolate the frequencies of the rest original quad points
  quad_tmp2 <- quadpt[quadpt > qstar_0 & quadpt < qstar_Q]
  if(length(quad_tmp2) > 0) {
    interval <- c(-Inf, quadpt_star, Inf)
    group <- as.numeric(cut(x=quad_tmp2, breaks=interval, dig.lab=8, labels=1:(length(interval)-1)))
    # freq_2 <-
    #   data.frame(Q_O=quad_tmp2,
    #              Q_S=interval[group],
    #              N_S=prior_freq[group - 1],
    #              N_S2=prior_freq[group]) %>%
    #   dplyr::transmute(freq_2=(((.data$Q_O - .data$Q_S) / delta) * (.data$N_S2 - .data$N_S)) + .data$N_S) %>%
    #   unlist() %>%
    #   unname()
    Q_O <- quad_tmp2
    Q_S <- interval[group]
    N_S <- prior_freq[group - 1]
    N_S2 <- prior_freq[group]
    freq_2 <- (((Q_O - Q_S) / delta) * (N_S2 - N_S)) + N_S

  } else {
    freq_2 <- NULL
  }

  # translated empirical histrogram at the original quad points
  prior_freq2 <- c(freq_1, freq_2, freq_3)

  # renormalize the empirical historgam
  prior_dense2 <- prior_freq2 / sum(prior_freq2)

  # return the prior density
  prior_dense2

}
