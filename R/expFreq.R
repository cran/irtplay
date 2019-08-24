# This function returns a contingency table of the expected frequencies
# to be used to compute S-X2 fit statistic
expFreq <- function(t.score, cats, prob.cat, lkhd_noitem, lkhd, wts, score.freq) {

  tmp1 <- array(0, c(t.score + 1, cats))
  tmp2 <- array(0, c(t.score + 1, cats))

  for(j in 1:cats) {
    tmp1[j:(t.score + 1 - cats + j), j] <- colSums(prob.cat[, j]  * lkhd_noitem * wts[, 2])
    tmp2[, j] <- tmp1[, j]  / colSums(lkhd * wts[, 2])
  }

  colnames(tmp2) <- paste0("score.", 0:(cats-1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  tmp2 <- score.freq * tmp2
  tmp2 <- tmp2[c(-1, -nrow(tmp2)), ]

  row.first <- purrr::map_dbl(1:cats, .f=function(i) sum(tmp2[1:(cats - 1), i]))
  row.end <- purrr::map_dbl(1:cats, .f=function(i) sum(tmp2[nrow(tmp2):(nrow(tmp2)- cats + 2), i]))

  first.name <- rownames(tmp2)[cats - 1]
  last.name <- rownames(tmp2)[nrow(tmp2)- cats + 2]

  tmp2 <- tmp2[-c(1:(cats - 1), nrow(tmp2):(nrow(tmp2)- cats + 2)), ]
  tmp3 <- rbind(row.first, tmp2, row.end)
  rownames(tmp3) <- c(first.name, rownames(tmp2), last.name)

  data.frame(tmp3)

}

# This function returns a contingency table of the observed frequencies
# to be used to compute S-X2 fit statistic
obsFreq <- function(rawscore, response, t.score, cats) {

  tmp <- data.frame(score=rawscore, response=response)
  tmp$score <- factor(tmp$score, levels=0:(t.score))
  tmp$response <- factor(tmp$response, levels=0:(cats-1))

  tmp2 <-
    table(tmp) %>%
    as.data.frame.matrix()

  colnames(tmp2) <- paste0("score.", 0:(cats-1))
  rownames(tmp2) <- paste0("score.", 0:t.score)

  tmp2 <- tmp2[c(-1, -nrow(tmp2)), ]

  row.first <- purrr::map_dbl(1:cats, .f=function(i) sum(tmp2[1:(cats - 1), i]))
  row.end <- purrr::map_dbl(1:cats, .f=function(i) sum(tmp2[nrow(tmp2):(nrow(tmp2)- cats + 2), i]))

  first.name <- rownames(tmp2)[cats - 1]
  last.name <- rownames(tmp2)[nrow(tmp2)- cats + 2]

  tmp2 <- tmp2[-c(1:(cats - 1), nrow(tmp2):(nrow(tmp2)- cats + 2)), ]
  tmp3 <- rbind(row.first, tmp2, row.end)
  rownames(tmp3) <- c(first.name, rownames(tmp2), last.name)

  data.frame(tmp3)

}
