# This function assigns numbers to the item parameters to be estimated
parloc <- function(x, loc_1p_const, loc_else, fix.a.1pl, fix.a.gpcm, fix.g) {

  par_order <- vector('list', nrow(x))
  for(i in 1:nrow(x)) {

    i.tmp <- x[i, ]
    cats.tmp <- i.tmp[, 2]
    mod.tmp <- i.tmp[, 3]

    if(mod.tmp == "1PLM") {
      if(fix.a.1pl) {
        n.par <- 1
        if(i == 1) {
          par_order[[i]] <- c(NA, seq(1, n.par, by = 1))
        } else {
          first.num <- tail(par_order[[i-1]], 1) + 1
          par_order[[i]] <- c(NA, seq(first.num, first.num + (n.par - 1), by = 1))
        }
      } else {
        n.par <- ifelse(i == loc_1p_const[1], 2, 1)
        if(i == 1) {
          par_order[[i]] <- seq(1, n.par, by = 1)
        } else {
          first.num <- tail(par_order[[i-1]], 1) + 1
          if(n.par == 1) {
            par_order[[i]] <- c(NA, seq(first.num, first.num + (n.par - 1), by = 1))
          } else {
            par_order[[i]] <- seq(first.num, first.num + (n.par - 1), by = 1)
          }
        }
      }
    }

    if(mod.tmp %in% c("2PLM", "3PLM", "GPCM", "GRM")) {

      if(mod.tmp == "2PLM") {
        n.par <- 2
      }
      if(mod.tmp == "3PLM") {
        n.par <- ifelse(fix.g, 2, 3)
      }
      if(mod.tmp == "GPCM") {
        n.par <- ifelse(fix.a.gpcm, cats.tmp - 1, cats.tmp)
      }
      if(mod.tmp == "GRM") {
        n.par <- cats.tmp
      }
      if(i == 1) {
        if(mod.tmp == "GPCM" & fix.a.gpcm) {
          par_order[[i]] <- c(NA, seq(1, n.par, by = 1))
        } else {
          par_order[[i]] <- seq(1, n.par, by = 1)
        }
      } else {
        first.num <- tail(par_order[[i-1]], 1) + 1
        if(mod.tmp == "GPCM" & fix.a.gpcm) {
          par_order[[i]] <- c(NA, seq(first.num, first.num + (n.par - 1), by = 1))
        } else {
          par_order[[i]] <- seq(first.num, first.num + (n.par - 1), by = 1)
        }
      }
    }

  }

  if(!fix.a.1pl) {
    reloc.par <- c(par_order[loc_1p_const], par_order[loc_else])
  } else {
    reloc.par <- par_order
  }
  loc.par <- bind.fill(par_order, type="rbind")
  if(ncol(loc.par) == 2L) loc.par <- cbind(loc.par, rep(NA, nrow(loc.par)))
  reloc.par <- unlist(reloc.par)
  reloc.par <- reloc.par[!is.na(reloc.par)]

  # return the results
  list(loc.par=loc.par, reloc.par=reloc.par)

}
