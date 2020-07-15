# This function computes the log of the prior distribution
logprior <- function(val, is.aprior=FALSE, D=NULL, dist=c("lnorm", "beta", "norm"), par.1, par.2) {

  dist <- match.arg(dist)
  if(is.aprior) {
    switch(dist,
           lnorm = stats::dlnorm(D * val, par.1, par.2, log=TRUE),
           beta = stats::dbeta(D * val, par.1, par.2, log=TRUE),
           norm = stats::dnorm(D * val, par.1, par.2, log=TRUE))
  } else {
    switch(dist,
           lnorm = stats::dlnorm(val, par.1, par.2, log=TRUE),
           beta = stats::dbeta(val, par.1, par.2, log=TRUE),
           norm = stats::dnorm(val, par.1, par.2, log=TRUE))
  }

}


# This function computes the first and second derivative of the prior distribution
logprior_deriv <- function(val, is.aprior=FALSE, D=NULL, dist=c("lnorm", "beta", "norm"), par.1, par.2) {

  if(is.aprior) {
    switch(dist,
           lnorm = aprior_lnorm(a=val, D=D, par.1=par.1, par.2=par.2),
           beta = aprior_beta(a=val, D=D, par.1=par.1, par.2=par.2),
           norm = aprior_norm(a=val, D=D, par.1=par.1, par.2=par.2))
  } else {
    switch(dist,
           lnorm = esprior_lnorm(es=val, par.1=par.1, par.2=par.2),
           beta = esprior_beta(es=val, par.1=par.1, par.2=par.2),
           norm = esprior_norm(es=val, par.1=par.1, par.2=par.2))
  }

}

aprior_norm <- function(a, D, par.1, par.2) {
  .expr4 <- 1/(par.2 * sqrt(2 * pi))
  .expr6 <- D * a - par.1
  .expr10 <- 2 * par.2^2
  .expr12 <- exp(-.expr6^2/.expr10)
  .expr13 <- .expr4 * .expr12
  .expr18 <- 2 * (D * .expr6)/.expr10
  .expr19 <- .expr12 * .expr18
  .expr20 <- .expr4 * .expr19
  .value <- -log(.expr13)
  .grad <- .expr20/.expr13
  .hessian <-
    .expr4 * (.expr12 * (2 * (D * D)/.expr10) - .expr19 * .expr18)/
    .expr13 + .expr20 * .expr20/.expr13^2
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}

aprior_lnorm <- function(a, D, par.1, par.2) {
  .expr1 <- D * a
  .expr4 <- sqrt(2 * pi)
  .expr5 <- .expr1 * par.2 * .expr4
  .expr6 <- 1/.expr5
  .expr8 <- log(.expr1) - par.1
  .expr12 <- 2 * par.2^2
  .expr14 <- exp(-.expr8^2/.expr12)
  .expr15 <- .expr6 * .expr14
  .expr18 <- D/.expr1
  .expr21 <- 2 * (.expr18 * .expr8)/.expr12
  .expr22 <- .expr14 * .expr21
  .expr25 <- D * par.2 * .expr4
  .expr26 <- .expr5^2
  .expr27 <- .expr25/.expr26
  .expr29 <- .expr6 * .expr22 + .expr27 * .expr14
  .expr43 <- .expr27 * .expr22
  .value <- -log(.expr15)
  .grad <- .expr29/.expr15
  .hessian <-
    (.expr6 * (.expr14 * (2 * (.expr18 * .expr18 - D * D/.expr1^2 * .expr8)/.expr12) - .expr22 * .expr21) -
       .expr43 - (.expr43 + .expr25 * (2 * (.expr25 * .expr5))/.expr26^2 * .expr14))/
    .expr15 + .expr29 * .expr29/.expr15^2
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}

aprior_beta <- function(a, D, par.1, par.2) {
  .expr7 <- gamma(par.1 + par.2)/(gamma(par.1) * gamma(par.2)) *
    D
  .expr8 <- par.1 - 1
  .expr10 <- .expr7 * a^.expr8
  .expr12 <- 1 - D * a
  .expr13 <- par.2 - 1
  .expr14 <- .expr12^.expr13
  .expr15 <- .expr10 * .expr14
  .expr18 <- .expr8 - 1
  .expr21 <- .expr7 * (a^.expr18 * .expr8)
  .expr23 <- .expr13 - 1
  .expr25 <- .expr13 * D
  .expr26 <- .expr12^.expr23 * .expr25
  .expr28 <- .expr21 * .expr14 - .expr10 * .expr26
  .expr37 <- .expr21 * .expr26
  .value <- -log(.expr15)
  .grad <- -(.expr28/.expr15)
  .hessian <-
    -((.expr7 * (a^(.expr18 - 1) * .expr18 * .expr8) * .expr14 - .expr37 -
         (.expr37 - .expr10 * (.expr12^(.expr23 - 1) * (.expr23 * D) * .expr25)))/
        .expr15 - .expr28 * .expr28/.expr15^2)
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}

esprior_norm <- function(es, par.1, par.2){
  .expr4 <- 1/(par.2 * sqrt(2 * pi))
  .expr5 <- es - par.1
  .expr9 <- 2 * par.2^2
  .expr11 <- exp(-.expr5^2/.expr9)
  .expr12 <- .expr4 * .expr11
  .expr16 <- 2 * .expr5/.expr9
  .expr17 <- .expr11 * .expr16
  .expr18 <- .expr4 * .expr17
  .value <- -log(.expr12)
  .grad <- .expr18/.expr12
  .hessian <-
    .expr4 * (.expr11 * (2/.expr9) - .expr17 * .expr16)/
    .expr12 + .expr18 * .expr18/.expr12^2
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}

esprior_lnorm <- function(es, par.1, par.2){
  .expr3 <- sqrt(2 * pi)
  .expr4 <- es * par.2 * .expr3
  .expr5 <- 1/.expr4
  .expr7 <- log(es) - par.1
  .expr11 <- 2 * par.2^2
  .expr13 <- exp(-.expr7^2/.expr11)
  .expr14 <- .expr5 * .expr13
  .expr17 <- 1/es
  .expr20 <- 2 * (.expr17 * .expr7)/.expr11
  .expr21 <- .expr13 * .expr20
  .expr23 <- par.2 * .expr3
  .expr24 <- .expr4^2
  .expr25 <- .expr23/.expr24
  .expr27 <- .expr5 * .expr21 + .expr25 * .expr13
  .expr40 <- .expr25 * .expr21
  .value <- -log(.expr14)
  .grad <- .expr27/.expr14
  .hessian <-
    (.expr5 * (.expr13 * (2 * (.expr17 * .expr17 - 1/es^2 * .expr7)/.expr11) -
                 .expr21 * .expr20) - .expr40 - (.expr40 + .expr23 * (2 * (.expr23 * .expr4))/
                                                   .expr24^2 * .expr13))/.expr14 + .expr27 * .expr27/.expr14^2
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}

esprior_beta <- function(es, par.1, par.2) {
  .expr6 <- gamma(par.1 + par.2)/(gamma(par.1) * gamma(par.2))
  .expr7 <- par.1 - 1
  .expr9 <- .expr6 * es^.expr7
  .expr10 <- 1 - es
  .expr11 <- par.2 - 1
  .expr12 <- .expr10^.expr11
  .expr13 <- .expr9 * .expr12
  .expr16 <- .expr7 - 1
  .expr19 <- .expr6 * (es^.expr16 * .expr7)
  .expr21 <- .expr11 - 1
  .expr23 <- .expr10^.expr21 * .expr11
  .expr25 <- .expr19 * .expr12 - .expr9 * .expr23
  .expr34 <- .expr19 * .expr23
  .value <- -log(.expr13)
  .grad <- -(.expr25/.expr13)
  .hessian <-
    -((.expr6 * (es^(.expr16 - 1) * .expr16 * .expr7) * .expr12 - .expr34 -
         (.expr34 - .expr9 * (.expr10^(.expr21 - 1) * .expr21 * .expr11)))/
        .expr13 - .expr25 * .expr25/.expr13^2)
  .hessian <- as.matrix(.hessian)
  attr(.value, "gradient") <- .grad
  attr(.value, "hessian") <- .hessian
  .value
}
