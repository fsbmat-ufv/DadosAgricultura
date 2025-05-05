loglik_HC <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Numero de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Numero de colunas de XO+1
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) return(NA)
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) return(NA)
  
  XS.g <- model.matrix(~XS) %*% g
  XO.b <- model.matrix(~XO) %*% b
  u2 <- YO - XO.b
  r <- sqrt( 1 - rho^2)
  B <- (XS.g + rho/sigma*u2)/r
  ll <- ifelse(YS == 0,
               (pnorm(-XS.g, log.p=TRUE)),
               dnorm(u2/sigma, log = TRUE) - log(sigma) +
                 (pnorm(B, log.p=TRUE))
  )
  sum(ll)
}

gradlik_HC <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Número de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Número de colunas de XO+1
  nObs <- length(YS)
  NO <- length(YS[YS > 0])
  nParam <- NXS + NXO + 2 #Total of parameters
  
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  N0 <- sum(YS==0)
  N1 <- sum(YS==1)
  
  w  <- rep(1,N0+N1 )
  w0 <- rep(1,N0)
  w1 <- rep(1,N1)
  NXS <- dim(model.matrix(~XS))[2]#Número de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Número de colunas de XO+1
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) return(matrix(NA, nObs, nParam))
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
  XS0.g <- as.numeric(model.matrix(~XS0) %*% g)
  XS1.g <- as.numeric(model.matrix(~XS1) %*% g)
  XO1.b <- as.numeric(model.matrix(~XO1) %*% b)
  #      u2 <- YO1 - XO1.b
  u2 <- YO1 - XO1.b
  r <- sqrt( 1 - rho^2)
  #      B <- (XS1.g + rho/sigma*u2)/r
  B <- (XS1.g + rho/sigma*u2)/r
  lambdaB <- exp( dnorm( B, log = TRUE ) - pnorm( B, log.p = TRUE ) )
  gradient <- matrix(0, nObs, nParam)
  gradient[YS == 0, ibetaS] <- - w0 * model.matrix(~XS0) *
    exp( dnorm( -XS0.g, log = TRUE ) - pnorm( -XS0.g, log.p = TRUE ) )
  gradient[YS == 1, ibetaS] <- w1 * model.matrix(~XS1) * lambdaB/r
  gradient[YS == 1, ibetaO] <- w1 * model.matrix(~XO1) * (u2/sigma^2 - lambdaB*rho/sigma/r)
  gradient[YS == 1, isigma] <- w1 * ( (u2^2/sigma^3 - lambdaB*rho*u2/sigma^2/r) - 1/sigma )
  gradient[YS == 1, irho] <- w1 * (lambdaB*(u2/sigma + rho*XS1.g))/r^3
  colSums(gradient)
}

hesslik <-function(beta){
  
  NXS <- dim(model.matrix(~XS))[2]#Número de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Número de colunas de XO+1
  nObs <- length(YS)
  NO <- length(YS[YS > 0])
  nParam <- NXS + NXO + 2 #Total of parameters
  
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  N0 <- sum(YS==0)
  N1 <- sum(YS==1)
  
  w  <- rep(1,N0+N1 )
  w0 <- rep(1,N0)
  w1 <- rep(1,N1)
  NXS <- dim(model.matrix(~XS))[2]#Número de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Número de colunas de XO+1
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) {
    return( matrix( NA, nrow = nParam, ncol = nParam ) )
  }
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) {
    return( matrix( NA, nrow = nParam, ncol = nParam ) )
  }
  XS0.g <- as.numeric(model.matrix(~XS0) %*% g)
  XS1.g <- as.numeric(model.matrix(~XS1) %*% g)
  XO1.b <- as.numeric(model.matrix(~XO1) %*% b)
  
  u2 <- YO1 - XO1.b
  r <- sqrt( 1 - rho^2)
  B <- (XS1.g + rho/sigma*u2)/r
  lambdaB <- exp( dnorm( B, log = TRUE ) - pnorm( B, log.p = TRUE ) )
  C <- ifelse(B > -500,
              -exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))*B -
                exp(2 * (dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))),
              -1)
  # recommended by Dimitrios Rizopoulos, KULeuven
  # This is a hack in order to avoid numerical problems.  How to do
  # it better?  How to prove the limit value?
  hess <- matrix(0, nParam, nParam)
  a <- ifelse( XS0.g < 500,
               -exp(dnorm(-XS0.g, log=TRUE) - pnorm(-XS0.g, log.p=TRUE))*XS0.g +
                 ( exp( dnorm(-XS0.g, log=TRUE) - pnorm(-XS0.g, log.p=TRUE)))^2, 1 )
  hess[ibetaS,ibetaS] <- -t(model.matrix(~XS0)) %*% (w0*model.matrix(~XS0)*a) + t(model.matrix(~XS1)) %*% (w1*model.matrix(~XS1)*C)/r^2
  hess[ibetaS,ibetaO] <- -t(model.matrix(~XS1)) %*% (w1*model.matrix(~XO1)*C)*rho/r^2/sigma
  hess[ibetaO,ibetaS] <- t(hess[ibetaS,ibetaO])
  hess[ibetaS,isigma] <- -rho/sigma^2/r^2*t(model.matrix(~XS1)) %*% (w1*C*u2)
  hess[isigma,ibetaS] <- t(hess[ibetaS,isigma])
  hess[ibetaS,irho] <- t(model.matrix(~XS1)) %*%
    (w1*(C*(u2/sigma + rho*XS1.g)/r^4 + lambdaB*rho/r^3))
  hess[irho,ibetaS] <- t(hess[ibetaS,irho])
  hess[ibetaO,ibetaO] <- t(model.matrix(~XO1)) %*%
    (w1*(model.matrix(~XO1) * ((rho/r)^2*C - 1)))/sigma^2
  hess[ibetaO,isigma] <- t(model.matrix(~XO1)) %*%
    (w1*(C*rho^2/sigma^3*u2/r^2 +
           rho/sigma^2*lambdaB/r - 2*u2/sigma^3))
  hess[isigma,ibetaO] <- t(hess[ibetaO,isigma])
  hess[ibetaO,irho] <- t(model.matrix(~XO1)) %*%
    (w1*(-C*(u2/sigma + rho*XS1.g)/r^4*rho -
           lambdaB/r^3))/sigma
  hess[irho,ibetaO] <- t(hess[ibetaO,irho])
  hess[isigma,isigma] <- sum( w1 *
                                ( -3*u2*u2/sigma^4
                                  +2*lambdaB* u2/r *rho/sigma^3
                                  +rho^2/sigma^4 *u2*u2/r^2 *C) ) +
    sum(w1) / sigma^2
  hess[isigma,irho] <- hess[irho,isigma] <-
    -sum(w1*(C*rho*(u2/sigma + rho*XS1.g)/r + lambdaB)*
           u2/sigma^2)/r^3
  hess[irho,irho] <-
    sum(w1*(C*((u2/sigma + rho*XS1.g)/r^3)^2 +
              lambdaB*(XS1.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 ))
  return(hess)
}
