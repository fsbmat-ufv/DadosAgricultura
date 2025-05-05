loglik_SK <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Numero de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Numero de colunas de XO+1
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  ilamb1 <- tail(irho, 1) + 1
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) return(NA)
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) return(NA)
  lamb1 <- beta[ilamb1]
  XS.g <- model.matrix(~XS) %*% g
  XO.b <- model.matrix(~XO) %*% b
  u2 <- YO - XO.b
  z <- u2/sigma 
  r <- sqrt( 1 - rho^2)
  u <- (1+(lamb1^2)-((lamb1*rho)^2))
  lstar <- -(lamb1*rho)/sqrt(u)
  # lstar <- (-lamb1*rho)/(sqrt(1+(lamb1^2)-((lamb1*rho)^2)))
  # nuc = function(t){2*dnorm(t)*pnorm(-lstar*t)}
  # h=Vectorize(function(ls) integrate(nuc,-Inf, -ls)$value)
  h <- sn::psn(-XS.g,0,1,-lstar)
  ll <- ifelse(YS == 0,
               log(h),
               log(2/sigma)+log(dnorm(z))+log(pnorm(lamb1*z))+
                 log(pnorm((XS.g+rho*z)/r))
  )
  sum(ll)
}

gradlik_SK <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Numero de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Numero de colunas de XO+1
  nObs <- length(YS)
  NO <- length(YS[YS > 0])
  nParam <- NXS + NXO + 3 #Total of parameters
  
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
  
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  ilamb1 <- tail(irho, 1) + 1
  
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) return(matrix(NA, nObs, nParam))
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
  lamb1 <- beta[ilamb1]
  XS0.g <- as.numeric(model.matrix(~XS0) %*% g)
  XS1.g <- as.numeric(model.matrix(~XS1) %*% g)
  XO1.b <- as.numeric(model.matrix(~XO1) %*% b)
  #      u2 <- YO1 - XO1.b
  u2 <- YO1 - XO1.b
  z <- u2/sigma
  r <- (1 - rho^2)
  u <- (1+(lamb1^2)-((lamb1*rho)^2))
  lstar <- -(lamb1*rho)/sqrt(u)
  dlrho=(-rho*(1+(lamb1^2)-(lamb1*rho)^2)+lamb1*rho*(lamb1-lamb1*(rho^2)))/((1+(lamb1^2)-(lamb1*rho)^2))^(3/2)
  dllam=(-rho/sqrt(u))+lamb1*rho*(lamb1-lamb1*(rho^2))/(sqrt(u)*u)
  omeg <- (XS1.g+rho*z)/sqrt(r)
  K1 <- dnorm(omeg)/pnorm(omeg)
  h <- function(t){sn::psn(-t,0,1,-lstar)}
  h2<- function(t){sn::psn(-t,0,1,lamb1*rho/sqrt(u))}
  K2 <- dnorm(-XS0.g)*pnorm(lstar*XS0.g)/h(XS0.g)
  eta <- lamb1*z
  K3 <- dnorm(eta)/pnorm(eta)
  K4 <- (dnorm((XS0.g*(sqrt(1+lamb1^2)))/(sqrt(u))))/h(XS0.g)
  gradient <- matrix(0, nObs, nParam)
  gradient[YS == 0, ibetaS] <-  w0 * model.matrix(~XS0)*(-2*K2)
  gradient[YS == 1, ibetaS] <-  w1 * model.matrix(~XS1) * (K1/sqrt(r))
  gradient[YS == 1, ibetaO] <-  w1 * model.matrix(~XO1) * ((z/sigma)-((lamb1/sigma)*K3)-((rho/(sigma*sqrt(r)))*K1))
  gradient[YS == 1, isigma] <-  w1 * ((-1/sigma)+((z^2)/sigma)-((lamb1*K3*z)/sigma)-((rho*K1*z)/(sigma*sqrt(r))))
  gradient[YS == 0, irho]   <-  w0 * ((-(2/sqrt(2*pi))*(lamb1*(1+lamb1^2)/(sqrt(u)*(u+(lamb1*rho)^2)))*dnorm(sqrt(1+(lamb1*rho/sqrt(u))^2)*XS0.g))/h2(XS0.g))
  gradient[YS == 1, irho]   <-  w1 * ((K1*(rho*XS1.g+z))/((sqrt(r))^3))
  gradient[YS == 0, ilamb1] <-  w0 * (dllam*sqrt(2/pi)*(1/(1+lstar^2))*dnorm(sqrt(1+lstar^2)*(XS0.g))/h(XS0.g))
  gradient[YS == 1, ilamb1] <-  w1 * (K3*z)
  colSums(gradient)
}
