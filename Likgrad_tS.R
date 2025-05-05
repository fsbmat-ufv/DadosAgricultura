loglik_tS <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Numero de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Numero de colunas de XO+1
  ## parameter indices
  ibetaS <- 1:NXS
  ibetaO <- seq(tail(ibetaS, 1)+1, length=NXO)
  isigma <- tail(ibetaO, 1) + 1
  irho <- tail(isigma, 1) + 1
  iv <- tail(irho, 1) + 1
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  sigma <- beta[isigma]
  if(sigma < 0) return(NA)
  rho <- beta[irho]
  if( ( rho < -1) || ( rho > 1)) return(NA)
  v <- beta[iv]
  XS.g <- model.matrix(~XS) %*% g
  XO.b <- model.matrix(~XO) %*% b
  u2 <- YO - XO.b
  z <- u2/sigma
  r <- sqrt( 1 - rho^2)
  Q <- ((v+1)/(v+(z^2)))^(1/2)
  eta <- Q*((rho*z+XS.g)/r)
  gam <- log(gamma((v+1)/2))-log(gamma(v/2))-0.5*log(pi)-0.5*log(v)-log(sigma)
  ll <- ifelse(YS == 0,
               (pt(-XS.g,v,log.p=TRUE)),
               (gam-((v+1)/2)*log(1+((z^2)/v))+pt(eta,v+1,log.p=TRUE))
  )
  sum(ll)
}

gradlik_tS <- function(beta) {
  NXS <- dim(model.matrix(~XS))[2]#Numero de colunas de XS+1
  NXO <- dim(model.matrix(~XO))[2]#Numero de colunas de XO+1
  nObs <- length(YS)
  NO <- length(YS[YS > 0])
  nParam <- NXS + NXO + 3 #Total of parameters
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
  iv <- tail(irho, 1) + 1
  g <- beta[ibetaS]
  b <- beta[ibetaO]
  
  chuteS=beta[isigma]
  lns=log(chuteS)
  sigma=exp(lns)
  chuteR <- beta[irho]
  tau=log((1+chuteR)/(1-chuteR))/2
  rho=(exp(2*tau)-1)/(exp(2*tau)+1)
  chuteV=beta[iv]
  lndf=log(chuteV)
  v=exp(lndf)
  
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  XS0.g <- as.numeric(model.matrix(~XS0) %*% g)
  XS1.g <- as.numeric(model.matrix(~XS1) %*% g)
  XO1.b <- as.numeric(model.matrix(~XO1) %*% b)
  #      u2 <- YO1 - XO1.b
  #u2S1 <- YO1 - XS1.g
  u2O <- YO1 - XO1.b
  #u2 <- YO - XO.b
  z0 <- u2O/sigma
  r <- sqrt( 1 - rho^2)
  Qv <- ((v+1)/(v+(z0)^2))^(1/2)
  Ar <- 1/sqrt(1-(rho^2))
  Arr <- rho*Ar
  QSI <- (Arr*z0)+Ar*XS1.g
  eta <- Qv*QSI
  dr=4*exp(2*tau)/((exp(2*tau)+1)^2)
  
  #tau <- XS1.g/r
  # myenv <- new.env()
  # assign("v", v, envir = myenv)
  # #assign("XS0.g",XS0.g,envir = myenv)
  # gv2 <-numericDeriv(quote(pt(-XS0.g, v,log.p=TRUE)), c("v"), myenv)
  
  f1 <- function(x){
    ff <- pt(-XS0.g,x,log.p=TRUE)
    return(ff)
  }
  
  gv2 <- numDeriv::grad(f1, rep(1,length(XS0.g))*v)
  #sum(gv2)
  # myenv2 <- new.env()
  # assign("v", v, envir = myenv2)
  # assign( "z0",z0, envir = myenv2)
  # assign( "QSI",QSI, envir = myenv2)
  # f <- quote(pt((((v+1)/(v+(z0)^2))^(1/2))*QSI,v+1,log.p=TRUE))
  # gv <- numericDeriv(f, c("v"), myenv2)
  
  f2 <- function(x){
  teste=(((x+1)/(x+(z0)^2))^(1/2))*QSI
  ff <- pt(teste,(x+1),log.p=TRUE)
  return(ff)
   }

  gv <- numDeriv::grad(f2, rep(1,length(z0))*v)
  
  gradient <- matrix(0, nObs, nParam)
  gradient[YS == 0, ibetaS] <- -w0 * model.matrix(~XS0)*(dt(-XS0.g,v)/pt(-XS0.g,v))
  gradient[YS == 1, ibetaS] <-  w1 * model.matrix(~XS1) * Ar * Qv*(dt(eta,v+1)/pt(eta,v+1))
  gradient[YS == 1, ibetaO] <-  w1 * model.matrix(~XO1) * (Qv/sigma)*
    ((Qv*z0)+(((QSI*((v+(z0^2))^(-1)))*z0-Arr)*(dt(eta,v+1)/pt(eta,v+1))))
  gradient[YS == 1, isigma] <- w1 * (-1+(Qv*z0)^2+(dt(eta,v+1)/pt(eta,v+1))*((Qv*(z0))*
                                                                               (QSI*((v+(z0^2))^(-1))*(z0)-Arr)))
  gradient[YS == 1, irho] <- w1 * (dt(eta,v+1)/pt(eta,v+1))*Qv*(Ar^(3))*(z0+rho*XS1.g)*dr
  # gradient[YS == 1, iv] <- w1 * (((1/2)*v*(digamma((v+1)/2)-digamma(v/2))-(1/2))-((1/2)*v*log(1+((z0^2)/v)))+
  #                                 (((Qv*(z0))^(2))/(2))+v*gv)
  # gradient[YS == 0, iv] <- w0*v*gv2
  # colSums(gradient)
  gradient[YS == 1, iv] <- w1 * ((((1/2)*digamma((v+1)/2)-(1/2)*digamma(v/2))-(1/(2*v)))-((1/2)*log(1+((z0^2)/v)))+
                                   (((Qv*(z0))^(2))/(2*v))+gv)
  gradient[YS == 0, iv] <- w0*gv2
  colSums(gradient)
}

# gradient[YS == 1, iv] <- w1 * ((((1/2)*digamma((v+1)/2)-(1/2)*digamma(v/2))-(1/(2*v)))-((1/2)*log(1+((z0^2)/v)))+
#                                  (((Qv*(z0))^(2))/(2*v))+gv)
# gradient[YS == 0, iv] <- w0*gv2
# colSums(gradient)

#digamma((100000+1)/2)
