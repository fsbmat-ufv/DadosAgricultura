loglik_BS <- function(par){
  n <- length(YO)
  NXS <- dim(model.matrix(~XS))[2]#N?mero de colunas da matriz de sele??o
  NXO <- dim(model.matrix(~XO))[2]#N?mero de colunas da matriz de regress?o prim?ria
  ## parameter indices
  igamma <- 1:NXS
  ibeta <- seq(tail(igamma, 1)+1, length=NXO)
  iphi1 <- tail(ibeta, 1) + 1
  irho <- tail(iphi1, 1) + 1
  gamma <- par[igamma]
  beta <- par[ibeta]
  phi1 <- par[iphi1]
  if(phi1 < 0) return(NA)
  rho <- par[irho]
  if( ( rho < -1) || ( rho > 1)) return(NA)
  phi2=1
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  N0 <- sum(YS==0)
  N1 <- sum(YS==1)
  
  XS0.g <- exp(as.numeric(model.matrix(~XS0) %*% gamma))
  XS1.g <- exp(as.numeric(model.matrix(~XS1) %*% gamma))
  XO1.b <- exp(as.numeric(model.matrix(~XO1) %*% beta))
  
  term0 <- ((YO1*(phi1+1)/(phi1*XO1.b))^(1/2)-((phi1*XO1.b)/(YO1*(phi1+1)))^(1/2))
  term1 <- exp((-phi1/4)*(term0)^2)
  term2 <- (((phi1+1)/(phi1*XO1.b*YO1))^(1/2)+((phi1*XO1.b)/((phi1+1)*(YO1^3)))^(1/2))
  term3 <- (1/(2*sqrt(2*pi)))*((phi1/2)^(1/2))
  term4 <- (((phi2+1))/(2*XS1.g*(1-rho^2)))^(1/2)
  term5 <- ((phi2*XS1.g)/(phi2+1))-1
  term6 <- rho*(phi1/(2*(1-rho^2)))^(1/2)
  integrand <- term4*term5+term6*term0
  term7 <- pnorm(integrand, log.p = TRUE)
  term8 <- ((phi2/2)^(1/2))*(((phi2+1)/(phi2*XS0.g))^(1/2)-((phi2*XS0.g)/(phi2+1))^(1/2))
  FT2 <- pnorm(term8, log.p = TRUE)
  term1[term1 < 1e-300] <- 1e-300
  term2[term2 < 1e-300] <- 1e-300
  
  logterm1 <- log(term1)
  logterm2 <- log(term2)
  
  if (any(!is.finite(logterm1)) || any(!is.finite(logterm2))) return(NA)
  
  ll <- sum(logterm2 + logterm1 + log(term3) + term7) + sum(FT2)
  (sum(ll))
}

gradlik_BS <- function(par) {
  n <- length(YO)
  NXS <- dim(model.matrix(~XS))[2]#N?mero de colunas da matriz de sele??o
  NXO <- dim(model.matrix(~XO))[2]#N?mero de colunas da matriz de regress?o prim?ria
  ## parameter indices
  igamma <- 1:NXS
  ibeta <- seq(tail(igamma, 1)+1, length=NXO)
  iphi1 <- tail(ibeta, 1) + 1
  irho <- tail(iphi1, 1) + 1
  gamma <- par[igamma]
  beta <- par[ibeta]
  phi1 <- par[iphi1]
  
  nObs <- length(YS)
  NO <- length(YS[YS > 0])
  nParam <- NXS + NXO + 2 #Total of parameters
  
 # if(phi1 < 0) return(matrix(NA, nObs, nParam))
  rho <- par[irho]
 # if( ( rho < -1) || ( rho > 1)) return(matrix(NA, nObs, nParam))
  phi2=1
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  N0 <- sum(YS==0)
  N1 <- sum(YS==1)
  
  XS0.g <- exp(as.numeric(model.matrix(~XS0) %*% gamma))
  XS1.g <- exp(as.numeric(model.matrix(~XS1) %*% gamma))
  XO1.b <- exp(as.numeric(model.matrix(~XO1) %*% beta))
  
  w  <- rep(1,N0+N1 )
  w0 <- rep(1,N0)
  w1 <- rep(1,N1)
  
  mu1 <- exp(as.numeric(model.matrix(~XO) %*% beta))
  mu2 <- exp(as.numeric(model.matrix(~XS) %*% gamma))
  
  term0 <- ((YO1*(phi1+1)/(phi1*XO1.b))^(1/2)-((phi1*XO1.b)/(YO1*(phi1+1)))^(1/2))
  term1 <- exp((-phi1/4)*(term0)^2)
  term2 <- (((phi1+1)/(phi1*XO1.b*YO1))^(1/2)+((phi1*XO1.b)/((phi1+1)*(YO1^3)))^(1/2))
  term3 <- (1/(2*sqrt(2*pi)))*((phi1/2)^(1/2))
  term4 <- (((phi2+1))/(2*XS1.g*(1-rho^2)))^(1/2)
  term5 <- ((phi2*XS1.g)/(phi2+1))-1
  term6 <- rho*(phi1/(2*(1-rho^2)))^(1/2)
  integrand <- term4*term5+term6*term0
  term7 <- pnorm(integrand, log.p = TRUE)
  term8 <- ((phi2/2)^(1/2))*(((phi2+1)/(phi2*XS0.g))^(1/2)-((phi2*XS0.g)/(phi2+1))^(1/2))
  FT2 <- pnorm(term8, log.p = TRUE)
  term9 <-  ((-1/2)*(((YO1*(phi1+1))/(phi1*(XO1.b)))^(1/2)+((phi1*XO1.b)/(YO1*(phi1+1)))^(1/2)))*model.matrix(~XO1) #Derivada de term0 em rela??o a beta
  term10 <-           (1/2)*(((XO1.b*phi1)/((YO1^3)*(phi1+1)))^(1/2)-((phi1+1)/(YO1*phi1*XO1.b))^(1/2))*model.matrix(~XO1) #Derivada de term2 em rela??o a beta
  term11 <- ((-1/2)*(sqrt(phi2/2))*(((phi2+1)/(phi2*(XS0.g)))^(1/2)+((phi2*XS0.g)/((phi2+1)))^(1/2)))*model.matrix(~XS0) #Derivada de term8 em rela??o a gamma
  term12 <- (-1/2)*((YO1/((phi1^3)*XO1.b*(phi1+1)))^(1/2)+(XO1.b/(YO1*((phi1+1)^3)*phi1))^(1/2)) #Derivada de term0 em rela??o a phi1
  term13 <- (1/2)*((XO1.b/(phi1*((phi1+1)^3)*(YO1^3)))^(1/2)-(1/((phi1^3)*(phi1+1)*XO1.b*YO1))^(1/2))#Derivada de term2 em rela??o a phi1
  term14 <- rho/(2*((2*phi1*(1-rho^2))^(1/2)))#Derivada de term6 em rela??o a phi1
  term15 <- (rho/(1-rho^2))*term4 #Derivada de term4 em rela??o a rho
  term16 <- (1/(1-rho^2))*((phi1/(2*(1-rho^2)))^(1/2))#Derivada de term6 em rela??o a rho
  #f1=log(term1), f2=log(term2), f3=log(term3), f4=term7, f5=FT2
  lambda_I <- exp( dnorm(integrand, log = TRUE ) - pnorm(integrand, log.p = TRUE ) )
  lambda_T8 <- exp( dnorm(term8, log = TRUE ) - pnorm(term8, log.p = TRUE ) )
  df1b <- (-phi1/2)*term0*term9 #Derivada de log(T1) em rela??o a beta
  df2b <- (term10/term2) #Derivada de log(T2) em rela??o a beta0
  df4b <- lambda_I*term9*term6 #Derivada de log(T7) em rela??o a beta

  df4g <- lambda_I*((1/4)*(((XS1.g+2)/sqrt(XS1.g*(1-rho^(2))))))*model.matrix(~XS1) #Derivada de f4 em rela??o a gamma

  df5g <- lambda_T8*term11 #Derivada de FT2 em rela??o a gamma
  df1phi1<- (-(term0^2)/4)-(phi1/2)*(term0*term12) #Derivada do log(T1) em rela??o a phi1
  df2phi1 <- (1/term2)*term13
  df3phi1 <- 1/(8*term3*((pi*phi1)^(1/2)))
  df4phi1 <- lambda_I*(term0*term14+term6*term12) 
 
  df4rho  <- lambda_I*(term5*term4*(rho/(1-rho^2))+(term0*term6/(rho*(1-rho^2))))
  gradient <- matrix(0, nObs, nParam)
  gradient[YS == 0, igamma] <- w0 *  df5g
  gradient[YS == 1, igamma] <- w1 * (df4g)
  gradient[YS == 1, ibeta]  <- w1 * (df1b+df2b+df4b)
  gradient[YS == 1, iphi1]  <- w1 * (df1phi1+df2phi1+df3phi1+df4phi1)
  gradient[YS == 1, irho]   <- w1 * (df4rho)
  colSums(gradient)
}




