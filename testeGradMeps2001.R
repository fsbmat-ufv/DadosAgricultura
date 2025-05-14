#Limpa todos os objetos na memoria
rm(list=ls())
#Limpa o console
cat("\014")
library(ssmodels)
data(MEPS2001)
data <- MEPS2001
attach(data)
selection <- dambexp ~ age + female + educ + blhisp + totchr + ins + income
outcome <- ambexp ~ age + female + educ + blhisp + totchr + ins
source("HeckmanBS_mod.R")
source("summary.HeckmanBS_mod.R")
theta_BS <- HeckmanBS_mod(selection, outcome, data = MEPS2001)
summary(theta_BS)
theta_BS <- theta_BS$coefficients


if (!inherits(selection, "formula")) stop("'selection' deve ser uma fórmula.")
if (!inherits(outcome, "formula")) stop("'outcome' deve ser uma fórmula.")

mfS <- model.frame(
  formula = selection,
  data = data,
  drop.unused.levels = TRUE,
  na.action = na.pass
)
mtS <- terms(mfS)
XS <- model.matrix(mtS, mfS)
NXS <- ncol(XS)
YS <- model.response(mfS)

mfO <- model.frame(
  formula = outcome,
  data = data,
  drop.unused.levels = TRUE,
  na.action = na.pass
)
mtO <- terms(mfO)
XO <- model.matrix(mtO, mfO)
NXO <- ncol(XO)
YO <- model.response(mfO)

loglik_BS <- function(par) {
  if (any(!is.finite(par))) return(NA)
  
  igamma <- 1:NXS
  ibeta <- seq(tail(igamma, 1) + 1, length.out = NXO)
  iphi1 <- tail(ibeta, 1) + 1
  irho <- iphi1 + 1
  gamma <- par[igamma]
  beta <- par[ibeta]
  phi1 <- par[iphi1]
  if (!is.finite(phi1) || phi1 < 1e-3 || phi1 > 100) return(NA)
  rho_star <- par[irho]
  rho <- 2 / (1 + exp(-rho_star)) - 1
  phi2 <- 1
  XS0 <- XS[YS == 0, , drop = FALSE]
  XS1 <- XS[YS == 1, , drop = FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS == 1]
  XO1 <- XO[YS == 1, , drop = FALSE]
  XS0.g <- exp(as.numeric((XS0) %*% gamma))
  XS1.g <- exp(as.numeric((XS1) %*% gamma))
  XO1.b <- exp(as.numeric((XO1) %*% beta))
  term0 <- sqrt(YO1 * (phi1 + 1)/(phi1 * XO1.b)) - sqrt((phi1 * XO1.b)/(YO1 * (phi1 + 1)))
  term1 <- exp((-phi1/4) * (term0)^2)
  term2 <- sqrt((phi1 + 1)/(phi1 * XO1.b * YO1)) + sqrt((phi1 * XO1.b)/((phi1 + 1) * (YO1^3)))
  term3 <- (1/(2 * sqrt(2 * pi))) * sqrt(phi1/2)
  term4 <- sqrt((phi2 + 1)/(2 * XS1.g * (1 - rho^2)))
  term5 <- ((phi2 * XS1.g)/(phi2 + 1)) - 1
  term6 <- rho * sqrt(phi1/(2 * (1 - rho^2)))
  integrand <- term4 * term5 + term6 * term0
  term7 <- pnorm(integrand, log.p = TRUE)
  term8 <- sqrt(phi2/2) * (sqrt((phi2 + 1)/(phi2 * XS0.g)) - sqrt((phi2 * XS0.g)/(phi2 + 1)))
  FT2 <- pnorm(term8, log.p = TRUE)
  term1[term1 < 1e-300] <- 1e-300
  logterm1 <- log(term1)
  logterm2 <- log(term2)
  if (any(!is.finite(logterm1)) || any(!is.finite(logterm2))) return(NA)
  ll <- sum(logterm1 + logterm2 + log(term3) + term7) + sum(FT2)
  return(sum(ll))
}

gradlik_BS <- function(par) {
  igamma <- 1:NXS
  ibeta <- seq(tail(igamma, 1)+1, length.out = NXO)
  iphi1 <- tail(ibeta, 1) + 1
  irho <- iphi1 + 1
  
  gamma <- par[igamma]
  beta <- par[ibeta]
  phi1 <- par[iphi1]
  rho_star <- par[irho]
  rho <- 2 / (1 + exp(-rho_star)) - 1
  phi2 <- 1
  
  XS0 <- XS[YS==0,,drop=FALSE]
  XS1 <- XS[YS==1,,drop=FALSE]
  YO[is.na(YO)] <- 0
  YO1 <- YO[YS==1]
  XO1 <- XO[YS==1,,drop=FALSE]
  
  XS0.g <- exp(as.numeric(XS0 %*% gamma))
  XS1.g <- exp(as.numeric(XS1 %*% gamma))
  XO1.b <- exp(as.numeric(XO1 %*% beta))
  
  temp_rho2 <- pmax(1 - rho^2, 1e-6)
  term0 <- sqrt(YO1*(phi1+1)/(phi1*XO1.b)) - sqrt((phi1*XO1.b)/(YO1*(phi1+1)))
  term1 <- exp((-phi1/4)*(term0)^2)
  term2 <- sqrt((phi1+1)/(phi1*XO1.b*YO1)) + sqrt((phi1*XO1.b)/((phi1+1)*YO1^3))
  term3 <- (1/(2*sqrt(2*pi)))*sqrt(phi1/2)
  term4 <- sqrt((phi2+1)/(2*XS1.g*temp_rho2))
  term5 <- (phi2*XS1.g)/(phi2+1) - 1
  term6 <- rho * sqrt(phi1 / (2 * temp_rho2))
  integrand <- term4*term5 + term6*term0
  term7 <- pnorm(integrand, log.p = TRUE)
  term8 <- sqrt(phi2 / 2) * (sqrt((phi2+1)/(phi2*XS0.g)) - sqrt((phi2*XS0.g)/(phi2+1)))
  FT2 <- pnorm(term8, log.p = TRUE)
  term9 <- (-1/2) * (sqrt(YO1*(phi1+1)/(phi1*XO1.b)) + sqrt((phi1*XO1.b)/(YO1*(phi1+1)))) * model.matrix(~XO1 - 1)
  term10 <- (1/2) * (sqrt((XO1.b*phi1)/(YO1^3*(phi1+1))) - sqrt((phi1+1)/(YO1*phi1*XO1.b))) * model.matrix(~XO1 - 1)
  term11 <- (-1/2) * sqrt(phi2/2) * (sqrt((phi2+1)/(phi2*XS0.g)) + sqrt((phi2*XS0.g)/(phi2+1))) * model.matrix(~XS0 - 1)
  term12 <- (-1/2) * (sqrt(YO1/(phi1^3*XO1.b*(phi1+1))) + sqrt(XO1.b/(YO1*(phi1+1)^3*phi1)))
  term13 <- (1/2) * (sqrt(XO1.b/(phi1*(phi1+1)^3*YO1^3)) - sqrt(1/(phi1^3*(phi1+1)*XO1.b*YO1)))
  term14 <- rho / (2 * sqrt(2 * phi1 * temp_rho2))
  lambda_I <- exp(dnorm(integrand, log = TRUE) - pnorm(integrand, log.p = TRUE))
  lambda_T8 <- exp(dnorm(term8, log = TRUE) - pnorm(term8, log.p = TRUE))
  df1b <- (-phi1/2) * term0 * term9
  df2b <- term10 / term2
  df4b <- lambda_I * term9 * term6
  df4g <- lambda_I * (1/4) * ((XS1.g + 2) / sqrt(XS1.g * temp_rho2)) * model.matrix(~XS1 - 1)
  df5g <- lambda_T8 * term11
  df1phi1 <- - (term0^2)/4 - (phi1/2) * term0 * term12
  df2phi1 <- term13 / term2
  df3phi1 <- 1 / (8 * term3 * sqrt(pi * phi1))
  df4phi1 <- lambda_I * (term0 * term14 + term6 * term12)
  dt5_drho <- if (abs(rho) < 1e-8) {
    sqrt(phi1 / 2)
  } else {
    sqrt(phi1 / (2 * temp_rho2)) + 
      rho^2 * sqrt(phi1) / ((temp_rho2)^(3/2) * sqrt(2))
  }
  d_rho <- (4 * exp(-rho_star)) / (1 + exp(-rho_star))^2
  df4rho_raw <- 0.5 * lambda_I * (term5 * term4 * (rho / temp_rho2) + dt5_drho * term0)
  gr_rho_star <- df4rho_raw * d_rho
  gradient <- matrix(0, nrow = length(YS), ncol = length(par))
  gradient[YS == 0, igamma] <- df5g
  gradient[YS == 1, igamma] <- df4g
  gradient[YS == 1, ibeta]  <- df1b + df2b + df4b
  gradient[YS == 1, iphi1]  <- df1phi1 + df2phi1 + df3phi1 + df4phi1
  gradient[YS == 1, irho]   <- gr_rho_star
  return(colSums(gradient))
}

library(numDeriv)

# Avaliação numérica do gradiente
grad_num <- grad(func = loglik_BS, x = theta_BS)

# Avaliação analítica do gradiente
grad_ana <- gradlik_BS(theta_BS)

# Comparação lado a lado
cbind(Gradiente_Analitico = round(grad_ana, 6),
      Gradiente_Numerico  = round(grad_num, 6),
      Diferenca           = round(grad_ana - grad_num, 6))
