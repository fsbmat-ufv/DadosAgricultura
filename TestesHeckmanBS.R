#########################################################################################
#PROGRAMA: DadosAGB.R                                                                  #
#                                                                                       #
# AUTOR: Fernando de Souza Bastos                                                       #                                                                                       #
# E-MAIL: fsbmat@gmail.com                                                              #                                                                                      #
# DATA: Maio/2025                                                                       #
#########################################################################################
#Limpa todos os objetos na memoria
rm(list=ls())
#Limpa o console
cat("\014")
#pacote para manipulacao de data.frames
library(tidyverse) 
#Pacote para construcao de tabelas
library(knitr) 
#Pacote para ajuste do modelo de selecao de Heckman e para ajuste de dois passos
#de Heckman
library("sampleSelection")
library(ssmodels)
#Pacote que contem o banco de dados MEPS2001
#library(ssmrob)
#Leitura do dados MEPS2001
data(MEPS2001)
df <- readRDS("AGB.Rds")
#tornando visiveis as colunas do data-frame
attach(MEPS2001)
attach(df)
#Local do arquivo
#setwd("~/GitHub/Tese/R/2020/BS/MEPS")
#Ajuste dos dois passos de Heckman
selectEq  <- sel ~ DAP + DENS + time_since_census
outcomeEq <- logAGB ~ DAP + DENS + b + AltCDano + Local
library(car)


# Suponha que você tenha um modelo de regressão da equação de resultado:
reg_outcome <- lm(logAGB ~ DAP + DENS + Prop + AltCDano + Local,
                  data = df)

vif(reg_outcome)
#summary(time_since_census)
set.seed(123)

base0 <- df[df$sel == 0, ]
base1 <- df[df$sel == 1, ]

# amostragem sobre as posições corretas (linha 1 até nrow(base1))
amostrados <- sample(1:nrow(base1), size = 5336)

df2 <- rbind(base0, base1[amostrados, ])

table(df2$sel, useNA = "ifany")  # deve dar apenas 0 e 1
anyNA(df2$sel)                   # deve ser FALSE


df2 <- df2 %>%
  mutate(logDAP          = log(DAP), 
    logDAP_z             = scale(log(DAP)),
    DENS_z               = scale(DENS),
    AltDAP_z             = scale(AltDAP),
    time_since_census_z  = scale(time_since_census),
    b_z                  = scale(PropVol))

selection_model <- glm(sel ~ logDAP_z + DENS_z + time_since_census_z,
                       family = binomial(link = "probit"),
                       data = df2)
hist(fitted(selection_model), breaks = 100, main = "Fitted values - Seleção (nova)", xlab = "P(seleção)")
table(round(fitted(selection_model), 2))


twopassos <- heckit(sel ~ logDAP_z + DENS_z + time_since_census_z,
                    logAGB ~ logDAP_z + DENS_z + AltCDano + Local, 
                    df2,
                    method = "2step")
#twopassos Local#twopassos <- heckit( dambexp ~ age+female+educ+blhisp+totchr+ins+income,
#                     lnambx ~ age+female+educ+blhisp+totchr+ins, 
#                     MEPS2001,
#                     method = "2step")
#Resultados
summary(twopassos)

sel ~ logDAP_z + DENS_z + time_since_census_z
logAGB ~ logDAP_z + DENS_z + AltCDano + Local


heckman_bal <- selection(
  selection = sel ~ logDAP_z + DENS_z + time_since_census_z,
  outcome   = logAGB ~ logDAP_z + DENS_z + AltCDano + Local,
  data = df2,
  method = "ml"
)

summary(heckman_bal)


#Ajuste do modelo de selecao de Heckman com estimacao via Maxima Verossimilhanca
fit<-selection(sel ~ logDAP_z + DENS_z + time_since_census_z,
               logAGB ~ logDAP_z + DENS_z + AltCDano + Local,
               df2)
#fit<-selection( dambexp ~ age+female+educ+blhisp+totchr+ins+income,
#                lnambx ~ age+female+educ+blhisp+totchr+ins, MEPS2001)
#Estimativas dos parametros 
beta<- fit$estimate
summary(fit)
#Pacote para ajuste do modelo probit
library(aod)
#Modelo probit ajustado a variavel Y2
sel <- df2$sel
sum(is.na(fit1$linear.predictors[sel == 1]))  # Agora funciona corretamente

fit1<-glm(sel ~ logDAP_z + DENS_z + time_since_census_z,
          family = binomial(link = "probit"),
          data=df2)
#Geracao de valores da covariavel razao inversa de Mills (IMR)
IMR <- dnorm(fit1$linear.predictors)/pnorm(fit1$linear.predictors)

#Acrescimo de IMR ao dataframe
dt <- data.frame(df2,IMR)

#Modelo lm(regressao multipla) ajustado a variavel Y1
fit2 <- lm(formula = logAGB ~ logDAP_z + DENS_z + AltCDano + IMR, data = dt[dt$sel==1, ])
#Resultado do ajuste da funcao lm
summary(fit2)
#Geracao de valores da nova covariavel delta 
delta <- (dt$IMR)*(dt$IMR+fit1$fitted.values)
sum(delta[sel==1])
sum(is.na(delta[sel == 1]))

#Acrescimo de delta ao dataframe
dt <- data.frame(dt,delta)
#Quantidade de valores u==1
q <- sum(sel)
#Calculo da variancia de Y1
Var <- (sum((fit2$residuals)^2)+(((coef(fit2)[5])^2)*sum(delta[sel==1])))/q
#Calculo da correlacao entre Y1 e Y2
cor <- (coef(fit2)[5])/sqrt(Var)
#Chute inicial para optim do modelo BS
par <- c(coef(fit1),coef(fit2)[-5],phi=sqrt(Var),cor)
names(par) <- c("(Intercept)",
                "logDAP_z",
                "DENS_z",
                "time_since_census_z",
                "(Intercept)",
                "logDAP_z", 
                "DENS_z",
                "AltCDano",
                "sigma", 
                "rho")
#Variavel indicadora u=1(Y2>1)
YS <- df2$sel
#YS <- u
#Variavel resposta para ajuste do modelo BS
YO <- df2$AGB
#YO <- y
#Matriz de covariaveis para calculo da media mu1 (Equacao primaria)
XO <- cbind(df2$logDAP_z, df2$DENS_z, df2$AltCDano)
#XO <- X
#sel ~ logDAP_z + DENS_z + time_since_census_z
#logAGB ~ logDAP_z + DENS_z + AltCDano + Local 
#Matriz de covariaveis para calculo da media mu2 (Equacao de selecao)
XS <- cbind(df2$logDAP_z, df2$DENS_z, df2$time_since_census)
#XS <- W
#Funcao log de verossimilhanca e gradiente dos modelos de selecao 
source("Likgrad_BS.R")
source("Likgrad_HC.R")
source("Likgrad_Skew.R")
source("Likgrad_tS.R")
#Funcao Optim para estimacao dos parametros do modelo
theta_BS <- optim(par, loglik_BS, gradlik_BS, method = "BFGS",hessian = T,control = list(fnscale=-1))
theta_BS$par
library(ssmodels)
selectEq  <- sel ~ logDAP_z + DENS_z + time_since_census_z
outcomeEq <- AGB ~ logDAP_z + DENS_z + AltCDano
source("HeckmanBS_mod.R")
source("summary.HeckmanBS_mod.R")
modelo_bs <- HeckmanBS_mod(
  selection = sel ~ logDAP_z + DENS_z + time_since_census_z,
  outcome   = AGB ~ logDAP_z + DENS_z + AltCDano+Local,
  data      = df2)
modelo_bs$coefficients
summary(modelo_bs)
#Variavel resposta para ajuste dos demais modelos
YO <- lnambx
theta_HC <- optim(par,loglik_HC,gradlik_HC,method = "BFGS",hessian = T,control = list(fnscale=-1))
par_SK <- c(coef(fit1),coef(fit2)[-8],phi=sqrt(Var),cor,1)
names(par_SK) <- c("(Intercept)","Age",
                   "Female","Educ","blhisp",
                   "totchr","ins","income","(Intercept)", "age","female",
                   "educ", "blhisp","totchr",
                   "ins","sigma", "rho","lambda")
theta_SK <- optim(par_SK,loglik_SK,gradlik_SK,method = "BFGS",hessian = T,control = list(fnscale=-1))
par_tS <- c(coef(fit1),coef(fit2)[-8],phi=sqrt(Var),cor,10)
names(par_tS) <- c("(Intercept)","Age",
                   "Female","Educ","blhisp",
                   "totchr","ins","income","(Intercept)", "age","female",
                   "educ", "blhisp","totchr",
                   "ins","sigma", "rho","gl")
theta_tS <- optim(par_tS,loglik_tS,gradlik_tS,method = "BFGS",hessian = T,control = list(fnscale=-1))
#Matriz de informa??o de Fisher
fisher_infoBS <- solve(-theta_BS$hessian)
fisher_infoHC <- solve(-theta_HC$hessian)
fisher_infoSK <- solve(-theta_SK$hessian)
fisher_infotS <- solve(-theta_tS$hessian)
#Raiz quadrada dos elementos da diagonal da matriz de informa??o de Fisher
prop_sigmaBS <- sqrt(diag(fisher_infoBS))
prop_sigmaHC <- sqrt(diag(fisher_infoHC))
prop_sigmaSK <- sqrt(diag(fisher_infoSK))
prop_sigmatS <- sqrt(diag(fisher_infotS))
#Intervalos de confian?a de 95%
upperBS<-theta_BS$par+ qnorm(0.975)*prop_sigmaBS
lowerBS<-theta_BS$par+ qnorm(0.025)*prop_sigmaBS
upperHC<-theta_HC$par+ qnorm(0.975)*prop_sigmaHC
lowerHC<-theta_HC$par+ qnorm(0.025)*prop_sigmaHC
upperSK<-theta_SK$par+ qnorm(0.975)*prop_sigmaSK
lowerSK<-theta_SK$par+ qnorm(0.025)*prop_sigmaSK
uppertS<-theta_tS$par+ qnorm(0.975)*prop_sigmatS
lowertS<-theta_tS$par+ qnorm(0.025)*prop_sigmatS
#Estimativas dos par?metros
coeffsBS <- theta_BS$par
coeffsHC <- theta_HC$par
coeffsSK <- theta_SK$par
coeffstS <- theta_tS$par
# fisher_info <- solve(-theta$hessian)
#stderr <- sqrt(diag(fisher_info))
#Valor-z
zscoreBS <- coeffsBS/prop_sigmaBS
zscoreHC <- coeffsHC/prop_sigmaHC
zscoreSK <- coeffsSK/prop_sigmaSK
zscoretS <- coeffstS/prop_sigmatS
#P-valor para uso no teste de Wald
#H0: theta_i=0, para todo i
#Ha: nao H0
pvalueBS <- 2*(1 - pnorm(abs(zscoreBS)))
pvalueHC <- 2*(1 - pnorm(abs(zscoreHC)))
pvalueSK <- 2*(1 - pnorm(abs(zscoreSK)))
pvaluetS <- 2*(1 - pnorm(abs(zscoretS)))
#Matriz de resultados
resultsBS <- cbind(coeffsBS,prop_sigmaBS,zscoreBS,pvalueBS)
resultsHC <- cbind(coeffsHC,prop_sigmaHC,zscoreHC,pvalueHC)
resultsSK <- cbind(coeffsSK,prop_sigmaSK,zscoreSK,pvalueSK)
resultstS <- cbind(coeffstS,prop_sigmatS,zscoretS,pvaluetS)
#P-valor do modelo de selecao de Heckman padrao
pvalor_Selection <- coef(summary(fit))[,4]
#Vetor de nome dos parametros
parameter <- c("(Intercepto)","Idade",
               "Fem","Educ","Blhisp",
               "Totcr","ins","Renda","(Intercepto)", "Idade","fem",
               "Educ", "Blhisp","Totcr",
               "ins","$\\sigma$", "$\\rho$")
parameterSK <- c("(Intercepto)","Idade",
                 "Fem","Educ","Blhisp",
                 "Totcr","ins","Renda","(Intercepto)", "Idade","fem",
                 "Educ", "Blhisp","Totcr",
                 "ins","$\\sigma$", "$\\rho$","$\\lambda$")
parametertS <- c("(Intercepto)","Idade",
                 "Fem","Educ","Blhisp",
                 "Totcr","ins","Renda","(Intercepto)", "Idade","fem",
                 "Educ", "Blhisp","Totcr",
                 "ins","$\\sigma$", "$\\rho$","$\\nu$")
#Data.Frame com todos os resultados
dfBS<-data.frame(Parametros=parameter,Heckman=beta,Pvalor=pvalor_Selection, BSmodel=theta_BS$par, z=zscoreBS, Pvalor=pvalueBS,Inferior=lowerBS, Superior=upperBS) 
dfHC<-data.frame(Parametros=parameter,Heckman=beta,Pvalor=pvalor_Selection, HCmodel=theta_HC$par, z=zscoreHC, Pvalor=pvalueHC,Inferior=lowerHC, Superior=upperHC) 
dfSK<-data.frame(Parametros=parameterSK,Heckman=c(beta,NA),Pvalor=c(pvalor_Selection,NA), SKmodel=theta_SK$par, z=zscoreSK, Pvalor=pvalueSK,Inferior=lowerSK, Superior=upperSK) 
dftS<-data.frame(Parametros=parametertS,Heckman=c(beta,NA),Pvalor=c(pvalor_Selection,NA), tSmodel=theta_tS$par, z=zscoretS, Pvalor=pvaluetS,Inferior=lowertS, Superior=uppertS)
#Matriz com todos os resultados
kable(dfBS, digits = 3, align = c('c','c','c','c','c','c','c','c','c'))
kable(dfHC, digits = 3, align = c('c','c','c','c','c','c','c','c','c'))
kable(dfSK, digits = 3, align = c('c','c','c','c','c','c','c','c','c'))
kable(dftS, digits = 3, align = c('c','c','c','c','c','c','c','c','c'))

selectEq  <- dambexp~ age + female + educ + blhisp + totchr + ins + income
outcomeEq <- lnambx ~ age + female + educ + blhisp + totchr + ins
meps.fit <- ssmrob(selectEq, outcomeEq, control = heckitrob.control(tcc = 3.2))
summary(meps.fit)
selectionSSMROB <- summary(meps.fit)$coefficients[1]
SelectionSSMROB <- selectionSSMROB$selection
SelectionSSMROB <- SelectionSSMROB[,4]

outcomeSSMROB <- summary(meps.fit)$coefficients[2]
OutcomeSSMROB <- outcomeSSMROB$outcome
OutcomeSSMROB <- OutcomeSSMROB[,4]


rhoSSMROB <- meps.fit$coefficients[16]/meps.fit$sigma

library(xtable)
parameter <- c("(Intercept)","age$~(x_{2})$",
               "fem$~(x_{3})$","educ$~(x_{4})$","blhisp$~(x_{5})$",
               "totcr$~(x_{6})$","ins$~(x_{7})$","inc$~(x_{8})$","(Intercept)", "age$~(x_{2})$","fem$~(x_{3})$",
               "Educ$~(x_{4})$", "blhisp$~(x_{5})$","totcr$~(x_{6})$",
               "ins$~(x_{7})$","$\\sigma~-~\\phi$", "$\\rho$")
#parameter <- c("","$\\gamma_{1}=$","","","$\\gamma_{2}=$","", "","$\\gamma_{3}=$","","","$\\gamma_{4}=$","","","$\\beta_{1}=$","","","$\\beta_{2}=$","","","$\\beta_{3}=$","","","$\\phi=$","", "","$\\rho=$","")

HCmodel <- c(theta_HC$par[1],theta_HC$par[2],theta_HC$par[3],theta_HC$par[4],
             theta_HC$par[5],theta_HC$par[6],theta_HC$par[7],theta_HC$par[8],
             theta_HC$par[9],theta_HC$par[10],theta_HC$par[11],theta_HC$par[12],
             theta_HC$par[13],theta_HC$par[14],theta_HC$par[15],theta_HC$par[16],theta_HC$par[17])
pvalue1 <- c(pvalueHC[1],pvalueHC[2],pvalueHC[3],pvalueHC[4],
             pvalueHC[5],pvalueHC[6],pvalueHC[7],pvalueHC[8],
             pvalueHC[9],pvalueHC[10],pvalueHC[11],pvalueHC[12],
             pvalueHC[13],pvalueHC[14],pvalueHC[15],pvalueHC[16],pvalueHC[17])

ssmrob1 <- c(meps.fit$coefficients[1],meps.fit$coefficients[2],meps.fit$coefficients[3],meps.fit$coefficients[4],
             meps.fit$coefficients[5],meps.fit$coefficients[6],meps.fit$coefficients[7],meps.fit$coefficients[8],
             meps.fit$coefficients[9],meps.fit$coefficients[10],meps.fit$coefficients[11],meps.fit$coefficients[12],
             meps.fit$coefficients[13],meps.fit$coefficients[14],meps.fit$coefficients[15],meps.fit$sigma,rhoSSMROB)
pvalue0 <- c(SelectionSSMROB[1],SelectionSSMROB[2],SelectionSSMROB[3],SelectionSSMROB[4],
             SelectionSSMROB[5],SelectionSSMROB[6],SelectionSSMROB[7],SelectionSSMROB[8],
             OutcomeSSMROB[1],OutcomeSSMROB[2],OutcomeSSMROB[3],OutcomeSSMROB[4],
             OutcomeSSMROB[5],OutcomeSSMROB[6],OutcomeSSMROB[7],NA,OutcomeSSMROB[8])

BSmodel <- c(theta_BS$par[1],theta_BS$par[2],theta_BS$par[3],theta_BS$par[4],
             theta_BS$par[5],theta_BS$par[6],theta_BS$par[7],theta_BS$par[8],
             theta_BS$par[9],theta_BS$par[10],theta_BS$par[11],theta_BS$par[12],
             theta_BS$par[13],theta_BS$par[14],theta_BS$par[15],theta_BS$par[16],theta_BS$par[17])
pvalue2 <- c(pvalueBS[1],pvalueBS[2],pvalueBS[3],pvalueBS[4],
             pvalueBS[5],pvalueBS[6],pvalueBS[7],pvalueBS[8],
             pvalueBS[9],pvalueBS[10],pvalueBS[11],pvalueBS[12],
             pvalueBS[13],pvalueBS[14],pvalueBS[15],pvalueBS[16],pvalueBS[17])

SKmodel <- c(theta_SK$par[1],theta_SK$par[2],theta_SK$par[3],theta_SK$par[4],
             theta_SK$par[5],theta_SK$par[6],theta_SK$par[7],theta_SK$par[8],
             theta_SK$par[9],theta_SK$par[10],theta_SK$par[11],theta_SK$par[12],
             theta_SK$par[13],theta_SK$par[14],theta_SK$par[15],theta_SK$par[16],theta_SK$par[17])
pvalue3 <- c(pvalueSK[1],pvalueSK[2],pvalueSK[3],pvalueSK[4],
             pvalueSK[5],pvalueSK[6],pvalueSK[7],pvalueSK[8],
             pvalueSK[9],pvalueSK[10],pvalueSK[11],pvalueSK[12],
             pvalueSK[13],pvalueSK[14],pvalueSK[15],pvalueSK[16],pvalueSK[17])

tSmodel <- c(theta_tS$par[1],theta_tS$par[2],theta_tS$par[3],theta_tS$par[4],
             theta_tS$par[5],theta_tS$par[6],theta_tS$par[7],theta_tS$par[8],
             theta_tS$par[9],theta_tS$par[10],theta_tS$par[11],theta_tS$par[12],
             theta_tS$par[13],theta_tS$par[14],theta_tS$par[15],theta_tS$par[16],theta_tS$par[17])
pvalue4 <- c(pvaluetS[1],pvaluetS[2],pvaluetS[3],pvaluetS[4],
             pvaluetS[5],pvaluetS[6],pvaluetS[7],pvaluetS[8],
             pvaluetS[9],pvaluetS[10],pvaluetS[11],pvaluetS[12],
             pvaluetS[13],pvaluetS[14],pvaluetS[15],pvaluetS[16],pvaluetS[17])

dt <- data.frame(Par?metros=parameter,BS=BSmodel,pvalue=pvalue2,HCmodel,pvalue=pvalue1,Hrobust=ssmrob1,pvalue=pvalue0,SK=SKmodel,pvalue=pvalue3, tS=tSmodel,pvalue=pvalue4)
head(dt)
#data.frame(parameter,truevalue,BBSB=mestBS,sd_optim=destBS,eqm_optim=eqmestBS) 
cat("\014")
print(xtable(dt,digits=c(1,3,3,3,3,3,3,3,3,3,3,3), caption = "Estimates of the parameters with their respectives 
             p-value under the classical Heckman, Birnbaum-Saunders (BS), Skew-Normal (SK) and Heckman-t (tS) sample 
             selection models.",align = c(rep("c",12))), 
      caption.placement = "top", include.rownames = FALSE,include.colnames = FALSE,
      type = "latex", 
      sanitize.text.function = function(x) {x},
      add.to.row = list(pos = list(0), 
                        command = c("Parameters&BS& p-value &HC& p-value & Robust& p-value & SK &p-value & tS &p-value\\\\"
                        )))



#Log de verossimilhan?a analisado em theta$par
theta_BS$value
YO=ambexp
loglik_BS(theta_BS$par)
YO=lnambx
loglik_HC(theta_HC$par)
theta_HC$value
loglik_SK(theta_SK$par)
theta_SK$value
loglik_tS(theta_tS$par)
theta_tS$value
#Pacote para geracao de tabelas para o latex
# library(xtable)
# #Tabela em formato latex
# print(xtable(df,digits=c(3,3,3,3,3,3,3,3,3), caption = "Estimativas do Modelo de Heckman padr?o via pacote sampleSelection e os respectivos p-valores juntamente com as
# estimativas do modelo de sele??o BS com os respectivos valores de desvio-padr?o (DP),
#              valor-z, p-valor e limites inferior e superior para o intervalo de confian?a de $95\\%.$",align = rep("c",9)),
#       caption.placement = "top", include.rownames = FALSE,include.colnames = TRUE,
#       type = "latex",
#       sanitize.text.function = function(x) {x})
# ################Sem a covariavel rendimento########################

# rm(list=ls())
# cat("\014")
# library(dplyr) # a useful data frame manipulation package
# library(knitr) # a useful package for displaying tables
# 
# library("sampleSelection")
# library(ssmrob)
# data(MEPS2001)
# attach(MEPS2001)
# # MEPS2001[984,21] <- mean(MEPS2001[,21])
# # attach(MEPS2001)
# 
# fit<-selection( dambexp ~ age+female+educ+blhisp+totchr+ins,
#                 lnambx ~ age+female+educ+blhisp+totchr+ins, MEPS2001)
# beta<- fit$estimate
# #Data frame para uso do modelo glm e do modelo lm
# dt=MEPS2001
# library(aod)
# #Modelo probit
# fit1<-glm(dambexp ~ age+female+educ+blhisp+totchr+ins, family = binomial(link = "probit"),data=dt)
# #Gera??o de valores da nova covari?vel raz?o inversa de Mills
# IMR <- dnorm(fit1$linear.predictors)/pnorm(fit1$linear.predictors)
# #Acrescimo de IMR ao dataframe
# dt <- data.frame(dt,IMR)
# #Modelo lm(regressao multipla)
# fit2 <- lm(formula = lnambx ~ age+female+educ+blhisp+totchr+ins+IMR, data = dt[dt$dambexp==1, ])
# summary(fit2)
# #Gera??o de valores da nova covari?vel delta 
# delta <- (dt$IMR)*(dt$IMR+fit1$fitted.values)
# #Acrescimo de delta ao dataframe
# dt <- data.frame(dt,delta)
# #Quantidade de valores u==1
# q <- sum(dambexp)
# #Calculo da vari?ncia
# Var <- (sum((fit2$residuals)^2)+(((coef(fit2)[8])^2)*sum(delta[dambexp==1])))/q
# #Calculo da correla??o
# cor <- (coef(fit2)[8])/sqrt(Var)
# #Chute inicial para optim do modelo BS
# par <- c(coef(fit1),coef(fit2)[-8],phi=sqrt(Var),cor)
# names(par) <- c("(Intercept)","Age",
#                 "Female","Educ","blhisp",
#                 "totchr","ins","(Intercept)", "age","female",
#                 "educ", "blhisp","totchr",
#                 "ins","sigma", "rho")
# #par<- coef(summary(fit))[,1]
# u <- dambexp
# y <- ambexp
# X <- cbind(age,female,educ,blhisp,totchr,ins)
# W <- cbind(age,female,educ,blhisp,totchr,ins)
# XS <- W
# XO <- X
# YO <- y
# YS <- u
# 
# source("Likgrad.R")
# #source("vero_grad_heckman.R")
# Log.lik(par,X,W,y)
# grad(par,X,W,y)
# 
# #beta<- c(fit$estimate[1:16],sin(fit$estimate[17]))
# #heckman<-optim(beta,loglikH,gradlikH,method = "BFGS",hessian = T,control = list(fnscale=-1)) 
# #heckman$par
# theta <- optim(par,Log.lik,grad,y=y,X=X,W=W,method = "BFGS",hessian = T,control = list(fnscale=-1))
# theta$par
# 
# theta$hessian[1]
# 
# fisher_info <- solve(-theta$hessian)
# prop_sigma <- sqrt(diag(fisher_info))
# upper<-theta$par+ qnorm(0.975)*prop_sigma
# lower<-theta$par+ qnorm(0.025)*prop_sigma
# coeffs <- theta$par
# fisher_info <- solve(-theta$hessian)
# stderr <- sqrt(diag(fisher_info))
# zscore <- coeffs/stderr
# pvalue <- 2*(1 - pnorm(abs(zscore)))
# results <- cbind(coeffs,stderr,zscore,pvalue)
# 
# pvalor_Selection <- coef(summary(fit))[,4]
# parameter <- c("(Intercepto)","Idade",
#                "Fem","Educ","Blhisp",
#                "Totcr","ins","(Intercepto)", "Idade","fem",
#                "Educ", "Blhisp","Totcr",
#                "ins","$\\sigma$", "$\\rho$")
# df<-data.frame(Parametros=parameter,Heckman=beta,Pvalor=pvalor_Selection, BSmodel=theta$par, z=zscore, Pvalor=pvalue,Inferior=lower, Superior=upper) 
# kable(df, digits = 3, align = c('c','c','c','c','c','c','c','c','c'))
# library(xtable)
# print(xtable(df,digits=c(3,3,3,3,3,3,3,3,3), caption = "Estimativas do Modelo de Heckman padr?o via pacote sampleSelection e os respectivos p-valores juntamente com as 
# estimativas do modelo de sele??o BS com os respectivos valores de desvio-padr?o (DP),
#              valor-z, p-valor e limites inferior e superior para o intervalo de confian?a de $95\\%.$",align = rep("c",9)), 
#       caption.placement = "top", include.rownames = FALSE,include.colnames = TRUE,
#       type = "latex", 
#       sanitize.text.function = function(x) {x})



# library(xtable)
# parameter <- c("","(Intercept)","","","age","","",
#                "fem","","","educ","","","blhisp","","",
#                "totcr","","","ins","","","inc","","","(Intercept)","","", "age","","","fem","","",
#                "Educ","","", "blhisp","","","totcr","","",
#                "ins","","","$\\sigma$","","", "$\\rho$","")
# #parameter <- c("","$\\gamma_{1}=$","","","$\\gamma_{2}=$","", "","$\\gamma_{3}=$","","","$\\gamma_{4}=$","","","$\\beta_{1}=$","","","$\\beta_{2}=$","","","$\\beta_{3}=$","","","$\\phi=$","", "","$\\rho=$","")
# 
# HCmodel <- c(NA,theta_HC$par[1],NA,NA,theta_HC$par[2],NA,NA,theta_HC$par[3],NA,NA,theta_HC$par[4],NA,NA,
#              theta_HC$par[5],NA,NA,theta_HC$par[6],NA,NA,theta_HC$par[7],NA,NA,theta_HC$par[8],NA,NA,
#              theta_HC$par[9],NA,NA,theta_HC$par[10],NA,NA,theta_HC$par[11],NA,NA,theta_HC$par[12],NA,NA,
#              theta_HC$par[13],NA,NA,theta_HC$par[14],NA,NA,theta_HC$par[15],NA,NA,theta_HC$par[16],NA,NA,theta_HC$par[17],NA)
# pvalue1 <- c(NA,pvalueHC[1],NA,NA,pvalueHC[2],NA,NA,pvalueHC[3],NA,NA,pvalueHC[4],NA,NA,
#              pvalueHC[5],NA,NA,pvalueHC[6],NA,NA,pvalueHC[7],NA,NA,pvalueHC[8],NA,NA,
#              pvalueHC[9],NA,NA,pvalueHC[10],NA,NA,pvalueHC[11],NA,NA,pvalueHC[12],NA,NA,
#              pvalueHC[13],NA,NA,pvalueHC[14],NA,NA,pvalueHC[15],NA,NA,pvalueHC[16],NA,NA,pvalueHC[17],NA)
# 
# BSmodel <- c(NA,theta_BS$par[1],NA,NA,theta_BS$par[2],NA,NA,theta_BS$par[3],NA,NA,theta_BS$par[4],NA,NA,
#              theta_BS$par[5],NA,NA,theta_BS$par[6],NA,NA,theta_BS$par[7],NA,NA,theta_BS$par[8],NA,NA,
#              theta_BS$par[9],NA,NA,theta_BS$par[10],NA,NA,theta_BS$par[11],NA,NA,theta_BS$par[12],NA,NA,
#              theta_BS$par[13],NA,NA,theta_BS$par[14],NA,NA,theta_BS$par[15],NA,NA,theta_BS$par[16],NA,NA,theta_BS$par[17],NA)
# pvalue2 <- c(NA,pvalueBS[1],NA,NA,pvalueBS[2],NA,NA,pvalueBS[3],NA,NA,pvalueBS[4],NA,NA,
#              pvalueBS[5],NA,NA,pvalueBS[6],NA,NA,pvalueBS[7],NA,NA,pvalueBS[8],NA,NA,
#              pvalueBS[9],NA,NA,pvalueBS[10],NA,NA,pvalueBS[11],NA,NA,pvalueBS[12],NA,NA,
#              pvalueBS[13],NA,NA,pvalueBS[14],NA,NA,pvalueBS[15],NA,NA,pvalueBS[16],NA,NA,pvalueBS[17],NA)
# 
# SKmodel <- c(NA,theta_SK$par[1],NA,NA,theta_SK$par[2],NA,NA,theta_SK$par[3],NA,NA,theta_SK$par[4],NA,NA,
#              theta_SK$par[5],NA,NA,theta_SK$par[6],NA,NA,theta_SK$par[7],NA,NA,theta_SK$par[8],NA,NA,
#              theta_SK$par[9],NA,NA,theta_SK$par[10],NA,NA,theta_SK$par[11],NA,NA,theta_SK$par[12],NA,NA,
#              theta_SK$par[13],NA,NA,theta_SK$par[14],NA,NA,theta_SK$par[15],NA,NA,theta_SK$par[16],NA,NA,theta_SK$par[17],NA)
# pvalue3 <- c(NA,pvalueSK[1],NA,NA,pvalueSK[2],NA,NA,pvalueSK[3],NA,NA,pvalueSK[4],NA,NA,
#              pvalueSK[5],NA,NA,pvalueSK[6],NA,NA,pvalueSK[7],NA,NA,pvalueSK[8],NA,NA,
#              pvalueSK[9],NA,NA,pvalueSK[10],NA,NA,pvalueSK[11],NA,NA,pvalueSK[12],NA,NA,
#              pvalueSK[13],NA,NA,pvalueSK[14],NA,NA,pvalueSK[15],NA,NA,pvalueSK[16],NA,NA,pvalueSK[17],NA)
# 
# tSmodel <- c(NA,theta_tS$par[1],NA,NA,theta_tS$par[2],NA,NA,theta_tS$par[3],NA,NA,theta_tS$par[4],NA,NA,
#              theta_tS$par[5],NA,NA,theta_tS$par[6],NA,NA,theta_tS$par[7],NA,NA,theta_tS$par[8],NA,NA,
#              theta_tS$par[9],NA,NA,theta_tS$par[10],NA,NA,theta_tS$par[11],NA,NA,theta_tS$par[12],NA,NA,
#              theta_tS$par[13],NA,NA,theta_tS$par[14],NA,NA,theta_tS$par[15],NA,NA,theta_tS$par[16],NA,NA,theta_tS$par[17],NA)
# pvalue4 <- c(NA,pvaluetS[1],NA,NA,pvaluetS[2],NA,NA,pvaluetS[3],NA,NA,pvaluetS[4],NA,NA,
#              pvaluetS[5],NA,NA,pvaluetS[6],NA,NA,pvaluetS[7],NA,NA,pvaluetS[8],NA,NA,
#              pvaluetS[9],NA,NA,pvaluetS[10],NA,NA,pvaluetS[11],NA,NA,pvaluetS[12],NA,NA,
#              pvaluetS[13],NA,NA,pvaluetS[14],NA,NA,pvaluetS[15],NA,NA,pvaluetS[16],NA,NA,pvaluetS[17],NA)
# 
# dt <- data.frame(Par?metros=parameter,HC=HCmodel,pvalue=pvalue1,Clear=rep("NA",51),BS=BSmodel,pvalue=pvalue2, Clear=rep("NA",51),SK=SKmodel,pvalue=pvalue3, Clear=rep("NA",51),tS=tSmodel,pvalue=pvalue4,Clear=rep("NA",51))
# head(dt)
# #data.frame(parameter,truevalue,BBSB=mestBS,sd_optim=destBS,eqm_optim=eqmestBS) 
# cat("\014")
# print(xtable(dt,digits=c(1,1,3,3,3,3,3,3,3,3,3,3,3,3), caption = "Estimates of the parameters with their respectives 
#              p-value under the classical Heckman, Birnbaum-Saunders (BS), Skew-Normal (SK) and Heckman-t (tS) sample 
#              selection models.",align = c(rep("c",14))), 
#       caption.placement = "top", include.rownames = FALSE,include.colnames = FALSE,
#       type = "latex", 
#       sanitize.text.function = function(x) {x},add.to.row = list(pos = list(0), command = c(
#         "Parameters&Heckman& p-value &Clear & BS& p-value &Clear & SK &p-value &Clear & tS &p-value&\\\\"
#       )))

