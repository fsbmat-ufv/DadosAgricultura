#Limpa todos os objetos na memoria
rm(list=ls())
#Limpa o console
cat("\014")
library(tidyverse)
df <- readRDS("AGB.Rds")
attach(df)
data <- df %>%
  mutate(logDAP          = log(DAP), 
         logDAP_z             = scale(log(DAP)),
         DENS_z               = scale(DENS),
         AltDAP_z             = scale(AltDAP),
         time_since_census_z  = scale(time_since_census),
         b_z                  = scale(PropVol))
source("MelhoriasHeckmanBS.R")
source("MelhoriasSummary.R")
theta_BS <- HeckmanBS_mod(
  selection = sel ~ logDAP + DENS + time_since_census,
  outcome   = AGB ~ logDAP + DENS + AltCDano + Local,
  data      = data)
summary.HeckmanBS_mod(theta_BS)

theta_HC <- sampleSelection::selection(sel ~ logDAP + DENS_z + time_since_census,
                 logAGB ~ logDAP + DENS + AltCDano + Local,
                 data)
summary(theta_HC)

