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
selection = sel ~ logDAP_z + DENS_z + time_since_census_z
outcome   = AGB ~ logDAP_z + DENS_z + AltCDano+Local
source("HeckmanBS_mod.R")
source("summary.HeckmanBS_mod.R")
theta_BS <- HeckmanBS_mod(
  selection = sel ~ logDAP_z + DENS_z + time_since_census_z,
  outcome   = AGB ~ logDAP_z + DENS_z + AltCDano+Local,
  data      = data)
summary.HeckmanBS_mod(theta_BS)
