rm(list = ls())
cat("\014")
library(tidyverse)
library(sampleSelection)
library(ssmodels)
# Lê o CSV com atenção para strings vazias ou lixo
df <- read.csv("ess/data/ForestGEO.csv",
               na.strings = c("", "NA", "?", "NULL", "null", "–", " ", "-", "N/A"))

# Pre-processamento dos dados
df <- df %>%
  mutate(
    dbh.full.census.num = suppressWarnings(as.numeric(as.character(dbh.full.census))),
    meanWD              = as.numeric(meanWD),
    hom                 = as.numeric(hom),
    H_considering_damage = as.numeric(H_considering_damage),
    weights.ind         = as.numeric(weights.ind),
    site                = factor(site),
    sel                     = as.integer(status == "A"),
    log_dbh                 = log(dbh.full.census.num),
    date.ams                = as.Date(date.ams),
    date.full.census        = as.Date(date.full.census),
    time_since_census       = as.numeric(date.ams - date.full.census)
  ) %>%
  filter(
    !is.na(sel),
    !is.na(log_dbh),
    !is.na(meanWD),
    !is.na(hom),
    !is.na(site),
    !is.na(time_since_census),
    !is.na(weights.ind),
    is.finite(log_dbh),
    is.finite(meanWD),
    is.finite(time_since_census)
  ) %>%
  droplevels()

boxplot(df$H_considering_damage)
table(df$H_considering_damage)

df %>%
  mutate(ano = year(date.ams)) %>%
  group_by(ano) %>%
  summarise(
    total_arvores = n(),
    mortas = sum(sel == 0),
    percentual_mortalidade = 100 * mortas / total_arvores
  ) %>%
  pull(percentual_mortalidade) %>%
  mean(na.rm = TRUE)


# Criar logH somente onde H_considering_damage está disponível (ou zero se for censurado)
df <- df %>%
  mutate(
    H_considering_damage = ifelse(sel, H_considering_damage, 0),
    logH = log1p(H_considering_damage),
    b = as.numeric(b),
    b = ifelse(is.na(b) & sel == 0, 0, b)
  ) %>%
  filter(!is.na(b), is.finite(b))

df <- df %>% mutate(b = b / 100,
                    log_dbh_z            = scale(log_dbh),
                    meanWD_z             = scale(meanWD),
                    hom_z                = scale(hom),
                    time_since_census_z  = scale(time_since_census),
                    b_z                  = scale(b)
)

# Verificar a censura

table(df$sel)
# Espera-se algo como:
# FALSE   TRUE 
#  5000  140000

selectEq  <- sel ~ log_dbh_z + meanWD_z + time_since_census_z
outcomeEq <- logH ~ log_dbh_z + b_z + site

modelo_heckit <- selection(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df)
summary(modelo_heckit)

modelo_cl <- HeckmanCL(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df)
summary(modelo_cl)

#start_ts <- c(start_manual, 5)  # df = 5 (graus de liberdade)
modelo_ts <- HeckmantS(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df,
  df        = 5)
summary(modelo_ts)

#start_sk <- c(start_manual, 0.5)  # lambda = 0.5 (assimetria moderada)
modelo_sk <- HeckmanSK(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df,
  lambda    = 0.5)
summary(modelo_sk)

start_manual <- coef(modelo_heckit, part = "full")
start_manual <- unname(start_manual)  # remove nomes para o ssmodels


start_bs <- c(start_manual, 1)  # phi1 = 1
start_bs <- pmax(start_bs, -0.99)  # evita valores fora do domínio
start_bs <- pmin(start_bs, 20)

df_bs <- df %>%
  mutate(
    sel = as.integer(sel),
    H_considering_damage = ifelse(sel == 1, H_considering_damage, NA),  # Censuradas ficam como NA
    log_dbh_z = scale(log_dbh),
    meanWD_z  = scale(meanWD),
    b_z       = scale(b)
  ) %>%
  filter(
    is.finite(H_considering_damage) | sel == 0,  # Mantém as censuradas (com NA em Y)
    is.finite(log_dbh_z),
    is.finite(meanWD_z),
    is.finite(b_z)
  )

table(df_bs$sel)

selectEq  <- sel ~ log_dbh_z
outcomeEq <- H_considering_damage ~ log_dbh_z

start_bs <- c(0.1, 0.1,    # seleção
              1.0,  0.1,   # resultado
              0.8,         # phi1
              0.5)         # rho


modelo_bs <- HeckmanBS(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df_bs,
  start     = start_bs
)


table(df$site)
teste1 <- df %>% filter(site == "bci")
selectEq  <- sel ~ log_dbh_z + meanWD_z + time_since_census_z
outcomeEq <- logH ~ log_dbh_z + b_z 
modelo_heckit <- selection(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = teste1)
summary(modelo_heckit)
table(teste1$sel)
