rm(list = ls())
cat("\014")
library(tidyverse)
library(sampleSelection)
library(ssmodels)
library(BIOMASS)
# Lê o CSV com atenção para strings vazias ou lixo
df <- read.csv("ess/data/ForestGEO.csv",
               na.strings = c("", "NA", "?", "NULL", "null", "–", " ", "-", "N/A"))
names(df) <- c("Local", "IDMort", "ID", "DataCenso", "DAP", "Alt", "Dens", "DataPesquisaAnual", "Status", "AltArv", "ComplVol", "Freq")
length(table(df$IDMort))
length(table(df$ID))
length(table(df$Local))
df %>%
  group_by(Local) %>%
  summarise(numero_IDs_distintos = n_distinct(IDMort))

df %>%
  group_by(Local) %>%
  summarise(numero_IDs_distintos = n_distinct(IDMort)) %>%
  pull(numero_IDs_distintos) %>%
  sum(na.rm = TRUE)

teste <- df %>%
  group_by(DataPesquisaAnual) %>%
  summarise(numero_IDs_distintos = n_distinct(ID))
# Define as coordenadas
latitude <- 24.7614
longitude <- 121.555
df$agb1 <- computeAGB(D = df$DAP, WD = df$Dens, H = df$Alt)
df$AltSD <- retrieveH(df$DAP, coord = c(longitude, latitude))

# Exemplo: definir latitude e longitude do local
coord <- c(24.7614, 121.555)  # (latitude, longitude)

# Estimar a altura usando retrieveH
altura_estimadas <- retrieveH(D = df$DAP, coord = coord)

# Ver as alturas estimadas
head(altura_estimadas)

library(BIOMASS)
library(httr2)

# Definindo coordenadas
coord <- c(-3.8091,-70.2678)

# Estimando alturas
altura_estimadas <- retrieveH(D = df$DAP, coord = coord)

# Adicionando ao seu banco de dados
df$AlturaEstimadas <- altura_estimadas
library(ggplot2)

ggplot(df, aes(x = DAP, y = AlturaEstimadas)) +
  geom_point(alpha = 0.6) +         # Pontinhos
  geom_smooth(method = "loess", se = FALSE, color = "blue") +  # Linha suavizada
  labs(
    title = "Relação entre DAP e Altura Estimada",
    x = "DAP (cm)",
    y = "Altura Estimada (m)"
  ) +
  theme_minimal()
