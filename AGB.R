rm(list = ls())
cat("\014")
library(tidyverse)
library(sampleSelection)
library(ssmodels)
library(BIOMASS)
library(httr2)
# Lê o CSV com atenção para strings vazias ou lixo
df <- read.csv("ess/data/ForestGEO.csv",
               na.strings = c("", "NA", "?", "NULL", "null", "–", " ", "-", "N/A"))

# Pre-processamento dos dados
df <- df %>%
  mutate(
    dbh.full.census     = suppressWarnings(as.numeric(as.character(dbh.full.census)))/10,#Transformar DAP em cm
    meanWD              = as.numeric(meanWD),
    hom                 = as.numeric(hom),
    H_considering_damage= as.numeric(H_considering_damage),
    weights.ind         = as.numeric(weights.ind),
    site                = factor(site),
    sel                 = as.integer(status == "A"),
    date.ams            = as.Date(date.ams),
    date.full.census    = as.Date(date.full.census),
    time_since_census   = as.numeric(date.ams - date.full.census),
    b                   = as.numeric(b)/100,
    b                   = ifelse(is.na(b), 0, b),
    b                   = ifelse(sel == 0, 0, b),
    H_considering_damage= ifelse(sel, H_considering_damage, 0),
    logH = ifelse(sel, log(H_considering_damage), 0)
  ) %>%
  filter(
    !is.na(sel),
    !is.na(dbh.full.census),
    !is.na(H_considering_damage),
    !is.na(meanWD),
    !is.na(hom),
    !is.na(site),
    !is.na(time_since_census),
    !is.na(weights.ind),
    is.finite(dbh.full.census),
    is.finite(weights.ind),
    is.finite(time_since_census)
  ) %>%
  droplevels()

names(df) <- c("Local",
              "IDMort",
              "ID",
              "dateCensus",
              "DAP",
              "AltDAP",
              "DENS",
              "dateAms",
              "Status",
              "AltCDano",
              "PropVol",
              "Prop",
              "sel", 
              "time_since_census",
              "logH")

df <-  df %>% 
  filter(!(sel == 1 & AltCDano == 0) & !(sel == 1 & Prop==0))
#head(df)
#df$DAP <- df$DAP/10
# Definir coordenadas dos locais
coord_local <- data.frame(
  Local = c("amacayacu", "yasuni", "pasoh", "bci", "hkk", "kc", "fushan"),
  long = c(-70.2678, -76.397, 102.313, -79.87, 99.217, 99.798, 121.555),
  lat = c(-3.8091, -0.6859, 2.982, 9.177125, 15.6324, 7.5435, 24.7614)
)

# Calcular E para cada local
coord <- cbind(coord_local$long, coord_local$lat)
#E <- computeE(coord)

# Juntar E com os locais
#IndexE <- data.frame(Local = coord_local$Local, E = as.numeric(E))
#saveRDS(IndexE, "IndexE.Rds")
IndexE <- readRDS("IndexE.Rds")
# Juntar o índice E no seu dataframe
df <- merge(df, IndexE, by = "Local")
# Aplicar a fórmula de altura diretamente:
df$H1 <- exp(0.893 - df$E + 0.760 * log(df$DAP) - 0.0340 * (log(df$DAP))^2)

# Aplicar a fórmula de McEwan et al. (2011) apenas para o Local 'fushan'
# Para todas as árvores de Fushan, calcular altura com fórmula local
df$H1[df$Local == "fushan"] <- 36.1 * (1 - exp(-0.013 * df$DAP[df$Local == "fushan"]))
summary(df$H1[df$Local == "fushan"])

# Corrigir altura máxima para cada Local
df$H1[df$Local == "amacayacu" & df$H1 > 50] <- 50
df$H1[df$Local == "yasuni"    & df$H1 > 50] <- 50
df$H1[df$Local == "pasoh"     & df$H1 > 60] <- 60
df$H1[df$Local == "bci"       & df$H1 > 40] <- 40
df$H1[df$Local == "hkk"       & df$H1 > 45] <- 45
df$H1[df$Local == "kc"        & df$H1 > 45] <- 45
df$H1[df$Local == "fushan"    & df$H1 > 23] <- 23
#df$AGB2 <- exp(-2.024-
#                 0.896*df$E+
#                 0.920*log(df$DENS)+
#                 2.795*log(df$DAP)-
#                 0.0461*((log(df$DAP))^2))
# Fórmula de Chave et al. (2014) com DAP, H e densidade
#computeAGB retorna medidas em Mg: 1 megagrama (Mg) = 1.000 kg
df$AGB <- (1000*computeAGB(df$DAP, df$DENS, H=df$AltCDano))*df$Prop #Valor em Toneladas
#df$AGB2 <- (0.0673 * (df$DENS * (df$DAP^2) * df$AltCDano)^0.976)*df$Prop #Valor em Toneladas

head(df)
summary(df$DAP) # DAP em cm
summary(df$DENS) # DENS em g/cm^3
summary(df$AltCDano) # AltCDano em metros
summary(df$AGB)
#summary(df$AGB2)


head(df[df$AGB > 1, c("DAP", "H1", "DENS", "AGB")])
# Histograma da AGB2
hist(df$AGB, breaks = 100, main = "Distribuição da Biomassa (AGB)", 
     xlab = "AGB (Toneladas)", col = "lightblue", xlim = c(0, quantile(df$AGB, 0.99, na.rm = TRUE)))

# Quantos indivíduos tem biomassa acima de 1000 kg (1 tonelada)
cat("Número de árvores com AGB > 1 kg:", sum(df$AGB > 1, na.rm = TRUE), "\n")

# Listar as 10 árvores com maior AGB
top_agb <- df[order(-df$AGB), c("DAP", "H1", "DENS", "AGB")]
head(top_agb, 10)

# Criar logH somente onde H_considering_damage está disponível (ou zero se for censurado)
df <- df %>%
  mutate(
    logAGB               = ifelse(sel == 1 & is.finite(AGB) & AGB > 0, log(AGB), NA),
    AGB                  = ifelse(sel == 1 & is.finite(AGB) & AGB > 0, AGB, NA),
    logDAP               = log(DAP), 
    logDAP_z             = scale(log(DAP)),
    DENS_z               = scale(DENS),
    AltDAP_z             = scale(AltDAP),
    time_since_census_z  = scale(time_since_census),
    b_z                  = scale(PropVol))

#saveRDS(df, "AGB.Rds")

selectEq  <- sel ~ logDAP_z + DENS_z + time_since_census_z
outcomeEq <- logAGB ~ logDAP_z + DENS_z + b_z + AltCDano + Local

modelo_heckit <- selection(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df)
summary(modelo_heckit)

start_manual <- coef(modelo_heckit, part = "full")
start_manual <- unname(start_manual)  # remove nomes para o ssmodels



# Criar logH somente onde H_considering_damage está disponível (ou zero se for censurado)
df2 <- df %>%
  mutate(
    logAGB = ifelse(sel == 1 & is.finite(AGB) & AGB > 0, log(AGB), 0))

modeloCL <- HeckmanCL(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df2,
  start = start_manual)
summary(modeloCL)


start_manual <- c(start_manual, 5) 

modelo_ts <- HeckmantS(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df2,
  start = start_manual)
summary(modelo_ts)

# Criar logH somente onde H_considering_damage está disponível (ou zero se for censurado)
df2 <- df %>%
  mutate(
    logAGB = ifelse(sel == 1 & is.finite(AGB) & AGB > 0, log(AGB), NA),
    AGB2   = ifelse(sel == 1 & is.finite(AGB) & AGB > 0, AGB, 0))

selectEq  <- sel ~ logDAP_z + DENS_z + time_since_census_z
outcomeEq <- AGB2 ~ logDAP_z + DENS_z + b_z + AltCDano + Local

outcomeCL <- logAGB ~ logDAP_z + DENS_z + b_z + AltCDano + Local

modelo_heckit <- selection(
  selection = selectEq,
  outcome   = outcomeCL,
  data      = df2)
summary(modelo_heckit)

start_manual <- coef(modelo_heckit, part = "full")
start_bs <- c(unname(start_manual), 1)  # adiciona phi1 = 1

modelo_bs <- HeckmanBS(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df2,
  start = start_bs)
