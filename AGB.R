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
    DAP                 = suppressWarnings(as.numeric(as.character(dbh.full.census)))/10,
    DENS                = as.numeric(meanWD),
    AltDAP              = as.numeric(hom),
    AltCDano            = as.numeric(H_considering_damage),
    Prop                = as.numeric(weights.ind),
    Local               = factor(site),
    sel                 = as.integer(status == "A"),
    date.ams            = as.Date(date.ams),
    date.full.census    = as.Date(date.full.census),
    time_since_census   = as.numeric(date.ams - date.full.census)
  ) %>%
  filter(
    !is.na(sel),
    !is.na(DAP),
    !is.na(DENS),
    !is.na(AltDAP),
    !is.na(Local),
    !is.na(time_since_census),
    !is.na(Prop),
    is.finite(DAP),
    is.finite(Prop),
    is.finite(time_since_census)
  ) %>%
  droplevels()

#names(df) <- c("Local", 
#               "IDMort", 
#               "ID", 
#               "DataCenso", 
#               "DAP", 
#               "AltDAP", 
#               "DENS", 
#               "DataPesq",
#               "Status",
#               "Alt",
#               "PropVol",
#               "Freq")
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
df$AGB <- (1000*computeAGB(df$DAP, df$DENS, H=df$H1))*df$Prop
df$AGB2 <- (0.0673 * (df$DENS * (df$DAP^2) * df$H1)^0.976)*df$Prop

summary(df$AGB)
summary(df$AGB2)
summary(df$DAP)
hist(df$Prop)

head(df[df$AGB > 10000, c("DAP", "H1", "DENS", "AGB", "AGB2")])
# Histograma da AGB2
hist(df$AGB, breaks = 100, main = "Distribuição da Biomassa (AGB2)", 
     xlab = "AGB (kg)", col = "lightblue", xlim = c(0, quantile(df$AGB2, 0.99, na.rm = TRUE)))

# Resumo estatístico detalhado
summary(df$AGB)

# Estatísticas adicionais
cat("Média:", mean(df$AGB, na.rm = TRUE), "\n")
cat("Mediana:", median(df$AGB, na.rm = TRUE), "\n")
cat("Máximo:", max(df$AGB, na.rm = TRUE), "\n")
cat("Mínimo:", min(df$AGB, na.rm = TRUE), "\n")
cat("Desvio padrão:", sd(df$AGB, na.rm = TRUE), "\n")

# Quantos indivíduos tem biomassa acima de 1000 kg (1 tonelada)
cat("Número de árvores com AGB > 1000 kg:", sum(df$AGB > 1000, na.rm = TRUE), "\n")

# Listar as 10 árvores com maior AGB
top_agb <- df[order(-df$AGB), c("DAP", "H1", "DENS", "AGB")]
head(top_agb, 10)

# Marcar gigantes (opcional)
df$Gigante <- ifelse(df$AGB > 10000, "Sim", "Não")

# Quantos gigantes
table(df$Gigante)

subset(df, weights.ind > 0.1)
head(df)



# Criar logH somente onde H_considering_damage está disponível (ou zero se for censurado)
df <- df %>%
  mutate(
    AltCDano = ifelse(sel, AltCDano, 0),
    logH = log1p(AltCDano),
    AGB = ifelse(sel, AGB, 0),
    logAGB = log1p(AGB),
    b = as.numeric(b),
    b = ifelse(is.na(b) & sel == 0, 0, b)
  ) %>%
  filter(!is.na(b), is.finite(b))

df <- df %>% mutate(b = b / 100,
                    log_DAP              = scale(log(DAP)),
                    DENS                 = scale(DENS),
                    AltDAP                = scale(AltDAP),
                    time_since_census_z  = scale(time_since_census),
                    b_z                  = scale(b)
)


selectEq  <- sel ~ log_DAP + DENS + time_since_census_z
outcomeEq <- logAGB ~ log_DAP + b_z + Local
modelo_heckit <- selection(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df)
summary(modelo_heckit)

start_manual <- coef(modelo_heckit, part = "full")
start_manual <- unname(start_manual)  # remove nomes para o ssmodels


start_manual <- c(start_manual, 5) 

modelo_ts <- HeckmantS(
  selection = selectEq,
  outcome   = outcomeEq,
  data      = df,
  df        = 5,
  start     = start_manual)
summary(modelo_ts)
