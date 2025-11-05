##### Penguins Monte Leon ####
rm(list = ls())
library(tidyverse)
library(dplyr)
library(ggplot2)

##############----- Penguins- Monte Leon -----##############
##Éxito reproductivo + data ambiental

data_PM <- read.csv("ER 2010-2024.csv")


data_amb <- read_csv("data_ambiental_PM_monte_leon.csv")
str(data_amb)

data_full <- data_PM %>% 
  inner_join(data_amb, by = c("season"))

head(data_full)

#write.csv(data_full, "PM_Monte_Leon_full.csv")

################## Análisis Multinomial  ####################
library(glmmTMB)
library(lme4)
library(MASS)
library(effects)

# Nueva Tabla con variables ambientales
data_PM_full <- read.csv("PM_Monte_Leon_full.csv")
str(data_PM_full)

# Columna pichon, analisis multinomial numero de pichones por nido#
# Ver que pasa por season primero:

# Convertir variable pichon a un factor ordenado
data_PM_full$pichon <- factor(data_PM_full$pichon, ordered = TRUE)

# Pasar season a factor (lo estaba tomando como chr por default)
data_PM_full$season <- as.factor(data_PM_full$season)

# Modelo ordinal para season
modelo_ordinal <- polr(pichon ~ season, data = data_PM_full, Hess = TRUE)
summary(modelo_ordinal)

# Obtener p-valores 
ctable <- coef(summary(modelo_ordinal))
p_values <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
cbind(ctable, "p value" = p_values)

# Para visualización de probabilidades de pichon 0,1, y 2 paquete effects

plot(allEffects(modelo_ordinal)) 

table(data_PM_full$pichon,data_PM_full$season)
table(data_PM_full$season)

#### Éxito + variables ambientales Con no linealidad 
library(ordinal)
library(splines)
library(ggeffects)
library(MuMIn)
library(sjPlot)

#pichon como factor ordenado
data_PM_full$pichon <- factor(data_PM_full$pichon, ordered = TRUE)
str(data_PM_full$pichon)


# Modelo ordinal con efecto no lineal 
# 1) Analisis de correlaciones entre variables ambientales

cor.test(data_full$SST,data_full$SSTA)#r=0.9, correlacionadas
cor.test(data_full$SST,data_full$Chla)# r=0.5
cor.test(data_full$SSTA,data_full$Chla)#0.7, correlacionadas
cor.test(data_full$SAM, data_full$SST)# 0.34
cor.test(data_full$SAM, data_full$SSTA)#0.27
cor.test(data_full$SAM, data_full$Chla)#0.15

#2) Modelo
# Con SST + Chla + SAM (excluimos SSTA por colinealidad con SST y Chla)
modelo_no_lineal <- clm(pichon ~ ns(SST, df = 2)+ns(Chla, df = 2)+ns(SAM, df = 2), 
                        data = data_PM_full,na.action = na.fail)
summary(modelo_no_lineal) #SST no significativa!--> decidimos sacarla del modelo

dredge(modelo_no_lineal)
tab_model(modelo_no_lineal)


#Sacando SST
modelo_no_lineal_2 <- clm(pichon ~ ns(Chla, df = 2)+ns(SAM, df = 2), 
                          data = data_PM_full,na.action = na.fail)


summary(modelo_no_lineal_2) 
dredge(modelo_no_lineal_2)
tab_model(modelo_no_lineal_2)

#Con SSTA (excluyendo SST y Chla por colinealidad) + SAM
modelo_no_lineal_3 <- clm(pichon ~ ns(SSTA, df = 1)+ns(SAM, df = 2), 
                          data = data_PM_full,na.action = na.fail)


summary(modelo_no_lineal_3) 
dredge(modelo_no_lineal_3)
tab_model(modelo_no_lineal_3)




# GRAFICAR Predicciones
plot(ggpredict(modelo_no_lineal_3))

# Renombrando paneles
# Obtener predicciones respetando los splines
pred_SSTA <- ggpredict(modelo_no_lineal_3, terms = "SSTA [all]")
pred_sam  <- ggpredict(modelo_no_lineal_3, terms = "SAM [all]")

# Modificar los niveles del factor que controla las facetas/categorias
pred_SSTA$response.level <- factor(pred_SSTA$response.level,
                                   levels = c("1", "2", "3"),
                                   labels = c("0 pichón", "1 pichón", "2 pichones"))
pred_sam$response.level  <- factor(pred_sam$response.level,
                                   levels = c("1", "2", "3"),
                                   labels = c("0 pichón", "1 pichón", "2 pichones"))

# Graficar 
plot(pred_SSTA) + ggtitle("Efecto de las Anomalias de Temperatura (SSTA)")
plot(pred_sam)  + ggtitle("Efecto del índice SAM ")
