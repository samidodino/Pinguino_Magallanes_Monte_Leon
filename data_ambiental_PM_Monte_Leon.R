############ Environmental data - Monte Leon ##################
### Agrupar rasters Chla, SST y SSTA por temporadas
### Extraer buffer al rededor de la colonia de Monte Leon
### Maxima distancia de alimentacion de Magallanicos adultos en cuidado temprano de pichones:
### De Rosciano et al. 2016: 21 a 35km
### Raya Rey et al. 2010: 14 to 45 km

################### Chla, SST, SSTA ################################
# Datos de Chla y SST fueron bajados de MODIS, de la pagina Oceandata: https://oceandata.sci.gsfc.nasa.gov/api/file_search
# Datos de SSTA bajados del NOAA: https://pae-paha.pacioos.hawaii.edu/erddap/griddap/dhw_5km.html

rm(list = ls())
library(terra)
library(rstatix)
library(dplyr)
library(sf)
library(tools)
library(tidyr)
library(purrr)
library(stringr)
library(readxl)

#### ----- Cuidado de pichones -----###
#### Agrupamiento de rasters por temporada ####

#Directorio con los rasters
base_dir <- "E:/SSTA"                    
output_dir <- "E:/SSTA/season_means_PM"     
dir.create(output_dir, showWarnings = FALSE) 

# Pingüino Magallames - etapas
seasons <- list(
  early_chick = c("NOVIEMBRE","DICIEMBRE","ENERO")# usamos solo esta etapa para este analisis pero se pueden sacar todas
  #pre_molt = c("FEBRERO", "MARZO"),
  #winter = c("ABRIL", "MAYO", "JUNIO", "JULIO", "AGOSTO", "SEPTIEMBRE"),
  #last_1_month_winter=c("SEPTIEMBRE"),
  #last_2months_winter= c("AGOSTO", "SEPTIEMBRE")
)

# Función para obtener los archivos de una temporada 

get_season_files <- function(year, season, season_name) {
  season_files <- c()
  for (month in season) {
    if (season_name == "early_chick") {
      # Noviembre y diciembre vienen del año anterior
      if (month %in% c("NOVIEMBRE", "DICIEMBRE")) {
        folder_year <- year - 1
      } else if (month == "ENERO") {
        folder_year <- year  # Enero del año actual
      } else {
        folder_year <- year
      }
    } else {
      folder_year <- year  # Para otras estaciones
    }
    
    folder <- file.path(base_dir, as.character(folder_year))
    month_files <- list.files(folder, pattern = paste0(month, folder_year, "_SSTA\\.nc$"), full.names = TRUE)
    season_files <- c(season_files, month_files)
  }
  
  return(season_files)
}


# Función para calcular la media de rásters
calculate_mean_raster <- function(raster_files) {
  if (length(raster_files) == 0) {
    return(NULL)  # Si no hay archivos, devolver NULL
  }
  rasters <- rast(raster_files)  
  mean_raster <- mean(rasters)  
  return(mean_raster)
}

# Procesar cada año y temporada

for (year in 2010:2024) {
  for (season_name in names(seasons)) {
    
    # Obtener los archivos correspondientes a la temporada
    season_files <- get_season_files(year, seasons[[season_name]], season_name)
    
    # Media de los rásters de la temporada
    mean_season_raster <- calculate_mean_raster(season_files)
    
    # Guardar ráster medio 
    if (!is.null(mean_season_raster)) {
      
      # Etiquetado del archivo de salida según la temporada
      if (season_name == "early_chick") {
        # Correcto: año de nov/dic es anterior (year - 1), año de enero es el actual (year)
        season_label <- paste0("s", substr(year - 1, 3, 4), "-", substr(year, 3, 4))
      } else {
        # Otros casos, solo el año
        season_label <- as.character(year)
      }
      
      output_file <- file.path(output_dir, paste0(season_label, "_", season_name, "_SSTA.tif"))
      
      writeRaster(mean_season_raster, output_file, overwrite = TRUE)
    }
  }
}


##Extraer datos de SST y Chla haciendo un buffer de 50km desde la colonia
# Definir colonias
sitios <- data.frame(
  lon = c(-68.95792),
  lat = c(-50.26419),
  colony_id = c("MLEON")
)

sitios_vect <- vect(sitios, geom = c("lon", "lat"), crs = "EPSG:4326")
sitios_proj <- project(sitios_vect, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

buffers_50km <- buffer(sitios_proj, width = 50000)
buffers_50km_geo <- project(buffers_50km, "+proj=longlat +datum=WGS84")

# Rutas a las carpetas de rásters
data_dirs <- list(
  PM = list(
    SST = "E:/SST_MODIS/season_means_PM",
    SSTA = "E:/SSTA/season_means_PM",
    Chla = "E:/Chla/season_means_PM"
  )
)



# Listar archivos ráster de ambas carpetas
raster_files <- lapply(data_dirs, function(paths) {
  list(
    SST = list.files(paths$SST, pattern = "\\_early_chick_SST.tif$", full.names = TRUE),
    SSTA = list.files(paths$SSTA, pattern = "\\_early_chick_SSTA.tif$", full.names = TRUE),
    Chla = list.files(paths$Chla, pattern = "\\_early_chick_CHLa4km.tif$", full.names = TRUE)
  )
})

# Función para extraer y promediar valores ráster dentro del buffer
extraer_valores <- function(raster_files, buffers, sitios, variable_name) {
  resultados <- list()
  for (file in raster_files) {
    raster <- rast(file)
    
    if (crs(raster) != crs(buffers)) {
      buffers <- project(buffers, crs(raster))
    }
    valores <- terra::extract(raster, buffers, fun = mean, na.rm = TRUE)
    year <- gsub("-.*", "", strsplit(basename(file), "_")[[1]][1])
    
    df <- data.frame(
      colname3 = sitios$colony_id,
      year = year,
      variable = variable_name,
      value = valores[, 2]
    )
    
    resultados[[length(resultados) + 1]] <- df
  }
  bind_rows(resultados)
}


# Extraer para PM 
datos_SST_PM <- extraer_valores(raster_files$PM$SST, buffers_50km_geo, sitios, "SST")
datos_SSTA_PM <- extraer_valores(raster_files$PM$SSTA, buffers_50km_geo, sitios, "SSTA")
datos_Chla_PM <- extraer_valores(raster_files$PM$Chla, buffers_50km_geo, sitios, "Chla")


# Unir todos los datos
datos_early_chick <- bind_rows(
  datos_SST_PM, datos_SSTA_PM, datos_Chla_PM
)

# Convertir a formato ancho
datos_early_chick <- datos_early_chick %>%
  group_by(colname3, year, variable) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = variable, values_from = value)
head(datos_early_chick)

################### Datos SAM desde la URL ##################
library(dplyr)
library(ggplot2)

url <- "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii"
sam_data <- read.table(url, header = FALSE)

head(sam_data)

colnames(sam_data) <- c("Year", "Month", "SAM_value")
sam_data$date <- as.Date(paste(sam_data$Year, sam_data$Month, "01", sep = "-"), format="%Y-%m-%d")

sam_filtered <- sam_data %>%
  filter(Year >= 2010 & Year <= 2025)
head(sam_filtered)
tail(sam_filtered)

# Plot de la serie temporal de SAM entre 2010 y 2025
ggplot(sam_filtered, aes(x = date, y = SAM_value)) +
  geom_line() +
  labs(title = "Serie temporal de SAM (2010-2025)", x = "Fecha", y = "Valor SAM") +
  theme_minimal()


# Agregar mes a sam_filtered y asignar season
sam_filtered <- sam_filtered %>%
  mutate(
    season = case_when(
      Month %in% c(10, 11, 12) ~ paste0("s", str_sub(Year, 3, 4), "-", str_sub(Year + 1, 3, 4)),
      Month %in% c(1, 2) ~ paste0("s", str_sub(Year - 1, 3, 4), "-", str_sub(Year, 3, 4)),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(season))

# Promedio SAM por temporada
sam_by_season <- sam_filtered %>%
  filter(Month %in% c(10, 11, 12, 1, 2)) %>%
  group_by(season) %>%
  summarise(SAM = mean(SAM_value, na.rm = TRUE))

head(sam_by_season)

### Unir SST, SSTA, Chla y SAM:
# Crear columna 'season' con formato "s10-11" y eliminar 'year'
datos_early_chick <- datos_early_chick %>%
  mutate(season = paste0(year, "-", sprintf("%02d", as.numeric(substr(year, 2, 3)) + 1))) %>%
  select(-year)

# Unir con SAM (ya está en formato correcto)
datos_completos <- datos_early_chick %>%
  left_join(sam_by_season, by = "season") %>%
  select(colname3, season, everything())  # Reordenar columnas


write.csv(datos_completos, "data_ambiental_PM_monte_leon.csv", row.names = FALSE)


