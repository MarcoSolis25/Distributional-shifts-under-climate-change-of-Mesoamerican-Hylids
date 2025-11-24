#### =====================================
####   CORRELATIVE NICHE MODEL CONSTRUCTION
####   Procesamiento hasta recorte con M.asc
#### =====================================

### CARGAR LIBRERÍAS NECESARIAS ###
library(raster)
library(sp)
library(sf)
library(dismo)
library(biomod2)
library(usdm)
library(ENMeval)
library(kuenm)
library(dplyr)
library(corrplot)
library(terra)

### 1. CONFIGURAR DIRECTORIOS Y PARÁMETROS ###
especie <- "Boana_rosenbergi"  # Cambiar por la especie en cuestión

dir_raiz <- "/Users/hibraim/Desktop/Frontiers"
dir_presencias <- file.path(dir_raiz, "ModelConstruction", "Presencias")
dir_m_variables <- file.path(dir_raiz, "ModelConstruction", "VariablesAmbientales", "M_variables")
dir_g_variables <- file.path(dir_raiz, "ModelConstruction", "VariablesAmbientales", "G_variables")
dir_capas <- file.path(dir_raiz, "ModelConstruction", "Ms_shp")


### 2. CREAR ESTRUCTURA DE DIRECTORIOS ###
dir_modelo <- file.path(dir_raiz, "ModelConstruction", "Modelos", especie)
dir.create(dir_modelo, recursive = TRUE, showWarnings = FALSE)

dir_M_out <- file.path(dir_modelo, "M_variables")
dir.create(dir_M_out, recursive = TRUE, showWarnings = FALSE)

dir_G_out <- file.path(dir_modelo, "G_variables")
dir.create(dir_G_out, recursive = TRUE, showWarnings = FALSE)

message("Directorios creados en: ", dir_modelo)

### 3. LEER DATOS DE PRESENCIA ###
archivo_presencias <- file.path(dir_presencias, paste0(especie, ".csv"))
if (!file.exists(archivo_presencias)) {
  stop(paste("Archivo de presencias no encontrado:", archivo_presencias))
}

data <- read.csv(archivo_presencias)

# Verificar estructura esperada
if (!all(c("species", "lon", "lat") %in% colnames(data))) {
  stop("El archivo de presencias debe contener columnas: species, lon, lat")
}

# Crear objeto espacial y eliminar NA
datos <- data.frame(
  species = data$species,
  lon = as.numeric(data$lon),
  lat = as.numeric(data$lat)
) %>% na.omit()

datos_sp <- SpatialPointsDataFrame(datos[,2:3], datos)

# Eliminar duplicados en misma celda
archivos_capas <- list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE)
if (length(archivos_capas) == 0) stop("No se encontraron archivos raster en M_variables")

capa_referencia <- raster(archivos_capas[1])
celdas <- cellFromXY(capa_referencia, datos_sp)
indices_unicos <- !duplicated(celdas)
datos_sp <- datos_sp[indices_unicos, ]
datos <- datos[indices_unicos, ]

message("Registros únicos por celda: ", nrow(datos))
if (sum(!indices_unicos) > 0) message("Se eliminaron ", sum(!indices_unicos), " puntos duplicados")

### 4. ANÁLISIS DE VARIABLES AMBIENTALES ###
# Cargar variables ambientales
capas_presente <- stack(list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE))
message("Variables cargadas: ", paste(names(capas_presente), collapse = ", "))

# Extraer valores ambientales
presencias_clima <- extract(capas_presente, datos_sp) %>%
  data.frame() %>%
  cbind(datos, .) %>%
  na.omit()

# Guardar coordenadas para MaxEnt
write.csv(presencias_clima[,1:3],
          file.path(dir_modelo, paste0(especie, "_maxent.csv")),
          row.names = FALSE)

# Matriz de correlación
cormatriz <- cor(presencias_clima[,4:ncol(presencias_clima)], use = "complete.obs")
png(file.path(dir_modelo, paste0(especie, "_correlacion_var.png")),
    width = 8, height = 8, units = "in", res = 300)
corrplot(cormatriz, method = "circle", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = paste("Correlación para", especie))
dev.off()

# Análisis VIF
vif_result <- vifstep(presencias_clima[,4:ncol(presencias_clima)], th = 10)
message("\nVariables seleccionadas después de análisis VIF:")
print(vif_result)

write.csv(vif_result@results, file.path(dir_modelo, "vif_results.csv"))
vars_seleccionadas <- vif_result@results$Variables

message("Variables seleccionadas: ", paste(vars_seleccionadas, collapse = ", "))

# Filtrar solo variables seleccionadas
capas_presente_filtradas <- subset(capas_presente, vars_seleccionadas)

### DELIMITACIÓN DE M (ÁREA ACCESIBLE) ###
# 1. Cargar el shapefile correspondiente a la especie ----
archivo_shp <- file.path(dir_capas, paste0(especie, "_M.shp"))

if (!file.exists(archivo_shp)) {
  stop(paste("No se encontró el shapefile del área accesible para la especie:", archivo_shp))
}

M_sf <- st_read(archivo_shp)

# 2. Validar geometría ----
M_sf <- st_make_valid(M_sf)

# 3. Reproyectar si es necesario (al CRS de las capas ambientales) ----
crs_ref <- st_crs(capa_referencia)
if (st_crs(M_sf) != crs_ref) {
  M_sf <- st_transform(M_sf, crs_ref)
  message("El shapefile de la especie fue reproyectado al CRS de las capas ambientales.")
}

# 4. Guardar el área accesible final ----
st_write(M_sf, file.path(dir_modelo, paste0(especie, "_M.shp")), delete_dsn = TRUE)

# 5. Visualización rápida ----
plot(st_geometry(M_sf), main = paste("Área accesible (M) para", especie))
plot(datos_sp, col = "red", pch = 19, cex = 0.5, add = TRUE)


### ==============================================
###   PROCESAMIENTO DE CAPAS (M y G alineadas) - CORREGIDO
### ==============================================

library(terra)
library(dplyr)

# ==========================
# 1. CARGAR REFERENCIA BASE
# ==========================
ref_path <- "/Users/hibraim/Desktop/Frontiers/annual hours_pools.tif"
if (!file.exists(ref_path)) stop("No se encontró la capa de referencia.")

ref_raster <- rast(ref_path)
message("✅ Capa de referencia cargada: ", ref_path)

# ================================
# 2. CARGAR POLÍGONO DE ÁREA M
# ================================
archivo_M <- file.path(dir_modelo, paste0(especie, "_M.shp"))
if (!file.exists(archivo_M)) stop(paste("No se encontró el shapefile:", archivo_M))

poligono <- vect(archivo_M)
poligono <- project(poligono, crs(ref_raster))
message("✅ Polígono cargado y reproyectado al CRS de referencia.")

# =====================================
# 3. VARIABLES SELECCIONADAS (VIF) - CON MAPEO CORREGIDO
# =====================================
archivo_vif <- file.path(dir_modelo, "vif_results.csv")
if (!file.exists(archivo_vif)) stop("No se encontró el archivo de resultados VIF.")

vif_table <- read.csv(archivo_vif)
vars_vif <- vif_table$Variables %>% unique() %>% na.omit()
message("✅ Variables seleccionadas según VIF: ", paste(vars_vif, collapse = ", "))

# MAPEO CORREGIDO DE NOMBRES - VIF a G
mapeo_vif_a_G <- function(nombres_vif) {
  nombres_G <- c()
  
  for (nombre in nombres_vif) {
    # Caso: bio_1, bio_2, etc. -> convertir a bio01, bio02
    if (grepl("bio_(\\d+)", nombre)) {
      numero <- gsub("bio_(\\d+)", "\\1", nombre)
      # Asegurar 2 dígitos
      numero_formateado <- sprintf("%02d", as.numeric(numero))
      nombre_G <- paste0("bio", numero_formateado)
    }
    # Caso: bio1, bio2, etc. -> convertir a bio01, bio02
    else if (grepl("bio(\\d+)", nombre)) {
      numero <- gsub("bio(\\d+)", "\\1", nombre)
      numero_formateado <- sprintf("%02d", as.numeric(numero))
      nombre_G <- paste0("bio", numero_formateado)
    }
    # Caso: ya está en formato bio01, bio02
    else if (grepl("bio\\d{2}", nombre)) {
      nombre_G <- nombre
    }
    # Caso: otros nombres, dejarlos igual
    else {
      nombre_G <- nombre
    }
    nombres_G <- c(nombres_G, nombre_G)
  }
  
  return(nombres_G)
}

# Aplicar mapeo
vars_vif_G <- mapeo_vif_a_G(vars_vif)
message("✅ Variables VIF mapeadas para G: ", paste(vars_vif_G, collapse = ", "))

# También crear mapeo inverso para reportes
mapeo_G_a_vif <- function(nombres_G) {
  nombres_vif <- c()
  
  for (nombre in nombres_G) {
    # Convertir bio01, bio02 a bio_1, bio_2
    if (grepl("bio(\\d{2})", nombre)) {
      numero <- as.numeric(gsub("bio(\\d{2})", "\\1", nombre))
      nombre_vif <- paste0("bio_", numero)
    } else {
      nombre_vif <- nombre
    }
    nombres_vif <- c(nombres_vif, nombre_vif)
  }
  
  return(nombres_vif)
}

# =====================================
# 4. FUNCIÓN MEJORADA PARA PROCESAR CAPAS
# =====================================
procesar_capas_selectivamente <- function(rutas_tif, poligono, ref_raster, vars_seleccionadas, dir_salida, tipo = "M") {
  dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
  capas_procesadas <- c()
  
  for (ruta in rutas_tif) {
    message("\n🔍 Analizando archivo: ", basename(ruta))
    
    # Leer el archivo multicapa
    capas <- rast(ruta)
    nombres_originales <- names(capas)
    
    message("   Variables disponibles: ", paste(nombres_originales, collapse = ", "))
    message("   Buscando variables: ", paste(vars_seleccionadas, collapse = ", "))
    
    # Encontrar variables que coinciden
    vars_coincidentes <- intersect(nombres_originales, vars_seleccionadas)
    
    if (length(vars_coincidentes) == 0) {
      message("   ⚠️ No se encontraron variables coincidentes")
      next
    }
    
    message("   ✅ Variables encontradas: ", paste(vars_coincidentes, collapse = ", "))
    
    # Filtrar solo las variables seleccionadas
    capas_filtradas <- subset(capas, vars_coincidentes)
    
    # Alinear con raster de referencia (resolution & CRS only - NO extent)
    message("   🔄 Alineando resolución y CRS con referencia...")
    
    # Create a target template with ref_raster properties but original extent
    target_template <- rast()
    crs(target_template) <- crs(ref_raster)  # Same CRS
    res(target_template) <- res(ref_raster)  # Same resolution
    ext(target_template) <- ext(capas_filtradas)  # BUT original extent
    
    # Now resample to this custom template
    capas_alineadas <- terra::resample(capas_filtradas, target_template, method = "bilinear")
    
    # Recortar y enmascarar con el polígono
    message("   ✂️ Recortando y enmascarando...")
    capas_recortadas <- crop(capas_alineadas, poligono, snap = "out")
    capas_enmascaradas <- mask(capas_recortadas, poligono)
    
    # Guardar capas individuales
    for (i in 1:nlyr(capas_enmascaradas)) {
      nombre <- names(capas_enmascaradas)[i]
      archivo_salida <- file.path(dir_salida, paste0(nombre, ".asc"))
      
      writeRaster(
        capas_enmascaradas[[i]],
        filename = archivo_salida,
        filetype = "AAIGrid",
        NAflag = -9999,
        overwrite = TRUE
      )
      capas_procesadas <- c(capas_procesadas, archivo_salida)
      message("   💾 Guardado: ", nombre, ".asc")
    }
  }
  
  message("\n✅ Procesamiento completado para tipo: ", tipo)
  message("   Total capas guardadas: ", length(capas_procesadas))
  return(capas_procesadas)
}
# Listar archivos M
rutas_M <- list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE)
message("Archivos M encontrados: ", length(rutas_M))

if (length(rutas_M) == 0) {
  stop("❌ No se encontraron archivos .tif en el directorio M_variables")
}

# Verificar contenido de los archivos M
for (ruta in rutas_M) {
  capas <- rast(ruta)
  message("\n📁 Archivo M: ", basename(ruta))
  message("   Capas disponibles: ", paste(names(capas), collapse = ", "))
  coincidentes_M <- intersect(names(capas), vars_vif)
  message("   Coinciden con VIF: ", ifelse(length(coincidentes_M) > 0, paste(coincidentes_M, collapse = ", "), "NINGUNA"))
}

# Procesar variables M
dir_M_out <- file.path(dir_modelo, "M_variables_recortadas/Set_1")
capas_M_procesadas <- procesar_capas_selectivamente(rutas_M, poligono, ref_raster, vars_vif, dir_M_out, "M")

# =====================================
# 6. VERIFICAR Y PROCESAR VARIABLES G
# =====================================
subdirs_G <- list.dirs(dir_g_variables, recursive = FALSE, full.names = TRUE)
message("Escenarios G encontrados: ", length(subdirs_G))

if (length(subdirs_G) == 0) {
  stop("❌ No se encontraron subdirectorios en G_variables")
}

# Mostrar el mapeo completo para verificación
message("\n🔍 MAPEO COMPLETO DE VARIABLES:")
for (i in 1:length(vars_vif)) {
  message("   VIF: ", vars_vif[i], " → G: ", vars_vif_G[i])
}

dir_G_out <- file.path(dir_modelo, "G_variables_recortadas/Set_1")
total_capas_G <- 0

for (sub in subdirs_G) {
  nombre_sub <- basename(sub)
  message("\n📂 Procesando escenario G: ", nombre_sub)
  
  rutas_G <- list.files(sub, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(rutas_G) == 0) {
    message("   ⚠️ No se encontraron archivos .tif")
    next
  }
  
  dir_G_sub_out <- file.path(dir_G_out, nombre_sub)
  capas_G_procesadas <- procesar_capas_selectivamente(rutas_G, poligono, ref_raster, vars_vif_G, dir_G_sub_out, "G")
  total_capas_G <- total_capas_G + length(capas_G_procesadas)
}

# =====================================
# 7. VERIFICACIÓN FINAL Y REPORTE
# =====================================
# Reporte M
message("\n📊 VARIABLES M PROCESADAS:")
if (length(capas_M_procesadas) > 0) {
  message("   ✅ Total: ", length(capas_M_procesadas), " capas")
  for (capa in capas_M_procesadas) {
    message("      ✓ ", basename(capa))
  }
} else {
  message("   ❌ No se procesaron capas M")
}

# Reporte G
message("\n📊 VARIABLES G PROCESADAS:")
if (total_capas_G > 0) {
  message("   ✅ Total: ", total_capas_G, " capas")
  escenarios_G <- list.dirs(dir_G_out, recursive = FALSE, full.names = FALSE)
  for (escenario in escenarios_G) {
    capas_escenario <- list.files(file.path(dir_G_out, escenario), pattern = "\\.asc$", full.names = TRUE)
    if (length(capas_escenario) > 0) {
      message("   📂 ", escenario, " (", length(capas_escenario), " capas):")
      for (capa in capas_escenario) {
        message("      ✓ ", basename(capa))
      }
    }
  }
} else {
  message("   ❌ No se procesaron capas G")
}

# Verificar cobertura de variables VIF
message("\n🔍 VERIFICACIÓN DE COBERTURA VIF:")
message("   Variables VIF originales: ", paste(vars_vif, collapse = ", "))

todas_capas <- c(capas_M_procesadas, 
                 list.files(dir_G_out, pattern = "\\.asc$", recursive = TRUE, full.names = TRUE))

nombres_procesados <- tools::file_path_sans_ext(basename(todas_capas))

# Para verificación, convertir nombres G procesados a formato VIF
nombres_procesados_vif <- mapeo_G_a_vif(nombres_procesados)

# Verificar variables M y G juntas
vars_faltantes <- setdiff(vars_vif, nombres_procesados_vif)
if (length(vars_faltantes) == 0) {
  message("   ✅ ¡TODAS las variables VIF fueron procesadas correctamente!")
} else {
  message("   ⚠️ Variables VIF faltantes: ", paste(vars_faltantes, collapse = ", "))
}

message("\n🎯 PROCESAMIENTO COMPLETADO PARA: ", especie)

# =====================================
# 8. RENOMBRAR ARCHIVOS G PARA COINCIDIR CON NOMBRES M (TODOS LOS ARCHIVOS)
# =====================================

renombrar_archivos_G <- function(dir_G_base) {
  message("\n🔄 RENOMBRANDO ARCHIVOS G PARA COINCIDIR CON NOMBRES M...")
  
  # Función para convertir nombres bio01 -> bio_1, etc.
  convertir_nombre <- function(nombre_archivo) {
    nombre_sin_ext <- tools::file_path_sans_ext(nombre_archivo)
    extension <- tools::file_ext(nombre_archivo)
    
    # Si es formato bio01, bio02, etc., convertir a bio_1, bio_2
    if (grepl("^bio\\d{2}$", nombre_sin_ext)) {
      numero <- as.numeric(gsub("bio(\\d{2})", "\\1", nombre_sin_ext))
      nuevo_nombre <- paste0("bio_", numero)
      return(paste0(nuevo_nombre, ".", extension))
    }
    # Si ya está en formato bio_1, bio_2, mantener igual
    else if (grepl("^bio_\\d+$", nombre_sin_ext)) {
      return(nombre_archivo)
    }
    # Para otros nombres, mantener igual
    else {
      return(nombre_archivo)
    }
  }
  
  # Buscar todos los subdirectorios en G_variables_recortadas
  subdirs <- list.dirs(dir_G_base, recursive = TRUE, full.names = TRUE)
  total_renombrados <- 0
  
  for (subdir in subdirs) {
    # Buscar TODOS los archivos en este subdirectorio (no solo .asc)
    todos_archivos <- list.files(subdir, full.names = FALSE)
    
    if (length(todos_archivos) > 0) {
      message("📂 Procesando directorio: ", basename(subdir))
      
      # Agrupar archivos por nombre base (sin extensión)
      archivos_por_base <- list()
      
      for (archivo in todos_archivos) {
        nombre_base <- tools::file_path_sans_ext(archivo)
        if (!nombre_base %in% names(archivos_por_base)) {
          archivos_por_base[[nombre_base]] <- character()
        }
        archivos_por_base[[nombre_base]] <- c(archivos_por_base[[nombre_base]], archivo)
      }
      
      # Procesar cada grupo de archivos
      for (nombre_base in names(archivos_por_base)) {
        archivos_grupo <- archivos_por_base[[nombre_base]]
        nuevo_nombre_base <- tools::file_path_sans_ext(convertir_nombre(paste0(nombre_base, ".asc")))
        
        # Si el nombre base cambió, renombrar todos los archivos del grupo
        if (nuevo_nombre_base != nombre_base) {
          for (archivo in archivos_grupo) {
            extension <- tools::file_ext(archivo)
            archivo_viejo <- file.path(subdir, archivo)
            archivo_nuevo <- file.path(subdir, paste0(nuevo_nombre_base, ".", extension))
            
            # Verificar que el archivo nuevo no exista ya
            if (!file.exists(archivo_nuevo)) {
              file.rename(archivo_viejo, archivo_nuevo)
              message("   ✅ Renombrado: ", archivo, " -> ", paste0(nuevo_nombre_base, ".", extension))
              total_renombrados <- total_renombrados + 1
            } else {
              message("   ⚠️ El archivo ya existe, omitiendo: ", paste0(nuevo_nombre_base, ".", extension))
            }
          }
        }
      }
    }
  }
  
  message("🎯 Total de archivos renombrados: ", total_renombrados)
  return(total_renombrados)
}

# Ejecutar el renombrado después del procesamiento de G variables
dir_G_base <- file.path(dir_modelo, "G_variables_recortadas")
if (dir.exists(dir_G_base)) {
  total_renombrados <- renombrar_archivos_G(dir_G_base)
  
  # Verificación final
  message("\n🔍 VERIFICACIÓN FINAL DE NOMBRES:")
  
  # Obtener nombres base de todos los archivos en M (sin extensiones)
  archivos_M <- list.files(dir_M_out, full.names = FALSE)
  nombres_base_M <- unique(tools::file_path_sans_ext(archivos_M))
  
  # Obtener nombres base de todos los archivos en G (sin extensiones)
  archivos_G <- list.files(dir_G_base, recursive = TRUE, full.names = FALSE)
  nombres_base_G <- unique(tools::file_path_sans_ext(archivos_G))
  
  message("Nombres base en M: ", paste(sort(nombres_base_M), collapse = ", "))
  message("Nombres base en G: ", paste(sort(nombres_base_G), collapse = ", "))
  
  if (setequal(nombres_base_M, nombres_base_G)) {
    message("✅ ¡LOS NOMBRES BASE COINCIDEN PERFECTAMENTE!")
  } else {
    message("❌ Aún hay diferencias en los nombres base")
    message("   Faltantes en G: ", paste(setdiff(nombres_base_M, nombres_base_G), collapse = ", "))
    message("   Faltantes en M: ", paste(setdiff(nombres_base_G, nombres_base_M), collapse = ", "))
  }
  
  # Mostrar archivos por variable para verificación
  message("\n📋 ARCHIVOS POR VARIABLE (ejemplo para bio_1):")
  ejemplo_var <- "bio_1"
  archivos_ejemplo_M <- list.files(dir_M_out, pattern = paste0("^", ejemplo_var, "\\."), full.names = FALSE)
  archivos_ejemplo_G <- list.files(dir_G_base, pattern = paste0("^", ejemplo_var, "\\."), recursive = TRUE, full.names = FALSE)
  
  message("Archivos para ", ejemplo_var, " en M: ", paste(archivos_ejemplo_M, collapse = ", "))
  message("Archivos para ", ejemplo_var, " en G: ", paste(unique(archivos_ejemplo_G), collapse = ", "))
  
} else {
  message("❌ No se encontró el directorio G_variables_recortadas")
}

# =====================================
# Limpiar registros sin datos climáticos
# =====================================
# Function to check occurrence points against climatic layers
check_occurrence_points <- function(occ_file, env_layers_dir, output_clean_file = NULL) {
  
  # Read occurrence data
  occ_data <- read.csv(occ_file)
  cat("Original occurrence points:", nrow(occ_data), "\n")
  
  # Get all environmental layers (excluding hidden files)
  env_files <- list.files(env_layers_dir, pattern = "\\.asc$", full.names = TRUE)
  cat("Number of environmental layers found:", length(env_files), "\n")
  
  if (length(env_files) == 0) {
    stop("No environmental layers found in the specified directory")
  }
  
  # Create raster stack
  env_stack <- stack(env_files)
  cat("Environmental layers stack created\n")
  cat("Extent of environmental data:\n")
  print(extent(env_stack))
  cat("CRS of environmental data:", as.character(crs(env_stack)), "\n\n")  # FIXED: Convert CRS to character
  
  # Convert occurrence data to spatial points
  coordinates(occ_data) <- ~lon+lat
  crs(occ_data) <- crs(env_stack)  # Set same CRS as environmental data
  
  # Extract environmental values for all points
  env_values <- extract(env_stack, occ_data)
  
  # Check for NA values (points without climatic data)
  na_per_point <- apply(env_values, 1, function(x) sum(is.na(x)))
  points_with_na <- which(na_per_point > 0)
  points_with_all_data <- which(na_per_point == 0)
  
  cat("=== CHECK RESULTS ===\n")
  cat("Points with complete climatic data:", length(points_with_all_data), "\n")
  cat("Points missing some climatic data:", length(points_with_na), "\n")
  
  if (length(points_with_na) > 0) {
    cat("\nPoints missing climatic data:\n")
    problematic_points <- as.data.frame(occ_data[points_with_na, ])
    print(problematic_points)
    
    # Show which layers have NA for each problematic point
    cat("\nMissing data details:\n")
    for (i in points_with_na) {
      missing_layers <- which(is.na(env_values[i, ]))
      cat("Point", i, "(", occ_data$lon[i], ",", occ_data$lat[i], 
          ") - Missing layers:", length(missing_layers), "/", ncol(env_values), "\n")
    }
  } else {
    cat("All points have complete climatic data! ✓\n")
  }
  
  # Create clean dataset (only points with all climatic data)
  clean_occ <- as.data.frame(occ_data[points_with_all_data, ])
  cat("\nClean dataset has", nrow(clean_occ), "points\n")
  
  # Calculate percentage kept
  if (nrow(occ_data) > 0) {
    percent_kept <- round((nrow(clean_occ) / nrow(occ_data)) * 100, 1)
    cat("Percentage of points kept:", percent_kept, "%\n")
  }
  
  # Save clean data if output file is specified
  if (!is.null(output_clean_file)) {
    write.csv(clean_occ, output_clean_file, row.names = FALSE)
    cat("Clean data saved to:", output_clean_file, "\n")
  }
  
  # Return results
  return(list(
    original_points = nrow(occ_data),
    clean_points = nrow(clean_occ),
    removed_points = length(points_with_na),
    problematic_points = if (length(points_with_na) > 0) as.data.frame(occ_data[points_with_na, ]) else NULL,
    clean_data = clean_occ
  ))
}
# Usage example for your specific case:
# Set your working directory
setwd(file.path(dir_modelo))
datos2 <- paste0(especie, "_maxent.csv")
# Check your occurrence points against M_variables
results <- check_occurrence_points(
  occ_file = datos2, 
  env_layers_dir = "M_variables_recortadas/Set_1",  # Adjust if you have multiple sets
  output_clean_file = "datos_clean.csv"
)

# =====================================
# Selección de puntos (train / test / joint / ind)
# =====================================

dir_kuenm <- file.path(dir_modelo)
datos_clean<- read.csv2("datos_clean.csv", sep=",", header=TRUE)
# Función para dividir datos y guardar archivos con las columnas species, lon, lat
dividir_datos <- function(datos_clean,
                          dir_salida = dir_kuenm,
                          prefijo = especie,
                          seed = 42,
                          prop_test_externo = 0.2,   # proporción para test externo
                          prop_internal = 0.8) {     # proporción dentro de train para joint (el resto será ind)
  set.seed(seed)
  
  # Verificaciones básicas
  if (!all(c("species", "lon", "lat") %in% colnames(datos_clean))) {
    stop("El data.frame 'datos' debe contener columnas: species, lon, lat")
  }
  if (nrow(datos_clean) < 30) {
    stop("Se requieren al menos 30 registros para proceder con la partición. Registros actuales: ", nrow(datos_clean))
  }
  
  # Crear identificador de coordenada para evitar duplicados entre particiones
  datos_clean$lon <- as.numeric(datos_clean$lon)
  datos_clean$lat <- as.numeric(datos_clean$lat)
  datos_clean$lon_lat <- paste0(round(datos_clean$lon, 6), "_", round(datos_clean$lat, 6))
  
  # Asegurarnos de que no haya filas NA
  datos_clean <- datos_clean %>% dplyr::filter(!is.na(lon) & !is.na(lat))
  
  # Paso 1: seleccionar test externo (por coordenadas únicas)
  unique_coords <- unique(datos_clean$lon_lat)
  n_test_coords <- max(1, round(length(unique_coords) * prop_test_externo))
  test_coords <- sample(unique_coords, n_test_coords)
  
  test <- datos_clean %>% dplyr::filter(lon_lat %in% test_coords)
  train_all <- datos_clean %>% dplyr::filter(!lon_lat %in% test_coords)
  
  # Paso 2: dentro de train_all, separar joint e ind
  unique_train_coords <- unique(train_all$lon_lat)
  n_joint_coords <- max(1, round(length(unique_train_coords) * prop_internal))
  joint_coords <- sample(unique_train_coords, n_joint_coords)
  
  joint <- train_all %>% dplyr::filter(lon_lat %in% joint_coords)
  ind <- train_all %>% dplyr::filter(!lon_lat %in% joint_coords)
  
  # Preparar data.frames con columnas en orden species, lon, lat
  df_train <- train_all %>% dplyr::select(species, lon, lat)
  df_test  <- test %>% dplyr::select(species, lon, lat)
  df_joint <- joint %>% dplyr::select(species, lon, lat)
  df_ind   <- ind %>% dplyr::select(species, lon, lat)
  
  # Nombres de archivo (consistentes con kuenm/MaxEnt)
  file_train <- file.path(dir_salida, paste0(prefijo, "_train.csv"))
  file_test  <- file.path(dir_salida, paste0(prefijo, "_test.csv"))
  file_joint <- file.path(dir_salida, paste0(prefijo, "_joint.csv"))
  file_ind   <- file.path(dir_salida, paste0(prefijo, "_ind.csv"))
  
  # Guardar (sobrescribe si existen)
  write.csv(df_train, file_train, row.names = FALSE)
  write.csv(df_test,  file_test,  row.names = FALSE)
  write.csv(df_joint, file_joint, row.names = FALSE)
  write.csv(df_ind,   file_ind,   row.names = FALSE)
  
  # Mensajes resumen
  message("Partición completada para: ", prefijo)
  message(" - Train (total): ", nrow(df_train), " registros (", length(unique(df_train$lon)), " coordenadas únicas)")
  message(" - Test (externo): ", nrow(df_test), " registros (", length(unique(df_test$lon)), " coordenadas únicas)")
  message(" - Joint (dentro de train): ", nrow(df_joint), " registros")
  message(" - Ind (dentro de train): ", nrow(df_ind), " registros")
  
  invisible(list(train = file_train,
                 test = file_test,
                 joint = file_joint,
                 ind = file_ind))
}

# Ejecutar la división usando el objeto 'datos' (ya generado antes en el script)
particiones <- dividir_datos(datos = datos_clean,
                             dir_salida = dir_kuenm,
                             prefijo = especie,
                             seed = 123,
                             prop_test_externo = 0.2,
                             prop_internal = 0.8)

# Verificar archivos creados
message("\nArchivos escritos en dir_kuenm:")
print(list.files(dir_kuenm, pattern = paste0("^", especie), full.names = TRUE))

####CALIBRACIÓN DE MODELOS CON KUENM####
# 1. Configurar parámetros ----
setwd(file.path(dir_modelo))
dir()
maxent_path <- "/Users/hibraim/Desktop/Frontiers/maxent/"
occ_joint <- paste0(especie, "_joint.csv")
occ_tra <- paste0(especie, "_train.csv")
occ_test <- paste0(especie, "_test.csv")
M_var_dir <- "./M_variables_recortadas/" # Relativo al directorio actual
list.dirs(M_var_dir, recursive = FALSE)
list.files(file.path(M_var_dir, "Set_1"), pattern = "\\.asc$")
batch_cal <- "Candidate_Models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
wait <- FALSE
run <- TRUE

# 1.5 Remover archivos fantasma "SOLO EN MACOS"
m_var_path <- "./M_variables_recortadas/"

# Function to remove hidden files recursively
remove_hidden_files <- function(dir_path) {
  # Get all files recursively
  all_files <- list.files(dir_path, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  
  # Identify hidden files (starting with ._)
  hidden_files <- grep("/\\._", all_files, value = TRUE)
  
  # Remove hidden files
  if (length(hidden_files) > 0) {
    file.remove(hidden_files)
    cat("Removed", length(hidden_files), "hidden files\n")
  } else {
    cat("No hidden files found\n")
  }
}

# Remove hidden files from all subdirectories in M_variables
subdirs <- list.dirs(m_var_path, recursive = TRUE)
for (subdir in subdirs) {
  remove_hidden_files(subdir)
}

# 2. Ejecutar calibración ----
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, 
          batch = batch_cal, out.dir = out_dir, reg.mult = reg_mult, 
          f.clas = f_clas, maxent.path = maxent_path, wait = FALSE, run = run)

# Verificar resultados
message("\nArchivos generados en calibración:")
print(list.files(out_dir))

####EVALUACIÓN DE MODELOS####
# 1. Configurar parámetros ----
occ_ind <- paste0(especie, "_ind.csv")
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE

# 2. Evaluar modelos ----
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)

# 3. Ver resultados ----
message("\nMejores modelos seleccionados:")
best_models <- read.csv(file.path(out_eval, "selected_models.csv"))
print(best_models)

############################################
######### The model creation ###############
############################################
batch_fin <- "Final_models2"
mod_dir <- "Final_Models"
rep_n <- 10 #el numero de repeticiones que quieres que haga el modelo
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "cloglog" #formato que ocupa MaxEnt
project <- TRUE #en caso de proyectar para futuro o pasado se tiene que poner T, si no se desea, debe ser FALSE
G_var_dir <- "./G_variables_recortadas" #el nombre de la carpeta donde estan las capas para modelar en tiempo
ext_type <- "ext_clam" #es clamping o puede ser extrapolatitudee "ext" o bien "all"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions
kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


############################################
####### The model evaluation final #########
############################################
replicates <- TRUE
threshold <- 10
out_feval <- "Final_Models_evaluation"

# Ejecutar evaluación final
fin_eval <- kuenm_feval(
  path = mod_dir, 
  occ.joint = occ_joint, 
  occ.ind = occ_ind, 
  replicates = replicates,
  out.eval = out_feval, 
  threshold = threshold, 
  rand.percent = rand_percent,
  iterations = iterations, 
  parallel.proc = FALSE
)

best <- read.csv("Calibration_results/selected_models.csv")
knitr::kable(best, caption = "Models selected based on significance, omission rates, and AICc, in that order.")
