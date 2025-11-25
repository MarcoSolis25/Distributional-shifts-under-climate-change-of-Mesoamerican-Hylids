#### =====================================
####   CORRELATIVE NICHE MODEL CONSTRUCTION
#### =====================================
# Pérez-Mendoza et al.(2025). Distributional shifts under climate change of Mesoamerican Hylids,
#a physiological approach considering their biphasic life cycle

### LOAD LIBRARIES ###
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


#### =====================================
####  DATA CURATION
#### =====================================
### 1. Set up directories and parameters ###
especie <- "Boana_rosenbergi"  # Species name

dir_raiz <- "/Users/hibraim/Desktop/Frontiers"
dir_presencias <- file.path(dir_raiz, "ModelConstruction", "Presencias")
dir_m_variables <- file.path(dir_raiz, "ModelConstruction", "VariablesAmbientales", "M_variables")
dir_g_variables <- file.path(dir_raiz, "ModelConstruction", "VariablesAmbientales", "G_variables")
dir_capas <- file.path(dir_raiz, "ModelConstruction", "Ms_shp")


### 2. Create directories structure ###
dir_modelo <- file.path(dir_raiz, "ModelConstruction", "Modelos", especie)
dir.create(dir_modelo, recursive = TRUE, showWarnings = FALSE)

dir_M_out <- file.path(dir_modelo, "M_variables")
dir.create(dir_M_out, recursive = TRUE, showWarnings = FALSE)

dir_G_out <- file.path(dir_modelo, "G_variables")
dir.create(dir_G_out, recursive = TRUE, showWarnings = FALSE)

message("Directories created in: ", dir_modelo)

### 3. Read presence data ###
archivo_presencias <- file.path(dir_presencias, paste0(especie, ".csv"))
if (!file.exists(archivo_presencias)) {
  stop(paste("Presence file no finded:", archivo_presencias))
}

data <- read.csv(archivo_presencias)

# verify structure
if (!all(c("species", "lon", "lat") %in% colnames(data))) {
  stop("Presence file must contain the following columns: species, lon, lat")
}

# Create spatial object and eliminate NA values
datos <- data.frame(
  species = data$species,
  lon = as.numeric(data$lon),
  lat = as.numeric(data$lat)
) %>% na.omit()

datos_sp <- SpatialPointsDataFrame(datos[,2:3], datos)

# Eliminate duplicates on the same cell
archivos_capas <- list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE)
if (length(archivos_capas) == 0) stop("No files founded on M_variables")

capa_referencia <- raster(archivos_capas[1])
celdas <- cellFromXY(capa_referencia, datos_sp)
indices_unicos <- !duplicated(celdas)
datos_sp <- datos_sp[indices_unicos, ]
datos <- datos[indices_unicos, ]

message("Unique records per cell: ", nrow(datos))
if (sum(!indices_unicos) > 0) message("We delete ", sum(!indices_unicos), " duplicated points")

### 4. Environmental variables analysis ###
# Load environmental variables
capas_presente <- stack(list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE))
message("Variables loaded: ", paste(names(capas_presente), collapse = ", "))

# Environmental values extraction
presencias_clima <- extract(capas_presente, datos_sp) %>%
  data.frame() %>%
  cbind(datos, .) %>%
  na.omit()

# Save coordinates for Maxent
write.csv(presencias_clima[,1:3],
          file.path(dir_modelo, paste0(especie, "_maxent.csv")),
          row.names = FALSE)

# Correlation matrix
cormatriz <- cor(presencias_clima[,4:ncol(presencias_clima)], use = "complete.obs")
png(file.path(dir_modelo, paste0(especie, "_correlacion_var.png")),
    width = 8, height = 8, units = "in", res = 300)
corrplot(cormatriz, method = "circle", type = "upper",
         tl.col = "black", tl.srt = 45,
         title = paste("Correlation of", especie))
dev.off()

# VIF analysis
vif_result <- vifstep(presencias_clima[,4:ncol(presencias_clima)], th = 10)
message("\nSelected variables after VIF analysis:")
print(vif_result)

write.csv(vif_result@results, file.path(dir_modelo, "vif_results.csv"))
vars_seleccionadas <- vif_result@results$Variables

message("Selected variables: ", paste(vars_seleccionadas, collapse = ", "))

# Filter only selected variables
capas_presente_filtradas <- subset(capas_presente, vars_seleccionadas)

### Accesible area delimitation (M's) ###
# 1. Load each species shapefile----
archivo_shp <- file.path(dir_capas, paste0(especie, "_M.shp"))

if (!file.exists(archivo_shp)) {
  stop(paste("No accesible shapefile were found for the species:", archivo_shp))
}

M_sf <- st_read(archivo_shp)

# 2. Geometry validation ----
M_sf <- st_make_valid(M_sf)

# 3. Match CRS to the environmental data set (if necessary)----
crs_ref <- st_crs(capa_referencia)
if (st_crs(M_sf) != crs_ref) {
  M_sf <- st_transform(M_sf, crs_ref)
  message("Specie's shapefile was reprojected to the environmental layers CRS")
}

# 4. Save final accessible area ----
st_write(M_sf, file.path(dir_modelo, paste0(especie, "_M.shp")), delete_dsn = TRUE)

# 5. Quick visualization ----
plot(st_geometry(M_sf), main = paste("Accesible area (M) of", especie))
plot(datos_sp, col = "red", pch = 19, cex = 0.5, add = TRUE)


### ==============================================
###   LAYERS PROCESS (M and G variables)
### ==============================================

library(terra)
library(dplyr)

# 1. Load reference
ref_path <- "/Users/hibraim/Desktop/Frontiers/annual hours_pools.tif"
if (!file.exists(ref_path)) stop("No reference file was found")

ref_raster <- rast(ref_path)
message("✅ Reference file was loaded in:", ref_path)

# 2. Load M polygon
archivo_M <- file.path(dir_modelo, paste0(especie, "_M.shp"))
if (!file.exists(archivo_M)) stop(paste("No shapefile was founded:", archivo_M))

poligono <- vect(archivo_M)
poligono <- project(poligono, crs(ref_raster))
message("✅ Loaded polygon and reprojected to the CRS")

# 3. Selected variables by VIF
archivo_vif <- file.path(dir_modelo, "vif_results.csv")
if (!file.exists(archivo_vif)) stop("No VIF results file founded.")

vif_table <- read.csv(archivo_vif)
vars_vif <- vif_table$Variables %>% unique() %>% na.omit()
message("✅ Selected variables based on VIF: ", paste(vars_vif, collapse = ", "))

# Maps correction based on names- VIF to G
mapeo_vif_a_G <- function(nombres_vif) {
  nombres_G <- c()
  
  for (nombre in nombres_vif) {
    # Case: bio_1, bio_2, etc. -> transform to bio01, bio02
    if (grepl("bio_(\\d+)", nombre)) {
      numero <- gsub("bio_(\\d+)", "\\1", nombre)
      # Secure to digits
      numero_formateado <- sprintf("%02d", as.numeric(numero))
      nombre_G <- paste0("bio", numero_formateado)
    }
    # Case: bio1, bio2, etc. -> transform to bio01, bio02
    else if (grepl("bio(\\d+)", nombre)) {
      numero <- gsub("bio(\\d+)", "\\1", nombre)
      numero_formateado <- sprintf("%02d", as.numeric(numero))
      nombre_G <- paste0("bio", numero_formateado)
    }
    # Case: they already have the format: bio01, bio02
    else if (grepl("bio\\d{2}", nombre)) {
      nombre_G <- nombre
    }
    # Case: other names, stay the same
    else {
      nombre_G <- nombre
    }
    nombres_G <- c(nombres_G, nombre_G)
  }
  
  return(nombres_G)
}

# Aplicate mapping
vars_vif_G <- mapeo_vif_a_G(vars_vif)
message("✅ VIF variables mapped to G:", paste(vars_vif_G, collapse = ", "))

# Inverse mapping for report
mapeo_G_a_vif <- function(nombres_G) {
  nombres_vif <- c()
  
  for (nombre in nombres_G) {
    # Transform bio01, bio02 to bio_1, bio_2
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


# 4. Layer's process function
procesar_capas_selectivamente <- function(rutas_tif, poligono, ref_raster, vars_seleccionadas, dir_salida, tipo = "M") {
  dir.create(dir_salida, recursive = TRUE, showWarnings = FALSE)
  capas_procesadas <- c()
  
  for (ruta in rutas_tif) {
    message("\n🔍 Analysing file: ", basename(ruta))
    
    # Read multilayer file
    capas <- rast(ruta)
    nombres_originales <- names(capas)
    
    message("   Available files: ", paste(nombres_originales, collapse = ", "))
    message("   Searching variables: ", paste(vars_seleccionadas, collapse = ", "))
    
    # Find matching variables
    vars_coincidentes <- intersect(nombres_originales, vars_seleccionadas)
    
    if (length(vars_coincidentes) == 0) {
      message("   ⚠️ No matching variables were found")
      next
    }
    
    message("   ✅ Matching variables ", paste(vars_coincidentes, collapse = ", "))
    
    # Filter by selected variables
    capas_filtradas <- subset(capas, vars_coincidentes)
    
    # Aling with reference raster (resolution & CRS only - NO extent)
    message("   🔄 Aligning resolution and CRS with reference...")
    
    # Create a target template with ref_raster properties but original extent
    target_template <- rast()
    crs(target_template) <- crs(ref_raster)  # Same CRS
    res(target_template) <- res(ref_raster)  # Same resolution
    ext(target_template) <- ext(capas_filtradas)  # BUT original extent
    
    # Now resample to this custom template
    capas_alineadas <- terra::resample(capas_filtradas, target_template, method = "bilinear")
    
    # Cut and mask polygon
    message("   ✂️ Cutting and masking...")
    capas_recortadas <- crop(capas_alineadas, poligono, snap = "out")
    capas_enmascaradas <- mask(capas_recortadas, poligono)
    
    # Save individual layers
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
      message("   💾 Saved: ", nombre, ".asc")
    }
  }
  
  message("\n✅ Completed process to type: ", tipo)
  message("   Saved layers: ", length(capas_procesadas))
  return(capas_procesadas)
}
# List M's files
rutas_M <- list.files(dir_m_variables, pattern = "\\.tif$", full.names = TRUE)
message("M's files founded: ", length(rutas_M))

if (length(rutas_M) == 0) {
  stop("❌ No .tif files were found in the M_variables directory")
}

# Verification of M files
for (ruta in rutas_M) {
  capas <- rast(ruta)
  message("\n📁 M file: ", basename(ruta))
  message("   Availible files: ", paste(names(capas), collapse = ", "))
  coincidentes_M <- intersect(names(capas), vars_vif)
  message("   VIF match: ", ifelse(length(coincidentes_M) > 0, paste(coincidentes_M, collapse = ", "), "NINGUNA"))
}

# Process of M varaibles
dir_M_out <- file.path(dir_modelo, "M_variables_recortadas/Set_1")
capas_M_procesadas <- procesar_capas_selectivamente(rutas_M, poligono, ref_raster, vars_vif, dir_M_out, "M")


# 6. Verification and process of G variables
subdirs_G <- list.dirs(dir_g_variables, recursive = FALSE, full.names = TRUE)
message("G scenarios founded: ", length(subdirs_G))

if (length(subdirs_G) == 0) {
  stop("❌ No subdirectories founded in G_varaibles")
}

# Show complete mapping for verification
message("\n🔍 Complete varaible mapping:")
for (i in 1:length(vars_vif)) {
  message("   VIF: ", vars_vif[i], " → G: ", vars_vif_G[i])
}

dir_G_out <- file.path(dir_modelo, "G_variables_recortadas/Set_1")
total_capas_G <- 0

for (sub in subdirs_G) {
  nombre_sub <- basename(sub)
  message("\n📂 Processing G scenario: ", nombre_sub)
  
  rutas_G <- list.files(sub, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(rutas_G) == 0) {
    message("   ⚠️ No .tif files were found")
    next
  }
  
  dir_G_sub_out <- file.path(dir_G_out, nombre_sub)
  capas_G_procesadas <- procesar_capas_selectivamente(rutas_G, poligono, ref_raster, vars_vif_G, dir_G_sub_out, "G")
  total_capas_G <- total_capas_G + length(capas_G_procesadas)
}


# 7. Final verification and report
# Report of M
message("\n📊 M's variables processed :")
if (length(capas_M_procesadas) > 0) {
  message("   ✅ Total: ", length(capas_M_procesadas), " layers")
  for (capa in capas_M_procesadas) {
    message("      ✓ ", basename(capa))
  }
} else {
  message("   ❌ No M layers were processed")
}

# Report of G
message("\n📊 G's variables processed:")
if (total_capas_G > 0) {
  message("   ✅ Total: ", total_capas_G, " layers")
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
  message("   ❌ No G layers were processed")
}

# VIF's layers covering verification
message("\n🔍 VIF covering verification:")
message("   Original VIF variables: ", paste(vars_vif, collapse = ", "))

todas_capas <- c(capas_M_procesadas, 
                 list.files(dir_G_out, pattern = "\\.asc$", recursive = TRUE, full.names = TRUE))

nombres_procesados <- tools::file_path_sans_ext(basename(todas_capas))

# To verification, transform processed G names to VIF format
nombres_procesados_vif <- mapeo_G_a_vif(nombres_procesados)

# M and G double verification
vars_faltantes <- setdiff(vars_vif, nombres_procesados_vif)
if (length(vars_faltantes) == 0) {
  message("   ✅ ¡All VIF's variables were correctly processed!")
} else {
  message("   ⚠️ Missing VIF varaibles: ", paste(vars_faltantes, collapse = ", "))
}

message("\n🎯 Completed process for : ", especie)


# 8. Rename G files to match M names (all files)
renombrar_archivos_G <- function(dir_G_base) {
  message("\n🔄 Renaming G's files to match M's names...")
  
  # Function to convert names bio01 -> bio_1, etc.
  convertir_nombre <- function(nombre_archivo) {
    nombre_sin_ext <- tools::file_path_sans_ext(nombre_archivo)
    extension <- tools::file_ext(nombre_archivo)
    
    # If format is bio01, bio02, etc., turn to bio_1, bio_2
    if (grepl("^bio\\d{2}$", nombre_sin_ext)) {
      numero <- as.numeric(gsub("bio(\\d{2})", "\\1", nombre_sin_ext))
      nuevo_nombre <- paste0("bio_", numero)
      return(paste0(nuevo_nombre, ".", extension))
    }
    # If format is already bio_1, bio_2, keep it
    else if (grepl("^bio_\\d+$", nombre_sin_ext)) {
      return(nombre_archivo)
    }
    # For other names, keep it
    else {
      return(nombre_archivo)
    }
  }
  
  # Look for all subdirectories in G_varaibles_recortadas
  subdirs <- list.dirs(dir_G_base, recursive = TRUE, full.names = TRUE)
  total_renombrados <- 0
  
  for (subdir in subdirs) {
    # Look for all files in the subdirectory (not only .asc)
    todos_archivos <- list.files(subdir, full.names = FALSE)
    
    if (length(todos_archivos) > 0) {
      message("📂 Processing directory: ", basename(subdir))
      
      # Group files by name (without extension) 
      archivos_por_base <- list()
      
      for (archivo in todos_archivos) {
        nombre_base <- tools::file_path_sans_ext(archivo)
        if (!nombre_base %in% names(archivos_por_base)) {
          archivos_por_base[[nombre_base]] <- character()
        }
        archivos_por_base[[nombre_base]] <- c(archivos_por_base[[nombre_base]], archivo)
      }
      
      # Process each group of files
      for (nombre_base in names(archivos_por_base)) {
        archivos_grupo <- archivos_por_base[[nombre_base]]
        nuevo_nombre_base <- tools::file_path_sans_ext(convertir_nombre(paste0(nombre_base, ".asc")))
        
        # If base-name changed, rename all files within the group 
        if (nuevo_nombre_base != nombre_base) {
          for (archivo in archivos_grupo) {
            extension <- tools::file_ext(archivo)
            archivo_viejo <- file.path(subdir, archivo)
            archivo_nuevo <- file.path(subdir, paste0(nuevo_nombre_base, ".", extension))
            
            # Verify that the new files does not longer exists.
            if (!file.exists(archivo_nuevo)) {
              file.rename(archivo_viejo, archivo_nuevo)
              message("   ✅ Renamed: ", archivo, " -> ", paste0(nuevo_nombre_base, ".", extension))
              total_renombrados <- total_renombrados + 1
            } else {
              message("   ⚠️ The files is no longer, omitting: ", paste0(nuevo_nombre_base, ".", extension))
            }
          }
        }
      }
    }
  }
  
  message("🎯 Total renamed files: ", total_renombrados)
  return(total_renombrados)
}

# Execute the renamed  after G variables process
dir_G_base <- file.path(dir_modelo, "G_variables_recortadas")
if (dir.exists(dir_G_base)) {
  total_renombrados <- renombrar_archivos_G(dir_G_base)
  
  # Final verification
  message("\n🔍 Name's final verification:")
  
  # Obtaining base-names from all files within M (no extensions)
  archivos_M <- list.files(dir_M_out, full.names = FALSE)
  nombres_base_M <- unique(tools::file_path_sans_ext(archivos_M))
  
  # Obtaining base-names from all files within G (no extensions)
  archivos_G <- list.files(dir_G_base, recursive = TRUE, full.names = FALSE)
  nombres_base_G <- unique(tools::file_path_sans_ext(archivos_G))
  
  message("Base-names in M: ", paste(sort(nombres_base_M), collapse = ", "))
  message("Base-names in G: ", paste(sort(nombres_base_G), collapse = ", "))
  
  if (setequal(nombres_base_M, nombres_base_G)) {
    message("✅ ¡Names perfectly match!")
  } else {
    message("❌ There are differences in base-names")
    message("   Missing in G: ", paste(setdiff(nombres_base_M, nombres_base_G), collapse = ", "))
    message("   Missing in M: ", paste(setdiff(nombres_base_G, nombres_base_M), collapse = ", "))
  }
  
  # Showing each file per varaible to verification
  message("\n📋 Files per variable (example for bio_1):")
  ejemplo_var <- "bio_1"
  archivos_ejemplo_M <- list.files(dir_M_out, pattern = paste0("^", ejemplo_var, "\\."), full.names = FALSE)
  archivos_ejemplo_G <- list.files(dir_G_base, pattern = paste0("^", ejemplo_var, "\\."), recursive = TRUE, full.names = FALSE)
  
  message("Files for ", ejemplo_var, " in M: ", paste(archivos_ejemplo_M, collapse = ", "))
  message("Files for ", ejemplo_var, " in G: ", paste(unique(archivos_ejemplo_G), collapse = ", "))
  
} else {
  message("❌ G_variables_recortadas directory not founded")
}



###### CORRELATIVE MODEL CREATION WITH KUENM AND MAXENT ####

# =====================================
# Cleaning of records with no environmental data
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
# Record selection (train / test / joint / ind)
# =====================================

dir_kuenm <- file.path(dir_modelo)
datos_clean<- read.csv2("datos_clean.csv", sep=",", header=TRUE)
# Function to split data and save files with species, lon and lat columns
dividir_datos <- function(datos_clean,
                          dir_salida = dir_kuenm,
                          prefijo = especie,
                          seed = 42,
                          prop_test_externo = 0.2,   # External test proportion
                          prop_internal = 0.8) {     # Proportion for joint in train (the rest will be ind)
  set.seed(seed)
  
  # Basic verification
  if (!all(c("species", "lon", "lat") %in% colnames(datos_clean))) {
    stop("data.frame must contain the following columns: species, lon, lat")
  }
  if (nrow(datos_clean) < 30) {
    stop("At least 30 records are necessary for data split. Actual records: ", nrow(datos_clean))
  }
  
  # Creation of coordinate ID to avoid duplicates within partitioning. 
  datos_clean$lon <- as.numeric(datos_clean$lon)
  datos_clean$lat <- as.numeric(datos_clean$lat)
  datos_clean$lon_lat <- paste0(round(datos_clean$lon, 6), "_", round(datos_clean$lat, 6))
  
  # Verification of no NA rows
  datos_clean <- datos_clean %>% dplyr::filter(!is.na(lon) & !is.na(lat))
  
  # Step 1: selection of external test (by unique coordinates)
  unique_coords <- unique(datos_clean$lon_lat)
  n_test_coords <- max(1, round(length(unique_coords) * prop_test_externo))
  test_coords <- sample(unique_coords, n_test_coords)
  
  test <- datos_clean %>% dplyr::filter(lon_lat %in% test_coords)
  train_all <- datos_clean %>% dplyr::filter(!lon_lat %in% test_coords)
  
  # Step 2: within train_all, split joint & ind
  unique_train_coords <- unique(train_all$lon_lat)
  n_joint_coords <- max(1, round(length(unique_train_coords) * prop_internal))
  joint_coords <- sample(unique_train_coords, n_joint_coords)
  
  joint <- train_all %>% dplyr::filter(lon_lat %in% joint_coords)
  ind <- train_all %>% dplyr::filter(!lon_lat %in% joint_coords)
  
  # Create a data.frames with column order: species, lon, lat
  df_train <- train_all %>% dplyr::select(species, lon, lat)
  df_test  <- test %>% dplyr::select(species, lon, lat)
  df_joint <- joint %>% dplyr::select(species, lon, lat)
  df_ind   <- ind %>% dplyr::select(species, lon, lat)
  
  # Naming of files (consistent with kuenm/MaxEnt)
  file_train <- file.path(dir_salida, paste0(prefijo, "_train.csv"))
  file_test  <- file.path(dir_salida, paste0(prefijo, "_test.csv"))
  file_joint <- file.path(dir_salida, paste0(prefijo, "_joint.csv"))
  file_ind   <- file.path(dir_salida, paste0(prefijo, "_ind.csv"))
  
  # Save (overwrite if necessary)
  write.csv(df_train, file_train, row.names = FALSE)
  write.csv(df_test,  file_test,  row.names = FALSE)
  write.csv(df_joint, file_joint, row.names = FALSE)
  write.csv(df_ind,   file_ind,   row.names = FALSE)
  
  # Summary message
  message("Splitting complete for: ", prefijo)
  message(" - Train (total): ", nrow(df_train), " records (", length(unique(df_train$lon)), " unique coordinates)")
  message(" - Test (external): ", nrow(df_test), " records (", length(unique(df_test$lon)), " unique coordinates)")
  message(" - Joint (within train): ", nrow(df_joint), " records")
  message(" - Ind (within train): ", nrow(df_ind), " records")
  
  invisible(list(train = file_train,
                 test = file_test,
                 joint = file_joint,
                 ind = file_ind))
}

# Execute splitting using the object "datos" (already created in the above script)
particiones <- dividir_datos(datos = datos_clean,
                             dir_salida = dir_kuenm,
                             prefijo = especie,
                             seed = 123,
                             prop_test_externo = 0.2,
                             prop_internal = 0.8)

# Created files verification
message("\nWrited files in dir_kuenm:")
print(list.files(dir_kuenm, pattern = paste0("^", especie), full.names = TRUE))


# =====================================
# Models calibration with kuenm
# =====================================
# 1. Parameters setting ----
setwd(file.path(dir_modelo))
dir()
maxent_path <- "/Users/hibraim/Desktop/Frontiers/maxent/"
occ_joint <- paste0(especie, "_joint.csv")
occ_tra <- paste0(especie, "_train.csv")
occ_test <- paste0(especie, "_test.csv")
M_var_dir <- "./M_variables_recortadas/" # Relative to actual directory
list.dirs(M_var_dir, recursive = FALSE)
list.files(file.path(M_var_dir, "Set_1"), pattern = "\\.asc$")
batch_cal <- "Candidate_Models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
f_clas <- "all"
args <- NULL
wait <- FALSE
run <- TRUE

# 1.5 Remove ghost files "ONLY IN MACOS"
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

# 2. Execute calibration ----
kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, 
          batch = batch_cal, out.dir = out_dir, reg.mult = reg_mult, 
          f.clas = f_clas, maxent.path = maxent_path, wait = FALSE, run = run)

# Results verification
message("\nFiles generated in the calibration:")
print(list.files(out_dir))

#### MODEL EVALUATION ####
# 1. Parameters setting ----
occ_ind <- paste0(especie, "_ind.csv")
out_eval <- "Calibration_results"
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE

# 2. Model evaluation ----
cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)

# 3. See results ----
message("\nBest model selected:")
best_models <- read.csv(file.path(out_eval, "selected_models.csv"))
print(best_models)

############################################
######### The model creation ###############
############################################
batch_fin <- "Final_models2"
mod_dir <- "Final_Models"
rep_n <- 10 # number of repetitions for the final model
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "cloglog" # Maxent's format file 
project <- TRUE # In case of projection to future =TRUE (In this case must be TRUE)
G_var_dir <- "./G_variables_recortadas" # Folder's name were all time scenarios layers are present
ext_type <- "ext_clam" # Clamping option (available options are "ext" or "all")
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
####### The model final evaluation #########
############################################
replicates <- TRUE
threshold <- 10
out_feval <- "Final_Models_evaluation"

# Execute final model evaluation
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
