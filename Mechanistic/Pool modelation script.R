###############################################################################################
########## Pool modelation with NicheMapR for the manuscript by Pérez-Mendoza et al. ########## 
############# "Distributional shifts under climate change of Mesoamerican Hylids, ############# 
############### a physiological approach considering their biphasic life cycle" ###############
############### To be submitted to Frontiers in Amphibian and Reptile Science #################
###############################################################################################


#### Divide the area into a grid and extract the coordinates of each cell ####

library(terra)
# Load the reference shapefile 
shp <- vect("mesoamerica.shp")

# Create and rasterize an empty raster at the desired resolution (5', 0.0833°)
base_raster <- rast(ext(shp), resolution = 0.083333) 
crs(base_raster) <- crs(shp) # Asign the same reference system
mesoam <- rasterize(shp, base_raster)

# Replace NaN cells for NA, count the target cells and save the raster
values(mesoam)[is.nan(values(mesoam))] <- NA
cells_no_na <- sum(!is.na(values(mesoam)))
writeRaster(mesoam, "Mesoam_5min.tif", overwrite = TRUE)

# Convert the raster into a data frame and save the coordinates
mesoamdf <- as.data.frame(mesoam, xy = TRUE)
colnames(mesoamdf)[1:2] <- c("lon", "lat")
write.csv(mesoamdf, "mesoam_5min.csv", row.names = FALSE)

# Visualize the raster with the points coordinates
points_mesoamcsv <- read.csv("mesoam_5min.csv")
points_mesoamvect <- vect(points_mesoamcsv, geom = c("lon", "lat"), crs = "EPSG:32614") 
plot(mesoam)
points(points_mesoamvect, col = "red", pch = 10, cex = 0.1)


#### Pool modelation with NicheMapR ####

# This part of the code runs a loop of the micro_global() from NicheMapR
# A model is run for each coordinate
# The results are stored with the coordinates as the names of the files

library(NicheMapR)
library(dplyr)
library(terra)

# Load points and create directories of the model's target results.
points <- read.csv("mesoam_5min.csv")
if (!dir.exists("soil")) dir.create("soil")
if (!dir.exists("metout")) dir.create("metout")

# Define a clay soil parameters and the hydraulic properties of such soil
# (According to Enriquez-Urtzelai et al. (2019) and comments by Michael Kearney in the NicheMapR Google group

soilpro <- matrix(data = 0, nrow = 6, ncol = 5)
colnames(soilpro) <- c('depth', 'blkdens', 'clay', 'silt', 'sand')
soilpro[,1] <- c(2.5, 7.5, 22.5, 45, 80, 150)
soilpro[,2] <- 1.4 # bulk density (Mg/m3)
soilpro[,3] <- 90 # % clay
soilpro[,4] <- 0.0 # % silt
soilpro[,5] <- 10 # % sand
DEP <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200)

soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro), DEP = DEP)
PE <- soil.hydro$PE
BB <- soil.hydro$BB
BD <- soil.hydro$BD
KS <- soil.hydro$KS
BulkDensity <- BD[seq(1, 19, 2)]

# Other micro_global() parameters
rainmult <- 4 # Rainfall multiplier to reflect catchment (1 will result in unadjusted rainfall)
LAI <- 0 # leaf area index (zero so no transpiration)
cap <- 0 # do not simulate the organic cap layer
soilgrids <- 0 # don't use SoilGrids soil database for thermal properties
timeinterval <- 365
Usrhyt <- 0.1

#Define the function to run the microclimate model for a single point
run_point <- function(lat, lon, soil_csv, metout_csv,
                      PE, BB, BD, KS, rainmult, LAI, cap, soilgrids, timeinterval, Usrhyt) {
  
  micro <- micro_global(loc = c(lon, lat), runmoist = 1, rainfrac = 0.1, Usrhyt = Usrhyt,
                        timeinterval = timeinterval, PE = PE, BB = BB, BD = BD, KS = KS,
                        rainmult = rainmult, LAI = LAI, cap = cap, soilgrids = soilgrids)
  
  # Save both files
  write.csv(micro$soil, soil_csv, row.names = FALSE)
  write.csv(micro$metout, metout_csv, row.names = FALSE)
  
  return("ok")
}

# Run the model loop with retry and error logging
start_time <- Sys.time()
errors <- list()   # Store potential errors (coordinates where the model failed)

for (i in 1:nrow(points)) {
  lat <- points$lat[i]
  lon <- points$lon[i]
  
  soil_csv <- file.path("soil", paste0(lat, ", ", lon, "_soil.csv"))
  metout_csv <- file.path("metout", paste0(lat, ", ", lon, "_metout.csv"))
  
  # Check if soil and metout files already exist
  if (file.exists(soil_csv) && file.exists(metout_csv)) {
    cat("Existing files, skipping:", lat, lon, 
        "(", i, "of", nrow(points), ")\n")
    next
  }
  
  success <- FALSE
  attempts <- 0
  while (!success && attempts < 2) {
    attempts <- attempts + 1
    cat("Processing point", lat, lon, "- attempt", attempts, 
        "(", i, "of", nrow(points), ")\n")
    
    res <- tryCatch({
      run_point(lat = lat, lon = lon, 
                soil_csv = soil_csv,
                metout_csv = metout_csv,
                PE = PE, BB = BB, BD = BD, KS = KS,
                rainmult = rainmult, LAI = LAI, cap = cap,
                soilgrids = soilgrids, timeinterval = timeinterval,
                Usrhyt = Usrhyt)
    }, error = function(e) {
      cat("Error in point", lat, lon, ":", conditionMessage(e), "\n")
      return("fail")
    })
    
    if (identical(res, "ok")) {
      success <- TRUE
      cat("Point processed successfully:", lat, lon, 
          "(", i, "of", nrow(points), ")\n")
    }
  }
  
  if (!success) {
    cat("Failed twice, skipping point:", lat, lon, 
        "(", i, "of", nrow(points), ")\n")
    errors[[length(errors) + 1]] <- paste(lat, lon)
  }
}

# Save error log if there were errors
if (length(errors) > 0) {
  writeLines(unlist(errors), "errors_log.txt")
  cat("Saved log with", length(errors), "errors in 'errors_log.txt'\n")
} else {
  cat("No errors, no log created\n")
}

end_time <- Sys.time()
cat("Total execution time:", difftime(end_time, start_time, units = "mins"), "minutes\n")


#### Combine the model results for each point into a single dataframe ####
library(dplyr)

# Get the names of all .csv files and the target column, which will be the depth of the simulated pools
directory <- "/Volumes/Seagate Por/Pruebas/metout"
files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
target_column <- "POOLDEP"

# Create a list to store the extracted columns
extracted_columns <- list()

# Loop through each file and extract the desired column. 
# Each column of the new dataset will have the coordinates of each site and the modeled pool depths

for (file in files) {
  data <- read.csv(file)
  column <- data[[target_column]]
  column_name <- gsub(".*/(.*?)(_metout)?\\.csv", "\\1", file)
  extracted_columns[[column_name]] <- column
}

# Combine all extracted columns into a single data frame and export the combined dataset  
site_hours <- bind_cols(extracted_columns)
write.csv(site_hours, "hours_pooldepth.csv", row.names = FALSE)


#### Determine which modeled ponds meet the criteria for allowing larval development ####

# The criteria are: a minimum depth of 10 cm, maintained for a run of 17 to 41 days, during the months of reproduction (March to November)
# The final database describes th dates and lenght of a sample of runs 

# Load hourly pool-depth data
pool_depth <- data.table::fread("hours_pooldepth.csv")

# Binarize depth (>=10 cm → 1)
depth_binary <- apply(pool_depth, 2, function(x) ifelse(x >= 10, 1, 0))

# Collapse hourly data into daily presence/absence
# Groups every 24 hours; a day is 1 if depth >=10 cm at least once
depth_daily <- apply(depth_binary, 2, function(x) {
  day_groups <- (seq_along(x) - 1) %/% 24 + 1
  days <- tapply(x, day_groups, function(d) ifelse(sum(d) > 0, 1, 0))
  as.numeric(days)
})

# Extract run-length statistics (run-length encoding)
extract_runs <- function(x) {
  r <- rle(x)
  run_lengths <- r$lengths[r$values == 1]
  run_end <- cumsum(r$lengths)[r$values == 1]
  run_start <- run_end - run_lengths + 1
  
  data.frame(
    start_day = run_start,
    end_day   = run_end,
    duration  = run_lengths
  )
}

# Apply to each site
runs_per_site <- lapply(as.data.frame(depth_daily), extract_runs)

# Filter runs to March–November (day 60 to 334)
filter_march_november <- function(df) {
  df[df$start_day >= 60 & df$end_day <= 334, ]
}
runs_per_site <- lapply(runs_per_site, filter_march_november)

# Descriptive metrics
n_total_runs <- sapply(runs_per_site, nrow)
max_run_duration <- sapply(runs_per_site, function(df)
  ifelse(nrow(df) == 0, 0, max(df$duration))
)

# Identify valid runs (duration 17–41 days)
valid_runs <- lapply(runs_per_site, function(df)
  df[df$duration >= 17 & df$duration <= 41, ]
)

n_valid_runs <- sapply(valid_runs, nrow)

# Binary indicator: site has ≥1 valid run
valid_site <- as.numeric(n_valid_runs > 0)

# Select runs according to rules
# If a site has valid runs: keep them; fill with longest remaining runs up to 3
# Otherwise: keep the 3 longest runs
select_runs <- function(df, meets_criteria) {
  
  if (nrow(df) == 0) {
    return(data.frame(
      start_day = rep(NA, 3),
      end_day   = rep(NA, 3),
      duration  = rep(NA, 3)
    ))
  }
  
  df <- df[order(-df$duration), ]
  
  if (meets_criteria) {
    valid_df <- df[df$duration >= 17 & df$duration <= 41, ]
    
    if (nrow(valid_df) < 3) {
      filler <- head(df[!(df$duration >= 17 & df$duration <= 41), ], 
                     3 - nrow(valid_df))
      valid_df <- rbind(valid_df, filler)
    }
    
    return(valid_df)
    
  } else {
    return(head(df, 3))
  }
}

final_runs <- mapply(select_runs, runs_per_site, as.logical(valid_site),
                     SIMPLIFY = FALSE)

# Convert runs to matrices (to produce columns 1:3)
max_runs <- max(sapply(final_runs, nrow))

start_matrix <- matrix(NA, nrow = length(final_runs), ncol = max_runs)
end_matrix   <- matrix(NA, nrow = length(final_runs), ncol = max_runs)
dur_matrix   <- matrix(NA, nrow = length(final_runs), ncol = max_runs)

for (i in seq_along(final_runs)) {
  n_i <- nrow(final_runs[[i]])
  start_matrix[i, 1:n_i] <- final_runs[[i]]$start_day
  end_matrix[i, 1:n_i]   <- final_runs[[i]]$end_day
  dur_matrix[i, 1:n_i]   <- final_runs[[i]]$duration
}

# Convert day-of-year to dates
days <- 1:365
dates <- as.Date("2024-01-01") + days - 1

day_to_date <- function(d) {
  if (is.na(d)) return(NA)
  format(dates[d], "%d-%b")
}

start_date_matrix <- apply(start_matrix, 2, function(col) sapply(col, day_to_date))
end_date_matrix   <- apply(end_matrix,   2, function(col) sapply(col, day_to_date))

# Extract coordinates from column names
coords_raw <- do.call(rbind, strsplit(colnames(pool_depth), ",\\s*"))
coords_raw[, 2] <- gsub("_metout$", "", coords_raw[, 2])

lat <- as.numeric(coords_raw[, 1])
lon <- as.numeric(coords_raw[, 2])

# Interleave columns for final output
interleaved <- do.call(cbind, lapply(1:max_runs, function(i) {
  cbind(
    start_matrix[, i, drop = FALSE],
    end_matrix[, i, drop = FALSE],
    dur_matrix[, i, drop = FALSE],
    start_date_matrix[, i, drop = FALSE],
    end_date_matrix[, i, drop = FALSE]
  )
}))

colnames(interleaved) <- unlist(lapply(1:max_runs, function(i) {
  c(paste0("start_", i),
    paste0("end_", i),
    paste0("duration_", i),
    paste0("start_date_", i),
    paste0("end_date_", i))
}))

# Final output table
final_output <- data.frame(
  lat = lat,
  lon = lon,
  valid_site = valid_site,
  n_total_runs = n_total_runs,
  max_run_duration = max_run_duration,
  n_valid_runs = n_valid_runs,
  interleaved
)

length(which(valid_site == 1))

# Save
write.csv(final_output, "pool_depth_10cm_17-41days_mar-nov.csv", row.names = FALSE)


#### Determine which pixels surpass the mean CTmax of tadpoles ####
library(data.table)

# Base directory containing the *_soil.csv files
dir_soil <- "/Volumes/Seagate Por/Pruebas/soil"

# Create output folder if it does not exist
dir_out <- file.path(dir_soil, "above_41")
if (!dir.exists(dir_out)) dir.create(dir_out)

# Read the coordinates of interest (CSV with columns lat, lon) 
coords_interest <- fread("/Volumes/Seagate Por/Pruebas/pool_depth_10cm_17-41days_mar-nov.csv")

if (!all(c("lat", "lon") %in% names(coords_interest))) {
  stop("The coordinate CSV must contain columns 'lat' and 'lon'")
}

# Round and format to 7 decimals (as strings) for safe matching
coords_interest[, lat7 := sprintf("%.7f", round(as.numeric(lat), 7))]
coords_interest[, lon7 := sprintf("%.7f", round(as.numeric(lon), 7))]
coords_interest[, key := paste(lat7, lon7, sep = ", ")]

# List files in dir_soil and extract coordinates from their names 
all_files <- list.files(dir_soil, pattern = "_soil\\.csv$", full.names = FALSE)

# Remove suffix and extract lat/lon: "lat, lon_soil.csv"
names_no_suffix <- gsub("_soil\\.csv$", "", all_files)
coords_from_files_list <- strsplit(names_no_suffix, ",\\s*")

# Ensure all entries contain 2 parts
valid_names <- sapply(coords_from_files_list, length) == 2
if (any(!valid_names)) {
  warning("Some filenames do not match expected format 'lat, lon_soil.csv' and will be skipped.")
}

all_files_valid <- all_files[valid_names]
coords_from_files_list <- coords_from_files_list[valid_names]

coords_from_files <- data.frame(
  lat = sapply(coords_from_files_list, function(x) x[1]),
  lon = sapply(coords_from_files_list, function(x) x[2]),
  file = all_files_valid,
  stringsAsFactors = FALSE
)

# Convert to numeric, round, and format to 7 decimals
coords_from_files$lat7 <- sprintf("%.7f", round(as.numeric(coords_from_files$lat), 7))
coords_from_files$lon7 <- sprintf("%.7f", round(as.numeric(coords_from_files$lon), 7))
coords_from_files$key  <- paste(coords_from_files$lat7, coords_from_files$lon7, sep = ", ")

# Select only files that match coords_interest 
keys_interest <- unique(coords_interest$key)
match_indices <- which(coords_from_files$key %in% keys_interest)
filtered_files <- coords_from_files$file[match_indices]

cat("Found", length(filtered_files), "files matching the coordinates of interest (out of",
    length(all_files), "total files).\n")

# Process each selected file
threshold <- 40.65
days_mar_nov <- 60:334  # March–November (approx.)

# Counters for summary
n_with_events <- 0
n_processed   <- 0

for (file in filtered_files) {
  
  n_processed <- n_processed + 1
  file_path <- file.path(dir_soil, file)
  
  # Read the file
  data <- tryCatch({
    data.table::fread(file_path, sep = ",", quote = "\"",
                      encoding = "UTF-8", header = TRUE, fill = TRUE)
  }, error = function(e) {
    warning(sprintf("Could not read '%s' (%s). Skipping.", file, e$message))
    return(NULL)
  })
  
  if (is.null(data)) next
  
  # Check DOY and TIME columns
  if (!all(c("DOY", "TIME") %in% names(data))) {
    warning(sprintf("File '%s' does not contain DOY and TIME columns. Skipping.", file))
    next
  }
  
  # Add original row index for traceability
  data[, original_row := .I]
  
  # Convert DOY to month
  data[, DOY := as.numeric(DOY)]
  data[, month := as.numeric(format(as.Date(DOY, origin = "2020-01-01"), "%m"))]
  
  # Filter March (3) to November (11)
  data <- data[month >= 3 & month <= 11]
  if (nrow(data) == 0) next
  
  # Must have at least 11 columns
  if (ncol(data) < 11) {
    warning(sprintf("File '%s' has fewer than 11 columns. Skipping.", file))
    next
  }
  
  # Select depth columns 4–9 (columns 6:11)
  depth_cols <- data[, 6:11, with = FALSE]
  
  # Identify rows where any soil depth exceeds the threshold
  rows_above_41 <- apply(depth_cols, 1, function(x) any(as.numeric(x) > threshold))
  
  # Extract rows that meet the condition
  above_41 <- data[rows_above_41]
  
  # If matching rows exist, save to output folder
  if (nrow(above_41) > 0) {
    n_with_events <- n_with_events + 1
    out_path <- file.path(dir_out, file)
    data.table::fwrite(above_41, out_path)
  }
  
  # Show progress every 500 files
  if (n_processed %% 500 == 0) {
    cat("Processed:", n_processed, "files...\n")
  }
}

# Final summary 
cat("Files processed:", n_processed, "\n")
cat("Files with at least one hour > 40.65°C (Mar–Nov):", n_with_events, "\n")
cat("Results saved to:", dir_out, "\n")

# Create list of coordinates with events > 40.65°C 
if (n_with_events > 0) {
  
  event_files <- list.files(dir_out, pattern = "_soil\\.csv$", full.names = FALSE)
  
  # Extract lat/lon from filenames
  coords_events_list <- strsplit(gsub("_soil\\.csv$", "", event_files), ",\\s*")
  valid_event_names  <- sapply(coords_events_list, length) == 2
  
  coords_events <- data.table(
    lat = as.numeric(sapply(coords_events_list[valid_event_names], `[`, 1)),
    lon = as.numeric(sapply(coords_events_list[valid_event_names], `[`, 2))
  )
  
  coords_events <- unique(coords_events)
  out_coords_path <- file.path(dir_out, "coords_above_41.csv")
  fwrite(coords_events, out_coords_path)
  
  cat("Coordinate list saved to:", out_coords_path, "\n")
  
} else {
  cat("No coordinates exceeded the threshold; no list file created.\n")
}


#### Extract the pixels in which pool temperature exceeds Ctmax and create the final raster ####

library(sf)
library(dplyr)
library(readr)
library(terra)

data_present <- read.csv("pool_depth_10cm_17-41days_mar-nov")
coords_exceed <- read.csv("coords_above_41.csv")

# Create unique key to match coordinates
data_present$key <- paste(data_present$lat, data_present$lon)
coords_exceed$key <- paste(coords_exceed$lat, coords_exceed$lon)

# Copy dataset to modify
data_modified <- data_present

# Replace pool value with 0 where coordinates exceed the threshold
data_modified$valid_site[data_modified$key %in% coords_exceed$key] <- 0

# Remove auxiliary column
data_modified$key <- NULL

# Count pools = 1 before and after modification
count_original <- sum(data_present$valid_site == 1, na.rm = TRUE)
cat("Pools = 1 in original dataset:", count_original, "\n")

count_modified <- sum(data_modified$valid_site == 1, na.rm = TRUE)
cat("Pools = 1 in modified dataset:", count_modified, "\n")

# Save final dataset
setwd("/Volumes/Seagate Por/Frontiers/binarios/17-41/rasters menos CTmax")
write.csv(data_modified, "final_pool_layer_present.csv", row.names = FALSE)

# Load final CSV to build raster
df <- read_csv("final_pool_layer_present.csv")

# Clean data
df_clean <- df %>%
  mutate(
    lon = round(lon, 6),
    lat = round(lat, 6)
  ) %>%
  distinct(lon, lat, .keep_all = TRUE)

# Convert cleaned data to SpatVector
vect_points <- vect(df_clean, geom = c("lon", "lat"), crs = "EPSG:4326")

# Extract unique coordinates
lons <- sort(unique(df_clean$lon))
lats <- sort(unique(df_clean$lat))

# Fixed buffer (2.5 arcmin)
buf <- 0.041665

# Compute resolution
dx <- mean(diff(lons))
dy <- mean(diff(lats))

# Create empty raster
template_raster <- rast(
  xmin = min(lons) - buf,
  xmax = max(lons) + buf,
  ymin = min(lats) - buf,
  ymax = max(lats) + buf,
  resolution = c(dx, dy),
  crs = "EPSG:4326"
)

# Rasterize using the "valid_site" column
final_raster <- rasterize(vect_points, template_raster, field = "valid_site")

# Plot raster
plot(final_raster)

# Save GeoTIFF
writeRaster(final_raster, "final_pool_layer_present.tif", overwrite = TRUE)


#### Obtain the offset values for pool models in future scenarios ####
#This process implies downloading current and future layers of temperature (we used WorldClim)

library(terra)

# Merge present-day raster layers into a multiband .tif file
# Get the list of raster files in the directory
files_t <- list.files("your.file.path/wc2.1_5m_tmax",
                      pattern = "\\.tif$", full.names = TRUE)

# Load rasters into a single multiband object
tmax_present <- rast(files_t)

# Save the multiband raster as a new file
writeRaster(tmax_present, "wc2.1_5m_tmax.tif", overwrite = TRUE)

# Load the points file and convert it to a spatial object
points_file <- read.csv(file.choose())
points <- vect(points_file, geom = c("lon", "lat"), crs = "EPSG:4326")

# Extract temperature and precipitation values at each point for the present
tmax_present  <- extract(tmax_present,  points)[, -1]

# Calculate differences between future and present temperatures 

# IPSL.CM6A.LR 

t_ipsl_3_2050 <- rast("wc2.1_5m_tmax_IPSL-CM6A-LR_ssp370_2041-2060.tif")
t_ipsl_3_2070 <- rast("wc2.1_5m_tmax_IPSL-CM6A-LR_ssp370_2061-2080.tif")

# Extract precipitation and temperature values for each future scenario
t_ipsl_3_2050.pol <- extract(t_ipsl_3_2050, points)[, -1]
t_ipsl_3_2070.pol <- extract(t_ipsl_3_2070, points)[, -1]

# Compute the mean differences (average of the 12 months)
offset.t_ipsl_3_2050 <- rowMeans(t_ipsl_3_2050.pol - tmax_present,  na.rm = TRUE)
offset.t_ipsl_3_2070 <- rowMeans(t_ipsl_3_2070.pol - tmax_present,  na.rm = TRUE)

# Create and save the datasets with offsets
offsets.ipsl_3_2050 <- points_file
offsets.ipsl_3_2050$temp  <- offset.t_ipsl_3_2050
write.csv(offsets.ipsl_3_2050, "offsets.ipsl_3_2050.csv", row.names = FALSE)

offsets.ipsl_3_2070 <- points_file
offsets.ipsl_3_2070$temp  <- offset.t_ipsl_3_2070
write.csv(offsets.ipsl_3_2070, "offsets.ipsl_3_2070.csv", row.names = FALSE)

# GFDL.ESM4 

t_gfdl_3_2050 <- rast("wc2.1_5m_tmax_GFDL-ESM4_ssp370_2041-2060.tif")
t_gfdl_3_2070 <- rast("wc2.1_5m_tmax_GFDL-ESM4_ssp370_2061-2080.tif")

# Extract values for future scenarios
t_gfdl_3_2050.pol <- extract(t_gfdl_3_2050, points)[, -1]
t_gfdl_3_2070.pol <- extract(t_gfdl_3_2070, points)[, -1]

# Compute mean differences
offset.t_gfdl_3_2050 <- rowMeans(t_gfdl_3_2050.pol - tmax_present,  na.rm = TRUE)
offset.t_gfdl_3_2070 <- rowMeans(t_gfdl_3_2070.pol - tmax_present,  na.rm = TRUE)

# Create and save outputs
setwd("/Users/leecar_1/Documents/Ricardo/Frontiers/Futuro/offsets/GFDL")

offsets.gfdl_3_2050 <- points_file
offsets.gfdl_3_2050$temp  <- offset.t_gfdl_3_2050
write.csv(offsets.gfdl_3_2050, "offsets.gfdl_3_2050.csv", row.names = FALSE)

offsets.gfdl_3_2070 <- points_file
offsets.gfdl_3_2070$temp  <- offset.t_gfdl_3_2070
write.csv(offsets.gfdl_3_2070, "offsets.gfdl_3_2070.csv", row.names = FALSE)


#### Pool modelation in future scenarios ####
# For this procedure we ran the code exactly the same way as detailed above, except for one change:
# In the "Pool modelation with NicheMapR" section, we added the "warm" parameter to the micro_global() function.
# The "warm" value corresponds to the "temp" value from the offset files. 
# For the sake of brevity, we only put here the modified function and loop to run the microclimate model

# Function to run a single point
run_point <- function(lat, lon, warm,
                      PE, BB, BD, KS, rainmult, LAI, cap, soilgrids, timeinterval, Usrhyt) {
  
  micro <- micro_global(loc = c(lon, lat), runmoist = 1, rainfrac = 0.1, Usrhyt = Usrhyt,
                        timeinterval = timeinterval, PE = PE, BB = BB, BD = BD, KS = KS,
                        rainmult = rainmult, LAI = LAI, cap = cap, soilgrids = soilgrids,
                        warm = warm)
  
  return(micro$metout)
}

#### Run model loop with retry and error log
start_time <- Sys.time()
errors <- list()   # errors will be stored here

for (i in 1:nrow(points)) {
  lat <- points$lat[i]
  lon <- points$lon[i]
  warm <- points$temp[i]   # <-- use 'temp' column from the CSV
  
  metout_csv <- file.path("metout", paste0(lat, ", ", lon, "_metout.csv"))
  
  # Only check if the metout file exists
  if (file.exists(metout_csv)) {
    cat("Existing file, skipping:", lat, lon,
        "(", i, "of", nrow(points), ")\n")
    next
  }
  
  success <- FALSE
  attempts <- 0
  while (!success && attempts < 2) {
    attempts <- attempts + 1
    cat("Processing point", lat, lon, "- attempt", attempts,
        "(", i, "of", nrow(points), ") - warm =", warm, "\n")
    
    res <- tryCatch({
      metout <- run_point(lat = lat, lon = lon, warm = warm,
                          PE = PE, BB = BB, BD = BD, KS = KS,
                          rainmult = rainmult, LAI = LAI, cap = cap,
                          soilgrids = soilgrids, timeinterval = timeinterval,
                          Usrhyt = Usrhyt)
      
      write.csv(metout, metout_csv, row.names = FALSE)
      "ok"
      
    }, error = function(e) {
      cat("Error at point", lat, lon, ":", conditionMessage(e), "\n")
      "fail"
    })
    
    if (identical(res, "ok")) {
      success <- TRUE
      cat("Point successfully processed:", lat, lon,
          "(", i, "of", nrow(points), ") - warm =", warm, "\n")
    }
  }
  
  if (!success) {
    cat("Failed twice, skipping point:", lat, lon,
        "(", i, "of", nrow(points), ")\n")
    errors[[length(errors) + 1]] <- paste(lat, lon)
  }
}

# Save error log if needed
if (length(errors) > 0) {
  writeLines(unlist(errors), "error_log.txt")
  
  
  cat("A log with", length(errors), "errors was saved in 'error_log.txt'\n")
} else {
  cat("No errors, no log created\n")
}

