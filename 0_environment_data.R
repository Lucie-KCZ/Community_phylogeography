#########################################
###     Getting Environmental Data    ###
#########################################

# Authors:
# Lucie Kuczynski, Katharine Marske

# Last edit: Nov 2024

# Clear memory and reset garbage collection
rm(list = ls())
gc(reset = TRUE)

# Load required libraries
library(raster)       # For raster data manipulation
library(sf)           # For modern spatial data handling
library(progress)     # For progress bar

# Function to generate meaningful file names
generate_file_names <- function(files) {
  # Define a mapping from directory names to desired prefixes
  prefix_mapping <- list(
    "GM_elevation_v1"         = "Ele",
    "GM_landcover_v3"         = "LC",
    "GM_vegetation_v2"        = "Veg",
    "WC_bioclim_v2"           = "Bio",
    "WC_solarradiation_v2"    = "Solar",
    "WC_wind_v2"              = "Wind"
  )
  
  # Extract the immediate directory name for each file
  dir_names <- basename(dirname(files))
  
  # Map each directory name to its corresponding prefix
  prefixes <- unlist(prefix_mapping[dir_names])
  
  # Handle unrecognized directory names
  if (any(is.na(prefixes))) {
    warning("Some files have unrecognized categories and will be named as 'Unknown_<number>'.")
    prefixes[is.na(prefixes)] <- "Unknown"
  }
  
  # Assign incremental numbers within each prefix group
  prefix_counts <- ave(prefixes, prefixes, FUN = seq_along)
  
  # Combine prefix and number to create the new file names
  new_names <- paste0(prefixes, "_", prefix_counts)
  
  return(new_names)
}

# Function to process a single raster and extract values
process_single_raster <- function(file_path, desired_name, xy_sf) {
  tryCatch({
    # Print the name of the file being processed
    message("Processing file: ", desired_name)
    
    # Load the raster file
    raster_i <- raster(file_path)
    
    # Extract raster values at the specified coordinates
    extract_i <- raster::extract(raster_i, xy_sf)
    
    # Return the extracted values with the desired column name
    return(setNames(data.frame(extract_i), desired_name))
  }, error = function(e) {
    # Handle errors gracefully
    warning("Failed to process file: ", desired_name, "\nError message: ", e$message)
    # Return NA to keep the data frame alignment
    return(setNames(data.frame(extract_i = NA), desired_name))
  })
}

# Load site data
sites <- read.csv('../data/processed/America_Anura_site.csv')[, c(2:4, 14)]
# Alternatively, for Urodela species:
# sites <- read.csv('../data/processed/America_Urodela_site.csv')[, c(2:4, 14)]

# Convert site data to an sf object for spatial operations
# EPSG code 4326 corresponds to WGS84
sites_sf <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)

# List all .tif files in the directories
files <- list.files(
  path = '~/Documents/Work/Data/Environment/data',
  pattern = '\\.tif$',
  recursive = TRUE,
  full.names = TRUE)

# Generate meaningful file names
files_names <- generate_file_names(files)

# Initialize a progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent in :elapsed \n",
  total = length(files), clear = FALSE, width = 60)

# Initialize a list to store extracted data
extracted_data_list <- vector("list", length(files))

# Loop through each raster file and process it
for (i in seq_along(files)) {
  pb$tick()  # Update progress bar
  
  # Process the current raster file and extract data
  extracted_data_list[[i]] <- process_single_raster(
    file_path = files[i],
    desired_name = files_names[i],
    xy_sf = sites_sf
  )
}

# Combine all extracted data into a single data frame
extracted_data <- do.call(cbind, extracted_data_list)

# Combine extracted data with site identifiers
sites <- cbind(st_drop_geometry(sites_sf), extracted_data)

# Post-processing: Aggregate environmental variables
sites$elevation <- rowSums(sites[, grepl('Ele_', names(sites))], na.rm = TRUE)
sites$vegetation <- rowSums(sites[, grepl('Veg_', names(sites))], na.rm = TRUE)
sites$solar <- rowMeans(sites[, grepl('Solar_', names(sites))], na.rm = TRUE)
sites$wind <- rowMeans(sites[, grepl('Wind_', names(sites))], na.rm = TRUE)
sites$landcover <- rowMeans(sites[, grepl('LC_', names(sites))], na.rm = TRUE)

# Handle extreme values
sites$elevation[sites$elevation < -4000] <- NA
sites$vegetation[sites$vegetation > 200] <- NA

# Select specific columns from the 'sites' data frame that are relevant for analysis.
# Only columns matching the names listed in the vector are retained.
sites <- sites[, which(colnames(sites) %in% c('id', 'year', 'Bio_1', 'Bio_4', 
                                              'Bio_12', 'Bio_15', 'elevation', 
                                              'vegetation', 'solar', 'wind', 
                                              'landcover'))]

# Rename the selected columns in 'sites' to more concise, standardized names.
colnames(sites) <- c('id', 'year', 'Temp', 'TSeas', 'Rain', 'RSeas', 'Ele', 
                     'Veg', 'Solar', 'Wind', 'LC')

# Recode the levels of the 'LC' (land cover) column with descriptive category names.
# This makes the categories more interpretable for analysis and visualization.
levels(sites$LC) <- c('Broad_ever', 'Broad_deciduous', 'Needle_ever', 
                      'Needle_deciduous', 'Mixed_forest', 'Tree_open', 
                      'Shrub', 'Herbaceous', 'Sparse_veg', 'Crop', 
                      'Paddy_field', 'Veg_mosaic', 'Wetland', 'Urban', 
                      'Water_bodies')

# Read in a separate CSV file containing additional site information (coordinates).
# Only select columns 2 through 4 and column 14 for merging.
coords <- read.csv('../data/processed/America_Anura_site.csv')[, c(2:4, 14)]

# Merge the 'sites' and 'coords' data frames by 'id' and 'year' columns.
# Keep all rows from 'sites' (all.x = T) but only matching rows from 'coords' (all.y = F).
sites <- merge(x = sites, y = coords, by = c('id', 'year'), all.x = TRUE, all.y = FALSE)

