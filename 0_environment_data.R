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
library(dplyr)         # For data manipulation

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

# Select and rename relevant columns using dplyr
cols_mapping <- c('Code' = 'Code', 
                  'Lat' = 'Lat', 
                  'Long' = 'Long', 
                  'Bio_1' = 'Temp', 
                  'Bio_4' = 'TSeas', 
                  'Bio_12' = 'Rain', 
                  'Bio_15' = 'RSeas', 
                  'elevation' = 'Ele', 
                  'vegetation' = 'Veg', 
                  'solar' = 'Solar', 
                  'wind' = 'Wind', 
                  'landcover' = 'LC')

sites <- sites %>%
  select(all_of(names(cols_mapping))) %>%
  rename(all_of(cols_mapping))

# Recode factor levels for 'LC'
new_levels <- c('Broad_ever', 'Broad_deciduous', 'Needle_ever', 
                'Needle_deciduous', 'Mixed_forest', 'Tree_open', 
                'Shrudb', 'Herbaceous', 'Sparse_veg', 
                'Crop', 'Paddy_field', 'Veg_mosaic', 
                'Wetland', 'Urban', 'Water_bodies')

sites <- sites %>%
  mutate(
    LC = factor(LC, levels = new_levels))

# Post-processing: Aggregate environmental variables
sites <- sites %>%
  mutate(
    elevation = rowSums(select(., starts_with("Ele_")), na.rm = TRUE),
    vegetation = rowSums(select(., starts_with("Veg_")), na.rm = TRUE),
    solar = rowMeans(select(., starts_with("Solar_")), na.rm = TRUE),
    wind = rowMeans(select(., starts_with("Wind_")), na.rm = TRUE),
    landcover = rowMeans(select(., starts_with("LC_")), na.rm = TRUE)) %>%
  mutate(
    elevation = ifelse(elevation < -4000, NA, elevation),
    vegetation = ifelse(vegetation > 200, NA, vegetation))

# View the first few rows of the final data frame
head(sites)
