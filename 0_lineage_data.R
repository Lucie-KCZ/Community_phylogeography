#############################################
###     Cleaning up the community data    ###
#############################################

# Authors:
# Lucie Kuczynski, Katharine Marske

# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Last edit: Oct 2024

# The script reads genetic lineage files, processes them by converting the data into spatial polygons,
# and matches them with corresponding community data to clean and prepare for further analyses.

# Clear memory and reset garbage collection
gc(reset = TRUE) ; rm(list = ls())

# Set the path to the directory containing genetic lineage files
path_lineages <- '/Users/lucie/Documents/Work/Data/Data_Norman/Genetic data/Genetic final data sheets 29 May 2018'

# List all the files in the directory, excluding complex data files and non-CSV files
files_lineages <- list.files(path = list.dirs(path = path_lineages), full.names = TRUE)
files_lineages <- files_lineages[grep('.csv', files_lineages)]
files_lineages <- files_lineages[-grep('_complex_', files_lineages)]

# Function to format the lineage data by filtering necessary columns and removing duplicates
format_lineage <- function(lineage) {
  
  # Read the lineage CSV file into a dataframe
  # Assumes 'latin1' encoding, no factors, and ',' as a separator
  lineage <- read.csv(file = lineage, header = TRUE, 
                      sep = ',', fileEncoding = 'latin1', 
                      stringsAsFactors = FALSE)
  
  # Keep only the relevant columns: 'Species', 'Genus', 'lineage', 'lat', 'long'
  # Remove duplicate rows (i.e., unique combinations of these columns)
  lineage <- unique(lineage[, colnames(lineage) %in% c('Species', 'Genus', 'lineage', 'lat', 'long')])
  
  # Remove any lineage group that has fewer than 4 records
  lineage <- lineage[-which(lineage$lineage %in% names(which(table(lineage$lineage) < 4))), ]
  
  # Check if there's more than one unique lineage remaining
  # If not, return NULL (as further processing is not meaningful with fewer lineages)
  if(length(unique(lineage$lineage)) < 2) {
    return(NULL)  # No valid lineages, stop processing
  } else {
    return(lineage)  # Return the formatted lineage data
  }
}

# Function to convert the lineage data into a polygon shapefile format
to_polygon_lineage <- function(lineage){
  
  # Create a list of polygons for each lineage by calculating the convex hull
  # The 'by' function splits the data by lineage and calculates convex hull coordinates
  lineage_shp <- as.list(by(
    data = lineage, 
    INDICES = lineage$lineage, 
    FUN = function(x) x[c(chull(x[, c('lat', 'long')]), chull(x[, c('lat', 'long')])[1]), c('lat', 'long')]
  ))
  
  # Combine the list into a single data frame
  lineage_shp <- do.call(rbind, lineage_shp)
  
  # Add a column for the lineage names extracted from the row names (split by '.')
  lineage_shp <- data.frame(
    lineage = unlist(lapply(strsplit(rownames(lineage_shp), '.', fixed = TRUE), function(x) x[1])), 
    lineage_shp
  )
  
  # Remove the row names for cleaner data
  rownames(lineage_shp) <- NULL
  
  # Convert the lineage data frame into an 'sf' polygon object using latitude/longitude coordinates
  lineage_shp <- sfheaders::sf_polygon(lineage_shp, x = 'long', y = 'lat', polygon_id = "lineage")
  
  # Set the coordinate reference system (CRS) to 'WGS84' (standard for geographic coordinates)
  sf::st_set_crs(lineage_shp, 'WGS84')
  
  # Return the final polygon shapefile
  return(lineage_shp)
}

# Function to import the appropriate community data based on the species' genus
import_right_com_data <- function(species_name){
  
  # Import the taxonomic order for each genus
  taxonomic_order <- read.csv('/Users/lucie/Documents/Work/Data/Data_Norman/Order.csv')
  
  # Extract the genus from the species name (assumes genus is the first part before '_')
  genus <- unlist(strsplit(species_name, '_'))[1]
  
  # Find the corresponding order for the genus
  order <- taxonomic_order[taxonomic_order$Genus == genus, 'Order']
  
  # If the order is found, import the community data; otherwise, print a message and return NULL
  if(length(order) != 0){
    com_data <- read.csv(paste0('../data/processed/America_', order, '_site.csv'))[, -1]
    return(com_data)
  } else {
    print(paste('Failed to import', species_name, ': Species not listed in the corresponding Orders file.'))
    return(NULL)  # Stop further processing if no data is found
  }
}

# Function to match the lineage polygons with the community data
match_lineage_to_community <- function(lineage_shp, lineage_species, community_data) {
  
  # Load the taxonomic order for each genus
  taxonomic_order <- read.csv('/Users/lucie/Documents/Work/Data/Data_Norman/Order.csv')
  
  # Extract the genus from the lineage species name
  genus <- unlist(strsplit(lineage_species, '_'))[1]
  
  # Identify the corresponding taxonomic order for the genus
  order <- taxonomic_order[taxonomic_order$Genus == genus, 'Order']
  
  # Convert community data to spatial points using latitude and longitude
  community_points_sf <- sf::st_as_sf(community_data, coords = c("longitude", "latitude"))
  
  # Perform a spatial join between community points and lineage polygons, keeping only 'id' and 'lineage'
  community_lineage <- unique(as.data.frame(sf::st_join(community_points_sf, lineage_shp))[, c('id', 'lineage')])
  rm(community_points_sf)  # Clean up intermediate variable
  
  # Rename the lineage column to the specific lineage_species
  colnames(community_lineage)[2] <- lineage_species
  
  # Path to save the matched data file
  match_file_path <- paste0('../data/processed/America_', order, '_lineages.csv')
  
  # Append data if the file exists, otherwise create a new file
  if (file.exists(match_file_path)) {
    match_data <- read.csv(match_file_path)[, -1]
    match_data <- unique(merge(x = match_data, y = community_lineage, by = 'id', all = TRUE))
    write.csv(match_data, file = match_file_path)
    print(paste0('Data appended to the file for ', lineage_species, '.'))
  } else {
    write.csv(community_lineage, file = match_file_path)
    print(paste0('New file created for ', lineage_species, ' and data saved.'))
  }
}

# Main function to process lineage files and match them with community data
process_lineage <- function(lineage) {
  
  gc(reset = TRUE)
  
  # Extract the species name from the lineage file path
  lineage_species <- unlist(strsplit(lineage, '/'))
  lineage_species <- lineage_species[length(lineage_species)]
  lineage_species <- gsub('.csv', '', lineage_species)
  
  # Check progress
  print(paste0(lineage_species, ' is being processed.'))
  
  # Format the lineage data
  lineage <- format_lineage(lineage)
  
  if(is.null(lineage)){
    # Skip species with only one documented lineage
    return(NULL)
    warning(paste0('Only one lineage detected for ', lineage_species, '. The species will be skipped.'))
  } else {
    # Convert the lineage data to a polygon shapefile
    lineage <- to_polygon_lineage(lineage)
    
    # Import relevant community data for the species
    com_data <- import_right_com_data(lineage_species)
    
    if(!is.null(com_data)){
      # Match the processed lineage with community data and save the results
      match_lineage_to_community(lineage_shp = lineage, lineage_species, community_data = com_data)
    } else {
      # Display a warning if no community data is found
      warning(paste0("No community data found for ", lineage_species, ". The species will be skipped."))
    }
  }
}

# Run the process on all lineage files
sapply(files_lineages, process_lineage)
