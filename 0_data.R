#############################################
###     Cleaning up the community data    ###
#############################################

# Authors:
# Lucie Kuczynski, Katharine Marske

# L. Kuczynski
# lucie.kuczynski@hotmail.com
# Last edit: Oct 2024

# Clear memory and reset garbage collection
gc(reset = TRUE) ; rm(list = ls())

# Dataset1: all versions - sites ------------------------------------------
# Load sites data
sites <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/all versions/sites CSV/site description table RSL 24-11-2014.csv', 
  sep = ';', dec = ','
)

# Select relevant columns
sites <- sites[, c(1, 9, 10, 13, 23, 26, 30, 34, 39, 44, 49, 54, 59)]

# Convert 'Site_area_.km.2.' to numeric
# Replace commas with dots to match numeric format
sites$Site_area_.km.2. <- gsub(',', '.', sites$Site_area_.km.2., fixed = TRUE)

# Remove '>' and '<' symbols for areas that are reported as greater or less than values
sites$Site_area_.km.2. <- gsub('>', '', sites$Site_area_.km.2., fixed = TRUE)
sites$Site_area_.km.2. <- gsub('<', '', sites$Site_area_.km.2., fixed = TRUE)

# Calculate the mean for areas defined as ranges (e.g., "xx-xx")
sites$Site_area_.km.2.[grepl("^[^a-zA-Z]*-[^a-zA-Z]*$", sites$Site_area_.km.2.)] <- 
  unlist(lapply(strsplit(sites$Site_area_.km.2.[grepl("^[^a-zA-Z]*-[^a-zA-Z]*$", sites$Site_area_.km.2.)], '-'), 
                function(x) mean(as.numeric(x))))

sites$Site_area_.km.2.[grepl("^[^a-zA-Z]*–[^a-zA-Z]*$", sites$Site_area_.km.2.)] <- 
  unlist(lapply(strsplit(sites$Site_area_.km.2.[grepl("^[^a-zA-Z]*–[^a-zA-Z]*$", sites$Site_area_.km.2.)], '-'), 
                function(x) mean(as.numeric(x))))

# Set area values with alphabetic characters (except scientific notation) to NA
sites$Site_area_.km.2.[grepl("[a-zA-Z]", sites$Site_area_.km.2.) & !grepl("E-", sites$Site_area_.km.2.)] <- NA

# Convert the cleaned 'Site_area_.km.2.' column to numeric
sites$Site_area_.km.2. <- as.numeric(sites$Site_area_.km.2.)

# Extracting the first year of sampling from 'Year_sampling'
# Step 1: Split 'Year_of_survey' by '+' (for multiple years in the same entry)
# Step 2: Further split by '-' (for year ranges) and extract the minimum year
first_sampling_year <- unlist(
  lapply(lapply(strsplit(sites$Year_sampling, '+', fixed = TRUE), function(x) 
    unlist(strsplit(x, '-', fixed = TRUE))), 
    function(y) data.frame(min(y)))
)

# Remove names from the 'first_sampling_year' vector
names(first_sampling_year) <- NULL

# Manual corrections for certain years
manual_corrections <- data.frame(
  incorrect_position = c(27, 33, 34, 37, 39, 40, 41), 
  corrected_year = c('2002', '1998', '2004', '2008', '1995', '2011', '2010')
)

# Get unique years from the first sampling year vector
unique_years <- unique(first_sampling_year)

# Replace incorrect years with corrected values
for (i in 1:nrow(manual_corrections)) {
  first_sampling_year[which(first_sampling_year == unique_years[manual_corrections[i, 1]])] <- manual_corrections[i, 2]
}

# Clean up variables and convert 'first_sampling_year' to numeric
rm(i, manual_corrections, unique_years)
first_sampling_year <- as.numeric(first_sampling_year)

# Add the cleaned 'first_sampling_year' to the site data
sites$first_sampling_year <- first_sampling_year
sites$Year_sampling <- NULL ; rm(first_sampling_year)

# Clean trap method binary variables (e.g., '?', '(yes)', and empty entries)
# Columns to clean and replace '?' or empty values with NA
cols_to_clean <- c('Terrestrial_trap', 'Aquatic_trap', 'Call_survey', 
                   'Terrestrial_survey', 'Aquatic_survey..Dip.net.surveys.', 
                   'Day_Visual_surveys', 'Night_Visual_surveys', 'Dip_net_surveys')

# Replace '?' or empty values with NA in selected columns
sites[cols_to_clean] <- lapply(sites[cols_to_clean], function(x) {
  x[x %in% c('?', '')] <- NA
  return(x)
}) ; rm(cols_to_clean)

# Special case for 'Terrestrial_survey' where '(yes)' needs to be replaced with 'yes'
sites$Terrestrial_survey[sites$Terrestrial_survey == '(yes)'] <- 'yes'

# Convert 'yes' values to TRUE, others to FALSE for binary trap method variables
for (i in 5:12) sites[, i] <- ifelse(sites[, i] == 'yes', TRUE, FALSE)

# Calculate the number of methods used for each site
sites$num_sampling_methods <- apply(sites[, 5:12], 1, function(x) sum(x, na.rm = TRUE))

# Updates names for clarity
colnames(sites) <- c('id', 'latitude', 'longitude', 'area', 
                     'terrestrial_trap', 'aquatic_trap', 'call_survey', 
                     'terrestrial_survey', 'aquatic_survey', 'day_visual', 
                     'night_visual', 'dip_net_survey', 
                     'year', 'number_methods')

# Dataset1: all versions - community --------------------------------------
# Load the species occurrence table
community <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/all versions/species CSV/species occurrence table RSL 18-11-2014.csv', 
  sep = ';', dec = ','
)

# Select relevant columns and ensure uniqueness
community <- unique(community[, c(1, 13, 42)])

# Extracting the first year of sampling from 'Year_of_survey'
first_sampling_year <- unlist(
  lapply(lapply(strsplit(community$Year_of_survey, '+', fixed = TRUE), function(x) 
    unlist(strsplit(x, '-', fixed = TRUE))), 
    function(y) data.frame(min(y)))
)

# Remove any names from the 'first_sampling_year' vector for a cleaner result
names(first_sampling_year) <- NULL

# Manual corrections for specific sampling years (errors or outliers in the data)
# 'incorrect_position' represents the index of the incorrect year, and 'corrected_year' is the valid value
manual_year_corrections <- data.frame(
  incorrect_position = c(14, 15, 30, 36, 40, 50, 51, 52, 53), 
  corrected_year = c('2008', '1995', '2004', '1997', '2011', '1994', '1995', '1996', '1996')
)

# Get unique years from the first sampling year vector
unique_years <- unique(first_sampling_year)

# Replace incorrect years with the corrected values from 'manual_year_corrections'
for (i in 1:nrow(manual_year_corrections)) {
  first_sampling_year[which(first_sampling_year == unique_years[manual_year_corrections[i, 1]])] <- manual_year_corrections[i, 2]
}

# Clean up unnecessary variables and convert 'first_sampling_year' to numeric for further analysis
rm(i, manual_year_corrections, unique_years)
first_sampling_year <- as.numeric(first_sampling_year)

# Add the cleaned 'first_sampling_year' to the species data
community$Year_of_survey <- first_sampling_year ; rm(first_sampling_year)

# Remove species where identification to the species level wasn't properly assessed/documented
# The 'Manually_converte_name' column contains species names, 
# and certain entries need to be removed due to unclear identification
species_to_remove <- 
  sort(unique(community$Manually_converte_name))[c(
    1, 11, 15, 37, 41, 85, 86, 89, 117, 122, 125, 133, 154, 155, 163, 174, 176, 196, 197)]

# Remove species with unclear identification 
# (this affects 438 occurrences out of over 33,700, spread across 9,100+ locations)
community <- community[!(community$Manually_converte_name %in% species_to_remove), ]

# Clean up unnecessary variables
rm(species_to_remove)

# Updates names for clarity
colnames(community) <- c('id', 'year', 'species')

# Dataset2: biosurvey -------------------------------------------------------------------
# Load amphibian biosurvey data
biosurvey_data <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Biosurvey/amphibians.csv', 
  sep = ',', dec = '.'
)

# Keep unique combinations of relevant columns: species name, year of survey, latitude, and longitude
# Columns: species (sname), latitude (decimallatitude), longitude (decimallongitude), year (eventdate)
biosurvey_data <- unique(biosurvey_data[, c(6, 22, 23, 2)])

# Extract year from 'eventdate' (format: 'dd/mm/yyyy'), keeping only the year portion
biosurvey_data$event_year <- as.numeric(
  unlist(lapply(strsplit(biosurvey_data$eventdate, '/'), function(x) x[3])))

# Standardize species names by removing subspecies information (keep genus and species, separated by an underscore)
biosurvey_data$species_name <- unlist(
  lapply(strsplit(biosurvey_data$sname, ' ', fixed = TRUE), function(x) paste(x[1], x[2], sep = '_')))

# Remove rows with missing latitude, longitude, or year information (1122 occurrences)
biosurvey_data <- na.omit(biosurvey_data)

# Create a unique site ID for each unique pair of latitude and longitude
site_id <- as.factor(paste(biosurvey_data$decimallatitude, biosurvey_data$decimallongitude))
levels(site_id) <- paste0('BIOSURVEY', 1:length(levels(site_id)))
biosurvey_data$id <- as.character(site_id) ; rm(site_id)

# Prepare site-level data frame with relevant columns, including placeholders for area and survey methods
sites_biosurvey <- data.frame(
  id = biosurvey_data$id, 
  latitude = biosurvey_data$decimallatitude, 
  longitude = biosurvey_data$decimallongitude, 
  area = NA,  # Placeholder for site area (if available later)
  terrestrial_trap = NA,  # Placeholder for terrestrial trap data
  aquatic_trap = NA,  # Placeholder for aquatic trap data
  call_survey = NA,  # Placeholder for call survey data
  terrestrial_survey = NA,  # Placeholder for terrestrial survey data
  aquatic_survey = NA,  # Placeholder for aquatic survey data
  day_visual = NA,  # Placeholder for daytime visual survey data
  night_visual = NA,  # Placeholder for nighttime visual survey data
  dip_net_survey = NA,  # Placeholder for dip net survey data
  year = biosurvey_data$event_year,  # Year of the survey
  number_methods = NA  # Placeholder for the number of sampling methods used
)

# Prepare species occurrence data frame with site ID, year, and species name
community_biosurvey <- data.frame(
  id = biosurvey_data$id, 
  year = biosurvey_data$event_year, 
  species = biosurvey_data$species_name
)

# Adding the data to initial dataset
sites <- unique(rbind(sites, sites_biosurvey))
community <- unique(rbind(community, community_biosurvey))

# Cleaning up the environment
rm(sites_biosurvey, community_biosurvey, biosurvey_data)

# Dataset3: GBIF -------------------------------------------------------------------
# Load amphibian GBIF data
gbif_data <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/GBIF/0029645-240906103802322.csv', 
  sep = '\t', dec = '.'
)

# Select unique rows with relevant columns (species name, latitude, longitude, occurrence status, and year)
gbif_data <- unique(gbif_data[, c(10, 19, 22, 23, 33)])

# Keep only records where species is present
gbif_data <- gbif_data[gbif_data$occurrenceStatus == 'PRESENT', ]
gbif_data$occurrenceStatus <- NULL

# Clean species names by standardizing to genus and species, separated by an underscore
gbif_data <- gbif_data[-which(gbif_data$species %in% c('', '5')), ]
gbif_data$species <- unlist(
  lapply(strsplit(gbif_data$species, ' ', fixed = TRUE), function(x) paste(x[1], x[2], sep = '_')))
gbif_data <- unique(gbif_data)

# Create unique site IDs based on latitude and longitude
site_id <- as.factor(paste(gbif_data$decimalLatitude, gbif_data$decimalLongitude))
levels(site_id) <- paste0('GBIF', 1:length(levels(site_id)))
gbif_data$site_id <- as.character(site_id)
rm(site_id)

# Prepare site-level data frame
sites_gbif <- data.frame(
  id = gbif_data$site_id,
  latitude = gbif_data$decimalLatitude,
  longitude = gbif_data$decimalLongitude,
  area = NA,
  terrestrial_trap = NA,
  aquatic_trap = NA,
  call_survey = NA,
  terrestrial_survey = NA,
  aquatic_survey = NA,
  day_visual = NA,
  night_visual = NA,
  dip_net_survey = NA,
  year = gbif_data$year,
  number_methods = NA
)

# Prepare species occurrence data frame
community_gbif <- data.frame(
  id = gbif_data$site_id, 
  year = gbif_data$year, 
  species = gbif_data$species
)

# Merge GBIF data with initial dataset
sites <- unique(rbind(sites, sites_gbif))
community <- unique(rbind(community, community_gbif))

# Clean up environment
rm(sites_gbif, community_gbif, gbif_data)

# Dataset4: Grinnell2005 -------------------------------------------------------------------
# List all files in the Grinnell directory
grinnell_files <- list.files('/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Grinnell', full.names = TRUE)

# Filter for downloaded files
grinnell_files <- grinnell_files[grep('download', grinnell_files)]

# Initialize an empty data frame to store Grinnell data
grinnell_data <- NULL

# Loop through each Grinnell file and load its content
for (i in grinnell_files) {
  # Read the current file (assuming the same path and structure for all)
  grinnell_temp <- read.table('/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Grinnell/bm2_download.txt', sep = '\t', header = TRUE)
  
  # Keep only the unique rows with relevant columns (Scientific Name, Latitude, Longitude)
  grinnell_temp <- unique(grinnell_temp[, c(2, 6, 7)])
  
  # Append the data to the main data frame
  grinnell_data <- rbind(grinnell_data, grinnell_temp)
  
  # Clean up the temporary variable after each iteration
  rm(grinnell_temp)
}

# Clean up loop control variables
rm(i, grinnell_files)

# Replace spaces in scientific names with underscores
grinnell_data$Scientific.Name <- gsub(' ', '_', grinnell_data$Scientific.Name, fixed = TRUE)

# Create unique site IDs based on Latitude and Longitude
site_id <- as.factor(paste(grinnell_data$Latitude, grinnell_data$Longitude))

# Generate unique IDs for each site by concatenating 'GRIN' and an incremental number
levels(site_id) <- paste0('GRIN', 1:length(levels(site_id)))

# Add the site IDs to the Grinnell data
grinnell_data$id <- as.character(site_id)

# Clean up the temporary site_id variable
rm(site_id)

# Prepare a site-level data frame with placeholders for area and sampling methods
sites_grinnell <- data.frame(
  id = grinnell_data$id,
  latitude = grinnell_data$Latitude,
  longitude = grinnell_data$Longitude,
  area = NA,
  terrestrial_trap = NA,
  aquatic_trap = NA,
  call_survey = NA,
  terrestrial_survey = NA,
  aquatic_survey = NA,
  day_visual = NA,
  night_visual = NA,
  dip_net_survey = NA,
  year = 2005,
  number_methods = NA
)

# Prepare a species occurrence data frame with the site ID, year, and species name
community_grinnell <- data.frame(
  id = grinnell_data$id,
  year = 2005,
  species = grinnell_data$Scientific.Name
)

# Combine Grinnell site data with the initial dataset of sites
sites <- unique(rbind(sites, sites_grinnell))

# Combine Grinnell species data with the initial community dataset
community <- unique(rbind(community, community_grinnell))

# Clean up the environment by removing temporary variables
rm(sites_grinnell, community_grinnell, grinnell_data)


# Dataset5: UGSG -------------------------------------------------------------------
# Load UGSG coordinate data
ugsg_coords <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/UGSG/Coordinates.csv', 
  sep = ',', dec = '.'
)

# Load UGSG species count data
ugsg_counts <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/UGSG/Counts.csv', 
  sep = ',', dec = '.'
)

# Load UGSG survey run data
ugsg_run <- read.csv(
  '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/UGSG/Runs.csv', 
  sep = ',', dec = '.'
)

# Merge counts and run data based on 'RunID', keeping all records from counts
ugsg_data <- merge(x = ugsg_counts, y = ugsg_run, by = 'RunID', all.x = TRUE, all.y = FALSE)

# Merge the resulting data with coordinate data based on 'RouteNumber'
ugsg_data <- merge(x = ugsg_data, y = ugsg_coords, by = 'RouteNumber', all.x = TRUE, all.y = FALSE)

# Remove unnecessary data frames to free up memory
rm(ugsg_coords, ugsg_counts, ugsg_run)

# Retain only unique records and select relevant columns (e.g., lat/lon, species)
ugsg_data <- unique(ugsg_data[, c(4, 9, 26, 27, 28)])

# Remove records where the species is labeled as part of a species complex
ugsg_data <- ugsg_data[-grep('complex', ugsg_data$Species), ]

# Format species names by replacing spaces with underscores
ugsg_data$Species <- gsub(' ', '_', ugsg_data$Species, fixed = TRUE)

# Format the site ID by adding the 'UGSG' prefix
ugsg_data$SiteID <- paste0('UGSG', ugsg_data$SiteID)

# Create a site-level data frame with placeholders for area and sampling methods
sites_ugsg <- data.frame(
  id = ugsg_data$SiteID,
  latitude = ugsg_data$lat,
  longitude = ugsg_data$lon,
  area = NA,
  terrestrial_trap = NA,
  aquatic_trap = NA,
  call_survey = NA,
  terrestrial_survey = NA,
  aquatic_survey = NA,
  day_visual = NA,
  night_visual = NA,
  dip_net_survey = NA,
  year = ugsg_data$SurveyYear,
  number_methods = NA
)

# Create a species occurrence data frame with site ID, survey year, and species name
community_ugsg <- data.frame(
  id = ugsg_data$SiteID,
  year = ugsg_data$SurveyYear,
  species = ugsg_data$Species
)

# Combine UGSG site data with the initial dataset of sites
sites <- unique(rbind(sites, sites_ugsg))

# Combine UGSG species occurrence data with the initial community dataset
community <- unique(rbind(community, community_ugsg))

# Clean up the environment by removing temporary variables
rm(sites_ugsg, community_ugsg, ugsg_data)


# Save outputs ------------------------------------------------------------
# Save sites and community data to CSV
write.csv(sites, file = '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Processed data/all_sites.csv')
write.csv(community, file = '/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Processed data/occurrences.csv')

# Import the order for each genus
taxonomic_order <- read.csv('/Users/lucie/Documents/Work/Data/Data_Norman/Order.csv')

# Function to subset site and community data based on continent and taxonomic order
subset_continent_taxa <- function(site_data, com_data, continent = 'America', order = 'Anura', write = TRUE) {
  
  # Filter data for American or non-American sites based on longitude
  if(continent == 'America') {
    site_data <- site_data[site_data$longitude < -50, ]  # For American sites (longitude > -50)
    com_data <- com_data[com_data$id %in% site_data$id, ] # Subset community data to match filtered site IDs
  } else {
    site_data <- site_data[site_data$longitude > -50, ]  # For non-American sites (longitude < -50)
    com_data <- com_data[com_data$id %in% site_data$id, ] # Subset community data to match filtered site IDs
  }
  
  # Subset species based on the specified taxonomic order
  species <- taxonomic_order[taxonomic_order$Order == order, 1]  # Get species from the taxonomic order data matching the specified order
  
  # Filter community data based on species that belong to the specified taxonomic order
  com_data <- com_data[unlist(lapply(strsplit(com_data$species, '_'), function(x) x[1])) %in% species, ]
  
  # Subset site data to only include sites present in the filtered community data
  site_data <- site_data[site_data$id %in% com_data$id, ]
  
  # Identify sites with fewer than 5 documented species
  low_species_sites <- names(which(table(site_data$id) < 5))
  
  # Filter out these sites and corresponding community data
  site_data <- site_data[!site_data$id %in% low_species_sites, ]
  com_data <- com_data[!com_data$id %in% low_species_sites, ]
  
  # Remove intermediate variables to free memory
  rm(low_species_sites, species)
  
  # Combine site and community data into a list for output
  out <- list(site = site_data, community = com_data)
  
  # If the write argument is TRUE, save the filtered data as CSV files
  if(write) {
    write.csv(site_data, 
              file = paste0('/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Processed data/', 
                            continent, '_', order, '_site.csv'))  # Save site data
    write.csv(com_data, 
              file = paste0('/Users/lucie/Documents/Work/Data/Data_Norman/Community data/Processed data/', 
                            continent, '_', order, '_community.csv'))  # Save community data
  } else {
    return(out)  # If write is FALSE, return the filtered data as a list
  }
}

# Loop over each continent ('America' and 'Europe')
for(i in c('America', 'Europe')) {
  
  # For each continent, loop over each taxonomic order ('Anura' and 'Urodela')
  for(j in c('Anura', 'Urodela')) {
    
    # Call the subset_continent_taxa function with the current continent (i) and taxonomic order (j)
    # The 'sites' and 'community' datasets are passed as inputs to be filtered
    subset_continent_taxa(sites, community, i, j)
    
  }
}
