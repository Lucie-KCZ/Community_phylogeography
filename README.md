# Community Phylogeography

## Data Overview

### 1. Data Cleaning and Import

#### A. **Initial Database: "All Versions"**
This dataset was compiled by K. Marske and collaborators in 2014. Below is an outline of the key data elements retained for the project:

**Location Information**
- **Site ID**: Unique identifier for each sampling site.
- **Latitude and Longitude**: Geographic coordinates of the sampling sites.
- **Site Area (km²)**: The area of each site.
- **Sampling Methods**: Binary variables indicating whether specific traps were used (Terrestrial trap, Aquatic trap, Call survey, Terrestrial survey, Aquatic survey, Day visual surveys, Night visual surveys, Dip net surveys)
- **Number of Sampling Methods**: The total number of traps used at each site.
- **First Sampled Year**: The earliest year of sampling at each site.

**Community Data**
- **Site ID**: Corresponding site identification.
- **Year of Survey**: The year when the community was surveyed.
- **Species Name** (`Manually_converte_name`): Name of the species at each site.

Note: Some records were excluded from the community data where identification did not reach the species level. This exclusion accounts for 438 occurrences, out of a total of over 33,700 occurrences, spread across 9,100+ locations.

Note: All life stages have been kept, based on the assumption that the presence of eggs or tadpoles implies the presence of adults capable of reproduction.

Note: Additional databases will follow this structure, although some data might be missing.

#### B. **Additional Database: OK Biosurvey**
This dataset is from the OU biosurvey (to be checked) - Oklahoma National Heritage Inventory and cover, to some extent, amphibian occurrences in Oklahoma.

Note: Several occurrences were missing spatial/temporal informations and were thus removed (1122 occurrences).

#### C. **Additional Database: GBIF**
This dataset was obtained from [GBIF](GBIF.org) (accessed on 25 September 2024): https://doi.org/10.15468/dl.rqwsh3. The following filters were applied:
- Location: Europe or North America
- Records with spatial coordinates
- Species reported as present
- Taxonomic class: Amphibia
- Time range: 1950-2024

#### D. **Additional Datanase: Grinnell Resurvey Project 2005**

### 2. Final dataset
Eventually, we are ending up with 69,498 sites for Anura in America, 5,115 sites for Urodela in America, 96,822 sites for Anura in Europe, and 82,656 sites for Urodela in Europe.

### 3. Genetic Lineage Processing Script

This script processes and cleans genetic lineage data, converting it into spatial polygons and matching it with community data for downstream analysis. The main steps are as follows:

1. **Reading and filtering lineage files**:  
   The script imports genetic lineage CSV files and selects only the relevant columns (species, genus, lineage, latitude, and longitude). Duplicate entries are removed, and lineages with fewer than 3 records are excluded.

2. **Creating spatial polygons for lineages**:  
   Using the latitude and longitude data, the script calculates convex hulls to represent each lineage's geographic range. These ranges are converted into polygon shapefiles in the WGS84 coordinate reference system.

3. **Importing community data**:  
   For each species, the script imports community data from a corresponding dataset based on the species' taxonomic order. These data include site information such as location and species occurrences.

4. **Matching genetic lineages with community data**:  
   A spatial join is performed between the polygon shapefiles representing lineages and the community data. This step assigns lineages to community sites based on spatial overlap. The results are saved into new CSV files, either by appending to existing files or creating new ones if none exist.

5. **Processing pipeline**:  
   The script applies these steps to all available lineage files and saves the final cleaned and matched data for further phylogeographic analysis.

---

### 4. Environmental Data Extraction Script

**Script Name**: `0_environment_data.R`

**Authors**: Lucie Kuczynski, Katharine Marske  
**Last Edit**: November 2024

This script is used to process and extract environmental data from multiple raster files, matching these data to specific sampling sites in the community dataset. It performs spatial data extraction, aggregation, and finalizes the data for downstream ecological analysis.

**Key Steps and Components**:

1. **Setup and Preparation**  
   - The script starts by clearing memory and loading required libraries (`raster`, `sf`, `progress`, `dplyr`).
   - The function `generate_file_names()` generates meaningful names for each raster file based on its directory, assigning prefixes (e.g., "Ele" for elevation) and numbers for easier identification.

2. **Processing Raster Data**  
   - The `process_single_raster()` function takes a raster file, extracts values at specified coordinates (sampling sites), and returns these values in a structured format.
   - A list of raster files is created by recursively searching a specified directory, and a progress bar tracks the extraction process.
   - For each raster, environmental data is extracted for sampling sites based on latitude and longitude.

3. **Combining and Aggregating Environmental Data**  
   - The extracted data from each raster file is combined with site data into a final data frame.
   - Key environmental variables such as elevation, vegetation cover, solar radiation, wind, and land cover are aggregated by summing or averaging specific columns.
   - The script applies data cleaning to handle extreme values for certain variables (e.g., replacing elevation values below -4000 with NA).

4. **Data Formatting**  
   - The script renames columns to provide consistent, meaningful labels (e.g., "Bio_1" becomes "Temp" for temperature).
   - The land cover (LC) column is recoded with specific labels (e.g., "Broad_ever" for broad-leaved evergreen forests).
   - The final processed data frame includes site identifiers, spatial coordinates, and aggregated environmental variables ready for ecological and statistical analysis.

**Output**  
The resulting data frame provides environmental attributes matched to each sampling site, which can be used in further ecological analyses, such as examining the influence of environmental variables on community structure or species distributions.

---

### 5. Additional Information

This repository includes data, scripts, and documentation to support the analysis of community phylogeography, providing resources to clean, integrate, and analyze community, genetic, and environmental data for large-scale phylogeographic studies.
