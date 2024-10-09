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

Note: All life stages have been kept, based on the assumption that the presence of eggs or tadpoles implies the existence of adults capable of reproduction.

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
Eventually, we are ending up with 169,196 sites that encompass 638,012 amphibian recordings. These files have been saved as `sites_5sp.csv` and `occurrences_5sp.csv`.