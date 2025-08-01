################################ 1.1 Title #####################################

# Script Title: Creating covariates for sites from the Great Lakes Marsh Monitoring
# Program and Coastal Wetland Monitoring Program

# Script Author: Rory Macklin (rmacklin@birdscanada.org) for Birds Canada

# Date: July 31, 2025


################################ 1.2 Setup #####################################

## Install requisite packages

#install.packages("librarian")

librarian::shelf(tidyverse, sf, mapview, terra, svMisc, landscapemetrics)

sf_use_s2(FALSE)

proj <- "+proj=aea +lon_0=-82.7929688 +lat_1=38.9739639 +lat_2=49.6093187 +lat_0=44.2916413 +datum=WGS84 +units=m +no_defs"

# Read in point count sites, create spatial object

sites <- read_csv("./Data/Sites/Build siteCovs FINAL.csv") %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(proj) %>%
  mutate(stationID = paste0(routeOrSite, "-", stationOrPoint),
         region = NA, basin = NA)

# Read in BCR polygons, transform CRS to match sites.

bcr <- read_sf("./Data/BCR") %>%
  st_transform(crs(sites))

# Find the nearest BCR feature for each site. If the site falls in BCR12, we
# assign it to the northern region. Otherwise, we assign it to the southern region.
# We do this by getting the nearest BCR feature to each point as some points fall just
# outside polygons as they are slightly offshore.

for(i in sites$stationID) {
  
  suppressMessages(
    
    bcr.num <- bcr$BCRNumber[st_nearest_feature(sites[sites$stationID == i,], bcr)]
    
  )
  
  sites$region[sites$stationID == i] <- ifelse(bcr.num == 12, "northern", "southern")
  
  progress(which(sites$stationID == i), nrow(sites))
  
}

# Read in basin polygons, transform CRS to match sites.

basins <- read_sf("./Data/Basin") %>%
  st_transform(crs(sites))

# Find the nearest basin feature for each site.  We do this by getting the nearest 
# basin feature to each point as some points fall just outside polygons as they 
# are slightly offshore.

for(i in sites$stationID) {
  
  suppressMessages(
    
    basin.name <- basins$merge[st_nearest_feature(sites[sites$stationID == i,], basins)]
    
  )
  
  sites$basin[sites$stationID == i] <- gsub("lk_", "", basin.name)
  
  progress(which(sites$stationID == i), nrow(sites))
  
}

# Fix a few abbreviated names

sites <- mutate(sites, basin = case_when(basin == "ont" ~ "ontario",
                                         basin == "mich" ~ "michigan",
                                         basin == "stLaw" ~ "stLawrence",
                                         .default = basin))

# Buffer sites to a local area with 250m radius to extract landcover data.

sites.local <- st_buffer(sites, 250)

landcover <- rast("./Data/Landcover/glahf_land_cover_11_12_00_nlcd_solris_plo.gdb")

# Replace category "1 - Great Lakes Waters" with NA as not to bias landscape metrics

landcover <- classify(landcover, cbind(1, NA))

sites <- st_transform(sites, crs = crs(landcover))
sites.local <- st_transform(sites.local, crs = crs(landcover))

for(i in sites$stationID) {
  
  # Add a check for whether the site point falls within the raster. If it doesn't,
  # leave that site's landcover value as NA.
  
  # Two options for removing sites outside of coverage - this one checks whether
  # the value at the site point is NA, assigns NA if it is. Produces
  # 533 NA sites.
  #
  # check <- extract(landcover, vect(sites[sites$stationID == i,]))
  # 
  # if(is.na(check[1, 2])) {
  
  # This option checks whether the entire 250m radius is NA, assigns NA if it is.
  # More forgiving. Produces 151 NA sites.
  
  tmp <- mask(crop(landcover, sites.local[sites.local$stationID == i,]), sites.local[sites.local$stationID == i,])
  
  if(is.null(unique(tmp))) {

    sites$percMarsh250m[sites$stationID == i] <- NA
    
    
  } else {
    
    lsm.tbl <- calculate_lsm(tmp, level = "class", metric = "pland")
    
    sites$percMarsh250m[sites$stationID == i] <- ifelse(95 %in% lsm.tbl$class, lsm.tbl$value[lsm.tbl$class == 95], 0)
    
  }

  progress(which(sites$stationID == i), nrow(sites))
  
}

# Buffer sites to a neighbourhood area with 2500m radius to extract landcover data.

sites.neighbourhood <- st_buffer(sites, 2500) %>%
  st_transform(crs(landcover))


for(i in sites$stationID) {
  
  # Two options for removing sites outside of coverage - this one checks whether
  # the value at the site point is NA, assigns NA if it is. Produces
  # 533 NA sites.
  #
  # check <- extract(landcover, vect(sites[sites$stationID == i,]))
  # 
  # if(is.na(check[1, 2])) {
  
  # This option checks whether the entire 2500m radius is NA, assigns NA if it is.
  # More forgiving. Produces 8 NA sites.
  
  tmp <- mask(crop(landcover, sites.neighbourhood[sites.neighbourhood$stationID == i,]), sites.neighbourhood[sites.neighbourhood$stationID == i,])
  
  if(is.null(unique(tmp))) {
    
    sites$percMarsh2500m[sites$stationID == i] <- NA
    sites$percUrban2500m[sites$stationID == i] <- NA
    sites$percAg2500m[sites$stationID == i] <- NA
    
    
  } else {
    
    lsm.tbl <- calculate_lsm(tmp, level = "class", metric = "pland")
    
    sites$percMarsh2500m[sites$stationID == i] <- ifelse(95 %in% lsm.tbl$class, lsm.tbl$value[lsm.tbl$class == 95], 0)
    sites$percUrban2500m[sites$stationID == i] <- ifelse(2 %in% lsm.tbl$class, lsm.tbl$value[lsm.tbl$class == 2], 0)
    sites$percAg2500m[sites$stationID == i] <- ifelse(8 %in% lsm.tbl$class, lsm.tbl$value[lsm.tbl$class == 8], 0)
    
  }
  
  progress(which(sites$stationID == i), nrow(sites))
  
}

sites <- sites %>%
  st_transform(4326) %>%
  cbind(st_coordinates(.)) %>%
  st_drop_geometry() %>%
  select(program, routeOrSite, stationOrPoint, latitude = Y, longitude = X,
         region, basin, percMarsh250m, percMarsh2500m, percUrban2500m,
         percAg2500m)

write_csv(sites, "./Outputs/cwmp_glmmp_covdata.csv")


