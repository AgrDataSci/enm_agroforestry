####################################################
###### Calculate change in suitability

# .........................................
# .........................................
# Packages ####
library("tidyverse")
library("magrittr")
library("svglite")
library("gridExtra")
library("sp")
library("dismo")
library("raster")
library("rgdal")
library("rgeos")

# .........................................
# .........................................
# Read species names
sp <- "data/species_acronyms.csv"
sp %<>% 
  read_csv(.) %>% 
  filter(. , main_use != "Crop") %>%
  dplyr::select(species, acronym, main_use, cover_inventory)

# keep only species chosen for this analysis 
species_names <- as.vector(sp$acronym)

# .........................................
# .........................................
# Prepare data and files to process maps ####
# names of RCP scenarios
RCP <- c(45, 85)

output <- "processing/species_sets/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

#define extention of study area
ext <- raster::extent(-101, -77,  6.8, 22)
# define projection
proj <- "+proj=longlat +datum=WGS84"

### read shapefile of inland water
lakes <- "data/shapefiles/water_areas/mesoamerica_water_areas_dcw.shp"
lakes %<>%  
  readOGR(.) %>%
  subset(.$HYC_DESCRI == "Perennial/Permanent") %>%
  raster::crop(. , ext)


# Read country borders 
border <- "data/shapefiles/country_borders/Mesoamerica.shp"
border %<>% 
  readOGR(.) %>%
  raster::crop(. , ext)

# .........................................
# .........................................
# Run over species_names ####

# create NULL dataframe to keep data on area changes for all species 
changes_all <- NULL

for (i in seq_along(species_names)) {
  cat(i, "out of", length(species_names), "\n")
  
  #presence <- paste0("processing/enm/", species_names[i], "/ensembles/presence")
  presence <- paste0("E:/ensemble_modelling/", species_names[i], "/ensembles/presence/")
  
  #create directory to save rasters
  output2 <- paste0(output, species_names[i])
  
  dir.create(output2, showWarnings = FALSE, recursive = TRUE)
  
  # load rasters of current niches
  # raster of current presence-absence
  presence_current <- raster::stack(paste0(presence,species_names[i], "_bio_current.grd"))
  # crop raster using Mesoamerica extention
  presence_current <- raster::crop(presence_current, ext)
  # remove Caribbean islands
  presence_current <- raster::mask(presence_current, border, inverse=FALSE)
  # remove presence in lakes
  presence_current <- raster::mask(presence_current, lakes, inverse=TRUE)
  presence_current <- raster::stack(presence_current)

  # run over RCP models
  for(j in RCP) {
    # read presence-absence rasters under modelled RCP scenario##
    k <- paste(j, ".grd$", sep="") #insert ".grd" within "j" in order to stack only focal rasters
    presence_rcp <- raster::stack(list.files(presence, 
                                             pattern = k, 
                                             full.names = TRUE))
    #crop raster using Mesoamerica extention
    presence_rcp <- raster::crop(presence_rcp, ext)
    # remove Caribbean islands
    presence_rcp <- raster::mask(presence_rcp, border, inverse = FALSE)
    # remove presence in lakes
    presence_rcp <- raster::mask(presence_rcp, lakes, inverse = TRUE)
    
    # mean values of 1-0 raster (presence and absence)
    # this is also the measure agreement raster
    presence_rcp_mean <- raster::calc(presence_rcp, fun = mean)
    

    # Define likelihood RCP mask
    # more than 66% of the raster must agree on the presence
    # of the species in each grid cell
    # further information about likelihood at Mastrandrea el al (2010)
    # raster for future presence-absence
    presence_rcp <- presence_rcp_mean
    presence_rcp[presence_rcp[] < 0.659] <- 0 
    presence_rcp[presence_rcp[] >= 0.659] <- 1 
    #raster of threshold for future presence (defined by likelihood)
    thresh_presence_rcp <- presence_rcp
    thresh_presence_rcp[thresh_presence_rcp[] == 0] <- NA
    
    #identify change in suitability in RCP scenario
    #future minus current raster
    presence_rcp[presence_rcp[]==1] <- 2
    #change in suitability codes
    #  1 = always suitable
    #  0 = never suitable
    # -1 = no longer suitable
    #  2 = new habitat
    change_suit_rcp <- raster::overlay(presence_rcp, 
                                       presence_current, 
                                       fun = function(x,y) {(x-y)})
    
    # write tif files
    # current presence-absence layer 
    writeRaster(change_suit_rcp, 
                filename = paste0(output2,"/",species_names[i],"_",j,"_change.tif"), 
                format = "GTiff", 
                overwrite = TRUE)
    writeRaster(presence_rcp, 
                filename = paste0(output2,"/",species_names[i],"_",j,"_presence.tif"), 
                format = "GTiff", 
                overwrite = TRUE)
    writeRaster(presence_rcp_mean, 
                filename = paste0(output2,"/",species_names[i],"_",j,"_presence_agreement.tif"), 
                format = "GTiff", 
                overwrite = TRUE)
  }
  
  # current presence-absence
  writeRaster(presence_current, 
              filename = paste0(output2,"/",species_names[i],"_presence.tif"), 
              format="GTiff", 
              overwrite=TRUE)
}





