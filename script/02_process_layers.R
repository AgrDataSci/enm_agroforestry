####################################################
###### SCRIPT TO PROCESS LAYERS FROM ENSEMBLE MODELLING
# Updated 13Aug2018
####################################################

library("tidyverse")
library("magrittr")
library("sp")
library("dismo")
library("raster")
library("maptools")
library("rgdal")
library("rgeos")


# Define extention of study region
ext <- raster::extent(-101, -77,  6.8, 22)
proj <- "+proj=longlat +datum=WGS84" 

# Read species names and acronyms
sp <- read_csv("data/species_acronyms.csv")

# keep only species chosen for this analysis 
sp %<>% 
  filter(main_use != "Crop") %>%
  dplyr::select(acronym, main_use)

# names of RCP scenarios and crops 
crop <- c("COFFAR", "THEOCA")

RCP <- c(45, 85)

# Read layer of inland water and remove areas of raster in water (lakes, rivers, etc)
#source http://www.diva-gis.org/Data
lakes <- 
  readOGR(dsn = "data/shapefiles/water_areas/mesoamerica_water_areas_dcw.shp") %>%
  subset(.$HYC_DESCRI == "Perennial/Permanent") %>%
  raster::crop(. , ext)


# Read country borders 
border <- 
  readOGR(dsn = "data/shapefiles/country_borders/Mesoamerica.shp") %>%
  raster::crop(. , ext)

# Read rasters of tree species
# run over the raster of current presence
# add to a list then stack 
tree <- list()

pb <- txtProgressBar(min = 1, max = nrow(sp), style = 3)

spnames <- sort(unique(sp$acronym))

for (i in seq_along(spnames)) {
  
  r <- list.files(paste0("processing/enm/", spnames[i] , "/ensembles/presence"),
                  pattern = "_bio_current.gri$", 
                  full.names = TRUE)
  
  r %<>% 
    raster::stack( . ) %>% 
    raster::crop(. , ext) %>% 
    raster::mask(. , border, inverse = FALSE) %>% 
    raster::mask(. , lakes, inverse = TRUE)
  
  tree[[i]] <- r
     
  
  setTxtProgressBar(pb, i)
  
}

close(pb)

# convert this list into a raster stack 
tree %<>%
  raster::stack(.)

names(tree) <- gsub("_bio_current_presence", "", names(tree))

output <- "processing/layers_comb/"

dir.create(output,
           showWarnings = FALSE,
           recursive = TRUE)


# Export all layers combined as a single raster
writeRaster(tree, 
            filename = "./processing/layers_comb/trees_baseline.grd", 
            format = "raster", 
            overwrite = TRUE, 
            bylayer = FALSE)

# Export sum of all presences
tree %<>% 
  raster::calc(. , fun = sum)

writeRaster(tree, 
            filename = "./processing/layers_comb/sum_trees_baseline.tif", 
            format = "GTiff", 
            overwrite = TRUE)

# Read layers of future presence for all species
tree_rcp <- list()

pb <- txtProgressBar(min = 1, max = nrow(sp), style = 3)

for (i in seq_along(spnames)) {

  for (j in seq_along(RCP)){

    #read rasters of k rcp scenarios
    r <- list.files(paste0("processing/enm/", spnames[i] , "/ensembles/presence"),
                    pattern = paste0(RCP[j],".gri$"),
                    full.names = TRUE)
    r %<>%
      raster::stack(. ) %>%
      raster::crop(. , ext) %>% #crop rasters within defined extention
      raster::calc(. , fun = mean) %>% #calculate the mean of all k rcp scenarios
      raster::mask(. , border, inverse = FALSE) %>% 
      raster::mask(. , lakes, inverse = TRUE)
    
    # if less than 66% of models agree with the presence then 0
    r[r[] < 0.659 ] <- 0
    # if more than 66% of models agree with the presence then 1 
    r[r[] > 0.659 ] <- 1 
    
    rname <- paste0(sp$acronym[i],RCP[j])
    
    tree_rcp[[rname]] <- r
  
  }
  
  setTxtProgressBar(pb, i)
}

close(pb)

# Export layers 
for (i in seq_along(RCP)){
  
  r <- subset(tree_rcp, grepl(RCP[i] , names(tree_rcp)))
  
  r <- raster::stack(r)
  
  names(r) <- gsub("[0-9]+", "", names(r))
  
  writeRaster(r, 
              filename = paste0("processing/layers_comb/trees_rcp", RCP[i],".grd"), 
              format="raster",
              overwrite = TRUE, 
              bylayer = FALSE)
  
  # sum the future presence of all focal species in i RCP scenario
  r <- raster::calc(r, fun=sum)
  # add to list
  writeRaster(r, 
              filename = paste0("processing/layers_comb/sum_trees_rcp", RCP[i],".tif"), 
              format = "GTiff", 
              overwrite = TRUE)

}

# Generate rasters dividing species per main use
uses <- sort(unique(sp$main_use))

for(i in seq_along(uses)){

  sp_use <-  sp$acronym[sp$main_use %in% uses[i]]

  for(j in seq_along(RCP)){

    r <- subset(tree_rcp, names(tree_rcp) %in% paste0(sp_use , RCP[j]) )

    r <- raster::stack(r)

    names(r) <- gsub("[0-9]+", "", names(r))

    writeRaster(r,
                filename = paste0("processing/layers_comb/", uses[i] ,"_trees_rcp", RCP[j],".grd"),
                format = "raster",
                overwrite = TRUE,
                bylayer = FALSE)

    # sum the future presence of all focal species in i RCP scenario
    r <- raster::calc(r, fun = sum)
    # add to list
    writeRaster(r,
                filename = paste0("processing/layers_comb/sum_", uses[i] ,"_trees_rcp",RCP[j],".tif"),
                format = "GTiff",
                overwrite = TRUE)

  }

}



# CROP #####
# read rasters of i crop species
crop_r <- list()
#run over crop names and rcp scenarios to identify changes in suitability between scenarios
for (i in seq_along(crop)){

  r <- list.files(paste0("processing/enm/", crop[i] , "/ensembles/presence"),
                  pattern= paste0(crop[i],"_bio_current.gri$",sep=""),
                  full.names = TRUE)
  
  #read current raster of i species
  r %<>%
    raster::stack(.) %>%
    raster::crop(. , ext) %>% #crop rasters within defined extention
    raster::mask(. , border, inverse = FALSE) %>% 
    raster::mask(. , lakes, inverse = TRUE)
  
  writeRaster(r, 
              filename = paste0("processing/layers_comb/", crop[i], "_baseline.grd"), 
              format = "raster", 
              overwrite = TRUE)

  #add to crop_raster list
  crop_r[[crop[i]]] <- r
  # run over rcp scenarios for i crop species
  for (j in seq_along(RCP)){
    
    #read rasters of j rcp scenario
    rf <- list.files(paste0("processing/enm/", crop[i] , "/ensembles/presence"),
                     pattern = paste0(RCP[j],".gri$",sep=""),
                     full.names = TRUE)
      
    rf %<>%
      raster::stack(.) %>%
      raster::crop(. , ext) %>% #crop rasters within defined extention
      raster::mask(. , border, inverse=FALSE) %>% 
      raster::mask(. , lakes, inverse=TRUE)  %>%
      raster::calc(. , fun = mean) #calculate the mean of all k rcp scenarios
    
    #reclassify raster, values 
    rf[rf[] < 0.659 ] <- 0 #if less than 66% of models agree with the presence then 0
    rf[rf[] > 0.659 ] <- 2 #if more than 66% of models agree with the presence then 2 
    #add to list
    crop_r[[paste0(crop[i],RCP[j])]] <- rf
    #identify changes in suitability between current and j rcp scenario of each i crop species
    #change in suitability codes
    #  1 = always suitable
    #  0 = never suitable
    # -1 = no longer suitable
    #  2 = new habitat
    rf %<>%
      raster::overlay( . , r, fun = function(x,y) {(x-y)})
    
    writeRaster(rf, 
                filename = paste0("processing/layers_comb/", crop[i], "_rcp", RCP[j],"_change.tif"), 
                format="GTiff", 
                overwrite=TRUE)
    
    crop_r[[ paste0(crop[i],RCP[j],"_change",sep="") ]] <- rf
  
  }

}

# separate rasters of change in suitability for each code
changes <- tibble(code=c( 1, -1, 2, 0),
                  change=c("suitable","no_longer","new_habitat", "never"))

#run loop over codes of suitability change

for (i in seq_along(crop)){
  for (j in seq_along(RCP)){
    
    rf <- NULL
    
    for(k in 1:nrow(changes)) {
      
      r <- crop_r[[ paste0(crop[i], RCP[j], "_change") ]]
      
      r[r[] != as.integer(as.character(changes[k,1])) ] <- NA
      
      r[r[] == as.integer(as.character(changes[k,1])) ] <- 1
      
      names(r) <- changes[k,2]
      
      rf <- stack(r, rf)
      
      
    }
    
    names(rf)
    
    rf <- subset(rf , changes$change)
    
    writeRaster(rf, 
                filename = paste0("processing/layers_comb/",crop[i], "_rcp", RCP[j],".grd"), 
                format = "raster", 
                overwrite = TRUE, 
                bylayer = FALSE)
    
  }
}





