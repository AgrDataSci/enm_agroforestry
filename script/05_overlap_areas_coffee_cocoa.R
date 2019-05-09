####################################################
###### OVERLAP AREAS COFFEE COCOA
# This script generate Fig. 2 showing areas where 
# cocoa could replace coffee
####################################################

library("maptools")
library("rgeos")
library("dismo")
library("raster")
library("rgdal")
library("grid")
library("ggsn")
library("tidyverse")
library("svglite")

#define extention of MesoAmerica
ext <- raster::extent(-101, -77,  6.8, 22)
crs_proj <- "+proj=longlat +datum=WGS84" 

#run over layers of change in suitability of crop areas
crop_names <- c("THEOCA","COFFAR")
rcp <- c(45,85)
crop <- list()
for (i in crop_names){
  for( j in rcp) {
    file <- raster::raster(paste0("processing/layers_comb/", i, "_rcp", j ,"_change.tif"))
    # keep only 1 (maintain) and 2 (new areas) in cocoa rasters
    if(i=="THEOCA") file[file[] < 1] <- 0
    if(i=="THEOCA") file[file[] >= 1] <- 2 
    if(i=="COFFAR") file[file[] >= 1] <- 3
    crop[[paste(i,j,sep="")]] <- file
    rm(file)
  }
}

# codes for the next raster as follow
# -1 = No longer suitable for coffee / Never suitable for cocoa
# 0 = Never suitable for coffee or cocoa
# 1 = No longer suitable for coffee / Suitable for cocoa
# 2 = Never suitable for coffee / Suitable for cocoa
# 3 = Suitable for coffee / Never suitable for cocoa
# 5 = Suitable for coffee / Suitable for cocoa

# .....................................
# ..................................... 

# Read gadm data 
gadm <- "data/shapefiles/country_borders/Mesoamerica.shp"
gadm %<>%
  readOGR(.) %>%
  raster::crop(. , ext)

# create the breaks and xy (lon lat) label vectors
lon_brks <- seq(-100,-78,3)
lat_brks <- seq(7,22,3)
lon_lbls <- unlist(lapply(lon_brks, function(x) {
  ifelse(x < 0, paste(x*-1, " W", sep=""), 
         ifelse(x > 0, paste(x, " E",sep=""),x))}
))

lat_lbls <- unlist(lapply(lat_brks, function(x) {
  ifelse(x < 0, paste(x*-1, " S", sep=""), 
         ifelse(x > 0, paste(x, " N",sep=""),x))}
))

map_colours <- c("#2c7bb6","#abd9e9","#d7191c","#ffffbf","#FFFFFF")

# output directory
output <- "output/cocoa_replace_coffee/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

#generate maps 
rcp45 <- raster::stack(crop$COFFAR45,crop$THEOCA45)
rcp85 <- raster::stack(crop$COFFAR85,crop$THEOCA85)

scenario <- list(rcp45 = rcp45, rcp85 = rcp85)

supp_info <- NULL
overlap_count <- NULL
overlap <- list()

for (i in names(scenario)){
  area <- raster::calc(scenario[[i]], fun=sum)
  #add to list 
  overlap[[i]] <- area
  count_area <- data.frame(raster::freq(area, useNA="no"), scenario=i)
  overlap_count <- rbind(count_area, overlap_count)
  
  area[area[] == 1] <- 101 # 1 = No longer suitable for coffee / Suitable for cocoa
  area[area[] == 5] <- 102 # 5 = Suitable for coffee / Suitable for cocoa
  area[area[] == -1] <- 103 # -1 = No longer suitable for coffee / Never suitable for cocoa
  area[area[] == 3] <- 105 # 3 = Suitable for coffee / Never suitable for cocoa
  area[area[] == 2] <- 106 # 2 = Never suitable for coffee / Suitable for cocoa
  area[area[] == 0] <- 106 # 0 = Never suitable for coffee or cocoa
  
  area <- as.data.frame(area, xy=T)
  area <- subset(area, !is.na(area$layer))
  #add dataframe to supporting info
  area$rcp <- i
  supp_info <- rbind(supp_info, area)
  
  #export this as an geoTiff 
  coordinates(area) <- ~x+y
  gridded(area) <- TRUE
  sp::proj4string(area) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  area <- raster(area)
  spplot(area)
  
  writeRaster(area, 
              filename = paste0(output, "cocoa_replace_coffee_",i,".tif"), 
              format = "GTiff", 
              overwrite = TRUE)
  
}

head(overlap_count)
overlap_count <- cbind(overlap_count, 
                       data.frame(label= rep(c("No longer suitable for coffee / Never suitable for cocoa",
                                                 "Never suitable for coffee or cocoa",
                                                 "No longer suitable for coffee / Suitable for cocoa",
                                                 "Never suitable for coffee / Suitable for cocoa",
                                                 "Suitable for coffee / Never suitable for cocoa",
                                                 "Suitable for coffee / Suitable for cocoa"), 2)))

write_csv(overlap_count, 
          paste0(output, "cocoa_replace_coffee_summary.csv"))

# .....................................
# .....................................
# Get supporting info ####
supp_info$label <- with(supp_info, ifelse(layer==101, "Nolongersuitable4coffee/Suitable4cocoa",
                                          ifelse(layer==102,"Suitable4coffee/Suitable4cocoa",
                                                 ifelse(layer==103, "Nolongersuitable4coffee/Neversuitable4cocoa",
                                                        ifelse(layer==105, "Suitable4coffee/Neversuitable4cocoa",
                                                               ifelse(layer==106, NA,NA))))))


supp_info <- na.omit(supp_info)

# Read gadm data 
gadm <- "data/shapefiles/country_borders/Mesoamerica_adm3.shp"
gadm %<>%
  readOGR(.) %>%
  raster::crop(. , ext)


# Read TNC ecorregion data
eco <- "data/shapefiles/ecorregion/ecorregion_MA.shp"
eco %<>%
  readOGR(.) %>%
  raster::crop(. , ext)

# Read geografic adminstrative areas
# read elevation data
elev <- "data/elevation/MA_elev30as.tif"
elev %<>%
  raster::raster(.) %>%
  raster::crop(., ext) %>%
  raster::stack(.)

#get coordinates, convert in spatial points
coord <- sp::SpatialPoints(supp_info[1:2], proj4string = CRS(proj4string(eco)))
# add gadm info
supp_info <- cbind(supp_info, sp::over(coord, gadm))
# add ecoregion
supp_info <- cbind(supp_info, sp::over(coord, eco))
# add elevation
supp_info <- cbind(supp_info, raster::extract(elev, supp_info[,c(1:2)]))

supp_info <- supp_info[,c("x","y","layer","label","rcp","ADM0_NAME","ADM1_NAME","ADM2_NAME","ECO_NAME","WWF_MHTNAM","MA_elev30as")]

head(supp_info)


# export outputs
write_csv(supp_info, 
          paste0(output,"si_cocoa_replace_coffee.csv"))

supp_info <- data.frame(xtabs(~ label + rcp + ADM0_NAME + ADM1_NAME + WWF_MHTNAM, data = supp_info))

supp_info <- subset(supp_info, supp_info$Freq!=0)

write_csv(supp_info, 
          paste0(output,"si_cocoa_replace_coffee_frequencies.csv"))


