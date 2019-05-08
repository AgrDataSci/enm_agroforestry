####################################################
## Assess the number of available species in coffee 
## and cocoa areas under climate change 
# This script produces Fig. 1, Fig. 3 and Fig. S6
####################################################

# .............................................
# .............................................
#  Packages  ####
library("tidyverse")
library("magrittr")
library("svglite")
library("gridExtra")
library("sp")
library("dismo")
library("raster")
library("maptools")
library("rgdal")
library("rgeos")
library("grid")
library("grDevices")
library("ggsn")
library("svglite")

#...................................................
#...................................................
# Prepare input ####

# Define spatial parameters 
ext <- raster::extent(-101, -77, 7, 22)
proj <- "+proj=longlat +datum=WGS84"

crop_names <- c("COFFAR","THEOCA")

# Read country borders 
border <- "data/shapefiles/country_borders/Mesoamerica.shp"
border %<>%
  readOGR(.) %>%
  raster::crop(. , ext)

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

# Read species names and acronyms
sp <- "data/species_acronyms.csv"
sp %<>% 
  read_csv() %>%
  dplyr::select(acronym, main_use)

# Identify change in suitabily change in cocoa/coffee areas
output <- "output/change_suitability/"
dir.create(output, 
           showWarnings = FALSE, 
           recursive = TRUE)


# Read layers of change in suitability for crop species
crop_change <- list()
crop_change[[1]] <- stack("processing/layers_comb/COFFAR_rcp45_change.tif")
crop_change[[2]] <- stack("processing/layers_comb/COFFAR_rcp85_change.tif")
crop_change[[3]] <- stack("processing/layers_comb/THEOCA_rcp45_change.tif")
crop_change[[4]] <- stack("processing/layers_comb/THEOCA_rcp85_change.tif")

for(i in seq_along(crop_change)) names(crop_change[[i]]) <- "layer"

names(crop_change) <- c("COFFAR45","COFFAR85","THEOCA45","THEOCA85")

#...................................................
#...................................................
# Prepare tables with supporting data for Fig. 1 ####

# get information of adminstrative units, elevation and ecoregions from where coffee will lose and gain areas
# get adm info about coffee losing or gaining suitability
supp_info <- rbind(cbind(as.data.frame(crop_change[[1]], xy = TRUE), data.frame(crop="COFFAR45")),
                   cbind(as.data.frame(crop_change[[2]], xy = TRUE), data.frame(crop="COFFAR85")),
                   cbind(as.data.frame(crop_change[[3]], xy = TRUE), data.frame(crop="THEOCA45")),
                   cbind(as.data.frame(crop_change[[4]], xy = TRUE), data.frame(crop="THEOCA85")))

# remove NA's and never suitable areas
supp_info <- subset(supp_info, supp_info[,3] != 0 & !is.na(supp_info[,3]))

# get coordinates, convert in spatial points
coord <- sp::SpatialPoints(supp_info[1:2], proj4string = CRS(proj4string(gadm)))
summary(coord)

# add gadm info
supp_info <- cbind(supp_info, sp::over(coord, gadm))
# add ecoregion
supp_info <- cbind(supp_info, sp::over(coord, eco))
# add elevation
supp_info <- cbind(supp_info, raster::extract(elev, supp_info[,c(1:2)]))
head(supp_info)
supp_info <- supp_info[,c("x","y","layer","crop","ADM0_NAME","ADM1_NAME","ADM2_NAME","ECO_NAME","WWF_MHTNAM","MA_elev30as")]

# save it as csv file
write_csv(supp_info, 
          paste0(output, "si_change_suitability.csv"))

supp_info <- data.frame(xtabs(~ crop + layer + ADM0_NAME + ADM1_NAME + WWF_MHTNAM, data = supp_info))

supp_info <- supp_info[supp_info$Freq != 0, ]


write_csv(supp_info, 
          paste0(output, "si_change_suitability_frequencies.csv"))

#...................................................
#...................................................
# Prepare data for SI Fig. S6 ####

# # Define labels for map 
# # create the breaks and xy (lon lat) label vectors
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


# Vulnerability of cocoa/coffee by number of agroforestry options
#Read layers processed in script 2
tree <- list()
tree[[1]] <- stack("processing/layers_comb/trees_baseline.grd")
tree[[2]] <- stack("processing/layers_comb/trees_rcp45.grd")
tree[[3]] <- stack("processing/layers_comb/trees_rcp85.grd")
names(tree) <- c("curr","r45","r85")

cropr <- list()
cropr[[1]] <- stack("processing/layers_comb/COFFAR_rcp45.grd")
cropr[[2]] <- stack("processing/layers_comb/COFFAR_rcp85.grd")
cropr[[3]] <- stack("processing/layers_comb/THEOCA_rcp45.grd")
cropr[[4]] <- stack("processing/layers_comb/THEOCA_rcp85.grd")
names(cropr) <- c("COFFAR45","COFFAR85","THEOCA45","THEOCA85")

# Split rasters of change in suitability for each code
changes <- tibble(code = c( 1, -1, 2, 0),
                  change = c("suitable","no_longer","new_habitat", "never"))

# prepare maps and raster os scale of vulnerability over i crop species in rpc scenarios
# define labels and colours of maps
map_colours <- grDevices::colorRampPalette(c("#FFFF80", "#addd8e","#41ab5d", "#006837","#0C1078"))
map_colours2 <- grDevices::colorRampPalette(c("#FFFF80",'#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84'))

output <- "output/species_available"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# create dataframe for supporting info
supp_info <- NULL

# run over multiple i scenarios and gradient of vulnerability
for (i in seq_along(cropr)){

  layer <- stack(cropr[[i]])

  if(grepl(45, names(layer)[1]))  r <- tree[[2]] else r <- tree[[3]]

  r %<>%
    raster::calc( . , fun = sum)


  #use future crop distribution as mask for future vulnerability of agroforestry tree species
  masks <- list()
  for (k in seq_along(changes$change)){

    m <- subset(layer, k)

    masks[[ changes$change[k] ]] <- raster::mask(r, m, inverse = FALSE)

  }
  #create a dataframe for each class of change
  crop_changes <- rbind(cbind(as.data.frame(masks$suitable, xy=TRUE), change="suitable"),
                        cbind(as.data.frame(masks$no_longer, xy=TRUE), change="no_longer"),
                        cbind(as.data.frame(masks$new_habitat, xy=TRUE), change="suitable"),
                        cbind(as.data.frame(masks$never, xy=TRUE), change="never"))
  crop_changes <- na.omit(crop_changes)
  crop_changes <- subset(crop_changes, crop_changes$change!="never")

  names(crop_changes)[3] <- "layer"


  suitable <- crop_changes


  suitable$layer <- ifelse(suitable$change!="suitable", NA,suitable$layer)

  no_suitable <- crop_changes

  no_suitable$layer <- ifelse(no_suitable$change!="no_longer", NA,no_suitable$layer)

  # species available in suitable areas
  map1 <- ggplot() +
    #the raster info
    geom_raster(data=suitable, aes(y=y,x=x, fill=layer)) +
    #define colours and legend labels
    scale_fill_gradientn(name=NULL, colours=map_colours(50),
                         breaks=c(0, 20, 40, 60, 80, 100),
                         limits=c(0, 100), na.value = "grey80") +
    #add polygon of Mesoamerican administrative units (countries)
    geom_path(data = border, aes(x=long,y=lat,group=id), colour="grey30", size=0.3) +
    #labels for axis x and yi
    labs(x = NULL, y = NULL, title= NULL) +
    scale_x_continuous(breaks = lon_brks, labels = lon_lbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = lat_brks, labels = lat_lbls, expand = c(0, 0)) +
    #define theme, size and colours of axis and other elements
    theme(legend.text= element_text(size=13, colour="black"),
          axis.text.x = element_text(size=9, angle = 0, hjust=1, vjust=1, face="plain"),
          axis.text.y = element_text(size=9, angle = 0, hjust=1, vjust=1, face="plain"),
          legend.position = c(.94,.80),
          legend.key = element_rect(colour = "black",fill=NULL,size=0.5,linetype=1),
          legend.background = element_rect(fill = "transparent"),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_blank(),
          panel.background = element_blank())


  #export the ggplot 
  ggsave(paste0(output,"/", names(cropr[i]),"_species_available.svg"),
         plot = map1, dpi = map_res,
         width = 20, height = 13 , units = "cm")
  
  
  # species available in no longer suitable areas
  map2 <- ggplot() +
    #the raster info
    geom_raster(data=no_suitable, aes(y=y,x=x, fill=layer)) +
    #define colours and legend labels
    scale_fill_gradientn(name=NULL, colours=map_colours2(50),
                         breaks=c(0, 20, 40, 60, 80, 100),
                         limits=c(0, 100), na.value = "grey80") +
    #add polygon of Mesoamerican administrative units (countries)
    geom_path(data = border, aes(x=long,y=lat,group=id), colour="grey30", size=0.3) +
    #labels for axis x and yi
    labs(x = NULL, y = NULL, title= NULL) +
    scale_x_continuous(breaks = lon_brks, labels = lon_lbls, expand = c(0, 0)) +
    scale_y_continuous(breaks = lat_brks, labels = lat_lbls, expand = c(0, 0)) +
    #set borders into top and right sides of the map
    #theme_bw() +
    #define theme, size and colours of axis and other elements
    theme(legend.text= element_text(size=13, colour="black"),
          axis.text.x = element_text(size=9, angle = 0, hjust=1, vjust=1, face="plain"),
          axis.text.y = element_text(size=9, angle = 0, hjust=1, vjust=1, face="plain"),
          legend.position = c(.94,.80),
          legend.key = element_rect(colour = "black",fill=NULL,size=0.5,linetype=1),
          legend.background = element_rect(fill = "transparent"),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.background = element_blank(),
          panel.background = element_blank())
  
  # export the ggplot 
  ggsave(paste0(output,"/", names(cropr[i]),"_species_available_no_longer_suitable_areas.svg"),
         plot = map2, dpi = map_res,
         width = 20, height = 13 , units = "cm")


  # get supporting information
  # get coordinates, convert in spatial points
  coord <- sp::SpatialPoints(crop_changes[1:2], proj4string = CRS(proj4string(gadm)))

  #add gadm info
  crop_changes <- cbind(crop_changes, sp::over(coord, gadm))

  # add ecoregion
  crop_changes <- cbind(crop_changes, sp::over(coord, eco))

  #add elevation
  crop_changes <- cbind(crop_changes, raster::extract(elev, crop_changes[,c(1:2)]))

  crop_changes <- crop_changes[,c("x","y","layer","change","ADM0_NAME","ADM1_NAME","ADM2_NAME","ECO_NAME","WWF_MHTNAM","MA_elev30as")]

  crop_changes$MA_elev30as <- round(crop_changes$MA_elev30as, digits = -2)

  crop_changes$crop <- names(cropr[i])

  supp_info <- rbind(supp_info, crop_changes)
}


write_csv(supp_info, 
          paste0(output, "/si_species_available.csv"))


#...................................................
#...................................................
# Prepare data for Fig. 3 ####




# Check change in frequency of species per main use in current coffee areas 
# Current areas correspond to remain suitable and no longer suitable areas within the future maps
# combine both layers to get current distribution
cropcurrent <- list()
cropcurrent[[1]] <- merge(cropr[[1]][[1]], cropr[[1]][[2]])
cropcurrent[[2]]  <- merge(cropr[[3]][[1]], cropr[[3]][[2]])

uses <- sort(unique(sp$main_use))[-1]

sp %>% 
  dplyr::select(acronym,main_use) %>% 
  filter(main_use != "Crop") %>% 
  group_by(main_use) %>% 
  summarise(n=n()) ->
  n_uses

rcp <- c("current", 45,85)

supp_info <- NULL


for(i in seq_along(crop_names)){

  for(j in seq_along(uses)){
    
    for(k in seq_along(rcp)){
      
      if(rcp[k]=="current"){
        
        r <- tree[[1]][[sp$acronym[sp$main_use==uses[j]]]]
        r <- calc(r , sum)
      
      }else{
        r <- stack(paste0("processing/layers_comb/sum_",uses[j],"_trees_rcp",rcp[k],".tif"))
      }
      
      r %<>%
        mask(. , cropcurrent[[i]])
      
      r <- as.data.frame(r, xy=TRUE)
      
      r <- r[!is.na(r[,3]), ]
      
      #add elevation
      r <- cbind(r, raster::extract(elev, r[,c(1:2)]))
      
      r$MA_elev30as <- round(r$MA_elev30as, digits = -2)
      
      r$crop <- crop_names[i]
      
      r$use <- uses[j]
      
      r$rcp <- rcp[k]
      
      names(r)[3:4] <- c("freq","elev") 
      
      r$freq_r <- r$freq / as.integer(n_uses[j,2])
      
      supp_info <- rbind(supp_info, r)
      
      
      
    }
    
  }

}


#Summarise results 
head(supp_info)


# Filter data by altitudinal range
keep <- supp_info$crop=="COFFAR" & supp_info$elev >= 400 & supp_info$elev <= 2400 |
  supp_info$crop=="THEOCA" & supp_info$elev >= 0 & supp_info$elev <= 1200 

supp_info <- supp_info[keep, ]

write_csv(supp_info, 
          paste0(output,"./si_species_available_uses.csv"))


# # number of species per use per grid
# info <- supp_info
# info[1:2] <- lapply(info[1:2], function(X) round(X, 3))
# info$id <- paste0(info$x, info$y, "rcp", info$rcp)
# info %<>% 
#   dplyr::select(. , id, freq, use, crop) %>% 
#   as_tibble()
# 
# info <- reshape2::dcast(info, id ~ crop + use, value.var="freq")
# 
# write_csv(info, paste0(output,"./vulnerability_sp_uses_complete_dcast.csv"))

supp_info %<>%
  dplyr::select(. , elev, crop, use, rcp, freq_r) %>%
  group_by(., crop, use, elev, rcp) %>%
  summarise(. , freq = median(freq_r))

supp_info %<>%
  mutate(., rcp = ifelse(rcp=="current", "Current", ifelse(rcp=="45","RCP45","RCP85")))

write_csv(supp_info, paste0(output,"/sp_uses.csv"))

# Export charts 
colors <- c("#1b9e77","#d95f02","#757077")

plots <- list()

k <- 0

for(i in seq_along(crop_names)){
  
  for(j in seq_along(uses)){
    
    crop_i <- crop_names[i]
    
    X <- supp_info %>% 
      filter(crop == crop_i, 
             use == uses[j])
    
    p <- ggplot(X, aes(x = elev, y = freq)) +
      geom_line(aes(colour = rcp, linetype = rcp), 
                show.legend = TRUE, size = 1.5) + 
      labs(x = "Elevation (m)", 
           y = "Fraction of available tree species", 
           title = NULL) +
      scale_x_continuous(limits=c(min(X$elev), (max(X$elev)+100)), 
                         breaks = seq(min(X$elev),max(X$elev),by=200), 
                         expand = c(0,0)) +
      scale_y_continuous(limits = c(0, 1.05), 
                         breaks = as.numeric(quantile(0:1, probs=seq(0, 1, .2))),
                         expand = c(0.05,0)) +
      scale_colour_manual(values = colors, name = "") +
      guides(linetype = FALSE) +
      theme(axis.text.x = element_text(size=14, angle = 0,  
                                       hjust=0.5,vjust=1, face="plain", colour = "black"),
            axis.title.x = element_text(size=14, face="bold", colour = "black"),
            axis.text.y = element_text(size=14, angle = 0,  
                                       hjust=1,vjust=0.5, face="plain", colour = "black"),
            axis.title.y = element_text(size=14, face="bold", 
                                        colour = "black"),
            axis.line = element_line(colour = "black"),
            plot.background = element_blank(),
            panel.background = element_blank(), 
            legend.text = element_text(size=15, colour="black"),
            legend.position = "bottom",
            legend.key = element_rect(colour = "white",fill= "white" ,size=0.5,linetype=1),
            legend.background = element_rect(fill = "white"),
            plot.margin=grid::unit(c(10,10,10,10), "mm")) 
    
    k <- k+1
    
    plots[[k]] <- p
    
    # ggsave(paste0(output,"/",crop_names[i], "_freq_",tolower(uses[j]),"species" ,".png"),
    #        plot = p, dpi = 600,
    #        width = 26, height = 20 , units = "cm")
    
  }
  
}


p <- grid.arrange(
  plots[[1]],
  plots[[4]],
  plots[[2]],
  plots[[5]],
  plots[[3]],
  plots[[6]],
  nrow = 3)

ggsave(paste0(output,"/Fig3_freq_elevation.svg"),
       plot = p,
       width = 20,
       height = 30,
       units = "cm")

  







