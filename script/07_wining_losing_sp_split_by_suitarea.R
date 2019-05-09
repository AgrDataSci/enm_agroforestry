####################################################
# #### CALCULATE THE CHANGES IN SUITABILITY OF FOCAL 
# SPECIES IN FUTURE SUITABLE AREAS FOR COCOA/COFFEE
# This script generates  Fig. S3, S4 and S5
# Winning and losing agroforestry species
# Updated 11Mar2019

#..............................................
#..............................................
# Packages ####
library("tidyverse")
library("svglite")
library("magrittr")
library("raster")
library("gridExtra")
library("scales")


#..............................................
#..............................................
# Read tree data ####
sp <- "data/species_acronyms.csv"
sp %<>% 
  read_csv(.) %>% 
  filter(. , main_use != "Crop") %>%
  dplyr::select(species, acronym, main_use, cover_inventory)

species_names <- sort(as.character(sp$acronym))

rcp <- c("current","45","85")

# read coffee and cocoa future areas
# create two masks
# one for areas where crop will remain or maintain it suitability
# the other for areas where crop lost suitability
crop_names <- c("COFFAR","THEOCA")

crop <- list()

for (i in seq_along(crop_names)) {
  for(j in seq_along(rcp)) {
    
    if(rcp[j]=="current"){
      r <- raster(paste0("processing/layers_comb/", crop_names[i],"_baseline.grd"))
      }else{
      r <- raster(paste0("processing/layers_comb/", crop_names[i], "_rcp", rcp[j], "_change.tif"))
      }
      
      suit <- r
      suit[suit[] < 1] <- NA 
      suit[suit[] >= 1] <- 1
      crop[[paste0(crop_names[i],rcp[j],"suit")]] <- suit
      
      # no longer suitable areas
      not <- r
      not[not[] != -1] <- NA
      not[not[] == -1] <- 1
      crop[[paste0(crop_names[i],rcp[j],"not")]] <- not
    }
}

rm(not, suit, r)

names(crop)
changes <- NULL

#..............................................
#..............................................
# Read changes of focal species ####
# this loop calculate the frequencies of suitability of each 
# focal species per masks generated above 
# how is the change in suitability of each species under 
# new/remain crop areas and lost crop areas

for(i in species_names) {
  for (j in rcp) {
    if (j == "current") {
      r <-
        raster::raster(paste(
          "./processing/species_sets/",
          i,
          "/",
          i,
          "_presence.tif",
          sep = ""
        ))
    } else{
      r <-
        raster::raster(paste(
          "./processing/species_sets/",
          i,
          "/",
          i,
          "_",
          j ,
          "_change.tif",
          sep = ""
        ))
    }
    
    for (k in crop_names) {
      for (l in c("suit", "not")) {
        land <- mask(r, crop[[paste(k, j, l, sep = "")]])
        frequencies <- data.frame(array(dim = c(1, 11)))
        names(frequencies) <-
          c(
            "acronym",
            "system",
            "never_suitable(0)",
            "remains_suitable(1)",
            "no_longer_suitable(-1)",
            "new_habitat(2)",
            "area_km2",
            "remains_suitable_(1)_perc",
            "no_longer_suitable(-1)_perc",
            "new_habitat_(2)_per",
            "net_change_perc"
          )
        freq1 <- as.data.frame(raster::freq(land, useNA = "no"))
        freq1 <- freq1[is.na(freq1[, 1]) == F,]
        
        frequencies[1, 1] <- i
        frequencies[1, 2] <- paste(k, j, l, sep = "_")
        
        if (length(freq1[freq1[, 1] == 0, 2]) > 0) {
          frequencies[1, 3] <- freq1[freq1[, 1] == 0, 2]
        }
        if (length(freq1[freq1[, 1] == 1, 2]) > 0) {
          frequencies[1, 4] <- freq1[freq1[, 1] == 1, 2]
        }
        if (length(freq1[freq1[, 1] == -1, 2]) > 0) {
          frequencies[1, 5] <- freq1[freq1[, 1] == -1, 2]
        }
        if (length(freq1[freq1[, 1] == 2, 2]) > 0) {
          frequencies[1, 6] <- freq1[freq1[, 1] == 2, 2]
        }
        
        frequencies[is.na(frequencies)] = 0
        
        #calculate current areas in km2 - each grid cell has ~4.9km2
        frequencies[7] <-
          as.integer(apply(
            frequencies[, c(4, 5)],
            MARGIN = 1,
            FUN = sum,
            na.rm = TRUE
          )) * 4.9
        #percent of always suitable areas
        frequencies[8] <-
          frequencies[, (4)] / as.integer(apply(
            frequencies[, c(4, 5)],
            MARGIN = 1,
            FUN = sum,
            na.rm = TRUE
          )) * 100
        #percent of no longer suitable areas
        frequencies[9] <-
          frequencies[, (5)] / as.integer(apply(
            frequencies[, c(4, 5)],
            MARGIN = 1,
            FUN = sum,
            na.rm = TRUE
          )) * 100
        #percent of new suitable areas
        frequencies[10] <-
          frequencies[, (6)] / as.integer(apply(
            frequencies[, c(4, 5)],
            MARGIN = 1,
            FUN = sum,
            na.rm = TRUE
          )) * 100
        #percent of net change
        frequencies[11] <- frequencies[, 10] - frequencies[, 9]
        #add to changes dataframe
        changes <- rbind(changes, frequencies)
      }
    }
  }
}


changes %<>% 
  inner_join(., sp, by="acronym")

# export dataframe
output <- "output/wining_losing_sp/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

write_csv(changes, 
          paste(output, "wining_losing_sp_splited.csv"))

#..............................................
#..............................................
# Calculate changes ####
changes <- changes[,-c(7:11)]
names(changes) <- c("acronym","System","never","remain","no_longer","new","species","main_use","abundance")

changes %<>% 
  mutate(area = rep("suitable", nrow(changes)),
         area_raster = never + remain + no_longer + new, #number of cells (areas) the raster cover
         cover = ((remain + no_longer) * 100) / area_raster, #proportion of current area
         future_cover = ((remain + new) * 100) / area_raster, #proportion of future areas
         change = future_cover  - cover) %>%  #change between current and future suitable areas
  as_tibble()


# generate labels of systems (land use)
changes %<>%
  mutate(System_label = factor(ifelse(grepl("COFFAR_45",System) & area=="suitable", "Coffee RCP 4.5",
                                      ifelse(grepl("COFFAR_45",System) & area=="no_suitable", "Coffee RCP 4.5",
                                             ifelse(grepl("COFFAR_85",System) & area=="suitable", "Coffee RCP 8.5",
                                                    ifelse(grepl("COFFAR_85",System) & area=="no_suitable", "Coffee RCP 8.5",
                                                           ifelse(grepl("THEOCA_45",System) & area=="suitable", "Cocoa RCP 4.5",
                                                                  ifelse(grepl("THEOCA_45",System) & area=="no_suitable", "Cocoa RCP 4.5",
                                                                         ifelse(grepl("THEOCA_85",System) & area=="suitable", "Cocoa RCP 8.5",
                                                                                ifelse(grepl("THEOCA_85",System) & area=="no_suitable", "Cocoa RCP 8.5",
                                                                                       NA)))))))),
                               levels = c("Coffee RCP 4.5","Coffee RCP 8.5","Cocoa RCP 4.5","Cocoa RCP 8.5")))

# add names and labels of focal species 
uses <- sort(unique(changes$main_use))

#..............................................
#..............................................
# Make charts ####

for (i in uses) {
  for(j in c("45","85")){
  df <- subset(changes, changes$main_use== i & grepl(j, changes$System))
  df$area <- factor(ifelse(df$area=="suitable","New/Remaining suitable","No longer suitable"),
                    levels=c("New/Remaining suitable","No longer suitable"))

  df %<>% 
    arrange(. , abundance)
  
  df$id <- as.integer(factor(df$acronym, levels = unique(df$acronym)))

  if (i != "N-fixing"){ lab <- tolower(i) }else { lab <- i}

  p <- ggplot(data = df) +
    geom_point(aes(x=cover,y=id), pch  = 21, size = 2, fill="grey",colour="black") +
    geom_segment(data=subset(df, df$change <= 0), 
                 mapping = aes(x=cover, xend = future_cover,
                               y=id, yend = id),
                 arrow=arrow(length = unit(0.1, "cm")),
                 size = 1, 
                 color = "#ca0020") +
    geom_segment(data = subset(df, df$change > 0), 
                 mapping = aes(x=cover, xend=future_cover,
                               y = id, yend = id),
                 arrow=arrow(length = unit(0.1, "cm")),
                 size=1, 
                 color= "#0571b0") +
    scale_y_continuous(breaks = as.integer(df$id), labels = df$species) +
    scale_x_continuous(breaks = seq(0,100, by=25), limits = c(0,100)) +
    labs(x = "Suitability (%)" ,
         y=NULL) +
    facet_grid(area ~ System_label) +
    theme(axis.text.y = element_text(size=15, angle = 0, 
                                     hjust=1,vjust=0.5, face="italic", colour = "black"),
          axis.text.x = element_text(size=15, angle = 0, 
                                     hjust=0.5,vjust=1, face="plain", colour = "black"),
          axis.title = element_text(size=14, face="bold"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          strip.text.x = element_text(size = 14, colour = "black"),
          strip.text.y = element_text(size=14, colour = "black"),
          strip.background = element_rect(colour="black", fill="#FFFFFF")) 
  
  
  ggsave(paste0(output,"wining_losing_RCP", j, "_", lab,"_species.svg"),
         plot = p, 
         dpi = 600, 
         width = 30,
         height = 40, 
         units = "cm")
  }
}


