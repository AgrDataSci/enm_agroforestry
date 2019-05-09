####################################################
# #### CALCULATE THE CHANGES IN SUITABILITY OF FOCAL SPECIES 
# IN FUTURE SUITABLE AREAS FOR COFFEE AND COCOA
# This script generates the dataframe for Fig. 4

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
crop_names <- c("COFFAR","THEOCA")

crop <- list()

for (i in seq_along(crop_names)) {
      
  r <- raster(paste0("processing/layers_comb/",crop_names[i],"_baseline.grd"))
      
  r[r[] != 1] <- NA 
  crop[[crop_names[i]]] <- r

}

rm(r)

names(crop)
changes <- NULL


#..............................................
#..............................................
# Read changes of focal species ####
# this loop calculate the frequencies of suitability of each 
# focal species using the mask generated above.

for (i in species_names) {
  for (j in rcp) {
    
    if (j=="current") {
      r <- raster(paste("processing/species_sets/", i, "/", i, "_presence.tif",sep="" )) 
      }else{
        r <- raster(paste("processing/species_sets/", i, "/", i, "_", j ,"_change.tif",sep="" ))
        }
    
    for (k in crop_names) {
      
      land <- mask(r, crop[[k]])
      frequencies <- data.frame(array(dim = c(1, 11)))
      names(frequencies) <- c("acronym","system","never_suitable(0)", "remains_suitable(1)", 
                              "no_longer_suitable(-1)", "new_habitat(2)","area_km2","remains_suitable_(1)_perc",
                              "no_longer_suitable(-1)_perc", "new_habitat_(2)_per","net_change_perc")
      freq1 <- as.data.frame(raster::freq(land, useNA="no"))
      freq1 <- freq1[is.na(freq1[, 1]) == F, ]
      
      frequencies[1,1] <- i
      frequencies[1,2] <- paste(k,j,sep="_")
      
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
      
      #percent of always suitable areas
      frequencies[8] <- frequencies[,4] / sum(frequencies[,c(4,5)], na.rm = TRUE) * 100
      #percent of no longer suitable areas
      frequencies[9] <- frequencies[,5] / sum(frequencies[,c(4,5)], na.rm = TRUE) * 100
      #percent of new suitable areas
      frequencies[10] <- frequencies[,6] / sum(frequencies[,c(4,5)], na.rm = TRUE) * 100
      #percent of net change 
      frequencies[11] <- frequencies[,10] - frequencies[,9]
      #add to changes dataframe
      changes <- rbind(changes, frequencies)
      
    }
  }
}

changes %<>% 
  inner_join(., sp, by="acronym") %>% 
  filter(. , !grepl("current", system)) 

# export dataframe
output <- "output/wining_losing_sp/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

write_csv(changes, paste0(output, "wining_losing_sp.csv"))

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

# make labels of systems (land use)
changes %<>% 
  mutate(System_label = factor(ifelse(grepl("COFFAR_45", System), "Coffee RCP 4.5",
                                      ifelse(grepl("COFFAR_85", System), "Coffee RCP 8.5",
                                             ifelse(grepl("THEOCA_45", System), "Cocoa RCP 4.5",
                                                    ifelse(grepl("THEOCA_85", System), "Cocoa RCP 8.5",
                                                           NA)))),
                               levels = c("Coffee RCP 4.5","Coffee RCP 8.5","Cocoa RCP 4.5","Cocoa RCP 8.5")))

changes

# add names and labels of focal species 
uses <- sort(unique(changes$main_use))

#..............................................
#..............................................
# Make charts ####

for(j in c(1:2)){
  
  plots <- list()
  
  for (i in seq_along(uses)) {
  
    if(j == 1) {r <- 45} else {r <- 85}
    
    df <- subset(changes, changes$main_use== uses[i] & grepl(r, changes$System))
    df$area <- " "
    head(df)
    
    df %<>% arrange(. , abundance)
    
    df$id <- as.integer(factor(df$acronym, levels = unique(df$acronym)))
    
    if (uses[i] != "N-fixing"){ lab <- tolower(uses[i]) } else { lab <- uses[i]}
    
    p <- ggplot(data=df) +
      geom_point(aes(x=cover,y=id), pch  = 21, size=2, fill="grey",colour="black") +
      geom_segment(data=subset(df, df$change <= 0), 
                   mapping = aes(x=cover, xend=future_cover,
                                 y=id, yend=id),
                   arrow=arrow(length = unit(0.1, "cm")),
                   size=1, 
                   color= "#ca0020") +
      geom_segment(data=subset(df, df$change > 0), 
                   mapping = aes(x=cover, xend=future_cover,
                                 y=id, yend=id),
                   arrow=arrow(length = unit(0.1, "cm")),
                   size=1, 
                   color= "#0571b0") +
      scale_y_continuous(breaks = as.integer(df$id), labels = df$species) +
      scale_x_continuous(breaks = seq(0,100, by=25), limits = c(0,100)) +
      labs(x = paste0("Suitability (%) of ", lab, " trees"),
           y = NULL) +
      facet_grid(area ~ System_label) +
      theme(axis.text.y = element_text(size = 12, angle = 0, hjust = 1, 
                                       vjust = 0.5, face="italic", colour = "black"),
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, 
                                       vjust = 1, face="plain", colour = "black"),
            axis.title = element_text(size = 12, face="bold"),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            plot.background = element_blank(),
            panel.background = element_blank(),
            strip.text.x = element_text(size = 12, colour = "black"),
            strip.text.y = element_text(size=12, colour = "black"),
            strip.background = element_rect(colour="black", fill="#FFFFFF")) 
    
    
    plots[[i]] <- p
  
  }
  
  p <- grid.arrange(
    plots[[1]],
    plots[[2]],
    plots[[3]],
    nrow = 1)
  
  ggsave(paste0(output,"wining_losing_RCP", r, ".svg"),
         plot = p,
         width = 75,
         height = 22,
         units = "cm")

}



