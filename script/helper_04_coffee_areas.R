# Calculate changes in suitability of cocoa and coffee 
# over the elevation gradient

library("tidyverse")
library("magrittr")

df <- read_csv("output/change_suitability/si_change_suitability.csv")

df

# replace coded values in layer by their meaning
df %<>% 
  rename(elevation = MA_elev30as) %>% 
  filter(!is.na(elevation)) %>%
  mutate(layer = ifelse(layer == -1, "notsuitable",
                        ifelse(layer == 1, "suitable",
                               ifelse(layer == 2, "newarea", NA))))

# create elevation classes
df %<>%
  mutate(elev_strata = ifelse(elevation < 400, "0-400", 
                              ifelse(elevation >= 400 & elevation < 800, "400-800",
                                     ifelse(elevation >= 800 & elevation < 1200, "800-1200",
                                            ifelse(elevation >= 1200 & elevation < 1600, "1200-1600",
                                                   ifelse(elevation >= 1600, "1600-2500", "0-400"))))))

df %<>%
  group_by(crop, layer, elev_strata) %>%
  summarise(pixels = length(layer))


write_csv(df, 
          "output/change_suitability/count_areas.csv")
