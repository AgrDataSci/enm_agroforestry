# Calculate areas where cocoa replace coffee
# under the elevation gradient

library("tidyverse")
library("magrittr")

df <- read_csv("output/cocoa_replace_coffee/si_cocoa_replace_coffee.csv")

df

# replace coded values in layer by their meaning
df %<>% 
  rename(elevation = MA_elev30as) %>% 
  filter(!is.na(elevation))

# create elevation classes
df %<>%
  mutate(elev_strata = ifelse(elevation < 400, "0-400", 
                              ifelse(elevation >= 400 & elevation < 800, "400-800",
                                     ifelse(elevation >= 800 & elevation < 1200, "800-1200",
                                            ifelse(elevation >= 1200 & elevation < 1600, "1200-1600",
                                                   ifelse(elevation >= 1600, "1600-2500", "0-400"))))))


df

df %<>%
  group_by(label, rcp, elev_strata) %>%
  summarise(pixels = length(label))


df

write_csv(df, 
          "output/cocoa_replace_coffee/count_areas_cocoa_replace_coffee.csv")
