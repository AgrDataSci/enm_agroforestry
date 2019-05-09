# Make a boxplot of AUCs from model outputs
# this produces Fig. S7

library("tidyverse")
library("magrittr")
library("schoolmath")
library("svglite")


sp <- "data/species_acronyms.csv"
sp %<>%
  read_csv() %>% 
  filter( main_use != "Crop") %>% 
  select(1, 3) %>%
  rename(sp = acronym)

# AUC values
auc <- "processing/enm/AUC_values.csv" 
auc %<>% 
  read_csv(., na = c("0","NA")) %>% 
  select(1, 7) %>% 
  inner_join(., sp, by = "sp")

# models used 
mod <- read_csv("processing/enm/models_used.csv")

mod$sp <- as.factor(mod$sp)

mod <- split(mod, mod$sp)

mod <- lapply(mod, function(x){
  if (is.odd(nrow(x))) {
    x <- x[-nrow(x), ]
  }
  
  x <- unlist(x[2:ncol(x)])
  
  i <- seq(1, length(x), 2)
  
  x <- cbind(model = x[i], auc = x[-i])
  
  x <- x[as.numeric(x[,2]) > 0.05 & !is.na(x[,2]), 1]
  
  x <- as.vector(x)
  
})


auc <- split(auc, auc$sp)


for(i in seq_along(auc)){
  x <- auc[[i]]
  y <- mod[[i]]
  
  x <- x[x$model %in% y, ]
  
  x <- x[!is.na(x$MEAN.T), ]
  
  auc[[i]] <- x
  
}

auc <- bind_rows(auc)


p <- ggplot(auc, aes(species, MEAN.T)) +
  geom_boxplot() + 
  labs(x = NULL,
       y = "AUC") + 
  theme(axis.text.y = element_text(size = 12, angle = 0, 
                                   hjust = 1, vjust = 0.5, colour = "black"),
        axis.text.x = element_text(size = 12, angle = 70, 
                                   hjust = 1, vjust = 1, colour = "black"),
        axis.title = element_text(size = 12, face="bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 12, colour = "black"),
        strip.text.y = element_text(size=12, colour = "black"),
        strip.background = element_rect(colour="black", fill="#FFFFFF"))


ggsave("manuscript/FigS7.svg",
       plot = p,
       width = 35,
       height = 15,
       units = "cm")



