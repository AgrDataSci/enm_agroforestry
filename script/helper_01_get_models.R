# Retrieve model outputs 
# and take AUC values from callibration
library("tidyverse")
library("magrittr")

sp <- "data/species_acronyms.csv"
sp %<>%
  read_csv() %>% 
  filter( main_use != "Crop")


spnames <- unique(sp$acronym)

result <- NULL

# get models 
for (i in seq_along(spnames)) {
  
  isp <- spnames[i]
  
  
  file <- paste0("processing/enm/", isp, "/outputs/", isp, "_output.txt")
  
  
  x <- read.delim(file)
  
  index <-  which(grepl("Results of ensemble.calibrate.weights|Results of ensemble.test.splits", x[,1]))
  
  
  x <- x[index:(index+29), ]
  
  x <- x[c(6:30)]
  
  x <- as.character(x)
  
  x <- cbind(sp = isp, x)
  
  result <- rbind(result, x)

}


write.csv(result, 
          "processing/enm/AUC_values.csv", 
          row.names = FALSE)


# which models were used
result2 <- NULL

for (i in seq_along(spnames)) {
  
  isp <- spnames[i]
  
  
  file <- paste0("processing/enm/", isp, "/outputs/", isp, "_output.txt")
  
  
  x <- read.delim(file)
  
  index <-  which(grepl("Ensemble weights based directly on input weights scaled", x[,1]))
  
  index <- index[length(index)]
  
  x <- x[index:(index+8), ]
  
  x <- as.character(x)
  
  x <- cbind(sp = isp, x)
  
  result2 <- rbind(result2, x)
  
}


write.csv(result2, 
          "processing/enm/models_used.csv",
          row.names = FALSE)

# clean it by hand and use the next script helper_02_auc_boxplot




