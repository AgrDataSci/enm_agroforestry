####################################################
###### SCRIPT TO MODEL CURRENT AND FUTURE NICHES OF 
###### PLANT SPECIES USING ENSEMBLE MODELLING 
# k.desousa (at) cgiar (dot) org
# Updated 28Jun2018
####################################################

###  LOAD PACKAGES   ######
if(!require(tidyverse)) install.packages("tidyverse") else library(tidyverse)
if(!require(magrittr)) install.packages("magrittr") else library(magrittr)
if(!require(sp)) install.packages("sp") else library(sp)
if(!require(dismo)) install.packages("dismo") else library(dismo)
if(!require(raster)) install.packages("raster") else library(raster)
if(!require(maptools)) install.packages("maptools") else library(maptools)
if(!require(rgeos)) install.packages("rgeos") else library(rgeos)
if(!require(rJava)) install.packages("rJava") else library(rJava)
if(!require(rgdal)) install.packages("rgdal") else library(rgdal)
if(!require(geosphere)) install.packages("geosphere") else library(geosphere)
if(!require(maps)) install.packages("maps") else  library(maps)
if(!require(alphahull)) install.packages("alphahull") else library(alphahull)
if(!require(tm)) install.packages("tm") else  library(tm)
if(!require(maxent)) install.packages("maxent") else library(maxent)
if(!require(gbm)) install.packages("gbm") else library(gbm)
if(!require(gam)) install.packages("gam") else library(gam)
if(!require(earth)) install.packages("earth") else library(earth)
if(!require(vegan)) install.packages("vegan") else library(vegan)
if(!require(MASS)) install.packages("MASS") else library(MASS)
if(!require(mgcv)) install.packages("mgcv") else library(mgcv)
if(!require(cluster)) install.packages("cluster") else library(cluster)
if(!require(rpart)) install.packages("rpart") else  library(rpart)
if(!require(effects)) install.packages("effects") else library(effects)
if(!require(multcomp)) install.packages("multcomp") else library(multcomp)
if(!require(ellipse)) install.packages("ellipse") else library(ellipse)
if(!require(maptree)) install.packages("maptree") else library(maptree)
if(!require(splancs)) install.packages("splancs") else library(splancs)
if(!require(spatial)) install.packages("spatial") else library(spatial)
if(!require(akima)) install.packages("akima") else library(akima)
if(!require(nnet)) install.packages("nnet") else library(nnet)
if(!require(randomForest)) install.packages("randomForest") else library(randomForest)
if(!require(mda)) install.packages("mda") else library(mda)
if(!require(kernlab)) install.packages("kernlab") else library(kernlab)
if(!require(e1071)) install.packages("e1071") else library(e1071)
if(!require(sem)) install.packages("sem") else library(sem)
if(!require(car)) install.packages("car") else library(car)
if(!require(maxlike)) install.packages("maxlike") else library(maxlike)
if(!require(glmnet)) install.packages("glmnet") else library(glmnet)
if(!require(PresenceAbsence)) install.packages("PresenceAbsence") else library(PresenceAbsence)
if(!require(tcltk2)) install.packages("tcltk2") else library(tcltk2)
if(!require(BiodiversityR)) install.packages("BiodiversityR") else  library(BiodiversityR)

# ================================================
# ================================================
## This scrip requires an external directory with bioclimatic layers 
## from WorldClim, you must download it via http://www.worldclim.org/version1
## and update the local directory in your R project
## we recommend enm_agroforestry/data/worldclim

# check if the files are available in the directory
if (length(list.files("data/worldclim/current")) == 1) {
  stop("Please download WorldClim data and add it to data/worldclim \n")
}

files <- list.dirs("data/worldclim/gcm")[-1]
gcm <- do.call(rbind, strsplit(files, "/"))
gcm <- gcm[,ncol(gcm)]
gcm <- gsub("bi50","",gcm)

# Create a vector of bioclim variables 
bio_names <- gsub(".tif","",list.files(files[1], pattern = ".tif$"))

# Define projection and extension
myproj <- "+proj=longlat +datum=WGS84"
mesoam <- raster::extent(-105, -75, 1, 25) 

# ================================================
# ================================================
# Prepare data  ####
# read occurrence data
sp <- read_csv("data/species_acronyms.csv")
df <- read_csv("data/occurrence_data.csv")

# remove duplicates in lat lon and add acronyms from sp
df %<>% 
  filter(!is.na(lon) & !is.na(lat)) %>%
  distinct(lon, lat, species, .keep_all = TRUE) %>%
  dplyr::select(acronym, lon, lat) %>%
  rename(x = lon, y = lat)

# ================================================
# ================================================
# Read current climate predictors variables #
files <- list.files("data/worldclim/current", 
                    pattern = ".bil$", 
                    full.names = TRUE)

bio_current <- raster::stack(files)

crs(bio_current) <- myproj

# crop bioclim layers into the predefined extension
bio_current %<>%
  raster::crop(., mesoam) %>%
  raster::stack(.)

# ================================================
# ================================================
# Variable selection with VIF ####

# # run once for all presence points
# # create background points across the full extension
# xy <- df %>% 
#   dplyr::select(x:y) %>%
#   distinct(x, y, .keep_all = TRUE)
# 
# xy <- xy[sample(row.names(xy), nrow(xy)*0.15), ]
# 
# largedist <- xy %>% 
#   raster::pointDistance(., longlat = FALSE) %>%
#   max(., na.rm = TRUE)
# 
# hull <- convHull(xy,lonlat=TRUE)
# 
# ext.hull <- gBuffer(hull@polygons, width=0.1*largedist)
# 
# bg1 <- spsample(ext.hull,10000, type="random")
# 
# bg <- gridSample(bg1, bio_current, n=1)
# 
# xy <- df %>% 
#   dplyr::select(x:y) %>%
#   as.matrix()
# 
# VIF.result <- ensemble.calibrate.models(x = bio_current, p=xy, an=bg,
#                             VIF = TRUE,
#                             MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, 
#                             GAMSTEP=0, MGCV=0, MGCVFIX=0,EARTH=0, RPART=0, NNET=0, FDA=0, 
#                             SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, models.save = TRUE)
# 
# #Predictors with VIF values > 10 are
# #eliminated from the analysis, see Rogerson (2001) Statistical Methods for Geography
# var.drops <- c("bio07","bio05","bio06")
# VIF.max <- 10 
# for (j in 1:length(names(bio_current))) {
#   
#   VIF.result <- ensemble.calibrate.models(x=bio_current, p=xy, an=10000,
#                               layer.drops=var.drops,
#                               VIF=TRUE,
#                               MAXENT=0, GBM=0, GBMSTEP=0, RF=0, GLM=0, GLMSTEP=0, GAM=0, 
#                               GAMSTEP=0, MGCV=0, MGCVFIX=0,EARTH=0, RPART=0, NNET=0, FDA=0, 
#                               SVM=0, SVME=0, BIOCLIM=0, DOMAIN=0, MAHAL=0, models.save = TRUE)
#   print(VIF.result$VIF)
#   if (VIF.result$VIF[1] > VIF.max) {var.drops <- c(var.drops, names(VIF.result$VIF)[1])}
# }
# 
# bio_vif1 <- as.data.frame(VIF.result$VIF)
# bio_vif  <- rownames(bio_vif1)

# Final seletion of variables for modelling
bio_current <- subset(bio_current, subset = bio_vif)

# ================================================
# ================================================
# Run ensemble modelling ####

species <- sort(unique(df$acronym))


for (i in seq_along(species) ) {
  
  cat("\n######## \n Ensemble modelling for", species[i], "\n Time:", date())
  
  output <- paste0("processing/enm/", species[i], "/")
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  # BiodiversityR saves outputs on the current working directory
  # we set this to the intended output directory
  patentwd <- getwd()
  
  setwd(output)
  
  # sampling data
  coord <- df[df$acronym == species[i], ]
  
  coord <- coord[c("x","y")]
  
  if(nrow(coord) > 5000){
    # calculate largest distance
    largedist <- coord[ sample(row.names(coord), 3000) ,] %>%
      raster::pointDistance(., longlat = FALSE) %>%
      max(., na.rm = TRUE)

  } else{
    # calculate largest distance
    largedist <- coord %>%
      raster::pointDistance(., longlat = FALSE) %>%
      max(., na.rm = TRUE)
    
  }
  
  # make a convex hull and remove duplicated coordinates in the same grid-cell
  hull <- convHull(coord, lonlat = TRUE)
  # extent convex hull
  ext_hull <- gBuffer(hull@polygons, width = 0.1 * largedist)
  crs(ext_hull) <- myproj
  # define raster
  r <- raster(ext_hull)
  # set the resolution of the cells to 2-5 minutes see Ranjitkar et al (2016)
  res(r) <- res(bio_current)
  
  coord %<>%
    as.matrix() %>% 
    as.data.frame() %>%
    dismo::gridSample(., r, n=1)
  
  n <- nrow(coord)
  
  # define threshold 
  # see Liu et al (2013) doi:10.1111/jbi.12058
  if(n >  150) thres <- "MaxSens+Spec"
  
  # Run ensemble modelling
  # step 1: 4-fold cross-validation
  cat("\n Step 1: calibrating ENM algorithms \n")
  
  enm_step1 <- ensemble.calibrate.weights(x = bio_current, p = coord, an = 1000,
                                          k=4,layer.drops=NULL, excludep = TRUE,
                                          SINK=TRUE, species.name= species[i],
                                          MAXENT=1, GBM=1, GBMSTEP=1, RF=1, 
                                          GLM=1, GLMSTEP=1, GAM=1, 
                                          GAMSTEP=1, MGCV=1, MGCVFIX=1, 
                                          EARTH=1, RPART=1, NNET=1, FDA=1, 
                                          SVM=1, SVME=1, BIOCLIM=1, 
                                          DOMAIN=1, MAHAL=1,
                                          ENSEMBLE.tune=TRUE, PROBIT=TRUE,
                                          threshold.method =  thres,
                                          threshold.sensitivity = 0.9,
                                          threshold.PresenceAbsence = TRUE, 
                                          ENSEMBLE.best=0, 
                                          ENSEMBLE.exponent=c(1, 2, 4, 6, 8),
                                          ENSEMBLE.min=0.7,
                                          Yweights="BIOMOD", factors=NULL,
                                          PLOTS=FALSE, formulae.defaults=TRUE,
                                          GBMSTEP.learning.rate=0.002)
  
  # step 2: create models that will be used for the raster predictions
  # models with input.weights <0.05 are excluded
  output_weights <- enm_step1$output.weights
  output_weights[output_weights < 0.05] <- 0

  cat("Step 2: model species distribution with selected ENM algorithms \n")
  
  enm_step2 <- ensemble.calibrate.models(x=bio_current, p=coord, 
                                         a=enm_step1$a,
                                         k=4, layer.drops=NULL,
                                         SINK=TRUE, species.name = species[i],
                                         models.keep=TRUE,
                                         input.weights=output_weights,
                                         threshold.method =  thres,
                                         threshold.sensitivity = 0.9,
                                         threshold.PresenceAbsence = TRUE, 
                                         ENSEMBLE.tune=FALSE, PROBIT=TRUE,
                                         Yweights="BIOMOD", factors=NULL,
                                         PLOTS=FALSE, formulae.defaults=TRUE,
                                         GBMSTEP.learning.rate=0.002,
                                         models.save = TRUE)

  
  cat("Step 3.1: Generate map of current distribution \n")
  #step3: use previously calibrated models to construct consensus layers
  ensemble_current <- ensemble.raster(xn = bio_current ,
                                      models.list = enm_step2$models,
                                      input.weights=output_weights,
                                      thresholds = enm_step2$models$thresholds,
                                      SINK=TRUE,
                                      RASTER.species.name = species[i], 
                                      RASTER.stack.name="bio_current")
  
  ### write raster for each gcm model in RCP45 and RCP85 
  for (k in seq_along(gcm)){
    cat("Step 3.2: Predict future distribution, GCM", toupper(gcm[k]), "\n")
    
    #load GCM layers
    gcmfiles <- list.files(paste0(parentwd, "/data/worldclim/gcm/", gcm[k], "bi50/"), 
                           pattern = ".tif$",
                           full.names = TRUE)
    
    gcm_model <- raster::stack(gcmfiles)
    
    gcm_model <- subset(gcm_model, subset = bio_vif)
    
    crs(gcm_model) <- myproj
    
    gcm_model <- crop(gcm_model, mesoam)
    
    gcm_model <- stack(gcm_model)

    ensemble_gcm_model <- ensemble.raster(xn = gcm_model,
                                          models.list = enm_step2$models,
                                          input.weights = output_weights,
                                          thresholds = enm_step2$models$thresholds,
                                          SINK = TRUE,
                                          RASTER.species.name = species[i], 
                                          RASTER.stack.name = gcm[k])
    
  }
  
  #return to parent directory
  setwd(parentwd)
  
}

#END





