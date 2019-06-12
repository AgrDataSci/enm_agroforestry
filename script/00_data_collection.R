#Data collection for distribution analysis 
#Author: Maarten van Zonneveld
#Affiliation: World Vegetable Center

#PACKAGES####
library(rgbif)
library(maptools)
library(rgeos)
library(raster)
library(dismo)
library(geosphere)
library(rworldmap)
library(ggplot2)


output <- "data/raw/"
#DATA INPUT####

#species information
setwd("C:/Users/CR-C003/Dropbox/Agroforestry tree niche modelling")

#enter data from csv file, sometimes seperated by ";", sometimes by ","
basic1 <- read.csv("additional_agf_species.csv",header=T, sep=",")
basic2 <- read.csv("additional_agf_species.csv",header=T, sep=";")
ifelse(length(basic1) > length(basic2),
       basic <- basic1,
       basic <- basic2)
taxa <- as.character(basic[,"species"])
taxa <- unique(taxa)


#COUNTRIES
setwd("C:/Users/CR-C003/Dropbox/World")
ctries <- readOGR("Countries.shp", "Countries")
proj4string(ctries) <- "+proj=longlat +datum=WGS84"
ctries_buffer <- readOGR("ctries_buffer.shp","ctries_buffer")

#LATIN AMERICA
#selected species
setwd("C:/Users/CR-C003/Dropbox/Fruits megafauna")
#enter country layer 
ctries<-readOGR(dsn=input, layer="lac")

#points to assign coordinates to records in sea near coast
country_coord <- ctries_points[,1:2]
country_coord <- SpatialPoints(country_coord)


#ISO
setwd("C:/Users/CR-C003/Dropbox/World")
iso <- read.csv("ISO.csv",sep=";")

#BIOCLIM
setwd("F:/MAPFORGEN/input/climate/bioclim")
elayers <- stack(list.files(pattern = ".asc"))
crs(elayers) <- "+proj=longlat +datum=WGS84"


#FUNCTIONS####

#CHECK DECIMAL PLACES https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r/5173906
decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


#DATA COLLECTION####

# Now we can use that taxonconceptkey to get occurrence data for spp and its synonyms
coord <- data.frame()
archive <- list()
for (i in 1:length(taxa)){
  print(taxa[i])
  
#check for all scientific names that include the species names  
  name <- as.data.frame(name_suggest(q = taxa[i]))
  coord2 <- data.frame()
  spec2 <- data.frame()
  if (nrow(name) > 0 ){ print("GO")
  for (j in 1:length(name[,1])){  
    tck <- occ_search(taxonKey = name[j,1], 
                      hasCoordinate = TRUE, 
                      fields = c("name","scientificName","countryCode","decimalLongitude","decimalLatitude","key","basisOfRecord","publishingOrg"),limit = 30000)
    if(length(tck$data) == 0) {"NA"} else 
    {tck$data[,"name"] <- taxa[i]
    coord2 <- rbind(coord2,tck$data)
    archive[[i]] <- occ_search(taxonKey = name[j,1])}
  }}
  coord <- rbind(coord,coord2)
 }

#Make species matrix
species <- as.data.frame(unique(coord[,c("name")]))
species_list <- matrix(NA,nrow=length(species[,1]),ncol=2)

for (i in 1: length(species[,1])){
  coord1 <- subset(coord,coord[,c("name")] == as.character(species[i,]))
  coord1 <- as.data.frame(coord1)
  species_list[i,] <- c(as.character(species[i,]),length(coord1[,1]))
}

coord2 <- coord 

#DATA CLEANING####

#Basic clean up
coord2$decimalLongitude <- as.numeric(as.character(coord2$decimalLongitude))
coord2$decimalLatitude <- as.numeric(as.character(coord2$decimalLatitude))
basic <- coord2[!is.na(coord2$decimalLatitude),]
basic <- basic[!is.na(basic$decimalLongitude),]
basic <- subset(basic, !countryCode == "NONE")
basic <- as.data.frame(basic)

#select countries with location records
cntr2 <- as.character(unique(basic[,"countryCode"]))
cntr3 <- data.frame()
for (i in 1:length(cntr2)){
  cn4 <- subset(iso,iso[,2] == cntr2[i])
  cntr3 <- rbind(cntr3,cn4)}

#check per country points outside country buffer
alldata <- data.frame()
for (i in 1: length(cntr2)){
  #select polygon per country
  cntry_buffer  <- subset(ctries_buffer, ctries_buffer@data[,2] == as.character(cntr3[i,4]))
  if (length(cntry_buffer@polygons) == 0){"NA"} else 
  #select presence points per country
  {cntr_data <- subset(basic, basic[,"countryCode"] == as.character(cntr3[i,2]))
  if (nrow(cntr_data) == 0) {"NA"} else {
    coord <- SpatialPoints(cntr_data[,c("decimalLongitude","decimalLatitude")])
    cntr_data1 <- over(coord,cntry_buffer)
    cntr_data2 <- cbind(cntr_data,is.na(cntr_data1[,2]))
    cntr_data3 <- subset(cntr_data2, cntr_data2[,8] == FALSE)
    alldata <- rbind(alldata,cntr_data3[,1:7])}}
}

species <- as.data.frame(unique(alldata[,c("name")]))
species_list <- matrix(NA,nrow=length(species[,1]),ncol=2)

for (i in 1: length(species[,1])){
  coord1 <- subset(alldata,alldata[,c("name")] == as.character(species[i,]))
  coord1 <- as.data.frame(coord1)
  species_list[i,] <- c(as.character(species[i,]),length(coord1[,1]))
}

#REMOVE DUPLICATE COORDINATES
alldata_nondup <- data.frame()
for (i in 1:length(species[,1])){
  spp <- subset(alldata, name == as.character(species[i,]))
  dup <- cbind(spp,duplicated(spp[,c("decimalLongitude","decimalLatitude")]))
  non_dup <- subset(dup, dup[,8] == FALSE)
  alldata_nondup <- rbind(alldata_nondup,non_dup[,1:7])}

alldata <- alldata_nondup

#Remove coordinates with ZERO decimals in both LAT and LONG 
basic1 <- data.frame()
for (i in 1:length(alldata[,1])){
  print(i)
  ifelse(decimalplaces(alldata[i,"decimalLatitude"]) + 
           decimalplaces(alldata[i,"decimalLongitude"]) >0,
         basic1 <- rbind(basic1,alldata[i,]),
         print("NA"))
}

alldata <- basic1 

#REVISE POINTS IN SEA AND CORRECT POINTS WITHIN 10 ARC MIN FROM COASTAL BORDER
#///idenitfy data points outside polygon and update coordinates for points near polygon\\\\\ AUTHOR: Maarten van Zonneveld

#identify points outside continent polygon in ocean
coords <- SpatialPoints(alldata[,c("decimalLongitude","decimalLatitude")],proj4string =CRS("+proj=longlat +datum=WGS84"))
#points over countries
sea_true <- over(coords, ctries)
#combine passportdata with colomun of points in/out continent
alldata3 <- cbind(alldata,is.na(sea_true[,3]))
#select points in sea
alldata4 <- subset(alldata3, alldata3[,8] == TRUE)
#convert to spatial points
sea_points <- SpatialPoints(alldata4[,c("decimalLongitude","decimalLatitude")])

#create 10 minutes buffer around sea points
sea_buffer <- gBuffer(sea_points, byid=TRUE, width=1/6)
#project sea points
proj4string(sea_buffer) <- "+proj=longlat +datum=WGS84"
proj4string(country_coord) <- "+proj=longlat +datum=WGS84"
#la border points over buffer to select close by points
la_true <- over(country_coord, sea_buffer)
#remove la border points outside buffers
la_true2 <- na.omit(cbind(coordinates(country_coord),la_true))
#project la points
la_points2 <- SpatialPoints(la_true2[,1:2])

#identify distance to nearest vetrix points of polygon and coordinates of that vetrix
sea_points2 <- matrix(NA, ncol = 5, nrow = length(alldata4[,1]))
for(i in 1:length(alldata4[,1])){
  print(i)
  m <- pointDistance(sea_points[i], la_points2, lonlat=TRUE)
  la_points3 <- cbind(la_true2 ,m)
  dups <- duplicated(m)
  la_points3 <-na.omit(la_points3[!dups,])
  k <- as.numeric(subset(la_points3, la_points3[,4] == min(la_points3[,4])))
  sea_points2[i,] <- c(i,k)
}

#update coordinates of points outside polygons 
alldata4[,c("decimalLongitude","decimalLatitude")] <- sea_points2[,2:3]
#Select points near vetrix less then 10 minutes (18000m)
alldata5 <- alldata4[sea_points2[,5] < 18000,]

#update dataset
#subset initially correct data points
alldata6 <- subset(alldata3, alldata3[,8] == FALSE)
#bind initially correct data points with updated coastal points
alldata2 <- rbind(alldata6,alldata5)[1:7]
alldata <- data.frame(alldata2)


#REMOVE COUNTRY MIDDLE POINTS
#Brazil
check <- subset(alldata, !decimalLongitude == -52.8731 & !decimalLatitude == -10.8339)
check <- subset(check, !decimalLongitude == -53.08972 & !decimalLatitude == -10.77306)

#Paraguay
check <- subset(check, !decimalLongitude == -58.1689 & !decimalLatitude == -23.2025)

#El Salvador
check <- subset(check, !decimalLongitude == -88.866511 & !decimalLatitude == 13.736897)

#French guyana
check <- subset(check, !decimalLongitude == -52.9708 & !decimalLatitude == 3.8558)

#Surinam
check <- subset(check, !decimalLongitude == -55.625 & !decimalLatitude == 4.1006)

#Guyana
check <- subset(check, !decimalLongitude == -59 & !decimalLatitude == 5)

#Belize
check <- subset(check, !decimalLongitude ==  -88.684767 & !decimalLatitude == 17.217668)

#Peru
check <- subset(check, !decimalLongitude == -74.14 & !decimalLatitude == -9.1839)

#Colombia
check <- subset(check, !decimalLongitude == -72.8667 & !decimalLatitude == 3.8811)

#Bolivia
check <- subset(check, !decimalLongitude == -64.435 & !decimalLatitude == -16.7261)

#Honduras
check <- subset(check, !decimalLongitude == -86.619143 & !decimalLatitude == 14.819222)

#Costa Rica
check <- subset(check, !decimalLongitude == -83.9519 & !decimalLatitude == 10.0114)

#Guatemala
check <- subset(check, !decimalLongitude == -90.1792 & !decimalLatitude == 15.7422)

#Nicaragua
check <- subset(check, !decimalLongitude == -85.516667 & !decimalLatitude == 12.766667)

#Panama
check <- subset(check, !decimalLongitude == -79.83638 & !decimalLatitude == 9.16361)

#Venezuela
check <- subset(check, !decimalLongitude == -65.9119 & !decimalLatitude == 7.0758)

#Argentina
check <- subset(check, !decimalLongitude == -64.9208 & !decimalLatitude == -35.3869)
check <- subset(check, !decimalLongitude == -63.5848 & !decimalLatitude == -38.4213)

#Mexico
check <- subset(check, !decimalLongitude == -102 & !decimalLatitude == 23)

#Russia
check <- subset(check, !decimalLongitude == 108.9844 & !decimalLatitude == 61.6064)

#France
check <- subset(check, !decimalLongitude == 2.21954 & !decimalLatitude == 46.25869)

#Liberia
check <- subset(check, !decimalLongitude == -9.5 & !decimalLatitude == 6.4)

#Angola
check <- subset(check, !decimalLongitude == 18.5 & !decimalLatitude == -12.5)

#Nepal
check <- subset(check, !decimalLongitude == 84.122776 & !decimalLatitude == 28.39103)
check <- subset(check, !decimalLongitude == 84.1228 & !decimalLatitude == 28.391)
check <- subset(check, !decimalLongitude == 84.12 & !decimalLatitude == 28.39)

#China
check <- subset(check, !decimalLongitude == 104.165813 & !decimalLatitude == 35.860085)
check <- subset(check, !decimalLongitude == 104.1658 & !decimalLatitude == 35.8601)

#Spain
check <- subset(check, !decimalLongitude == -4.1 & !decimalLatitude == 40.68)

#Iraq
check <- subset(check, !decimalLongitude == 43.68561 & !decimalLatitude == 33.240547)
check <- subset(check, !decimalLongitude == 43.6856 & !decimalLatitude == 33.24054)

#Myanmar
check <- subset(check, !decimalLongitude == 96.0667 & !decimalLatitude == 21.95)

cleandata <- check


#ONLY POINTS FROM LATIN AMERICA
points <- SpatialPoints(cleandata[,c("decimalLongitude","decimalLatitude")],
                        proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")) 
cl <- over(points, ctries)
cntr_data2 <- cbind(cleandata,is.na(cl[,3]))
cleandata_lac <- subset(cntr_data2, cntr_data2[,8] == FALSE)[,1:7]

cleandata <- cleandata_lac

#WRITE ARCHIVES
setwd("C:/Users/CR-C003/Dropbox/Agroforestry tree niche modelling")
write.csv(as.data.frame(table(cleandata[,"name"])),"agf_spp_summary.csv",row.names = F)
write.csv(cleandata,"agf_spp_records.csv",row.names = F)




#


sp <- read_csv("./input/passport_data/species_acronyms.csv")
df <-read_csv("./input/passport_data/passport_data.csv")

#remove duplicates in lat lon
df %<>% 
  filter(!is.na(lon) & !is.na(lat)) %>%
  distinct(lon, lat, species, .keep_all = TRUE)

#count number of points per species and add to 'sp'
sp <- df %>% 
  group_by(species) %>%  
  count(species) %>%
  filter(n > 60) %>%
  inner_join(sp, . , by = "species", all.x = TRUE)

write_csv(sp, "./input/passport_data/species_acronyms.csv")
