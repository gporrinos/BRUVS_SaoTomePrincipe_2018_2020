
################################################################################
################################################################################
###                            VARIABLE DEFINITION                           ###
################################################################################
################################################################################












          #############################################################
          ###                       LOAD DATA                       ###
          #############################################################


### STEP 1: LOAD DATABASES
source(file.path(getwd(),"data_analysis/script/01_data_and_directories.R"))



### STEP 2: EXCLUDE DISCARDED DEPLOYMENTS
bruvs   <-  bruvs[!bruvs$sampling_time=="afternoon_repeat",]    ## Exclude afternoon repeats (they were tests and were not used for analysis)
bruvs   <-  bruvs[!bruvs$sampling_time=="not_valid",]








         #############################################################
         ###             LOAD/INSTALL REQUIRED PACKAGES            ###
         #############################################################


### INSTALL / UPDATE PACKAGES ###
if(!require("sf")) install.packages("sf",   repos = "https://cloud.r-project.org")
if(!require("geosphere")) install.packages("geosphere",   repos = "https://cloud.r-project.org")
if(!require("terra")) install.packages("terra",   repos = "https://cloud.r-project.org")







         #############################################################
         ###         CALCULATE ENVIRONMENTAL VARIABLES             ###
         #############################################################



                 #--------------------------------------------#
                 ####---------- Distance to shore ---------####
                 #--------------------------------------------#


#### STEP 1: EXTRACT COORDINATES OF DEPLOYMENTS
deployments              <- data.frame(Longitude= bruvs$lon,
                                       Latitude = bruvs$lat)

#### STEP 2: IMPORT STP'S SHAPEFILE
STP                      <- st_read(dsn   = file.path(getwd(), "geodata"),    layer = "STP boundaries")

#### STEP 3: CONVERT STP COASTLINE TO POINTS AND EXTRACT COORDINATES
STPcoords   <- as.data.frame(st_coordinates(st_sample(st_cast(STP, "LINESTRING"), size = 100000, type = "regular")))
STPcoords   <- data.frame(Longitude = STPcoords$X, Latitude = STPcoords$Y)

#### STEP 4: CALCULATE DISTANCE OFFSHORE
bruvs$disttoshore <- 
  unlist(lapply(c(1:nrow(deployments)), function(x){
    deployment   = deployments[x,]
    module       = ((STPcoords[,1]-deployment[,1])^2 + 
                      (STPcoords[,2]-deployment[,2])^2)^(1/2)
    nearestpoint = STPcoords[which(module==min(module)),][1,]
    disttoshore  = distVincentyEllipsoid (p1 = deployment, p2 = nearestpoint)
    return(disttoshore)
  }))




                 #--------------------------------------------#
                 ####---------------- Slope --------------#####
                 #--------------------------------------------#

deployments    <- data.frame(Longitude= bruvs$lon,
                             Latitude = bruvs$lat)

depth          <- rast(file.path(getwd(),"geodata","STP_depth.TIF",sep=""))
slope          <- terrain(depth, "slope",unit = "degrees", neighbors=8)
bruvs$slope    <- extract(slope,deployments)$slope





                 #--------------------------------------------#
                 #####-------------- Season --------------#####
                 #--------------------------------------------#

bruvs$season <- substr(bruvs$date,4,5)
bruvs$season[bruvs$season %in% c("12","01","02","03","04")] <- "Summer"
bruvs$season[bruvs$season %in% c("06","07","08","09")]      <- "Gravana"

bruvs$month <- substr(bruvs$date,4,5)





                 #--------------------------------------------#
                 #####----- Season 2 (ordinal date) ------#####
                 #--------------------------------------------#


startyear            <- c("2017-01-01","2018-01-01","2019-01-01","2020-01-01","2021-01-01","2022-01-01","2023-01-01")
startyear            <- as.integer(as.POSIXct(startyear,format = c("%Y-%m-%d"),tz="UTC"))/(24*60*60)
bruvs$dateoriginal <- bruvs$date

# Reclassify date as an integer
bruvs$date        <- as.integer(as.POSIXct(bruvs$date,format = c("%d/%m/%Y"),tz="UTC"))/(24*60*60)

# Calculate ordinal date
bruvs$Season      <- unlist(lapply(bruvs$date,function(x){x-max(startyear[!startyear>x])}))



rm(STP,STPcoords,deployments,depth,slope)





          #############################################################
          ###    CALCULATE DIVERSITY INDICATORS AND SPECIES MAXN    ###
          #############################################################



                 #--------------------------------------------#
                 ###        Install and load packages       ###
                 #--------------------------------------------#


if(!require("vegan")) install.packages("vegan",   repos = "https://cloud.r-project.org")




                 #--------------------------------------------#
                 ####          Clean MaxN database         ####
                 #--------------------------------------------#


### STEP 1: MAXN  |  REMOVE OBSERVATIONS OVER 60 MINUTES 
maxn <- maxn[!maxn$minutes>=60,]



### STEP 2: RECLASSIFY SPECIES NAMES + ADD CO-VARIATES AND FAMILY INFORMATION
sp <- match(maxn$species,species$species_original)
maxn$sp          <- unlist(lapply(sp,function(x){return(species$species[x])}))
maxn$speciessimple <- unlist(lapply(sp,function(x){return(species$speciessimple[x])}))
maxn$specieslist <- unlist(lapply(sp,function(x){return(species$specieslist[x])}))
maxn$type        <- unlist(lapply(sp,function(x){return(species$type[x])}))
maxn$family      <- unlist(lapply(sp,function(x){return(species$family[x])}))

video <- match(maxn$video,bruvs$video)
maxn$habitat <- unlist(lapply(video,function(x){return(bruvs$habitat[x])}))
maxn$island  <- unlist(lapply(video,function(x){return(bruvs$island[x])}))
maxn        <- maxn[!is.na(maxn$island),]



### STEP 3: KEEP ONLY FISH OBSERVATIONS
maxn <- rbind(maxn[maxn$type=="Teleost",],maxn[maxn$type=="Elasmobranch",])
maxn <- maxn[!is.na(maxn$maxn),]






                 #--------------------------------------------#
                 ###        Create max MaxN database        ###
                 #--------------------------------------------#

### STEP 1: CREATE NEW VARIABLES
mxMaxN                   <- maxn
mxMaxN$videospecies      <- paste(mxMaxN$video,mxMaxN$species,sep="")
mxMaxN$videosp           <- paste(mxMaxN$video,mxMaxN$sp,sep="") 
videospecies             <- levels(as.factor(mxMaxN$videospecies))
videosp                  <- levels(as.factor(mxMaxN$videosp))


### STEP 2: KEEP ONLY MAXIMUM OBSERVATION OF EACH SPECIES
temp  <-  mxMaxN[0,]
for(i in 1:length(videospecies)){
  x = mxMaxN[mxMaxN$videospecies==videospecies[i],]
  x = x[!is.na(x$maxn),]
  z = x[which(x$maxn == max(x$maxn)),]
  temp=rbind(temp,z)}
mxMaxN <- temp


### STEP 3: SUM RECORDS OF THE SAME SPECIES (juvenile/adult & male/female)
temp  <-  mxMaxN[0,]
for(i in 1:length(videosp)){
  x = mxMaxN[mxMaxN$videosp==videosp[i],]
  x = x[!is.na(x$maxn),]
  y = sum(x$maxn)
  x = x[1,]
  x$maxn[1] = y
  temp=rbind(temp,x)}
mxMaxN <- temp
mxMaxN <- mxMaxN[,-which(colnames(mxMaxN)=="videospecies")] 
mxMaxN <- mxMaxN[,-which(colnames(mxMaxN)=="videosp")] 


### STEP 3: REMOVE UNIDENTIFIED SPECIES
mxMaxN <- mxMaxN[-which(is.na(mxMaxN$sp)),] 





                 #--------------------------------------------#
                 ###    Calculate diversity indicators      ###
                 #--------------------------------------------#


#### STEP 1: SPECIES RICHNESS
sites <- bruvs$video
bruvs$richness <- unlist(lapply(sites,function(x){
  return(length(levels(as.factor(mxMaxN[mxMaxN$video==x,]$sp))))}))

#### STEP 2: ABUNDANCE (MaxN)
bruvs$abundance <- unlist(lapply(sites,function(x){
  return(sum(mxMaxN[mxMaxN$video==x,]$maxn))}))



#### STEP 3: CREATE SPECIES VS SITES MATRIX
divdata <- levels(as.factor(mxMaxN$sp))
temp <- data.frame(bruvs$video)
for(i in 1:length(divdata)){
  y = unlist(lapply(bruvs$video, function(x){
    site = mxMaxN[mxMaxN$video==x,]
    site = site[!is.na(site$maxn),]
    spec = site[site$sp==divdata[i],
                which(colnames(site)=="maxn")]
    spec = spec[!is.na(spec)]
    return(sum(spec))
  }))
  temp =data.frame(temp,y)}
colnames(temp)<-c("site",divdata)
rownames(temp)<-temp$site
divdata <- temp[,-1]


#### STEP 4: CALCULATE E1/D EVENNESS INDEX
bruvs$evenness <- as.numeric(vegan::diversity(divdata,index="invsimpson")) / bruvs$richness
bruvs$evenness[which(bruvs$richness==0)] <- as.numeric(NA)
bruvs$evenness[which(bruvs$richness==1)] <- as.numeric(NA)



rm(video,videosp,videospecies,y,sp,sites,i,x,z,temp)



