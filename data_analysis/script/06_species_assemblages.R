
################################################################################
################################################################################
###                           SPECIES ASSEMBLAGES                            ###
################################################################################
################################################################################








#          #############################################################
#          ###                 LOAD DATA AND PACKAGES                ###
#          #############################################################






### LOAD DATA
source(file.path(getwd(),"data_analysis/script/02_variable_definition.R"))


### Statistical analyses
if(!require("vegan"))   install.packages("vegan",   repos = "https://cloud.r-project.org")
if(!require("remotes")) install.packages("remotes",   repos = "https://cloud.r-project.org")
remotes::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

### Data management
if(!require("dplyr")) install.packages("dplyr",   repos = "https://cloud.r-project.org")

### Plot things
if(!require("ggplot2")) install.packages("ggplot2",   repos = "https://cloud.r-project.org")
if(!require("ggpubr")) install.packages("ggpubr",   repos = "https://cloud.r-project.org")
if(!require("ggh4x")) install.packages("ggh4x",   repos = "https://cloud.r-project.org")




#          #############################################################
#          ###                 CREATE DISTANCE MATRIX                ###
#          #############################################################



#                #--------------------------------------------#
#                ####               Functions              ####
#                #--------------------------------------------#

### FUNCTION TO TAG ROWS WITH NAs AND ALL-ZERO
keep <- function(divdata){
  # Tag rows that have NA with a 1, else 0)
  NArows  = unlist(lapply(c(1:nrow(divdata)),function(i){
    sum(unlist(lapply(divdata[i,],function(x){
      if(is.na(x)){1} else {0}})))
  }))
  divdata[is.na(divdata)] = 0
  # Tag all-zero rows
  allzero = unlist(lapply(c(1:nrow(divdata)),function(i){
    if(sum(divdata[i,])==0) {1} else {0}
  }))
  # Remove NAs and allzero
  keeps <- NArows + allzero == 0
  return(keeps)}



#--------------------------------------------#
############ Create species matrix ###########
#--------------------------------------------#

#### STEP 1: CREATE 1 COLUMN PER SPECIES
sp <- na.omit(mxMaxN[mxMaxN$type == "Teleost",])
sp <- levels(as.factor(sp$speciessimple))
temp <- data.frame(bruvs$video)
for(i in 1:length(sp)){
  y = unlist(lapply(bruvs$video, function(x){
    site = mxMaxN[mxMaxN$video==x,]
    site = site[!is.na(site$maxn),]
    spec = site[site$sp==sp[i],
                which(colnames(site)=="maxn")]
    spec = spec[!is.na(spec)]
    return(sum(spec))
  }))
  temp =data.frame(temp,y)}
colnames(temp)<-c("site",sp)
rownames(temp)<-temp$site
raw.spMaxN <- temp[,-1]
divdata <- temp[,-1]^(1/4)

#### STEP 2: REMOVE NAs AND ALL-ZERO ROWS (vegan package does not accept them)
keeps       <- keep(divdata)
divdata    <- divdata[keeps == TRUE,]
which(is.na(divdata))

### STEP 3: CREATE VECTORS WITH CO-VARIATES
habitat     <- as.factor(bruvs$habitat)[keeps == TRUE]
island      <- as.factor(bruvs$island)[keeps == TRUE]
depth       <- bruvs$depth[keeps == TRUE]
season      <- as.factor(bruvs$season[keeps == TRUE])
disttoshore <- bruvs$disttoshore[keeps == TRUE]
slope       <- bruvs$slope[keeps == TRUE]
Season      <- bruvs$Season[keeps == TRUE]






#          #############################################################
#          ###           DISTANCE-BASED REDUNDANCY ANALYSIS          ###
#          #############################################################



#                #--------------------------------------------#
#                ####               Run dbRDA              ####
#                #--------------------------------------------#

# Run dbRDA ordination
bruvs_dbrda <- dbrda(divdata ~ habitat + island + slope + disttoshore + depth + Season, distance = "bray")

# Calculate adjusted R2 (variance explained by environmental variables)
RsquareAdj(bruvs_dbrda)

# Test significance of the ordination and the axis
anova(bruvs_dbrda, step = 999)
(dbrdaANOVA <- anova.cca(bruvs_dbrda, step = 999, by = "axis"))

# Write table with dbRDA permutation test output
write.csv(data.frame(as.data.frame(dbrdaANOVA), 
                     importance = c(summary(bruvs_dbrda)$concont$importance[2,], NA)),
          paste(tab_out,"TABLE_S9_dbrdaANOVA.csv",sep="/"))

# Write table with dbRDA scores
(dbrdascrs <- rbind(summary(bruvs_dbrda)$biplot, summary(bruvs_dbrda)$centroids))
write.csv(dbrdascrs,
          paste(tab_out,"TABLE_S10_dbrdascores.csv",sep="/"))





#                #--------------------------------------------#
#                ####              Plot dbRDA              ####
#                #--------------------------------------------#


sc_si <- as.data.frame(scores(bruvs_dbrda, display="sites", choices=c(1,2,3,4), scaling=2))

Island = as.character(island)
for(i in 1:3) Island[which(island == c("PC", "ST","TI")[i])] = c("Príncipe", "São Tomé", "Tinhosas")[i]
Habitat <- factor(as.character(habitat), levels= c("Sand","Maerl","Rock")) 

themeobject <- theme_bw() + theme(
  panel.grid = element_blank(),
  legend.justification = c(1, 0), 
  legend.position = c(1, 0),
  legend.background = element_blank(),
  legend.key.height = unit(0.7,"line"))

dbRDAhabitat = ggplot(data = sc_si,aes(x = dbRDA1, y= dbRDA2, fill = Habitat)) + 
  geom_polygon(data = sc_si[which(habitat == "Sand"),] %>%  slice(chull(dbRDA1, dbRDA2)), fill = "white", colour = colourSand, alpha = 0, linewidth = 1) +
  geom_polygon(data = sc_si[which(habitat == "Rock"),] %>%  slice(chull(dbRDA1, dbRDA2)), fill = "white", colour = colourRock, alpha = 0, linewidth = 1) +
  geom_polygon(data = sc_si[which(habitat == "Maerl"),] %>%  slice(chull(dbRDA1, dbRDA2)), fill = "white", colour = colourMaerl, alpha = 0, linewidth = 1) +
  geom_point(shape = 21, size = 2) + scale_fill_manual(values = c(colourSand, colourMaerl, colourRock)) +
  xlab(paste0("dbRDA1 (",round(summary(bruvs_dbrda)$concont[[1]][2,1]*100,1),"%)")) +
  ylab(paste0("dbRDA2 (",round(summary(bruvs_dbrda)$concont[[1]][2,2]*100,1),"%)")) +
  themeobject

dbRDAisland = ggplot(data = sc_si,aes(x = dbRDA1, y= dbRDA3, fill = Island)) + 
  geom_polygon(data = sc_si[which(island == "ST"),] %>%  slice(chull(dbRDA1, dbRDA3)), colour = colourST, fill = "white",alpha = 0, linewidth = 1) +
  geom_polygon(data = sc_si[which(island == "TI"),] %>%  slice(chull(dbRDA1, dbRDA3)), colour = colourTI, fill = "white",alpha = 0, linewidth = 1) +
  geom_polygon(data = sc_si[which(island == "PC"),] %>%  slice(chull(dbRDA1, dbRDA3)), colour = colourPC, fill = "white", alpha = 0, linewidth = 1) +
  geom_point(shape = 21, size = 2) + 
  scale_fill_manual(values = c(colourPC, colourST, colourTI)) +
  xlab(paste0("dbRDA1 (",round(summary(bruvs_dbrda)$concont[[1]][2,1]*100,1),"%)")) +
  ylab(paste0("dbRDA3 (",round(summary(bruvs_dbrda)$concont[[1]][2,3]*100,1),"%)")) +
  themeobject



tiff(paste0(fig_out,"/Fig4.tiff"), unit = "cm", width = 20, height =10, res = 500)
ggarrange(dbRDAhabitat,dbRDAisland, labels = c("A","B"))
dev.off()


sc_bp <- 1.6 * as.data.frame(scores(bruvs_dbrda, display="bp", choices=c(1,2,3,4), scaling=2))[c("slope","disttoshore","depth","Season"),]
sc_cn <- summary(bruvs_dbrda)$centroid
fontsize = 3
dbRDA12 <- ggplot(data=sc_si,aes(x=dbRDA1,y=dbRDA2)) + geom_point(colour = "grey75", fill = "white", shape = 21) + 
  geom_segment(dat=sc_bp,aes(x=0,y=0,xend=dbRDA1,yend=dbRDA2), 
               arrow = arrow(length=unit(0.2, 'cm')),
               colour = "blue") +
  geom_segment(dat=sc_cn,aes(x=0,y=0,xend=dbRDA1,yend=dbRDA2), 
               colour = "blue", linetype = 2) +
  geom_point(dat=sc_cn, aes(x=dbRDA1,y=dbRDA2),shape = 3, colour = "blue", size = 3, stroke = 1.5) +
  geom_text(dat=sc_bp, aes(x=dbRDA1,y=dbRDA2+0.1, label = rownames(sc_bp)),colour = "blue", vjust = c(0.2,0.2,0.2,2.7), size = fontsize) +
  geom_text(dat=sc_cn, aes(x=dbRDA1,y=dbRDA2+0.1, label = rownames(sc_cn)),colour = "blue", vjust = 0, hjust = c(0,0,0,0,0,0), size = fontsize) +
  xlab(paste0("dbRDA1 (",round(summary(bruvs_dbrda)$concont[[1]][2,1]*100,1),"%)")) +
  ylab(paste0("dbRDA2 (",round(summary(bruvs_dbrda)$concont[[1]][2,2]*100,1),"%)")) +
  theme_bw()

dbRDA13 <- ggplot(data=sc_si,aes(x=dbRDA1,y=dbRDA3)) + geom_point(colour = "grey75", fill = "white", shape = 21) + 
  geom_segment(dat=sc_bp,aes(x=0,y=0,xend=dbRDA1,yend=dbRDA3), 
               arrow = arrow(length=unit(0.2, 'cm')),
               colour = "blue") +
  geom_segment(dat=sc_cn,aes(x=0,y=0,xend=dbRDA1,yend=dbRDA3), 
               colour = "blue", linetype = 2) +
  geom_point(dat=sc_cn, aes(x=dbRDA1,y=dbRDA3),shape = 3, colour = "blue", size = 3, stroke = 1.5) +
  geom_text(dat=sc_bp, aes(x=dbRDA1,y=dbRDA3+0.1, label = rownames(sc_bp)),colour = "blue", vjust = 0, size = fontsize) +
  geom_text(dat=sc_cn, aes(x=dbRDA1,y=dbRDA3+0.1, label = rownames(sc_cn)),colour = "blue", vjust = 0, hjust = c(0,0,1,0,0,0), size = fontsize) +
  xlab(paste0("dbRDA1 (",round(summary(bruvs_dbrda)$concont[[1]][2,1]*100,1),"%)")) +
  ylab(paste0("dbRDA3 (",round(summary(bruvs_dbrda)$concont[[1]][2,3]*100,1),"%)")) +
  theme_bw()


tiff(paste0(fig_out,"/FigS5_dbRDA.tiff"), unit = "cm", width = 10, height =20, res = 500)
ggarrange(dbRDA12,dbRDA13, ncol = 1, labels = c("A","B"))
dev.off()



#          #############################################################
#          ###                        PERMANOVA                      ###
#          #############################################################

perms = 9999
adonis2(divdata ~ depth + disttoshore + slope + habitat + island + season, permutations = perms)
islandposthoc  <- EcolUtils::adonis.pair(vegdist(divdata,method="bray"), island)
habitatposthoc <- EcolUtils::adonis.pair(vegdist(divdata,method="bray"), habitat)




#          #############################################################
#          ###                         SIMPER                        ###
#          #############################################################



###### STEP 0: SUMMARY FUNCTIONS
simperextract <- function(i,a,b){temp <- simperout[[i]]
colnames(temp) <- c("Average","sd","ratio",a,b,"cumsum","p")
return(temp)}

significant <- function(dat,threshold){
  dat = dat[dat$cumsum<threshold,]
  return(dat[which(dat$p < 0.05),])
}




#### STEP 1: SPECIES SIMPER
perms = 9999
simperout <- summary(simper(divdata,habitat,permutations = perms))
Sand_Rock  <- simperextract(1,"Sand","Rock")   
Sand_Maerl <- simperextract(2,"Sand","Maerl")   
Rock_Maerl <- simperextract(3,"Rock","Maerl") 
nrow(Sand_Rock[Sand_Rock$cumsum<0.7,])
nrow(Sand_Maerl[Sand_Maerl$cumsum<0.7,])
nrow(Rock_Maerl[Rock_Maerl$cumsum<0.7,])
nrow(significant(Sand_Rock,0.7))
nrow(significant(Sand_Maerl,0.7))
nrow(significant(Rock_Maerl,0.7))

simperout <- summary(simper(divdata[-which(island=="TI"),],island[-which(island=="TI")],permutations = perms))
PC_ST  <- simperextract(1,"PC","ST")
simperout <- summary(simper(divdata[which(habitat=="Rock"),],island[which(habitat=="Rock")],permutations = perms))
PC_TI  <- simperextract(1,"PC","TI")   
TI_ST <- simperextract(3,"TI","ST")
nrow(PC_ST[PC_ST$cumsum<0.7,])
nrow(PC_TI[PC_TI$cumsum<0.7,])
nrow(TI_ST[TI_ST$cumsum<0.7,])
nrow(significant(PC_TI,0.7))
nrow(significant(TI_ST,0.7))
nrow(significant(PC_ST,0.7))


simperout <- summary(simper(divdata[which(island=="PC"),],season[which(island=="PC")],permutations = perms))
Gravana_Summer_PC       <- simperextract(1,"Gravana","Summer")   
nrow(Gravana_Summer_PC[Gravana_Summer_PC$cumsum<0.7,])
nrow(significant(Gravana_Summer_PC,0.7))

simperout <- summary(simper(divdata[which(island=="ST"),],season[which(island=="ST")],permutations = perms))
Gravana_Summer_ST       <- simperextract(1,"Gravana","Summer")   
nrow(Gravana_Summer_ST[Gravana_Summer_ST$cumsum<0.7,])
nrow(significant(Gravana_Summer_ST,0.7))



###### STEP 2: FUNCTION TO PRINT BARPLOT
# Print SIMPER species' barplot
simperbarplot <- function(cat1,cat2,variable,nspecies,n){
  dat     = get(paste(cat1,cat2,sep="_"))
  average = dat$Average
  p       = dat$p[1:nspecies]
  p       = unlist(lapply(p,function(p){if(p<0.05){"bold.italic"} else {"italic"}}))
  species = rownames(dat)
  label   = data.frame(cat   = c("PC","ST","TI","Sand","Maerl","Rock"),     colour= c(colourPC,colourST,colourTI,colourSand,colourMaerl,colourRock))
  if("TI" %in% c(cat1,cat2)) {label   = cbind(label,label = c("Príncipe (Rock)","São Tomé (Rock)","Tinhosas","Sand","Maerl","Rock"))
  } else {label   = cbind(label,label = c("Príncipe","São Tomé","Tinhosas","Sand","Maerl","Rock"))}
  lab1    = label$label[which(label$cat==cat1)]
  lab2    = label$label[which(label$cat==cat2)]
  col1    = label$colour[label$cat == cat1]
  col2    = label$colour[label$cat == cat2]
  colour  = c(label$colour[which(label$label==sort(c(lab1,lab2))[1])],label$colour[which(label$label==sort(c(lab1,lab2))[2])])
  divdat = raw.spMaxN[keeps==TRUE,]
  if("TI" %in% c(cat1,cat2)) {divdat = divdat[which(habitat=="Rock"),]} else {divdat = divdat}
  dat     = data.frame(species = "",category = "",mean=as.numeric(0),se = as.numeric(0))[-1,]
  yaxis   = c(0,1,5,10,20,50,100,150)
  for(i in 1:nspecies) { x=species[i]
  dat1 = divdat[which(variable==cat1),which(colnames(divdat)==x)]
  dat2 = divdat[which(variable==cat2),which(colnames(divdat)==x)]
  dat  = rbind(dat,data.frame(species = x, category = lab1, mean = mean(dat1)^0.25,se = 0.03+(mean(dat1)+ 2*(sd(dat1)/length(dat1^0.5)))^0.25,average = average[i]*30))
  dat  = rbind(dat,data.frame(species = x, category = lab2, mean = mean(dat2)^0.25,se = 0.03+(mean(dat2)+ 2*(sd(dat2)/length(dat2^0.5)))^0.25,average = average[i]*30))}
  ggplot(data=dat,aes(x=species,y=mean,ymin=0,ymax=3.6,fill=category)) + 
    geom_bar(position="dodge", stat='identity',colour='black',width=0.7, alpha = 0.8)  +
    scale_fill_manual(values=colour)+ 
    geom_errorbar(aes(ymin=mean,ymax=se),position="dodge",width=0.7) +
    geom_line(aes(x=species,y=average,group=1),colour="#65BFFF",linewidth = 0.2*n) +
    # Legend
    geom_rect(aes(xmin = nspecies-0.5, xmax = nspecies+0.1, ymin = 3.25, ymax = 3.5),fill=col1,colour="black", alpha = 0.8) +
    annotate(geom="text", label = lab1, x = nspecies-0.7, y = 3.5,hjust = 1, vjust = 1,size = 0.65*n)      +
    geom_rect(aes(xmin = nspecies-0.5, xmax = nspecies+0.1, ymin = 2.85, ymax = 3.1),fill=col2,colour="black", alpha = 0.8) +
    annotate(geom="text", label = lab2, x = nspecies-0.7, y =3.08,hjust = 1, vjust = 1,size = 0.65*n)      +
    geom_segment(aes(x=nspecies-0.5, xend=nspecies+0.1, y=2.5, yend=2.5),colour ="#65BFFF", linewidth=0.2*n) +
    annotate(geom="text", label = "Contribution to dissimilarities", x = nspecies-0.7, y = 2.6,hjust = 1, vjust = 1,size = 0.65*n)      +
    scale_x_discrete(limits = species[1:nspecies], breaks = species[1:nspecies], labels=species[1:nspecies]) +
    scale_y_continuous(breaks = yaxis^0.25, labels=yaxis,sec.axis = sec_axis(~.*(1/30),name="Contribution to overall dissimilarity", breaks = seq(0,12/30,0.03))) +
    ylab("MaxN") + labs(fill = "MaxN") + theme_minimal() + 
    theme(  
      axis.line     = element_line(linewidth=0.1*n),
      axis.text.x   = element_text(size = 2.1*n, angle = 35,  hjust=1, colour="black",face=p),
      axis.text.y   = element_text(size = 2.1*n,   colour="black"),
      axis.title.x  = element_blank(),    
      axis.title.y  = element_text(size=2*n,face="bold",colour="black"),
      panel.grid = element_blank() ,
      plot.margin        = margin(t = 1*n, r = 1*n, b = 3*n, l = 6.5*n,"points"),
    ) + 
    guides(fill = "none")
}  

###### STEP 3: EXPORT BARPLOT

tiff(filename=paste(fig_out,"/fig5.tiff",sep=""),width = 18,height=20, res = 500, unit = "cm")
n=3
precision = 1

short_margin <- theme(plot.margin        = margin(t = 1*n, r = 1*n, b = 3*n, l = 3*n,"points"))
large_margin <- theme(plot.margin        = margin(t = 1*n, r = 1*n, b = 7*n, l = 3*n,"points"))

cowplot::plot_grid(
  simperbarplot("Rock","Maerl",habitat,12,n),
  simperbarplot("PC","ST",island,12,n) + short_margin,
  simperbarplot("Sand","Rock",habitat,12,n),
  simperbarplot("PC","TI",island[which(habitat=="Rock")],12,n),
  simperbarplot("Sand","Maerl",habitat,12,n),
  simperbarplot("TI","ST",island[which(habitat=="Rock")],12,n), 
  ncol = 2, align = "v", axis = "tblr",
  labels = c("A","D","B","E","C","F"), label_size = 4*n
)
dev.off()


###### STEP 3: TABLES
meansdfunction <- function(variable){
  for(j in 1:length(levels(as.factor(variable)))){
    category = levels(as.factor(variable))[j]
    temp = unlist(lapply(colnames(raw.spMaxN),function(species){
      dat = raw.spMaxN[which(variable==category),colnames(raw.spMaxN)==species]
      paste(round(mean(dat),2), " (",round(sd(dat),2),  ")",sep="")}))
    if(j == 1) {dat = data.frame(temp = temp)
    colnames(dat) = category
    } else {
      names = c(colnames(dat),category)
      dat   = cbind(dat,temp)
      colnames(dat) = names}}
  return(dat)}
simpers <-  function(dat){unlist(lapply(colnames(divdata),function(species){
  temp = dat[which(rownames(dat)==species),]
  p = temp$p
  if(p<0.01) {p = "p<0.01"} else {p = paste("p=",round(p,2),sep="")}
  paste(round(temp$Average,3),p,sep=", ")}))}


habitatsimper <- cbind(meansdfunction(habitat),Sand_Maerl = simpers(Sand_Maerl),Rock_Maerl = simpers(Rock_Maerl),Sand_Rock = simpers(Sand_Rock))
rownames(habitatsimper) <- colnames(divdata)
write.csv(habitatsimper,paste(tab_out,"TABLE_S11_simperhabitat.csv",sep="/"))

islandsimper <- cbind(meansdfunction(island),PC_ST = simpers(PC_ST),TI_ST = simpers(TI_ST),PC_TI = simpers(PC_TI))
rownames(islandsimper) <- colnames(divdata)
write.csv(islandsimper,paste(tab_out,"TABLE_S12_simperisland.csv",sep="/"))




