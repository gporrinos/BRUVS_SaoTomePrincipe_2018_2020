
################################################################################
################################################################################
###                             OBSERVED SPECIES                             ###
################################################################################
################################################################################





### STEP 1: LOAD DATABASES
source(file.path(getwd(),"data_analysis/script/02_variable_definition.R"))


### step 2: LOAD PACKAGES
if(!require("ggplot2")) install.packages("ggplot2",   repos = "https://cloud.r-project.org")
if(!require("dplyr")) install.packages("dplyr",   repos = "https://cloud.r-project.org")
if(!require("ggpubr")) install.packages("ggpubr",   repos = "https://cloud.r-project.org")









#                #--------------------------------------------#
#                ####           Stacked barplot            ####
#                #--------------------------------------------#





### STEP 1: FUNCTION THAT BUILDS A TABLE OF MAXN AND OCCURRENCES BY FAMILY (ORGANISED BY NUMBER OF OCCURRENCES)
createdat <- function(island,variable,numberoffamilies){
  ## Count occurrences of each family, and organise families by number of occurrences
  dat = mxMaxN[mxMaxN$island==island,]            # Filter by island
  dat = dat[dat$type=="Teleost",]                 # Keep teleosts only
  nBRUVS = length(levels(as.factor(dat$video)))   # Count number of deployments
  familylist = data.frame(familylist = levels(as.factor(dat$family)),
                          occurr.all = unlist(lapply(levels(as.factor(dat$family)),function(x){
                            length(levels(as.factor(dat[which(dat$family==x),]$video)))
                          })))
  familylist = arrange(familylist,-occurr.all)    # Organise families by number of occurrences
  ## Create family mean MaxN table, organised by higher number of family occurrences
  if(variable == "maxn"){
    otherfamilies = familylist$familylist[(numberoffamilies+1):length(familylist$familylist)]
    dat$family[which(dat$family %in% otherfamilies)] = "Other"
    familylist$familylist[(numberoffamilies+1):length(familylist$familylist)] = "Other"
    
    meanmaxnperfamily = lapply(c(familylist$familylist[1:numberoffamilies],"Other"),function(x){
      maxnpervideo = unlist(lapply(levels(as.factor(dat$video)),function(y){
        temp = dat[which(dat$family==x),]
        return(sum(temp[which(temp$video==y),]$maxn))}))
      return(
        data.frame(mean = mean(maxnpervideo),sd = sd(maxnpervideo)/(length(maxnpervideo)^0.5))
      )})
    for(i in 1:length(meanmaxnperfamily)){
      if(i == 1) {outputtable = meanmaxnperfamily[[i]]} else {outputtable = rbind(outputtable,meanmaxnperfamily[[i]])}
    }
    outputtable$family <- c(familylist$familylist[1:numberoffamilies],"Other")
  }      # Close "if variable == maxn" statement
  ## Create family ocurrences table, disaggregated  by habitat
  if(variable == "occurrences"){
    for(j in 1:length(familylist$familylist)){
      temp = dat[which(dat$family==familylist$familylist[j]),]
      temp = data.frame(Rock  = length(levels(as.factor(temp[which(temp$habitat=="Rock"),]$video))),
                        Maerl = length(levels(as.factor(temp[which(temp$habitat=="Maerl"),]$video))),
                        Sand  = length(levels(as.factor(temp[which(temp$habitat=="Sand"),]$video))))
      if(j==1) {occurrences = temp} else {occurrences = rbind(occurrences,temp)}          
    }
    temp = occurrences[(numberoffamilies+1):length(familylist$familylist),]
    occurrences = rbind(occurrences[1:numberoffamilies,],
                        data.frame(Rock = sum(temp[,1]),Maerl = sum(temp[,2]), Sand = sum(temp[,3])))
    occurrences$family <- c(familylist$familylist[1:numberoffamilies],"Other")
    outputtable = expand.grid(habitat = c("Rock","Maerl","Sand"),family = occurrences$family)
    outputtable$occurrences = 100*unlist(lapply(1:nrow(outputtable),function(i){
      occurrences[which(occurrences$family==outputtable$family[i]),
                  which(colnames(occurrences)==outputtable$habitat[i])]
    }))/nBRUVS
  }
  ## Return outputtable
  return(outputtable)}




### STEP 2: CREATE TABLES FOR THE PLOTS
numberoffamilies = 25         # Number of families that will be plotted (the rest will be aggregated into "Other")
PCoccurrences <- createdat("PC","occurrences",numberoffamilies)
SToccurrences <- createdat("ST","occurrences",numberoffamilies)
PCmaxn        <- createdat("PC","maxn",numberoffamilies)
STmaxn        <- createdat("ST","maxn",numberoffamilies)




### STEP 3: SET SCALE FACTOR
n=3.2      # Set scale factor





### STEP 4: CREATE THEME OBJECTS TO USE ON THE PLOTS
ggtheme <-  theme_minimal() + theme(
  axis.text.x   = element_text(size = 2.3*n,colour="black"),
  axis.text.y   = element_text(size = 2.3*n, colour="black"),
  axis.title.x  = element_text(size=2.4*n,face="bold", colour="black"),
  axis.title.y  = element_text(size=2.5*n,face="bold", colour="black"),
  panel.grid.major.y = element_blank() ,
  panel.grid.major.x = element_line(color = "darkgrey",linewidth = 0.01*n,linetype = 1) ,
  panel.grid.minor.y = element_blank() ,
  panel.grid.minor.x = element_line(color = "darkgrey",linewidth = 0.01*n,linetype = 1) )

legendtheme <- theme(
  legend.text        = element_text(size = 2.15*n),
  legend.title       = element_text(size=2.2*n,face = "bold"),
  legend.position    = c(0.79,0.13),
  legend.background  = element_rect(fill="white",linewidth=0.5, linetype="solid", colour ="white"),
  legend.key.width = unit(1.1,"line"),
  legend.key.height = unit(0.7,"line"))

short_margin <- theme(plot.margin        = margin(t = 1*n, r = 1*n, b = 3*n, l = 3*n,"points"))
large_margin <- theme(plot.margin        = margin(t = 1*n, r = 1*n, b = 7*n, l = 3*n,"points"))




### STEP 5: CREATE PLOT OBJECTS
PCoccurrences.plot <- ggplot(PCoccurrences, aes(x = factor(family,PCmaxn$family[(nrow(PCmaxn)+1):1]),  y = occurrences, fill = habitat)) + 
  geom_bar(position='stack', stat='identity',fill='white', colour='white') + 
  geom_bar(position='stack', stat='identity',colour='black', alpha = 0.75) + ggtheme  + legendtheme + coord_flip() + ylim(0,100) +
  scale_fill_manual(values=c(colourRock,colourMaerl,colourSand)) + 
  xlab("") + ylab("Occurrences (% deployments, Príncipe)") + labs(fill = "Habitat")

SToccurrences.plot <- ggplot(SToccurrences, aes(x = factor(family,STmaxn$family[(nrow(STmaxn)+1):1]),  y = occurrences, fill = habitat)) + 
  geom_bar(position='stack', stat='identity',fill='white', colour='white') + 
  geom_bar(position='stack', stat='identity',colour='black', alpha = 0.75) + ggtheme  + coord_flip() + legendtheme + ylim(0,100) +
  scale_fill_manual(values=c(colourRock,colourMaerl,colourSand)) + 
  xlab("") + ylab("Occurrences (% deployments, São Tomé)") + labs(fill = "Habitat")

PCmaxn.plot <- ggplot(PCmaxn, aes(x = factor(family,family[(nrow(PCmaxn)+1):1]),  y = mean)) + 
  geom_bar(stat='identity',colour='black',fill='grey50') + ggtheme  + coord_flip() + ylim(0,40) +
  geom_errorbar(aes(ymin=mean,ymax=mean+sd)) +
  xlab("") + ylab("Mean MaxN (Príncipe)")

STmaxn.plot <- ggplot(STmaxn, aes(x = factor(family,family[(nrow(STmaxn)+1):1]),  y = mean)) + 
  geom_bar(stat='identity',colour='black',fill='grey50') + ggtheme  + coord_flip() + ylim(0,40) +
  geom_errorbar(aes(ymin=mean,ymax=mean+sd)) +
  xlab("") + ylab("Mean MaxN (São Tomé)")




tiff(filename=paste(fig_out,"/Fig2.tiff",sep=""),unit = "cm", width = 18,height=19, res = 500)
ggarrange(PCoccurrences.plot + large_margin, PCmaxn.plot + large_margin, SToccurrences.plot + short_margin, STmaxn.plot  + short_margin,     
          labels = c("A", "B","C","D"),       
          font.label = list(size = 3*n,face="bold"),        ncol = 2, nrow = 2)
dev.off()





rm(n,numberoffamilies,PCmaxn,PCmaxn.plot,PCoccurrences,PCoccurrences.plot,
   STmaxn,STmaxn.plot,SToccurrences,SToccurrences.plot,short_margin,large_margin,
   legendtheme, ggtheme, createdat)



#                #--------------------------------------------#
#                ######### Species and family list  ###########
#                #--------------------------------------------#


### STEP 0: CREATE FUNCTION TO SUMMARISE TAXA
taxasummary <- function(dat,taxalist,taxa){
  dat$taxa   = dat[,which(colnames(dat)==taxa)]
  temp     = data.frame()
  for(i in 1:length(taxalist)){
    taxa     = taxalist[i]
    taxadata = dat[dat$taxa==taxalist[i],]
    taxadata = taxadata[!is.na(taxadata$maxn),]
    taxadata = taxadata[!is.na(taxadata$video),]
    Rock     = taxadata[taxadata$habitat == "Rock",]
    Maerl    = taxadata[taxadata$habitat == "Maerl",]
    Sand     = taxadata[taxadata$habitat == "Sand",]
    x = data.frame(
      taxa        = taxa,
      sp.richness = length(levels(as.factor(taxadata$specieslist))),
      occurr      = length(levels(as.factor(taxadata$video))),
      oc.perc     = (1/length(levels(as.factor(dat$video))))
      * length(levels(as.factor(taxadata$video))),
      sp.occurr   = nrow(taxadata),
      MaxN        = sum(taxadata$maxn),
      avMaxN      = (1/length(levels(as.factor(dat$video))))
      * sum(taxadata$maxn),
      avMaxnperOc = sum(taxadata$maxn)/length(levels(as.factor(taxadata$video)))
    )
    temp = rbind(temp,x)}
  return(temp)}




### STEP 1: SPECIES LIST
teleostsp             <- levels(as.factor(species[species$type=="Teleost",]$specieslist))
elasmobsp             <- levels(as.factor(species[species$type=="Elasmobranch",]$specieslist))
specieslist           <- c(teleostsp,elasmobsp)
specieslist           <- data.frame(taxasummary(mxMaxN[mxMaxN$island=="PC",],specieslist,"sp")[,c(1,3,4,6,7,8)],
                                    taxasummary(mxMaxN[mxMaxN$island=="ST",],specieslist,"sp")[,c(3,4,6,7,8)],
                                    taxasummary(mxMaxN[mxMaxN$island=="TI",],specieslist,"sp")[,c(3,4,6,7,8)])
colnames(specieslist) <- c("species","occurr.PC","oc.perc.PC","MaxN.PC","MeanMaxNPC","MaxN.per.Oc.PC","occurr.ST","oc.perc.ST","MaxN.ST","MeanMaxNST","MaxN.per.Oc.ST","occurr.TI","oc.perc.TI","MaxN.TI","MeanMaxNTI","MaxN.per.Oc.TI")
family                <- unlist(lapply(specieslist$species,function(x){species[match(x,species$specieslist),]$family}))
type                  <- unlist(lapply(specieslist$species,function(x){species[match(x,species$specieslist),]$type}))
specieslist           <- data.frame(specieslist,family,type)[,c(17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
specieslist           <- arrange(specieslist,family)
specieslist           <- arrange(specieslist,type)
write.csv(specieslist,paste(tab_out,"TABLE_S3_species_list.csv",sep="/"))



### STEP 2: COUNT SPECIES AND FAMILIES OF TELEOSTS, SHARKS AND RAYS
elasm <- function(variable,island,functionalgroup){
  isl   = mxMaxN[which(mxMaxN$island==island),]
  elasm = isl[which(isl$functionalgroup == functionalgroup),]
  if(variable == "cpue") {output = sum(elasm$maxn)/nrow(bruvs[which(bruvs$island==island),])}
  if(variable == "occurrences") {output = length(elasm$maxn)/nrow(bruvs[which(bruvs$island==island),])}
  return(output)}

countfamilies <- function(island){
  elasm = specieslist[specieslist$type == "Elasmobranch",]
  elasm = elasm[!elasm[,which(colnames(elasm) == paste("occurr",island,sep="."))] ==0, ]
  teleo = specieslist[specieslist$type == "Teleost",]
  teleo = teleo[!teleo[,which(colnames(teleo) == paste("occurr",island,sep="."))] ==0, ]
  return(
    data.frame(n.teleost.sp = nrow(teleo),n.teleost.fam = length(levels(as.factor(teleo$family))),
               n.elasmobranch.sp = nrow(elasm),n.elasmobranch.fam = length(levels(as.factor(elasm$family))),
               occurr.sharks = elasm("occurrences",island,"large_sharks"),cpue.sharks = elasm("cpue",island,"large_sharks"),
               occurr.rays = elasm("occurrences",island,"medium-large_rays"),cpue.rays = elasm("cpue",island,"medium-large_rays"))
  )}

nsp.families <- rbind(countfamilies("PC"),countfamilies("ST"),countfamilies("TI"))
rownames(nsp.families) <- c("PC","ST","TI")
nsp.families    # Table with number of species and families of teleosts, sharks and rays



### STEP 3: FAMILY TABLE
teleostfam            <- levels(as.factor(species[species$type=="Teleost",]$family))
elasmobfam            <- levels(as.factor(species[species$type=="Elasmobranch",]$family))
familylist            <- c(teleostfam,elasmobfam)
familylist            <- data.frame(taxasummary(mxMaxN[mxMaxN$island=="PC",],familylist,"family")[,c(1,2,3,4,6,7)],
                                    taxasummary(mxMaxN[mxMaxN$island=="ST",],familylist,"family")[,c(2,3,4,6,7)],
                                    taxasummary(mxMaxN[mxMaxN$island=="TI",],familylist,"family")[,c(2,3,4,6,7)])
colnames(familylist)  <- c("family","sp.richness.PC","occurr.PC","oc.perc.PC","MaxN.PC","MaxN.av.PC","sp.richness.ST","occurr.ST","oc.perc.ST","MaxN.ST","MaxN.av.ST",
                           "sp.richness.TI","occurr.TI","oc.perc.TI","MaxN.TI","MaxN.av.TI")
type                  <- unlist(lapply(familylist$family,function(x){species[match(x,species$family),]$type}))
familylist            <- data.frame(familylist,type)[,c(17,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
familylist            <- arrange(familylist,type)
write.csv(familylist,paste(tab_out,"TABLE_S4_family_list.csv",sep="/"))




rm(elasmobfam,elasmobsp,teleostfam,teleostsp,family,type,countfamilies,elasm,taxasummary)
