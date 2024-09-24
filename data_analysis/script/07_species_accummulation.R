
################################################################################
################################################################################
###                           SPECIES ACCUMULATION                          ###
################################################################################
################################################################################








#          #############################################################
#          ###                 LOAD DATA AND PACKAGES                ###
#          #############################################################






### LOAD DATA
source(file.path(getwd(),"data_analysis/script/02_variable_definition.R"))




                  #-------------------------------------------#
                  ############ Species accumulation ###########
                  #-------------------------------------------#

#### STEP 1: CREATE 1 COLUMN PER SPECIES AND DIVIDE DATABASES BY ISLAND
speciesmatrix <- divdata

ST <-speciesmatrix[which(bruvs$island =="ST"),]
PC <-speciesmatrix[which(bruvs$island =="PC"),]


#### STEP 2: PLOT SPECIES ACCUMULATION CURVES
plotspecacum <- function(dat,label){
       spac <- specaccum(dat,method="random",permutations=1000)
       spac <- data.frame(sites =spac$sites,richness=spac$richness,sd=spac$sd)
       plot(c(0,270),c(0,150),col="white",
            ylab="Species richness",xlab="number of sites",main=label)
           polygon(x=c(spac$sites,rev(spac$sites)),
                   y=c(spac$richness+spac$sd,rev(c(spac$richness-spac$sd))),
                   col="lightgrey")
           lines(spac$sites,spac$richness,col="black",lwd=2)
           lines(spac$sites,spac$richness+spac$sd,col="black",lwd=1.5)
           lines(spac$sites,spac$richness-spac$sd,col="black",lwd=1.5)}








                #----------------------------------------------#
                ######## SPECIES ACCUMMULATION CURVES 1 ########
                #----------------------------------------------#



## STEP 1: CREATE "videospecies" AND "videospeciestime" COLUMNS
earliestoccurrence <- na.omit(maxn)
earliestoccurrence$videospecies <- paste(earliestoccurrence$video,earliestoccurrence$sp,sep="") 
earliestoccurrence$videospeciestime <- paste(earliestoccurrence$video,earliestoccurrence$sp,earliestoccurrence$minutes,sep="") 
earliestoccurrence <- na.omit(earliestoccurrence)


## STEP 2: KEEP ONLY FIRST OBSERVATION OF EACH SPECIES
spvideo <- levels(as.factor(earliestoccurrence$videospecies))
temp  <-  earliestoccurrence[0,]
for(i in 1:length(spvideo)){
    x = earliestoccurrence[earliestoccurrence$videospecies==spvideo[i],]
    y = paste(spvideo[i],min(x$minutes),sep="")
    z = earliestoccurrence[earliestoccurrence$videospeciestime == y,]
    temp=rbind(temp,z)}
earliestoccurrence <- temp


## STEP 3: REMOVE INVALID VIDEOS
bruvsspac <- bruvs[!bruvs$viewer == "Jemima Dimbleby",]


## STEP 4: BUILD SPECIES ACCUMMULATION MATRIX
temp <- data.frame()
for(i in 1:length(bruvsspac$video)){
      x = bruvsspac$video[i]
      time = as.numeric(c(1,5,10,15,20,25,30,35,40,45,50,55,60))
      for(j in 1:length(time)){
           k = earliestoccurrence[earliestoccurrence$video==x,]
           k = na.omit(k)
           k = k[k$minutes<=time[j],]
           y = c(x,time[j],length(levels(as.factor(k$sp))))            
           temp = rbind(temp,y)}}
colnames(temp) <- c("video","time","richness")
spac<- temp
spac$habitat <- bruvs$habitat[match(spac$video,bruvs$video)]
spac$island <- bruvs$island[match(spac$video,bruvs$video)]
spac$time <- as.numeric(spac$time)
spac$richness <- as.numeric(spac$richness)


## STEP 5: FUNCTION TO SUMMARISE AVERAGE RICHNESS AND SD (TIME SLOT & HABITAT)
summarise.mean.sd <- function(habitat){
      time = as.numeric(c(1,5,10,15,20,25,30,35,40,45,50,55,60))
      temp = data.frame()
      for(i in 1:length(time)) {
             x       = time[[i]]
             dat     = habitat[habitat$time==x,]
             dat =na.omit(dat)
             average = mean(as.numeric(dat$richness))
             sd      = sqrt(var(dat$richness) / length(dat$richness))
             result  = c(x,average,sd)
             temp    = rbind(result,temp)}
       colnames(temp) = c("time","richness_mean","richness_sd")
       return(temp)}

# Create function to plot points and error bars
ploterrorbars <- function(x,y){
                 points(x$time,x$richness_mean,pch=19,col=y,cex=n*1.1)
                 arrows(x$time,x$richness_mean-x$richness_sd,
                        x$time,x$richness_mean+x$richness_sd,
                        length=0.025, angle=90, code=3,col=y)}

richnessvstime <- function(dat,main,col){
     rock  = dat[which(dat$habitat=="Rock"),]
     maerl = dat[which(dat$habitat=="Maerl"),]
     sand = dat[which(dat$habitat=="Sand"),]
     plot(c(0,60),c(0,22),ylim=c(0, 22),pch=19, 
          xlab="Time", ylab="Mean species richness",main=main,
          col="white")
     for(i in c(1,2,3)){
        hab = list(rock,maerl,sand)[[i]]
        hab = summarise.mean.sd(hab)
        ploterrorbars(hab,col[i])
          lines(data.frame(rev(hab$time),rev(hab$richness_mean)),col = col[i])}
     legend("topleft",c("Rock","Maerl","Sand"),
            fill=col,xjust=1,bty = "n",cex=0.8*n)}

### STEP 6: PRINT FIGURE S6
tiff(filename=paste(fig_out,"/figS6_specacum.tiff",sep=""),width = 18,height=18, res = 500, unit = "cm")
      n=1.2
      par(mar=c(5,4,2,2),mgp=c(2,0.5,0),family="serif",font.axis=1,font.lab=2,
          cex.lab=1.1*n,cex.axis=1.1*n,cex.main=1.2*n,col="black",mfrow=c(2,2),lwd=1.8)
richnessvstime(spac[spac$island=="PC",],main="Príncipe",col=c("gray0","gray47","gray60"))
             mtext('A', side=2,line=2,cex=n*1.4,font=2,at=24,las=2)
plotspecacum(speciesmatrix[which(bruvs$island =="PC"),],"Príncipe")
             mtext('B', side=2,line=2,cex=n*1.4,font=2,at=162,las=2)
richnessvstime(spac[spac$island=="ST",],main="São Tomé",col=c("gray0","gray47","gray60"))
plotspecacum(speciesmatrix[which(bruvs$island =="ST"),],"São Tomé")

dev.off()






