
################################################################################
################################################################################
###                           DIVERSITY INDICATORS                           ###
################################################################################
################################################################################








#          #############################################################
#          ###                 LOAD DATA AND PACKAGES                ###
#          #############################################################






### LOAD DATA (LAST SCRIPT)
source(file.path(getwd(),"data_analysis/script/02_variable_definition.R"))

### Statistical packages ###
if(!require("MASS")) install.packages("MASS",   repos = "https://cloud.r-project.org")
if(!require("betareg")) install.packages("betareg",   repos = "https://cloud.r-project.org")
if(!require("car")) install.packages("betareg",   repos = "https://cloud.r-project.org")
if(!require("mgcv")) install.packages("mgcv",   repos = "https://cloud.r-project.org")
if(!require("MuMIn")) install.packages("MuMIn",   repos = "https://cloud.r-project.org")
if(!require("gam.hp")) install.packages("gam.hp",   repos = "https://cloud.r-project.org")


### Data management
if(!require("dplyr")) install.packages("dplyr",   repos = "https://cloud.r-project.org")

### Plot things ###
if(!require("ggplot2")) install.packages("ggplot2",   repos = "https://cloud.r-project.org")
if(!require("gridExtra")) install.packages("gridExtra",   repos = "https://cloud.r-project.org")
if(!require("gratia"))  install.packages("gratia", repos = "https://cloud.r-project.org")








#          #############################################################
#          ###            OBSERVED RICHNESS AND ABUNDANCE            ###
#          #############################################################




# Create function to summarise indicators
by.habitat <- function(data,variable){
  x           = c(mean(na.omit(data[,which(colnames(data)==variable)])),
                  sd(na.omit(data[,which(colnames(data)==variable)])),
                  length(na.omit(data[,which(colnames(data)==variable)])))
  for (i in 1:length(levels(as.factor(data$habitat)))){
    temp=data[which(data$habitat==levels(as.factor(data$habitat))[i]),]
    x = c(x,
          mean(na.omit(temp[,which(colnames(temp)==variable)])),
          sd(na.omit(temp[,which(colnames(temp)==variable)])),
          length(na.omit(temp[,which(colnames(temp)==variable)])))}
  empty       = data.frame()
  x           = rbind(x,empty)
  colnames(x) = c("all.mean","all.sd","all.n","maerl.mean","maerl.sd","maerl.n",
                  "rock.mean","rock.sd","rock.n",
                  "sand.mean","sand.sd","sand.n") 
  return(x)}

summary.table <- function(variable){
  temp = rbind(by.habitat(bruvs,variable),
               by.habitat(bruvs[which(bruvs$island=="ST"),],variable),
               by.habitat(bruvs[which(bruvs$island=="PC"),],variable))
  temp$variable = variable
  temp = temp[,c(13,1,2,3,4,5,6,7,8,9,10,11,12)]
  rownames(temp) <- c("ST+PC","ST","PC")
  return(temp)}

# Summarise indicators
(summary.table <- rbind(summary.table("richness"),summary.table("abundance"),summary.table("evenness")))

write.csv(summary.table,paste(tab_out,"TABLE_S5_diversity_indicators_by_hab_and_isl.csv",sep="/"))




#          #############################################################
#          ###            OBSERVED RICHNESS AND ABUNDANCE            ###
#          #############################################################



#                             #----------------------#
#                             ###### Run models ######
#                             #----------------------#



### STEP 1: RICHNESS
# Prepare data
dat <- bruvs[,c("richness","depth","disttoshore","slope","habitat","island","Season")]
dat <- na.omit(dat[-which(dat$island == "TI"),])
dat$island <- as.factor(dat$island)

# Check Variable Inflation Factor (excluding seasonality)
mod <-     MASS::glm.nb(richness ~ depth + disttoshore + slope + habitat + island, 
                        data = dat)
car::vif(mod)

# Run model
mod <-     gam(richness ~ s(depth) + s(disttoshore) + s(slope) + habitat + island + 
                 s(Season, bs = 'cc', k=-1) +  s(Season, by = island, bs = 'cc', k=-1) , 
               data = dat,   family = nb(link = "log"), method = "REML")


options(na.action = "na.fail")
bestmod <- dredge(global.model = mod, rank="AIC",extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))
bestmod <- subset(bestmod,delta<6)
bestmod <- subset(bestmod, !nested(.))
richness.bestmod <- as.data.frame(bestmod)
richness = model.avg(bestmod, fit = TRUE)
sw(richness)




### STEP 2: ABUNDANCE
# Prepare data
dat <- bruvs[c("abundance","depth","disttoshore","slope","habitat","island","Season")]
dat <- na.omit(dat[-which(dat$island == "TI"),])
dat$island <- as.factor(dat$island)

# Check Variable Inflation Factor (excluding seasonality)
mod <-     MASS::glm.nb(abundance ~ depth + disttoshore + slope + habitat + island, 
                        data = dat)
car::vif(mod)

# Run model
mod <-     gam(abundance ~ s(depth) + s(disttoshore) + s(slope) + habitat + island + 
                 s(Season, bs = 'cc', k=-1) +  s(Season, by = island, bs = 'cc', k=-1), 
               data = dat,   family = nb(link = "log"), method = "REML")
options(na.action = "na.fail")
bestmod <- dredge(global.model = mod, rank="AIC",extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))
bestmod <- subset(bestmod,delta<6)
bestmod <- subset(bestmod, !nested(.))
abundance.bestmod <- as.data.frame(bestmod)
abundance <- model.avg(bestmod, fit = TRUE)
sw(abundance)




### STEP 3: EVENNESS
# Prepare data
dat <- bruvs[c("evenness","abundance", "depth","disttoshore","slope","habitat","island","Season")]
dat <- na.omit(dat[-which(dat$island == "TI"),])
dat$island <- as.factor(dat$island)
dat$evenness <- (((dat$evenness/dat$abundance)*(dat$abundance-1))+0.5)/dat$abundance

# Check Variable Inflation Factor (excluding seasonality)
mod <-     betareg::betareg(evenness ~ depth + disttoshore + slope + habitat + island, 
                            data = dat)
car::vif(mod)



# Run model
mod <-     gam(evenness ~ s(depth) + s(disttoshore) + s(slope) + habitat + island + 
                 s(Season, bs = 'cc', k=-1) +  s(Season, by = island, bs = 'cc', k=-1), 
               data = dat,   family = betar(link = "log"), method = "REML")
options(na.action = "na.fail")
bestmod <- dredge(global.model = mod, rank="AIC", extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))
bestmod <- subset(bestmod,delta<6)
bestmod <- subset(bestmod, !nested(.))
evenness.bestmod <- as.data.frame(bestmod)
evenness <- model.avg(bestmod, fit = TRUE)
sw(evenness)


### STEP 4: MODEL SELECTION TABLE
write.csv(rbind(data.frame(var = "Richness", richness.bestmod),
                data.frame(var = "Abundance",abundance.bestmod),
                data.frame(var = "Evenness", evenness.bestmod)),
          paste(tab_out,"TABLE_S6_modelselection.csv",sep="/"))






### STEP 5: EXPORT COEFFICIENTS AND RELATIVE IMPORTANCE
calculatecoefs <- function(indicator, names = "columns"){
  effect = c("(Intercept)", "depth","slope","disttoshore","habitatRock","habitatSand","islandST","s(Season)","s(Season, by = island)")
  vars   = c("(Intercept)", "s(depth)","s(slope)","s(disttoshore)","habitat","habitat","island","s(Season, bs = \"cc\", k = -1)","s(Season, by = island, bs = \"cc\", k = -1)")
  mod = get(indicator)
  relativeimportance = unlist(lapply(vars,function(x)    if(x %in% names(sw(mod)))   {sw(mod)[which(names(sw(mod)) == x)]} else {if(grepl("Intercept",x)) {NA} else {0}}))
  names(relativeimportance) = effect
  coefs = as.data.frame((summary(mod))$coefmat.full[,c(1,2)])
  estimate = unlist(lapply(effect,   function(x)    if(x %in% rownames(coefs))   {
    coefs[which(rownames(coefs) == x),1]
  } else {
    if(grepl("Season",x) | x %in% c("depth","slope","disttoshore")) {NA} else {0}}))
  stderr   = unlist(lapply(effect,   function(x)    if(x %in% rownames(coefs))   {
    coefs[which(rownames(coefs) == x),2]
  } else {
    if(grepl("Season",x) | x %in% c("depth","slope","disttoshore")) {NA} else {0}}))
  names(estimate)  = effect
  names(stderr)    = effect
  output           = data.frame(vars,effect,estimate,stderr,relativeimportance)
  rownames(output) = 1:nrow(output)
  if(names == "columns") colnames(output) = paste(indicator,colnames(output))
  if(names == "rows") output$indicator = indicator
  return(output)}

coefs <- rbind(calculatecoefs("richness", "rows"),calculatecoefs("abundance", "rows"),calculatecoefs("evenness", "rows"))
coefs$vars <- gsub(', bs = "cc", k = -1', '', coefs$vars)
write.csv(cbind(calculatecoefs("richness"),calculatecoefs("abundance"),calculatecoefs("evenness")),
          paste(tab_out,"TABLE_S7_modeloutputs.csv",sep="/"))






### STEP 6: CONTRIBUTION OF INDIVIDUAL PREDITORS TO R2
predictors <- labels(terms(as.formula("richness ~ s(depth) + s(disttoshore) + s(slope) + habitat + island + 
  s(Season, bs = 'cc', k=-1) +  s(Season, by = island, bs = 'cc', k=-1)")))

calculate.contributions <- function(divindicator){
  mods = get.models(get(divindicator), subset = delta < 6)
  for(i in 1:length(mods)) {
    hier.part = gam.hp(mods[[i]])$hierarchical.partitioning
    positions = which(substr(rownames(hier.part),1,1)==" ")
    rownames(hier.part)[positions] = substr(rownames(hier.part)[positions],
                                            2, nchar(rownames(hier.part)[positions]))
    nchar.rwnames = nchar(rownames(hier.part))
    positions = which(substr(rownames(hier.part),nchar.rwnames,nchar.rwnames)==" ")
    rownames(hier.part)[positions] = substr(rownames(hier.part)[positions],
                                            1, nchar.rwnames-1)
    contribution = unlist(lapply(predictors,function(predictor){
      temp = which(rownames(hier.part) == predictor)
      if(length(temp) == 0) NA else hier.part[temp,4]/100
    }))
    contribution = c(contribution,summary(mods[[i]])$dev.expl)
    if(i == 1) output = contribution else output = rbind(output,contribution)
  }
  colnames(output) = c(predictors,"dev.expl")
  output = as.data.frame(output)
  output$weight = model.sel(get(divindicator))$weight
  output$indicator = divindicator
  output$names = names(mods)
  return(output)}


write.csv(rbind(calculate.contributions("richness"),
                calculate.contributions("abundance"),
                calculate.contributions("evenness")),
          paste(tab_out,"TABLE_S8_relativeimportance.csv",sep="/"))









#                       #---------------------------------#
#                       ###### Calculate predictions ######
#                       #---------------------------------#






### STEP 1: CALCULATE EFFECTS USING GRATIA
### When an effect was not present, we substituted by zero
getpredictions <- function(diversityindicator, output = "partial.effect"){
  # Get names and ranges of smooth terms and their variables
  dat = bruvs[!bruvs$island == "TI",]
  dat$island <- as.factor(dat$island)
  if(diversityindicator == "evenness") dat = dat[!is.na((((dat$evenness/dat$abundance)*(dat$abundance-1))+0.5)/dat$abundance),]
  globalmodel <- gam(richness ~ 
                       s(depth) + 
                       s(disttoshore) + 
                       s(slope) + 
                       habitat + 
                       island + 
                       s(Season, bs = 'cc', k=-1) +  
                       s(Season, by = island, bs = 'cc', k=-1),
                     data = dat,
                     family = poisson(link = "log"),
                     method = "REML")
  smoothterms <- smooths(globalmodel)
  ranges = lapply(smoothterms, function(sm){
    temp = as.data.frame(smooth_estimates(globalmodel,smooth = sm))
    if(!"island" %in% names(temp)) temp$island = NA
    temp$variable = colnames(temp)[6]
    temp = temp[,c(6:8)]
    colnames(temp) = c("x","island", "variable")
    return(temp)
    })
  names(ranges) = smoothterms
  mods = get.models(get(diversityindicator), subset = delta < 6)  
  for(i in 1:length(mods)){
   mod = mods[[i]] 
   for(j in 1:length(smoothterms)){
     if(smoothterms[j] %in% smooths(mod)) {
       fit = as.data.frame(smooth_estimates(mod,smooth = smoothterms[j]))
       if(!"island" %in% names(fit)) fit$island = NA
       fit$variable = colnames(fit)[6]
       colnames(fit)[6] = "x"
       fit = fit[,c(6,7,8,4,5)]
       resids = as.data.frame(partial_residuals(mod,select = smoothterms[j]))
     } else {
       fit = data.frame(ranges[[which(names(ranges) == smoothterms[j])]], a = 0, b = 0)
       names(fit)[c(4,5)] = c(".estimate", ".se")
       resids = data.frame(a = rep(0,length(predict(mod))))
     }
     colnames(resids) = "residuals"
     variable = ranges[[which(names(ranges) == smoothterms[j])]]$variable[1]
     resids$variable = variable
     resids$island   = ranges[[which(names(ranges) == smoothterms[j])]]$island[1]
     resids$x = unlist(globalmodel$model[which(names(globalmodel$model) == variable)])
     if(j == 1 & i ==1) fitoutput = fit else fitoutput = rbind(fitoutput,fit)
     if(j == 1 & i ==1) residsoutput = resids else residsoutput = rbind(residsoutput,resids)
   }
  }
  colnames(fitoutput) = c("x", "island", "variable", "fit", "se")
  newdat = expand.grid(habitat = c("Rock", "Maerl", "Sand"),
                       island = c("PC","ST"),
                       depth  = median(bruvs$depth),
                       slope  = median(bruvs$slope),
                       disttoshore = median(bruvs$disttoshore),
                       Season = c(1,30,60,90,120,150,180,210,240,270,310,330,360))
  preds = predict(object = get(diversityindicator), newdat, se.fit = TRUE)
  fitoutput = rbind(fitoutput, 
                    data.frame(x = newdat$island, 
                               island = NA, 
                               variable = "island",
                               fit = preds$fit,
                               se  = preds$se.fit))
  fitoutput = rbind(fitoutput, 
                    data.frame(x = newdat$habitat, 
                               island = NA, 
                               variable = "habitat",
                               fit = preds$fit,
                               se = preds$se.fit))
  fitoutput$variable
  fitoutput    = as.data.frame(fitoutput %>% group_by(variable, island, x)  %>% summarise(fit = mean(fit), se = mean(se)))
  residsoutput = as.data.frame(residsoutput %>% group_by(variable, island, x)  %>% summarise(residuals = mean(residuals)))
  habitatpoints = fitoutput$fit[which(fitoutput$variable == "habitat")]
  fitoutput$fit[which(fitoutput$variable == "habitat")] = 
    fitoutput$fit[which(fitoutput$variable == "habitat")] - (max(habitatpoints) + min(habitatpoints))/2
  islandpoints = fitoutput$fit[which(fitoutput$variable == "island")]
  fitoutput$fit[which(fitoutput$variable == "island")] = 
    fitoutput$fit[which(fitoutput$variable == "island")] - (max(islandpoints) + min(islandpoints))/2
  divindincators = data.frame(name = c("Richness", "Abundance", "Evenness"), 
                              div  = c("richness", "abundance", "evenness"))
  fitoutput$diversityindicator =  divindincators$name[divindincators$div == diversityindicator]
  residsoutput$diversityindicator = divindincators$name[divindincators$div == diversityindicator]
  fitoutput$ymin = fitoutput$fit - fitoutput$se
  fitoutput$ymax = fitoutput$fit + fitoutput$se
  if(output == "partial.effect") fitoutput else if(output == "residuals") residsoutput
}

pred_full = rbind(getpredictions("richness"),
             getpredictions("evenness"),
             getpredictions("abundance"))
pred = pred_full[!(is.na(pred_full$island) & pred_full$variable == "Season"),]

diversityindicator = c("Richness","Abundance","Evenness")
pred$diversityindicator <- factor(pred$diversityindicator,diversityindicator)




### STEP 2: PLOT PREDICTIONS
set_ylims = data.frame(x = c(1), fit =c(1),diversityindicator = factor("Richness",levels=diversityindicator))[-1,]
for(i in 1:length(diversityindicator))   set_ylims <- rbind(set_ylims,data.frame(
  fit   = c(min(pred$ymin[which(pred$diversityindicator == diversityindicator[i])])+
              (0.1*min(pred$ymin[which(pred$diversityindicator == diversityindicator[i])])),
            max(pred$ymax[which(pred$diversityindicator == diversityindicator[i])])+
              (0.1*max(pred$ymax[which(pred$diversityindicator == diversityindicator[i])]))),
  diversityindicator = factor(diversityindicator[i],levels=diversityindicator)
))



ploteffects <- function(variable,nostriptext = TRUE, noyaxis = TRUE){
  dat = pred[which(pred$variable == variable),]
  # CREATE VECTOR WITH LABELS
  variables     = c("depth",      "slope",    "disttoshore",            "habitat",     "Season",         "island")
  variablenames = c("Depth (m)",  "Slope (º)", "Dist. to shore (m)", "Habitat",     "Season (ord. date)",    "Island")
  dat$variablename = dat$variable
  for(i in 1:6) dat$variablename[which(dat$variable == variables[i])] = variablenames[i]
  
  # ADJUST X AXIS
  if(variable == "habitat")   dat$x = factor(dat$x, levels = c("Rock","Maerl","Sand"))
  if(variable == "island")    {dat$x[dat$x == "PC"] = "Príncipe"
  dat$x[dat$x == "ST"] = "São Tomé"
  dat$x = factor(dat$x, levels = c("Príncipe","São Tomé"))}
  if(variable == "Season") {dat$island = as.character(dat$island)
  dat$island[which(dat$island == "PC")] = "Príncipe"
  dat$island[which(dat$island == "ST")] = "São Tomé" 
  dat$island = factor(dat$island, levels = c("Príncipe", "São Tomé"))      
  }
  if(variable %in% c("depth","slope","disttoshore","Season")) dat$x = as.numeric(dat$x)
  
  # EXTRACT RELATIVE IMPORTANCE OF EACH EFFECT
  rimp = which(coefs$vars == variable)
  if(variable == "Season") rimp = which(grepl("Season",coefs$vars))
  if(variable == "depth") rimp = which(grepl("depth",coefs$vars))
  if(variable == "slope") rimp = which(grepl("slope",coefs$vars))
  if(variable == "disttoshore") rimp = which(grepl("disttoshore",coefs$vars))
  if(variable == "habitat") rimp = rimp[c(2,4,6)]
  rimp = data.frame(rimp = coefs$relativeimportance[rimp], diversityindicator = coefs$indicator[rimp])
  for(i in 1:3) rimp$diversityindicator[rimp$diversityindicator == c("richness","abundance","evenness")[i]] = c("Richness","Abundance","Evenness")[i]
  rimp$diversityindicator = factor(rimp$diversityindicator,levels = c("Richness","Abundance","Evenness"))
  rimp$rimp = round(rimp$rimp,2)
  for(i in 1:nrow(rimp)) if(nchar(rimp$rimp[i])==1) rimp$rimp[i] = paste0(rimp$rimp[i],".00") else 
    if(nchar(rimp$rimp[i])==3)    rimp$rimp[i] = paste0(rimp$rimp[i],"0")
  rimp$variablename = dat$variablename[1]
  rimp$x = 0.5
  if(variable %in% c("depth","slope","disttoshore","Season")) rimp$x = max(dat$x)
  if(variable == "island") rimp$x = 2.5
  if(variable == "habitat") rimp$x[which(!rimp$diversityindicator == "Evenness")] = 3.5
  rimp$ymax = 0
  rimp$ymin = 0
  rimp$fit = unlist(lapply(rimp$diversityindicator,function(x) max(set_ylims$fit[set_ylims$diversityindicator == x])-
                             (max(set_ylims$fit[set_ylims$diversityindicator == x])*0.07)))
  if(variable == "habitat") hjust = c(1,1,0) else hjust = 1
  if(variable == "Season") rimp$fit[c(2,4,6)] = rimp$fit[c(2,4,6)]- (rimp$fit[c(2,4,6)]*0.195)
  if(variable == "Season") {rimp$rimp[c(1,3,5)] = paste("SW =",rimp$rimp[c(1,3,5)])
  rimp$rimp[c(2,4,6)] = paste("SW (by isl.) =",rimp$rimp[c(2,4,6)])} else 
  {rimp$rimp = paste("SW =",rimp$rimp)}
  rimp$island = "Príncipe"
  
  # SET Y LIMITS
  set_ylims$x = rep(dat$x[1],nrow(set_ylims))
  set_ylims$ymax = 0
  set_ylims$ymin = 0
  set_ylims$island = "Príncipe"
  # CREATE PLOTS
  if(variable %in% c("habitat","island"))
    output = ggplot(data = dat, aes(x = x, y = fit, ymin = ymin, ymax = ymax)) + 
    geom_point() + geom_errorbar(width = 0) 
  if(variable %in% c("depth","slope","disttoshore"))
    output = ggplot(data = dat, aes(x = x, y = fit, ymin = ymin, ymax = ymax)) + 
    geom_ribbon(alpha = 0.6) + geom_line() 
  if(variable == "Season")
    output = ggplot(data = dat, aes(x = x, y = fit, ymin = ymin, ymax = ymax, fill = island, colour = island)) + 
    geom_line() + geom_ribbon(alpha = 0.6 ) +
    scale_x_continuous(breaks = c(15, 46, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350), 
                       labels = c("J",  "F",  "M",  "A",  "M",   "J",   "J",   "A",   "S",   "O",   "N", "D")) + 
    theme(axis.ticks.major.x = element_blank(),
          axis.ticks.minor.x = element_blank())
  output = output + facet_grid(diversityindicator~variablename, scales = "free") + 
    geom_point(data = set_ylims,aes(x = x, y = fit),alpha = 0) + 
    geom_text(data = rimp,aes(x=x,y=fit, label = rimp), hjust = hjust, size =3, colour = "black") +
    xlab(NULL) + ylab("Effect") +
    theme_bw() + theme(
      panel.grid = element_blank(),
      legend.position = c(0,0),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.justification.inside = c(0, 0),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.35, "cm"),
      plot.margin = unit(c(0.1,0.05,0.1,0.05), "cm") )
  if(nostriptext) output = output + theme(strip.text.y = element_blank())
  if(noyaxis)     output = output + ylab(NULL) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  if(!variable == "Season") output = output + theme(panel.grid.minor = element_blank())
  return(output)}




wdth = 21
hgth = 12.5
tiff(paste0(fig_out,"/Fig3.tiff"), unit = "cm", width = wdth, height =hgth, res = 500)
grid.arrange(ploteffects("depth", noyaxis = FALSE), 
             ploteffects("slope"),
             ploteffects("disttoshore"),
             ploteffects("habitat"),
             ploteffects("island"),
             ploteffects("Season", nostriptext = FALSE),ncol=6, widths = c(1.35,1,1,1,1,1.19))
dev.off()




### STEP 3: PLOT PREDICTIONS WITH RESIDUALS
resids <- rbind(getpredictions("richness", output = "residuals"),
                getpredictions("evenness", output = "residuals"),
                getpredictions("abundance", output = "residuals"))



dat = pred_full
dat$x = as.numeric(dat$x)
dat <- dat[-which(dat$variable %in% c("habitat", "island")),]
dat$island[is.na(dat$island)] <- ""
dat$variable = paste0(dat$variable,dat$island)



resids$x = as.numeric(resids$x)
resids$island[is.na(resids$island)] <- ""
resids$variable = paste0(resids$variable,resids$island)
resids$ymax = 0
resids$ymin = 0
resids <- resids[!(resids$residuals > 10 | resids$residuals < -10),]



variables <- c("depth","disttoshore","Season","SeasonPC","SeasonST","slope")
variablenames <- c("Depth (m)", "Dist. to shore (m)", "Season\n(ordinal date)", 
                  "Season:Príncipe\n(ordinal date)", "Season:SãoTomé\n(ordinal date)", "Slope (º)")

for(i in 1:length(variables))
  resids$variable[resids$variable == variables[i]] <- variablenames[i]

for(i in 1:length(variables))
  dat$variable[dat$variable == variables[i]] <- variablenames[i]





wdth = 21
hgth = 12.5
tiff(paste0(fig_out,"/FigS4_GAMdetail.tiff"), unit = "cm", width = wdth, height =hgth, res = 500)

ggplot(data = dat, 
       aes (x = x, y = fit, ymin = ymin, ymax = ymax)) +
  geom_point(dat = resids, aes(x=x, y = residuals), size = 1, alpha = 0.7, colour = "darkgrey") +
  geom_ribbon() + 
  geom_line() +
  facet_grid(diversityindicator~variable, scales = "free") + xlab("") + ylab("Effect")
dev.off()






