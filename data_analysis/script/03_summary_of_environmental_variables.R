
################################################################################
################################################################################
###                    SUMMARY OF ENVIRONMENTAL VARIABLES                    ###
################################################################################
################################################################################





### STEP 1: LOAD DATABASES
source(file.path(getwd(),"data_analysis/script/02_variable_definition.R"))

### step 2: LOAD PACKAGES
if(!require("ggplot2")) install.packages("ggplot2",   repos = "https://cloud.r-project.org")






#                #--------------------------------------------#
#                ##### Summary of environmental variables #####
#                #--------------------------------------------#



#### STEP 1: FUNCTION TO SUMMARISE ENVIRONMENTAL VARIABLES
summary.env <- function(variable){
  temp = lapply(c(1:3),function(i){
    dat = bruvs[bruvs$island == c("PC","TI","ST")[i],      which(colnames(bruvs)==variable)]
    return(data.frame(variable = variable,     island   = c("Principe","Tinhosas","São Tomé")[[i]],
                      max = max(dat),     min = min(dat),     mean = mean(dat),     sd = sd(dat)))
  })
  output = temp[[1]]
  for(i in c(1,2)){   output = rbind(output,temp[[i+1]])   }
  return(output)}





#### STEP 2: EXPORT SUMMARY OF ENVIRONMENTAL VARIABLES
(environmental.variables <- rbind(summary.env("depth"),summary.env("disttoshore"),summary.env("slope")))
write.csv(environmental.variables,paste(tab_out,"TABLE_S2_summary_environmental_variables.csv",sep="/"))




### STEP 3: DEPTH AND SLOPE HISTOGRAMS
bruvs$Habitat <- factor(bruvs$habitat, levels = c("Sand", "Maerl", "Rock"))
bruvs$island  <- bruvs$island
bruvs$Island[bruvs$island == "PC"] <- "Príncipe"
bruvs$Island[bruvs$island == "ST"] <- "São Tomé"
bruvs$Island[bruvs$island == "TI"] <- "Tinhosas"


ggtheme <- theme_bw() + theme(panel.grid = element_blank())

p2 <- ggplot(data = data.frame(bruvs, var = "Slope (º)"), aes(x = Habitat, y = slope, fill = Habitat)) + 
  geom_boxplot(alpha = 0.6) +  facet_grid(Island ~ var) + ylab("Slope (º)") + 
  scale_fill_manual(values = c(colourSand, colourMaerl, colourRock)) +
  ggtheme

p1 <- ggplot(data = data.frame(bruvs, var = "Depth (m)"), aes(x = Habitat, y = depth, fill = Habitat)) + 
  geom_boxplot(alpha = 0.6) +  facet_grid(Island ~ var) + ylab("Depth (m)") + 
  scale_fill_manual(values = c(colourSand, colourMaerl, colourRock)) +
  ggtheme



wdth = 15
hgth = 15
tiff(paste0(fig_out,"/FigS3_depth_slope.tiff"), unit = "cm", width = wdth, height =hgth, res = 500)
gridExtra::grid.arrange(p1 + theme(legend.position = "none"),p2,ncol=2, widths = c(1,1.4))
dev.off()
rm(summary.env,p1,p2)



