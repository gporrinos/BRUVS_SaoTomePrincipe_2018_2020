
################################################################################
################################################################################
###                   LOAD DATABASES AND CREATE DIRECTORIES                  ###
################################################################################
################################################################################







# Author: Guillermo Prieto Porrinos
# Database: BRUVS 
#            PC2018a, PC2018b
#                      Collected by: University of Exeter, Fundacao Principe
#                      Funded by   : Darwin Initiative, Forever Príncipe, Halpin Trust
#            PC2019, PC2020, ST2019, ST2020
#                      Collected by: Fauna & Flora International, Fundacao Príncipe, Oikos Coop. & Desenv., MARAPA
#                      Funded by   : Blue Action Fund, Arcadia Fund







#         #############################################################
#         ###                      DIRECTORIES                      ###
#         #############################################################



## Directories exist? If not, create them
fig_out    <- file.path(getwd(), "data_analysis", "outputs", "figures")
if(!file.exists(fig_out)) dir.create(fig_out)

tab_out    <- file.path(getwd(), "data_analysis", "outputs", "tables")
if(!file.exists(tab_out)) dir.create(tab_out)






         #############################################################
         ###                       DATABASES                       ###
         #############################################################


bruvs   <-  read.table(file.path(getwd(),"data", "videos.csv"),header=TRUE,sep=",")
maxn    <- read.table(file.path(getwd(),"data","maxn.csv"),header=TRUE,sep=",")
species <- read.table(file.path(getwd(),"data","species.csv"),header=TRUE,sep=",")






         #############################################################
         ###                  COLOURS FOR FIGURES                  ###
         #############################################################


colourPC <- "#F8766D"
colourST <- "#00BFC4"
colourTI <- "gray20"

colourSand <- "goldenrod3"
colourMaerl <- "dodgerblue3"
colourRock <- "gray15"

