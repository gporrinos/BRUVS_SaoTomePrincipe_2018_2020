
################################################################################
################################################################################
###                            DEPLOYMENT SUMMARY                            ###
################################################################################
################################################################################








         #############################################################
         ###                       LOAD DATA                       ###
         #############################################################


source(file.path(getwd(),"data_analysis/script/01_data_and_directories.R"))
## Exclude afternoon repeats (they were tests and were not used for analysis)
bruvs   <-  bruvs[!bruvs$sampling_time=="afternoon_repeat",]    







         #############################################################
         ###                    DEPLOYMENT SUMMARY                 ###
         #############################################################



### STEP 1: SUBSET DATABASES
ST <- bruvs[bruvs$island == "ST",]
PC <- bruvs[bruvs$island == "PC",]
TI <- bruvs[bruvs$island == "TI",]




### STEP 2: FUNCTION TO SUMMARISE DEPLOYMENTS
deployments <- function(data){ # Write function to summarise deployments
  date = paste(substr(data$date,7,19),substr(data$date,4,5),substr(data$date,1,2),sep="")      
  retained = data[-which(data$sampling_type=="not_valid"),]  
  return(data.frame(totaldeployments = nrow(data),
                    discarded   = nrow(data[which(data$sampling_type=="not_valid"),]),
                    retained    = nrow(data[-which(data$sampling_type=="not_valid"),]),
                    random      = nrow(data[which(data$sampling_type=="random"),]),
                    non_random  = nrow(data[which(data$sampling_type=="non-random"),]),
                    Rock        = nrow(retained[which(retained$habitat=="Rock"),]),
                    Maerl       = nrow(retained[which(retained$habitat=="Maerl"),]),
                    Sand        = nrow(retained[which(retained$habitat=="Sand"),]),
                    start       = data[which(date==min(date)),which(colnames(data)=="date")][[1]],
                    end         = data[which(date==max(date)),which(colnames(data)=="date")][[1]],
                    nsamplingdays = length(levels(as.factor(date)))))}





### STEP 3: SUMMARY TABLE OF DEPLOYMENTS, PER SEASON
deployment.sum <- deployments(PC)[-1,]
for(i in 1:3) {
  island         = c("PC","TI","ST")[i]
  dat            = get(island)
  dat$sampling_round = paste(dat$sampling_round,dat$season)
  temp           = deployments(dat)
  rownames(temp) <- c(island)
  deployment.sum <- rbind(deployment.sum,temp)
  for(i in 1:length(levels(as.factor(dat$sampling_round)))){
    sampling_round = levels(as.factor(dat$sampling_round))[i]
    temp           = deployments(dat[which(dat$sampling_round==sampling_round),])
    rownames(temp) = c(sampling_round)
    deployment.sum = rbind(deployment.sum,temp)}}




### STEP 4: WRITE OUTPUT TABLE AS A CSV
deployment.sum # Table that summarises all deployments
write.csv(deployment.sum,paste(tab_out,"TABLE_S1_deployment_summary.csv",sep="/"))




### STEP 5: CAUSES OF DEPLOYMENT FAILURE
fail <- levels(as.factor(bruvs[bruvs$sampling_type =="not_valid",
                               which(colnames(bruvs)=="failure")]))
cause_of_failure <- data.frame()
for(i in 1:length(fail)){cause_of_failure = rbind(cause_of_failure,
                                                  unlist(lapply(list(ST,PC,TI),function(dat){
                                                    length(dat$failure[which(dat$failure == fail[i])])})))}
rownames(cause_of_failure) <- fail
colnames(cause_of_failure) <- c("ST","PC","TI")
cause_of_failure

rm(dat,PC,ST,TI,temp,i,deployments,island,fail,sampling_round)
