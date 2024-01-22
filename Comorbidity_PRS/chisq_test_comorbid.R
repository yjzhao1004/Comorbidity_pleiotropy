library(dplyr)
setwd("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/Sur_data/")
rm(list = ls())

brain_disorder_collect<-c("AD","ADHD","ALD","ALS","ANX","ASD","BP","EAT","MDD","MS","OCD","PD","PTSD","SCZ","TS")
chronicphy_collect<-c("MCP","ASM","CHD","CLD","HTN","PUD","RA","Stroke","T2D")

library(data.table)
n=1
resultslist = list()
for ( i in c(1:length(brain_disorder_collect))){
  current_BD<-fread(paste("Sur_",brain_disorder_collect[i],"_20231030.csv",sep = ""))
  current_BD<-current_BD[,c("eid","baselinedate","status")]
  names(current_BD)[3]<-"BD_status"
  for(j in c(1:length(chronicphy_collect))){
    current_CP<-fread(paste("Sur_",chronicphy_collect[j],"_20231030.csv",sep = ""))
    current_CP<-current_CP[,c("eid","baselinedate","status")]
    names(current_CP)[3]<-"CP_status"
    intersect_eid<-intersect(current_BD$eid,current_CP$eid)
    
    current_BD<-subset(current_BD,eid%in%intersect_eid)
    current_CP<-subset(current_CP,eid%in%intersect_eid)
    
    current_people<-merge(current_BD[,c(1,3)],current_CP[,c(1,3)],by = "eid")
    
    results<-xtabs(~current_people$BD_status+current_people$CP_status,data=current_people)%>%chisq.test()
    table<-xtabs(~current_people$BD_status+current_people$CP_status,data=current_people)
    
   
    resultslist$Brain_disorder[n]<-brain_disorder_collect[i]
    resultslist$Chronic_physical_condition[n]<-chronicphy_collect[j]
    
    resultslist$comorbid_sample_0_0[n]<-table[1,1]
    resultslist$comorbid_sample_0_1[n]<-table[1,2]
    resultslist$comorbid_sample_1_0[n]<-table[2,1]
    resultslist$comorbid_sample_1_1[n]<-table[2,2]
    resultslist$X_squared[n]<-results$statistic
    resultslist$p_value[n]<-results$p.value
    n=n+1
    rm(results,table)
  }
}
resultslist<-as.data.frame(resultslist)
write.csv(resultslist,file = "/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/chisquared/resultslist_20231030.csv",row.names = FALSE)


### fdr correction
resultslist$id<-c(1:135)
Data = resultslist[order(resultslist[,8]),]
Data$fdr =
  p.adjust(Data$p_value,
           method = "fdr")

Data = Data[order(Data[,9]),]
write.csv(Data,file = "/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/chisquared/resultslist_fdr_20231030.csv",row.names = FALSE)


