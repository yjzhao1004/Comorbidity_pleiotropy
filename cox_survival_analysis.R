library(data.table)
library(dplyr)
library(survival)
rm(list = ls())
#load disease name
files<-list.files("/home1/yjzhao/COMMO_new/PRS_CS/outcome_event/outcome_finaluse_csv/")
disease_kind<-as.numeric()
for(i in c(1:length(files))){
  disease_kind[i]<-unlist(strsplit(files[i],"_"))[1]
}

somatic_condition<-c("Cirrhosis","CKD","COPD","Diabetes","HTN","ISH","Osteoarthritis","Osteoporosis")
psy_disorder<-c("ADHD","ANX","ASD","BP","CUD","EAT","MDD","OCD","PTSD","SCZ")

#load genomic ancestry
genomic_ancestry<-fread("/home1/yjzhao/COMMO_new/PRS_CS/genetic_ancestry.csv")%>%
  filter(!is.na(caucasian))%>%
  select(eid,caucasian)

#load cov
load("/home1/yjzhao/lifestyle_depression/SEM_1108/SEM_1108/cov_sem_impute_1109.RData")

###from somatic to psy
for(i in c(1:length(somatic_condition))){
  print(somatic_condition[i])
  somatic_c<-somatic_condition[i]
  
  somaticprs<-fread(paste0("/home1/yjzhao/COMMO_new/BOLT_LMM/PGS_somatic_score/",somatic_c,".sscore"))%>%
    rename(eid = IID,prs = BETA_SUM)%>%
    select(eid,prs)
  somaticprs$zscoreprs<-scale(somaticprs$prs)
  
  Quan<-quantile(somaticprs$zscoreprs)
  somaticprs$class<-NA
  somaticprs$class[which(somaticprs$zscoreprs<Quan[2])]<-1
  somaticprs$class[which(somaticprs$zscoreprs>Quan[4])]<-3
  somaticprs$class[which(somaticprs$zscoreprs>=Quan[2]&somaticprs$zscoreprs<=Quan[4])]<-2
  
  somaticprs<-subset(somaticprs,eid %in% genomic_ancestry$eid)
  for(j in c(1:length(psy_disorder))){
    psy_d<-psy_disorder[j]
    
    outcomedata<-fread(paste0("/home1/yjzhao/COMMO_new/PRS_CS/outcome_event/outcome_finaluse_csv/",psy_d,"_outcome.csv"))
    
    outcomedata<-merge(outcomedata,cov_sem,by = "eid")
    outcomedata<-merge(outcomedata,somaticprs,by = "eid")
    
    
    res.cox <- coxph(Surv(days,status) ~ factor(class, level = c(1:3)) + Sex + Townsend_index + BMI  + Age + Qualifications, data = outcomedata)
    cox.zph(res.cox)
    X<-summary(res.cox)
    
    ###create collection csv
    collection<-list()
    #reference
    t = 1
    collection$prsdisease[t]<-somatic_c
    collection$outcome[t]<-psy_d
    collection$prsclass[t]<-"lowPRS"
    ss<-subset(outcomedata,class==1)
    collection$NcaNcon[t]<-paste0(sum(ss$status),"/",length(which(ss$status==0)))
    collection$HR[t]<-1
    collection$LowerCI[t]<-"/"
    collection$HigherCI[t]<-"/"
    collection$pval[t]<-"/"
    
    t = 2
    collection$prsdisease[t]<-somatic_c
    collection$outcome[t]<-psy_d
    collection$prsclass[t]<-"MedPRS"
    ss<-subset(outcomedata,class==2)
    collection$NcaNcon[t]<-paste0(sum(ss$status),"/",length(which(ss$status==0)))
    collection$HR[t]<-X$conf.int[1,1]
    collection$LowerCI[t]<-X$conf.int[1,3]
    collection$HigherCI[t]<-X$conf.int[1,4]
    collection$pval[t]<-X$coefficients[1,5]
    
    t = 3
    collection$prsdisease[t]<-somatic_c
    collection$outcome[t]<-psy_d
    collection$prsclass[t]<-"HighPRS"
    ss<-subset(outcomedata,class==3)
    collection$NcaNcon[t]<-paste0(sum(ss$status),"/",length(which(ss$status==0)))
    collection$HR[t]<-X$conf.int[2,1]
    collection$LowerCI[t]<-X$conf.int[2,3]
    collection$HigherCI[t]<-X$conf.int[2,4]
    collection$pval[t]<-X$coefficients[2,5]
    
    collection<-as.data.frame(collection)
    
    write.csv(collection,file = paste0("/home1/yjzhao/COMMO_new/PRS_CS/outcome_event/Results_s_p/",somatic_c,"_PRS_pred_",psy_d,".csv"),row.names = F)
    
    rm(collection)
  }
  
}

#combind csv
setwd("/home1/yjzhao/COMMO_new/PRS_CS/outcome_event/Results_s_p/")
folder_path<-"/home1/yjzhao/COMMO_new/PRS_CS/outcome_event/Results_s_p//"
file_suffix <- ".csv" 
file_list <- list.files(folder_path, pattern = file_suffix, full.names = TRUE)
combined_data <- data.frame()
for (file in file_list) {
  data <- read.csv(file, header = TRUE)  
  combined_data <- rbind(combined_data, data)
}
write.csv(combined_data, "combined_somatic_pred_psy.csv", row.names = FALSE)


