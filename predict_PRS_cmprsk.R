rm(list = ls())
library(rms)
library(survival)
library(mice)
library(mitools)
library(Hmisc)
library(dplyr)
library(data.table)
library("cmprsk")
load("/Users/yujiezhao/Desktop/Lifestyle_Depression/SEM_1108/cov_sem_impute_1109.RData")
setwd("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/cmprsk_results/")


####ADD PRS

diseasename<-"MDD"
PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/ISS_PRS.csv")

#PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/lungfunction/lungfunction.all_score")
PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/chronicpain.all_score")
PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/depression_PRS.all_score")
#names(PRS)[1]<-"eid"
PRS<-PRS[,c(2,3)]
names(PRS)[2]<-"PRS"

#PRS<-PRS[,c(1,13)]
#names(PRS)[2]<-"PRS"


PRS<-na.omit(PRS)
#10 quan
Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))

PRS$class<-NA
PRS$class[which(PRS$PRS<Quan[2])]<-1
PRS$class[which(PRS$PRS>Quan[9])]<-3
PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2

length(which(PRS$class==2))
length(which(PRS$class==1))
length(which(PRS$class==3))

brain_disorder_collection<-c("AD","ADHD","ALD","ALS","ANX","ASD","BP","EAT",
                             "MDD","MS","OCD","PD","PTSD","SCZ","TS")

brain_disorder_collection<-c("AD","BP","MS","PD","SCZ")

brain_disorder_collection<-c("ALD","EAT","MDD","PTSD","OCD","ANX","ALS","ASD","ADHD","TS")


chronicphy_collection<-c("ASM","T2D","MCP","HTN","CHD","Stroke","CLD","RA","PUD")

chronicphy_collection<-c("JOINT","SPINE","LIMB","BACK","DENTAL_TMD","PUD")
chronicphy_collection<-c("JOINT","SPINE","LIMB","BACK","DENTAL_TMD")

chronicphy_collection<-c("CLD","MCP")

##cox ####
for(i in c(1:15)){
  braindisorder<-brain_disorder_collection[i]
  print(braindisorder)
  
  current_BD<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/sur_data/Sur_",brain_disorder_collection[i],"_20231010.csv",sep = ""))
  names(current_BD)[5]<-"BD_status"
  names(current_BD)[6]<-"BD_dia_days"
  names(current_BD)[4]<-"BD_group"
  current_CP<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/sur_data/Sur_",diseasename,"_20231010.csv",sep = ""))
  names(current_CP)[5]<-"CP_status"
  names(current_CP)[6]<-"CP_dia_days"
  names(current_CP)[4]<-"CP_group"
  
  intersect_eid<-intersect(current_BD$eid,current_CP$eid)
  
  current_BD<-subset(current_BD,eid%in%intersect_eid)
  current_CP<-subset(current_CP,eid%in%intersect_eid)
  
  current_people<-merge(current_BD[,c(1,4,5,6)],current_CP[,c(1,4,5,6)],by = "eid")
  
  current_people<-merge(current_people,PRS,by = "eid")
  current_people$comorbid<-0
  current_people$comorbid[which(current_people$BD_status==1&current_people$CP_status==1)]=1
  #current_people$comorbid[which(current_people$BD_status==1&current_people$CP_status==1)]=1
  #current_people$comorbid[which(current_people$BD_status==1&current_people$CP_status==0)]=2
  #current_people$comorbid[which(current_people$BD_status==0&current_people$CP_status==1)]=2
  #current_people$comorbid[which(current_people$BD_status==0&current_people$CP_status==0)]=0
  #current_people$comorbid[which(current_people$BD==0&current_people$CP_status==0)]=0
  
  
  current_people$comorbid_days<-NA
  current_people$comorbid_days[which(current_people$BD_status==1&current_people$CP_status==1)]=pmax(current_people$BD_dia_days[which(current_people$BD_status==1&current_people$CP_status==1)],current_people$CP_dia_days[which(current_people$BD_status==1&current_people$CP_status==1)])
  current_people$comorbid_days[which(current_people$BD_status==1&current_people$CP_status==0)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==1&current_people$CP_status==0)],current_people$CP_dia_days[which(current_people$BD_status==1&current_people$CP_status==0)])
  current_people$comorbid_days[which(current_people$BD_status==0&current_people$CP_status==1)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==0&current_people$CP_status==1)],current_people$CP_dia_days[which(current_people$BD_status==0&current_people$CP_status==1)])
  current_people$comorbid_days[which(current_people$BD_status==0&current_people$CP_status==0)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==0&current_people$CP_status==0)],current_people$CP_dia_days[which(current_people$BD_status==0&current_people$CP_status==0)])
  
  brain_sur_PRS<-merge(current_people,cov_sem,by ="eid")
  
  res.cox <- coxph(Surv(comorbid_days,comorbid) ~ factor(class, level = c(1,2,3)) + Sex + Townsend_index + BMI  + Age +Qualifications, data = brain_sur_PRS)
  #cox.zph(res.cox)
  X<-summary(res.cox)
  X
  
  braind_cox = list()
  t=1
  
  braind_cox$variable[t]<-"low_PRS"
  braind_cox$outcome[t]<-braindisorder
  
  braind_cox$HR[t]<-1
  braind_cox$lowerCI[t]<-" "
  braind_cox$upperCI[t]<-" "
  #braind_cox$num[t]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
  ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==1)]
  braind_cox$num2[t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
  braind_cox$hrci[t]<-" "
  braind_cox$pvalue[t]<-" "
  
  
  zz=1
  z=1
  braind_cox$variable[z+t]<-"intermediate_PRS"
  braind_cox$outcome[z+t]<-braindisorder
  braind_cox$HR[z+t]<-X$conf.int[zz,1]
  braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
  braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
  #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
  ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==2)]
  braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
  braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
  braind_cox$pvalue[z+t]<-X$coefficients[zz,5]
  
  
  zz = zz+1
  z = z+1    
  braind_cox$variable[z+t]<-"High_PRS"
  braind_cox$outcome[z+t]<-braindisorder
  braind_cox$HR[z+t]<-X$conf.int[zz,1]
  braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
  braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
  #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
  ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==3)]
  braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
  braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
  braind_cox$pvalue[z+t]<-X$coefficients[zz,5]
  
  
  braind_cox<-as.data.frame(braind_cox)
  filename<-paste(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_predict/",diseasename,"_",braindisorder,"_cox_PRS_20231010.csv",sep = ""))
  write.csv(braind_cox,file = filename,row.names = FALSE)
  
}




dirpath<-"/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_predict/"
files<-list.files(dirpath)
collection = list()
for(j in c(33:47)){
  collection_one<-read.csv(paste(dirpath,"/",files[j],sep = ""))
  collection_one$name<-files[j]
  collection<-rbind(collection,collection_one)
}
write.csv(collection,file = "CLD_PRS_results_allcomorbid1_20231010.csv")


####cmprsk ####

####ADD PRS


for (j in c(1:length(chronicphy_collection))){
  diseasename<-chronicphy_collection[j]
  
  print(diseasename)
  if (diseasename=="ASM"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/ASM_PRS.csv")
  }else if (diseasename=="CHD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/CAD_PRS.csv")
  }else if (diseasename=="HTN"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/HT_PRS.csv")
  }else if (diseasename=="Stroke"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/ISS_PRS.csv")
  }else if (diseasename=="RA"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/RA_PRS.csv")
  }else if (diseasename=="MCP"){
      PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/MCP.all_score")
  }else if (diseasename=="T2D"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/T2D_PRS.csv")
  }else if (diseasename=="PUD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PUD/PUDfinn.all_score")
  }else if (diseasename=="SCZ"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/SCZ_PRS.csv")
  }else if (diseasename=="AD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/AD_PRS.csv")
  }else if (diseasename=="BP"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/BP_PRS.csv")
  }else if (diseasename=="MS"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/MS_PRS.csv")
  }else if (diseasename=="PD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/PD_PRS.csv")
  }else if (diseasename=="EAT"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/EAT/EAT.all_score")
  }else if (diseasename=="ANX"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/ANX/ANX.all_score")
  }else if (diseasename=="PTSD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PTSD/PTSD.all_score")
  }else if (diseasename=="ADHD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/ADHD/ADHD.all_score")
  }else if (diseasename=="TS"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/TS/TS.all_score")
  }else if (diseasename=="MDD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/MDD/MDD.all_score")
  }else if (diseasename=="ALS"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/ALS/ALS.all_score")
  }else if (diseasename=="ASD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/ASD/ASD.all_score")
  }else if (diseasename=="OCD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/OCD/OCD.all_score")
  }else if (diseasename=="ALD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/ALD.all_score")
  }else if (diseasename=="CLD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/lungfunction/CLD_noukb.all_score")
  }else if (diseasename=="JOINT"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/FinnGen_gwas/JOINT.all_score")
  }else if (diseasename=="SPINE"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/FinnGen_gwas/SPINE.all_score")
  }else if (diseasename=="LIMB"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/FinnGen_gwas/LIMB.all_score")
  }else if (diseasename=="BACK"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/FinnGen_gwas/BACK.all_score")
  }else if (diseasename=="DENTAL_TMD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/FinnGen_gwas/DENTAL_TMD.all_score")
  }
  
  
  names(PRS)[1]<-"eid"
  PRS<-PRS[,c(1,4)]
  names(PRS)[2]<-"PRS"
  
  PRS<-na.omit(PRS)
  #10 quan
  Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
  
  PRS$class<-NA
  PRS$class[which(PRS$PRS<Quan[2])]<-1
  PRS$class[which(PRS$PRS>Quan[9])]<-3
  PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
  
  
  for(i in c(1:15)){
    braindisorder<-brain_disorder_collection[i]
    print(braindisorder)
    #chronicphy<-chronicphy_collection[i]
    #print(chronicphy)
   
    current_BD<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/Sur_data/Sur_",braindisorder,"_20231030.csv",sep = ""))
    
    names(current_BD)[5]<-"BD_status"
    names(current_BD)[6]<-"BD_dia_days"
    names(current_BD)[4]<-"BD_group"
    #current_CP<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/Sur_data/Sur_",diseasename,"_20231030.csv",sep = ""))
    current_CP<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/Sur_data/Sur_MCP_20231030.csv",sep = ""))
    
    names(current_CP)[5]<-"CP_status"
    names(current_CP)[6]<-"CP_dia_days"
    names(current_CP)[4]<-"CP_group"
    
    intersect_eid<-intersect(current_BD$eid,current_CP$eid)
    
    current_BD<-subset(current_BD,eid%in%intersect_eid)
    current_CP<-subset(current_CP,eid%in%intersect_eid)
    
    current_people<-merge(current_BD[,c(1,4,5,6)],current_CP[,c(1,4,5,6)],by = "eid")
    
    current_people<-merge(current_people,PRS,by = "eid")
    current_people$comorbid<-NA
    current_people$comorbid[which(current_people$BD_group=="control"&current_people$CP_group=="control")]=0
    current_people$comorbid[which(current_people$BD_group=="Death_Control"&current_people$CP_status==1)]=2
    current_people$comorbid[which(current_people$BD_group=="Death_Control"&current_people$CP_group=="Death_Control")]=2
    current_people$comorbid[which(current_people$BD_group=="control"&current_people$CP_status==1)]=0
    current_people$comorbid[which(current_people$BD_status==1&current_people$CP_group=="control")]=0
    current_people$comorbid[which(current_people$BD_status==1&current_people$CP_group=="Death_Control")]=2
    current_people$comorbid[which(current_people$BD_status==1&current_people$CP_status==1)]=1
    
    current_people$comorbid_days<-NA
    current_people$comorbid_days[which(current_people$BD_status==1&current_people$CP_status==1)]=pmax(current_people$BD_dia_days[which(current_people$BD_status==1&current_people$CP_status==1)],current_people$CP_dia_days[which(current_people$BD_status==1&current_people$CP_status==1)])
    current_people$comorbid_days[which(current_people$BD_status==1&current_people$CP_status==0)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==1&current_people$CP_status==0)],current_people$CP_dia_days[which(current_people$BD_status==1&current_people$CP_status==0)])
    current_people$comorbid_days[which(current_people$BD_status==0&current_people$CP_status==1)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==0&current_people$CP_status==1)],current_people$CP_dia_days[which(current_people$BD_status==0&current_people$CP_status==1)])
    current_people$comorbid_days[which(current_people$BD_status==0&current_people$CP_status==0)]<-pmax(current_people$BD_dia_days[which(current_people$BD_status==0&current_people$CP_status==0)],current_people$CP_dia_days[which(current_people$BD_status==0&current_people$CP_status==0)])
    
    brain_sur_PRS<-merge(current_people,cov_sem,by ="eid")
    
    attach(brain_sur_PRS)
    cov<-model.matrix(~factor(class, level = c(1,2,3))+Sex + Townsend_index + BMI  + Age +Qualifications)[,-1]
    fit3 <- crr(brain_sur_PRS$comorbid_days, brain_sur_PRS$comorbid, cov1=cov, failcode=1, cencode=0)
    X<-summary(fit3)
    
    braind_cox = list()
    t=1
    
    braind_cox$variable[t]<-"low_PRS"
    braind_cox$outcome[t]<-braindisorder
    
    braind_cox$HR[t]<-1
    braind_cox$lowerCI[t]<-" "
    braind_cox$upperCI[t]<-" "
    #braind_cox$num[t]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==1)]
    braind_cox$num2[t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[t]<-" "
    braind_cox$pvalue[t]<-" "
    
    
    zz=1
    z=1
    braind_cox$variable[z+t]<-"intermediate_PRS"
    braind_cox$outcome[z+t]<-braindisorder
    braind_cox$HR[z+t]<-X$conf.int[zz,1]
    braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
    braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
    #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==2)]
    braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
    braind_cox$pvalue[z+t]<-X$coef[zz,5]
    
    
    zz = zz+1
    z = z+1    
    braind_cox$variable[z+t]<-"High_PRS"
    braind_cox$outcome[z+t]<-braindisorder
    braind_cox$HR[z+t]<-X$conf.int[zz,1]
    braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
    braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
    #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==3)]
    braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
    braind_cox$pvalue[z+t]<-X$coef[zz,5]
    
    
    braind_cox<-as.data.frame(braind_cox)
    filename<-paste(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/cmprsk_results/cmprsk_1103/",diseasename,"pain_",braindisorder,"_cox_PRS_20231103.csv",sep = ""))
    #filename<-paste(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_cmprsk_20231015/",diseasename,"_",chronicphy_collection[i],"_cox_PRS_20231015.csv",sep = ""))
    
    write.csv(braind_cox,file = filename,row.names = FALSE)
    
  }
  
}
rm(list = ls())
dirpath<-"/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/cmprsk_results/cmprsk_1103"
files<-list.files(dirpath)
collection = list()
for(j in c(1:105)){
  collection_one<-read.csv(paste(dirpath,"/",files[j],sep = ""))
  collection_one$name<-files[j]
  collection<-rbind(collection,collection_one)
}
write.csv(collection,file = "/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/new_PRS_results_20231103.csv")


### PRS brain outcome ####


for (j in c(1:length(chronicphy_collection))){
  diseasename<-chronicphy_collection[j]
  #diseasename<-"RA"
  print(diseasename)
  if (diseasename=="CHD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/CAD_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="HTN"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/HT_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="ASM"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/ASM_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="Stroke"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/ISS_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="RA"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/RA_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="T2D"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/T2D_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="PUD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_load_compare/PRS_UKB/UC_PRS.csv")
    PRS<-PRS[,c(2,3)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }else if (diseasename=="CLD"){
    PRS<-fread("/Volumes/ZHAOYUJIE/COMMO/lungfunction/CLD.all_score")
    names(PRS)[1]<-"eid"
    PRS<-PRS[,c(1,4)]
    names(PRS)[2]<-"PRS"
    
    PRS<-na.omit(PRS)
    #10 quan
    Quan<-quantile(PRS$PRS,seq(0.1,1,0.1))
    
    PRS$class<-NA
    PRS$class[which(PRS$PRS<Quan[2])]<-1
    PRS$class[which(PRS$PRS>Quan[9])]<-3
    PRS$class[which(PRS$PRS>=Quan[2]&PRS$PRS<=Quan[9])]<-2
    
  }
  
  
  
  for(i in c(1:15)){
    braindisorder<-brain_disorder_collection[i]
    print(braindisorder)
    
    current_BD<-fread(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231027/sur_data/Sur_",brain_disorder_collection[i],"_20231010.csv",sep = ""))
    names(current_BD)[5]<-"BD_status"
    names(current_BD)[6]<-"BD_dia_days"
    names(current_BD)[4]<-"BD_group"
    
    
    brain_sur_PRS<-merge(current_BD,PRS,by ="eid")
    brain_sur_PRS<-merge(brain_sur_PRS,cov_sem,by ="eid")
    brain_sur_PRS<-na.omit(brain_sur_PRS)
    res.cox <- coxph(Surv(BD_dia_days,BD_status) ~ factor(class, level = c(1,2,3)) + Sex + Townsend_index + BMI  + Age +Qualifications, data = brain_sur_PRS)
    #cox.zph(res.cox)
    X<-summary(res.cox)
    X
    
    braind_cox = list()
    t=1
    
    braind_cox$variable[t]<-"low_PRS"
    braind_cox$outcome[t]<-braindisorder
    
    braind_cox$HR[t]<-1
    braind_cox$lowerCI[t]<-" "
    braind_cox$upperCI[t]<-" "
    #braind_cox$num[t]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==1)]
    braind_cox$num2[t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[t]<-" "
    braind_cox$pvalue[t]<-" "
    
    
    zz=1
    z=1
    braind_cox$variable[z+t]<-"intermediate_PRS"
    braind_cox$outcome[z+t]<-braindisorder
    braind_cox$HR[z+t]<-X$conf.int[zz,1]
    braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
    braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
    #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==2)]
    braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
    braind_cox$pvalue[z+t]<-X$coefficients[zz,5]
    
    
    zz = zz+1
    z = z+1    
    braind_cox$variable[z+t]<-"High_PRS"
    braind_cox$outcome[z+t]<-braindisorder
    braind_cox$HR[z+t]<-X$conf.int[zz,1]
    braind_cox$lowerCI[z+t]<-X$conf.int[zz,3]
    braind_cox$upperCI[z+t]<-X$conf.int[zz,4]
    #braind_cox$num[t+z]<-paste(X$nevent,"/",X$n,sep = "",collapse = "")
    ss<-brain_sur_PRS$comorbid[which(brain_sur_PRS$class==3)]
    braind_cox$num2[z+t]<-paste(length(ss[which(ss==1)]),"/",length(ss),sep = "",collapse = "")
    braind_cox$hrci[z+t]<-paste(round(X$conf.int[zz,1],2),"(",round(X$conf.int[zz,3],2),",",round(X$conf.int[zz,4],2),")",sep = "",collapse = "")
    braind_cox$pvalue[z+t]<-X$coefficients[zz,5]
    
    
    braind_cox<-as.data.frame(braind_cox)
    filename<-paste(paste("/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_predict_brain/",diseasename,"_",braindisorder,"_cox_PRS_20231010.csv",sep = ""))
    write.csv(braind_cox,file = filename,row.names = FALSE)
    
  }
  
}

dirpath<-"/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_cmprsk_20231019"
files<-list.files(dirpath)
collection = list()
for(j in c(1:length(files))){
  collection_one<-read.csv(paste(dirpath,"/",files[j],sep = ""))
  collection_one$name<-files[j]
  collection<-rbind(collection,collection_one)
}
write.csv(collection,file = "/Volumes/ZHAOYUJIE/COMMO/PRS_20231010/PRS_cmprsk_results/BD_new_PRS_20231019.csv")
