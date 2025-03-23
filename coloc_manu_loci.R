
library(coloc)
library(dplyr)
library(data.table)
rm(list = ls())
###calculate coloc loci in each pleiotropy loci 
combined_loci_fuma<-fread("/Volumes/yjzhao/COMMO_new/loci_fuma/combined_conjfdr_snp_FUMA.csv")
unique(combined_loci_fuma$Somatic_condition)
brain_d<-unique(combined_loci_fuma$Brain_disorder)
gwas_filepath<-"/Volumes/yjzhao/COMMO_new/new_gwas_qc/"
for(D in c(1:length(brain_d))){
  Brain_current<-brain_d[D]
  print(Brain_current)
  brain_gwas<-fread(paste0(gwas_filepath,Brain_current,"_gwas.txt"))
  brain_gwas$varbeta<-brain_gwas$SE^2
  loci_current<-subset(combined_loci_fuma,Brain_disorder==Brain_current)
  somatic_current<-unique(loci_current$Somatic_condition)
  print(somatic_current)
  for ( j in c(1:length(somatic_current))){
    Somatic_current<-somatic_current[j]
    print(Somatic_current)
    loci_current_current<-subset(loci_current,Somatic_condition==Somatic_current)
    loci_current_current$V1<-c(1:nrow(loci_current_current))
    somatic_gwas<-fread(paste0(gwas_filepath,Somatic_current,"_gwas.txt"))
    somatic_gwas$varbeta<-somatic_gwas$SE^2
    #merge two diseases
    input <- merge(brain_gwas, somatic_gwas, by="SNP", all=FALSE, suffixes=c("_brain","_somatic"))%>%
      filter(P_brain>0,P_somatic>0)
    input<-input[!duplicated(input$SNP),]
    
    Results_pleio=list()
    
    somatic_pre<-fread(paste0("/Users/yujiezhao/Desktop/COMMO_new/Disease_selection/pheno_15duseases_british_boltlmm/",Somatic_current,"_pheno.txt"))
    genotypeid<-fread("/Users/yujiezhao/Desktop/COMMO_new/BOLT-LMM/ukb_gpPopulation_500K_variants_chr1.fam")
    somatic_pre<-subset(somatic_pre,FID%in%genotypeid$V2)
    print(paste0(sum(somatic_pre$pheno),"/",nrow(somatic_pre)))
    
    for(i in c(1:length(loci_current_current$GenomicLocus))){
      
      inputcurrent<-subset(input,input$BP_brain<=loci_current_current$end[i]&input$BP_brain>=loci_current_current$start[i])
      inputcurrent<-subset(inputcurrent,inputcurrent$BP_somatic<=loci_current_current$end[i]&inputcurrent$BP_somatic>=loci_current_current$start[i])
      
      ####You cannot simply execute the loop; it needs to be done step by step because the ratios are different.
      result <- coloc.abf(dataset1=list(pvalues=inputcurrent$P_brain, type="cc",s =53386/130644, N=130644,snp = inputcurrent$SNP), 
                          dataset2=list(pvalues=inputcurrent$P_somatic,  type="cc",s = 110345/337026,N=337026,snp=inputcurrent$SNP),
                          MAF=inputcurrent$MAF)
      
      Results_pleio$V1[i]<-i
      Results_pleio$nsnps[i]<-result$summary[1]
      Results_pleio$PP.H0.abf[i]<-result$summary[2]
      Results_pleio$PP.H1.abf[i]<-result$summary[3]
      Results_pleio$PP.H2.abf[i]<-result$summary[4]
      Results_pleio$PP.H3.abf[i]<-result$summary[5]
      Results_pleio$PP.H4.abf[i]<-result$summary[6]
      
      Results_pleio$Best_causal_variant[i]<- result$results$snp[which(result$results$SNP.PP.H4==max(result$results$SNP.PP.H4))]
      Results_pleio$SNP.PP.H4[i]<- max(result$results$SNP.PP.H4)
      rm(result,inputcurrent)
      
    }
    Results_pleio<-as.data.frame(Results_pleio)
    Results_pleio$brain_d<-Brain_current
    Results_pleio$somatic_c<-Somatic_current
    
    write.csv(Results_pleio,file = paste0("/Volumes/yjzhao/COMMO_new/coloc/results/",Brain_current,"_",Somatic_current,"_coloc_loci.csv"))
    pleioloci_anno<-merge(loci_current_current,Results_pleio,by = "V1")
    write.csv(pleioloci_anno,file = paste0("/Volumes/yjzhao/COMMO_new/coloc/results/",Brain_current,"_",Somatic_current,"_coloc_loci_anno.csv"))
    
  }
}

##combined data
rm(list = ls())
setwd("/Volumes/yjzhao/COMMO_new/coloc/")
folder_path<-"/Volumes/yjzhao/COMMO_new/coloc/results/"
file_suffix <- "_coloc_loci_anno.csv" 
file_list <- list.files(folder_path, pattern = file_suffix, full.names = TRUE)
combined_data <- data.frame()
for (file in file_list) {
  data <- read.csv(file, header = TRUE)  
  combined_data <- rbind(combined_data, data)
}
write.csv(combined_data, "combined_coloc_loci_anno.csv", row.names = FALSE)


###filter pph4>0.8
rm(list = ls())
pleiocoloc_loci<-fread("/Volumes/yjzhao/COMMO_new/coloc/combined_coloc_loci_anno.csv")
pleiocoloc_loci<-pleiocoloc_loci%>%
  filter(PP.H4.abf>0.8)%>%
  arrange()

pleiocoloc_loci$bestSNP_discriminant = 0
pleiocoloc_loci$bestSNP_discriminant[which(pleiocoloc_loci$SNP==pleiocoloc_loci$Best_causal_variant)]<-1

write.csv(pleiocoloc_loci,file = "/Volumes/yjzhao/COMMO_new/coloc/combined_coloc_loci_anno_filter.csv")


coloc_loci<-fread("/Volumes/yjzhao/COMMO_new/coloc/combined_coloc_loci_anno_filter.csv")
coloc_loci_best<-subset(coloc_loci,bestSNP_discriminant==1)
coloc_loci_best<-as.data.frame(coloc_loci_best)
coloc_loci_best<-coloc_loci_best[,c(-1,-2)]
coloc_loci_best<-coloc_loci_best%>%
  mutate(disease_pair = paste0(brain_d,"_",somatic_c),loci_boundary = paste0(start,"-",end))%>%
  select(uniqID,SNP,MAF,loci_boundary,conjfdr,NearestGene,function_annotation,CADD,PP.H4.abf,SNP.PP.H4,disease_pair)
write.csv(coloc_loci_best,file = "/Volumes/yjzhao/COMMO_new/coloc/Best_causal_variants.csv",row.names = F)

frequency<-table(coloc_loci$SNP)
mulisnp<-names(frequency[frequency>1])

coloc_loci_multi<-subset(coloc_loci,SNP%in%mulisnp)

write.csv(coloc_loci_multi,file = "/Volumes/yjzhao/COMMO_new/coloc/multids_coloc_loci_anno.csv",row.names = F)
