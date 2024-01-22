rm(list = ls())
library(data.table)
placo_snp<-fread("/Users/yujiezhao/Desktop/COMMO/PLACO/all_sig_placo_bon_20230607.csv")

#select placo-significant SNP for this trait-pair
placo_snp_T2D_ALD<-subset(placo_snp,name=="T2D_ALD")


#select annotated SNP from placo SNPs
snp_ALD<-fread("/Volumes/ZHAOYUJIE/COMMO/ALD/FUMA_jobALD/snps.txt")
snp_ALD_placo<-subset(snp_ALD,rsID%in%placo_snp_T2D_ALD$SNP)
snp_T2D<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/snps.txt")
snp_T2D_placo<-subset(snp_T2D,rsID%in%placo_snp_T2D_ALD$SNP)
rm(snp_ALD,snp_T2D)


#find intersected SNP of the two trait-annotated-placo SNP
inter_snp<-intersect(snp_ALD_placo$rsID,snp_T2D_placo$rsID)
snp_ALD_placo<-subset(snp_ALD_placo,rsID%in%inter_snp)
snp_T2D_placo<-subset(snp_T2D_placo,rsID%in%inter_snp)


#pleiotropic snp
pleiotropic_loci = list()
pleiotropic_loci<-snp_ALD_placo[,c("uniqID","rsID","chr","pos","non_effect_allele",
                                        "effect_allele","MAF","nearestGene","func","CADD","RDB")]

pleiotropic_loci$gwasP_ALD<-snp_ALD_placo$gwasP
pleiotropic_loci$beta_ALD<-snp_ALD_placo$beta
pleiotropic_loci$beta_ALD<-log(snp_ALD_placo$or)

pleiotropic_loci$se_ALD<-snp_ALD_placo$se

pleiotropic_loci$gwasP_T2D<-snp_T2D_placo$gwasP
pleiotropic_loci$beta_T2D<-snp_T2D_placo$beta
pleiotropic_loci$se_T2D<-snp_T2D_placo$se

pleiotropic_loci$p.placo<-placo_snp_T2D_ALD$p.placo[match(snp_ALD_placo$rsID,placo_snp_T2D_ALD$SNP)]

write.csv(pleiotropic_loci,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_ALD_pleiosnp.csv")


#pleiotropic loci
uniquepleioloci_ALD<-unique(snp_ALD_placo$GenomicLocus)
uniquepleioloci_T2D<-unique(snp_T2D_placo$GenomicLocus)


snp_ALD_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/ALD/FUMA_jobALD/GenomicRiskLoci.txt")
snp_T2D_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/GenomicRiskLoci.txt")
snp_ALD_loci_pleio<-subset(snp_ALD_loci,GenomicLocus%in%uniquepleioloci_ALD)
snp_T2D_loci_pleio<-subset(snp_T2D_loci,GenomicLocus%in%uniquepleioloci_T2D)

pleiolocus = list()
i=1
for (j in uniquepleioloci_ALD){
  current1<-subset(snp_ALD_placo,GenomicLocus==j)
  current2<-subset(pleiotropic_loci,rsID%in%current1$rsID)
  
  pleiolocus$TopSNP[i]<-pleiotropic_loci$rsID[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$CHR[i]<-pleiotropic_loci$chr[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$BP[i]<-pleiotropic_loci$pos[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$A1[i]<-pleiotropic_loci$effect_allele[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$A2[i]<-pleiotropic_loci$non_effect_allele[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$MAF[i]<-pleiotropic_loci$MAF[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$Nearest_Gene[i]<-pleiotropic_loci$nearestGene[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$Functional_annotation[i]<-pleiotropic_loci$func[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$CADD[i]<-pleiotropic_loci$CADD[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$RDB[i]<-pleiotropic_loci$RDB[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  pleiolocus$P.placo[i]<-pleiotropic_loci$p.placo[which(pleiotropic_loci$p.placo==min(current2$p.placo))]
  start<-c(snp_T2D_loci_pleio$start[i],snp_ALD_loci_pleio$start[i])
  end<-c(snp_T2D_loci_pleio$end[i],snp_ALD_loci_pleio$end[i])
  pleiolocus$Locus_boundary[i]<-paste(as.character(max(start)),"-",as.character(min(end)))
  i = i+1
 }

pleiolocus<-as.data.frame(pleiolocus)

pleiolocus$Trait_pair<-"T2D_ALD"
pleiolocus
write.csv(pleiolocus,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_ALD_pleioloci.csv")



##pleiotropy loci collection
rm(list = ls())
fileorder_disease1<-c("MCP","MCP","MCP","MCP","MCP","MCP",
                      "MCP","MCP","MCP","CHD","CHD","CHD","T2D","T2D","T2D","T2D","T2D",
                      "T2D","T2D","T2D","T2D","RA","RA","CLD","PUD","PUD","PUD","Stroke","Stroke")

fileorder_disease2<-c("ADHD","ANX","ASD","BP","MDD","MS","PTSD","SCZ","TS","ADHD","ALD","MDD","ADHD","ANX","ALD",
                      "EAT","MDD","MS","OCD","PTSD","SCZ","BP","MS","SCZ","ADHD","ALD","SCZ","ADHD","BP")

###fileorder_disease1 and fileorder_disease2 are diseases name in the disease pair with pleiotropic loci



library(data.table)
j = 1
rm(qq_disease_file1,qq_disease_file)
qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/',fileorder_disease1[j],"_",fileorder_disease2[j],'_pleioloci.csv',sep = '')
qq_disease_file1<-fread(qq_disease,header=T)

for( j in c(2:length(fileorder_disease1))){
  qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/',fileorder_disease1[j],"_",fileorder_disease2[j],'_pleioloci.csv',sep = '')
  qq_disease_file<-fread(qq_disease,header=T)
  
  qq_disease_file1<-rbind(qq_disease_file1,qq_disease_file)
  
}

write.csv(qq_disease_file1,file = '/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci20230619.csv')



#### functional annotation in placo loci ####
rm(list = ls())
pleioloci<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci20230619.csv")

trait_unique_pair<-unique(pleioloci$Trait_pair)

setwd("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/")
pleioloci$beta_somatic<-0
pleioloci$beta_brain<-0
for (i in c(1:length(trait_unique_pair))){
  
  filename<-paste(trait_unique_pair[i],"_pleiosnp.csv",sep = "")
  pleiosnp<-read.csv(filename)
  
  pleioloci$beta_somatic[which(pleioloci$Trait_pair==trait_unique_pair[i])]<-pleiosnp[match(pleioloci[which(pleioloci$Trait_pair==trait_unique_pair[i]),]$TopSNP,pleiosnp$rsID),17]
  pleioloci$beta_brain[which(pleioloci$Trait_pair==trait_unique_pair[i])]<-pleiosnp[match(pleioloci[which(pleioloci$Trait_pair==trait_unique_pair[i]),]$TopSNP,pleiosnp$rsID),14]
  rm(pleiosnp,filename)
  
}
pleioloci$beta2<-pleioloci$beta_somatic*pleioloci$beta_brain
pleioloci$direction<-0
pleioloci$direction[pleioloci$beta2<0]="-"
pleioloci$direction[pleioloci$beta2>0]="+"
write.csv(pleioloci,file = "all_fuma_pleioloci_direction_20230717.csv")

###pleiotropy loci with multiple disease pairs
multi_pleioloci<-pleioloci[duplicated(pleioloci$Locus_boundary),]
all_multipleioloci = list()
for (i in c(1:length(multi_pleioloci$X))){
  current<-subset(pleioloci,Locus_boundary==multi_pleioloci[i,14])
  all_multipleioloci<-rbind(all_multipleioloci,current)
}
write.csv(all_multipleioloci,file = "all_multi_fuma_pleioloci_direction_20230717.csv")







