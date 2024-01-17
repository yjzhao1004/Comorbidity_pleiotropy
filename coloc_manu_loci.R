#####
#calculate coloc loci in each pleiotropic loci
input <- merge(T2D, Alcohol, by="SNP", all=FALSE, suffixes=c("_T2D","_Alcohol"))
head(input)
#input<-subset(input,MAF_T2D<1)
input<-subset(input,MAF<1)
input<-input[!duplicated(input$SNP),]
names(input)[9]<-"pval_nominal_T2D"
names(input)[14]<-"pval_nominal_Alcohol"

input<-subset(input,input$pval_nominal_T2D>0)
input<-subset(input,input$pval_nominal_Alcohol>0)



placo_snp<-fread("/Users/yujiezhao/Desktop/COMMO/PLACO/all_sig_placo_bon_20230607.csv")

#select placo-significant SNP for this trait-pair
placo_snp_T2D_Alcohol<-subset(placo_snp,name=="T2D_Alcohol")

#select annotated SNP from placo SNPs
snp_Alcohol<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/FUMA_jobAlcohol/snps.txt")
snp_Alcohol_placo<-subset(snp_Alcohol,rsID%in%placo_snp_T2D_Alcohol$SNP)
snp_T2D<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/snps.txt")
snp_T2D_placo<-subset(snp_T2D,rsID%in%placo_snp_T2D_Alcohol$SNP)
rm(snp_Alcohol,snp_T2D)


#find intersected SNP of the two trait-annotated-placo SNP
inter_snp<-intersect(snp_Alcohol_placo$rsID,snp_T2D_placo$rsID)
snp_Alcohol_placo<-subset(snp_Alcohol_placo,rsID%in%inter_snp)
snp_T2D_placo<-subset(snp_T2D_placo,rsID%in%inter_snp)


#pleiotropic loci
#snp_T2D_placo$GenomicLocus[which(snp_T2D_placo$GenomicLocus==75)]<-74

uniquepleioloci_Alcohol<-unique(snp_Alcohol_placo$GenomicLocus)
uniquepleioloci_T2D<-unique(snp_T2D_placo$GenomicLocus)


snp_Alcohol_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/FUMA_jobAlcohol/GenomicRiskLoci.txt")
snp_T2D_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/GenomicRiskLoci.txt")

snp_Alcohol_loci_pleio<-subset(snp_Alcohol_loci,GenomicLocus%in%uniquepleioloci_Alcohol)
snp_T2D_loci_pleio<-subset(snp_T2D_loci,GenomicLocus%in%uniquepleioloci_T2D)

T2D_Alcohol_pleioloci<-fread("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_Alcohol_pleioloci.csv")

Results_pleio=list()


a<-strsplit(T2D_Alcohol_pleioloci$Locus_boundary, " - ")
a<-as.matrix(a)

for(i in c(1:length(snp_Alcohol_loci_pleio$GenomicLocus))){
 # inputcurrent<-subset(input,input$Alcohol_T2D<snp_T2D_loci_pleio$end[i]&input$Alcohol_T2D>snp_T2D_loci_pleio$start[i])
 # inputcurrent<-subset(inputcurrent,inputcurrent$Alcohol_Alcohol<snp_Alcohol_loci_pleio$end[i]&inputcurrent$Alcohol_Alcohol>snp_Alcohol_loci_pleio$start[i])
  
 inputcurrent<-subset(input,input$BP_T2D<as.numeric(a[[i]][2])&input$BP_T2D>as.numeric(a[[i]][1]))
 inputcurrent<-subset(inputcurrent,inputcurrent$BP_Alcohol<as.numeric(a[[i]][2])&inputcurrent$BP_Alcohol>as.numeric(a[[i]][1]))
  
  #inputcurrent<-subset(input,input$Alcohol<as.numeric(a[[i]][2])&input$Alcohol>as.numeric(a[[i]][1]))
  
  result <- coloc.abf(dataset1=list(pvalues=inputcurrent$pval_nominal_Alcohol, type="cc", s =8485/28757, N=28757,snp = inputcurrent$SNP), dataset2=list(pvalues=inputcurrent$pval_nominal_T2D,  type="cc",s = 80154/933970,N=933970,snp=inputcurrent$SNP), MAF=inputcurrent$MAF)
  
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
Results_pleio$trait_pair="T2D_Alcohol"

write.csv(Results_pleio,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/T2D_Alcohol_coloc_loci.csv")
T2D_Alcohol_pleioloci<-merge(T2D_Alcohol_pleioloci,Results_pleio[,c(1,6:10)],by = "V1")
write.csv(T2D_Alcohol_pleioloci,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/T2D_Alcohol_pleio_coloc_loci.csv")


rm(a,T2D_Alcohol_pleioloci,input,placo_snp,placo_snp_T2D_Alcohol,Results_pleio,snp_T2D_loci,snp_T2D_loci_pleio,snp_T2D_placo,snp_Alcohol_loci,snp_Alcohol_loci_pleio,snp_Alcohol_placo)

