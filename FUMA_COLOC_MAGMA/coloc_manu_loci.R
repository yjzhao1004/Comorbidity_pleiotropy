if(!require("remotes"))
  install.packages("remotes")
library(remotes)
install_github("chr1swallace/coloc",build_vignettes=TRUE)

library("coloc")
library(dplyr)
library(data.table)

###calculate coloc loci in each pleiotropy loci 

#input T2D and ALD gwas data
#T2D
T2D<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/DIAMANTE-EUR.sumstat.txt",header=T)
names(T2D)[1]<-"CHR"
names(T2D)[2]<-"BP"
names(T2D)[4]<-"SNP"
names(T2D)[5]<-'A1'
names(T2D)[6]<-'A2'
names(T2D)[7]<-'MAF'
names(T2D)[8]<-"beta"
names(T2D)[9]<-"SE"
names(T2D)[10]<-"pval_nominal"
T2D$varbeta<-T2D$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
T2D2<-T2D[T2D$CHR==6,]
T2D3<-T2D2[T2D2$BP>25000000,]
T2D3<-T2D3[T2D3$BP<35000000,]
T2D_qc1<-subset(T2D,SNP %in% setdiff(T2D$SNP,T2D3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
T2D_qc2<-subset(T2D_qc1,MAF>0.01)
rm(T2D2,T2D3,T2D_qc1)
T2D_qc3<-T2D_qc2[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(T2D,T2D_qc2)
T2D<-T2D_qc3
T2D$pval_nominal<-as.numeric(T2D$pval_nominal)
rm(T2D_qc3)

#ALD
ALD<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/Alcohol_dependency_pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt",header=T)
ALD$beta<-log(ALD$OR)
names(ALD)[9]<-"pval_nominal"
ALD$varbeta<-ALD$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
ALD2<-ALD[ALD$CHR==6,]
ALD3<-ALD2[ALD2$BP>25000000,]
ALD3<-ALD3[ALD3$BP<35000000,]
ALD_qc1<-subset(ALD,SNP %in% setdiff(ALD$SNP,ALD3$SNP))
ALD_qc3<-ALD_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(ALD2,ALD3,ALD_qc1,ALD)
ALD<-ALD_qc3
rm(ALD_qc3)

#merge two diseases
input <- merge(T2D, ALD, by="SNP", all=FALSE, suffixes=c("_T2D","_ALD"))
head(input)
input<-subset(input,MAF<1)
input<-input[!duplicated(input$SNP),]
names(input)[9]<-"pval_nominal_T2D"
names(input)[14]<-"pval_nominal_ALD"

input<-subset(input,input$pval_nominal_T2D>0)
input<-subset(input,input$pval_nominal_ALD>0)


#select placo-significant SNP for this trait-pair
placo_snp<-fread("/Users/yujiezhao/Desktop/COMMO/PLACO/all_sig_placo_bon_20230607.csv")
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

#pleiotropic loci
uniquepleioloci_ALD<-unique(snp_ALD_placo$GenomicLocus)
uniquepleioloci_T2D<-unique(snp_T2D_placo$GenomicLocus)

snp_ALD_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/ALD/FUMA_jobALD/GenomicRiskLoci.txt")
snp_T2D_loci<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/GenomicRiskLoci.txt")

snp_ALD_loci_pleio<-subset(snp_ALD_loci,GenomicLocus%in%uniquepleioloci_ALD)
snp_T2D_loci_pleio<-subset(snp_T2D_loci,GenomicLocus%in%uniquepleioloci_T2D)

T2D_ALD_pleioloci<-fread("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_ALD_pleioloci.csv")

Results_pleio=list()

a<-strsplit(T2D_ALD_pleioloci$Locus_boundary, " - ")
a<-as.matrix(a)

for(i in c(1:length(snp_ALD_loci_pleio$GenomicLocus))){

 inputcurrent<-subset(input,input$BP_T2D<as.numeric(a[[i]][2])&input$BP_T2D>as.numeric(a[[i]][1]))
 inputcurrent<-subset(inputcurrent,inputcurrent$BP_ALD<as.numeric(a[[i]][2])&inputcurrent$BP_ALD>as.numeric(a[[i]][1]))
  
  #inputcurrent<-subset(input,input$ALD<as.numeric(a[[i]][2])&input$ALD>as.numeric(a[[i]][1]))
  
  result <- coloc.abf(dataset1=list(pvalues=inputcurrent$pval_nominal_ALD, type="cc", s =8485/28757, N=28757,snp = inputcurrent$SNP), dataset2=list(pvalues=inputcurrent$pval_nominal_T2D,  type="cc",s = 80154/933970,N=933970,snp=inputcurrent$SNP), MAF=inputcurrent$MAF)
  
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
Results_pleio$trait_pair="T2D_ALD"

write.csv(Results_pleio,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/T2D_ALD_coloc_loci.csv")
T2D_ALD_pleioloci<-merge(T2D_ALD_pleioloci,Results_pleio[,c(1,6:10)],by = "V1")
write.csv(T2D_ALD_pleioloci,file = "/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/T2D_ALD_pleio_coloc_loci.csv")



#### coloc loci collection
rm(list = ls())
fileorder_disease1<-c("MCP","MCP","MCP","MCP","MCP","MCP",
                      "MCP","MCP","MCP","CAD","CAD","CAD","T2D","T2D","T2D","T2D","T2D",
                      "T2D","T2D","T2D","T2D","RA","RA","CLD","PUD","PUD","PUD","Stroke","Stroke")

fileorder_disease2<-c("ADHD","ANX","ASD","BP","MDD","MS","PTSD","SCZ","TS","ADHD","ALD","MDD","ADHD","ANX","ALD",
                      "EAT","MDD","MS","OCD","PTSD","SCZ","BP","MS","SCZ","ADHD","ALD","SCZ","ADHD","BP")


library(data.table)
j = 1
rm(qq_disease_file1,qq_disease_file)
qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/',fileorder_disease1[j],"_",fileorder_disease2[j],'_pleio_coloc_loci.csv',sep = '')
qq_disease_file1<-fread(qq_disease,header=T)

for( j in c(2:length(fileorder_disease1))){
  qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/',fileorder_disease1[j],"_",fileorder_disease2[j],'_pleio_coloc_loci.csv',sep = '')
  qq_disease_file<-fread(qq_disease,header=T)
  qq_disease_file1<-rbind(qq_disease_file1,qq_disease_file)
}

write.csv(qq_disease_file1,file = '/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/all_fuma_coloc_pleioloci20230619.csv')

need_result= qq_disease_file1[,c(2:20)] %>% filter(PP.H4.abf > 0.7)

need_result$shifou = 0
need_result$shifou[which(need_result$TopSNP==need_result$Best_causal_variant)]<-1

write.csv(need_result,file = '/Users/yujiezhao/Desktop/COMMO/FUMA/COLOC/coloc_pleioloci20230619.csv')




