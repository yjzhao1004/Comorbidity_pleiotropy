
if(!require("remotes"))
  install.packages("remotes")
install.packages("dplyr")
library(remotes)
install_github("chr1swallace/coloc",build_vignettes=TRUE)

rm(list = ls())
library("coloc")
library(dplyr)
library(data.table)

#chronicpain
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/chronic_pain-bgen.stats",header=T)
names(Depression)[5]<-'A1'
names(Depression)[6]<-'A2'
names(Depression)[7]<-'MAF'
names(Depression)[11]<-"beta"
names(Depression)[16]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,MAF>0.01)
rm(Depression2,Depression3,Depression_qc1)
Depression_qc2$N<-387649 
Depression_qc3<-Depression_qc2[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression,Depression_qc2)
chronicpain<-Depression_qc3
rm(Depression_qc3)

#T2D
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/DIAMANTE-EUR.sumstat.txt",header=T)
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[4]<-"SNP"
names(Depression)[5]<-'A1'
names(Depression)[6]<-'A2'
names(Depression)[7]<-'MAF'
names(Depression)[8]<-"beta"
names(Depression)[9]<-"SE"
names(Depression)[10]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,MAF>0.01)
rm(Depression2,Depression3,Depression_qc1)
Depression_qc3<-Depression_qc2[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression,Depression_qc2)
T2D<-Depression_qc3
T2D$pval_nominal<-as.numeric(T2D$pval_nominal)
rm(Depression_qc3)


##CHD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/CHD/second/cad.additive.Oct2015.pub/cad.add.160614.website.txt",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"CHR"
names(Depression)[3]<-"BP"
names(Depression)[4]<-'A1'
names(Depression)[5]<-'A2'
names(Depression)[6]<-'MAF'
names(Depression)[10]<-"SE"
names(Depression)[11]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,MAF>0.01)
rm(Depression2,Depression3,Depression_qc1)
Depression_qc3<-Depression_qc2[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression,Depression_qc2)
CHD<-Depression_qc3
rm(Depression_qc3)


##lungfunction
rm(list=ls())
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/lungfunction/Shrine_30804560_FEV1_to_FVC_RATIO_meta-analysis.txt",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"CHR"
names(Depression)[3]<-"BP"
names(Depression)[4]<-'A1'
names(Depression)[5]<-'A2'
names(Depression)[7]<-'MAF'

Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,MAF>0.01)
rm(Depression2,Depression3,Depression_qc1)
Depression_qc3<-Depression_qc2[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","P","varbeta")]
rm(Depression,Depression_qc2)
lungfunction<-Depression_qc3
rm(Depression_qc3)
#write.table(lungfunction,file = "/Users/yujiezhao/Desktop/COMMO/LDSC/data/QC_GWAS_new/lungfunction_gwas_qc_meta.txt",quote = FALSE,row.names = FALSE)

#Stroke
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/Stroke/MEGASTROKE.1.AS.EUR.out",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"A1"
names(Depression)[3]<-"A2"
names(Depression)[4]<-"MAF"
names(Depression)[5]<-"beta"
names(Depression)[6]<-"SE"
names(Depression)[7]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression,MAF>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
Stroke<-Depression_qc3
rm(Depression_qc3)


#PUD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PUD/PUD_summary.txt",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"A1"
names(Depression)[3]<-"A2"
names(Depression)[4]<-"MAF"
names(Depression)[5]<-"beta"
names(Depression)[6]<-"SE"
names(Depression)[7]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression,MAF>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","A1","A2","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PUD<-Depression_qc3
rm(Depression_qc3)


#arthritis
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/arthritis/Rheumatoid_arthritis_stahl_2010_20453842_ra_efo0000685_1_gwas.sumstats.tsv",header=T)
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
names(Depression)[4]<-"A1"
names(Depression)[5]<-"A2"
names(Depression)[6]<-"pval_nominal"
names(Depression)[7]<-"beta"
names(Depression)[8]<-"SE"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal","varbeta")]
rm(Depression,Depression_qc2)
arthritis<-Depression_qc3
rm(Depression_qc3)


###MS
rm(MS)
Depression<-fread('/Volumes/ZHAOYUJIE/COMMO/MS/MS2019lnOR.csv',header=T)
names(Depression)[3]<-"beta"
names(Depression)[4]<-"SE"
Depression$varbeta<-Depression$SE^2
names(Depression)[9]<-"pval_nominal"
names(Depression)[1]<-"SNP"
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
#MS was unable to do this step due to the lost of CHR info
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression,ms.Freq>0.01)
names(Depression_qc2)[7]<-"MAF"
Depression_qc3<-Depression_qc2[,c("SNP","pval_nominal","MAF","beta","varbeta")]
rm(Depression,Depression_qc2,Depression2,Depression3)
MS<-Depression_qc3
rm(Depression_qc3)


#ADHD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ADHD/ADHD_daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,FRQ_A_19099>0.01)
Depression_qc2<-subset(Depression_qc2,FRQ_U_34194>0.01)
rm(Depression2,Depression3,Depression_qc1)
Depression_qc3<-Depression_qc2[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression,Depression_qc2)
ADHD<-Depression_qc3
rm(Depression_qc3)
  


#ASD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ASD/ASD_iPSYCH-PGC_ASD_Nov2017.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression)
ASD<-Depression_qc3
rm(Depression_qc3)


#BP
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/BP/BP_pgc-bip2021-all.txt",header=T)
names(Depression)[6]<-"beta"
names(Depression)[8]<-"pval_nominal"
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,FCAS>0.01)
Depression_qc2<-subset(Depression_qc2,FCON>0.01)
names(Depression_qc2)[10]<-"MAF"
Depression_qc3<-Depression_qc2[,c("SNP","BP","MAF","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
BP<-Depression_qc3
rm(Depression_qc3)



#MDD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/MDD/MDD_daner_pgc_mdd_meta_w2_no23andMe_rmUKBB_N.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,FRQ_A_45396>0.01)
Depression_qc2<-subset(Depression_qc2,FRQ_U_97250>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
MDD<-Depression_qc3
rm(Depression_qc3)


#PTSD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PTSD/pts_eur_freeze2_overall.results",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,FRQ_A_23212>0.01)
Depression_qc2<-subset(Depression_qc2,FRQ_U_151447>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PTSD<-Depression_qc3
rm(Depression_qc3)



#SCZ
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",header=T)
names(Depression)[9]<-"beta"
names(Depression)[11]<-"pval_nominal"
names(Depression)[1]<-"CHR"
names(Depression)[3]<-"BP"
names(Depression)[2]<-"SNP"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,FCAS>0.01)
Depression_qc2<-subset(Depression_qc2,FCON>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
SCZ<-Depression_qc3
rm(Depression_qc3)




#TS
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/TS/Tourette_syndrome_Oct2018.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
TS<-Depression_qc3
rm(Depression_qc3)



#ANX
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ANX/Anxietydisorder.meta.full.cc.txt",header=T)
names(Depression)[7]<-"beta"
names(Depression)[9]<-"pval_nominal"
names(Depression)[1]<-"SNP"
names(Depression)[8]<-"SE"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
Depression_qc2<-subset(Depression_qc1,Freq1>0.01)
Depression_qc3<-Depression_qc2[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
ANX<-Depression_qc3
rm(Depression_qc3)



#OCD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/OCD/OCD_aug2017.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression)
OCD<-Depression_qc3
rm(Depression_qc3)



#Alcohol
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/Alcohol_dependency_pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression)
Alcohol<-Depression_qc3
rm(Depression_qc3)


#Eat
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/EAT/pgcAN2.2019-07_EAT.txt",header=T)
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
names(Depression)[6]<-"beta"
names(Depression)[8]<-"pval_nominal"
Depression$varbeta<-Depression$SE^2
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
Depression2<-Depression[Depression$CHR==6,]
Depression3<-Depression2[Depression2$BP>25000000,]
Depression3<-Depression3[Depression3$BP<35000000,]
Depression_qc1<-subset(Depression,SNP %in% setdiff(Depression$SNP,Depression3$SNP))
Depression_qc3<-Depression_qc1[,c("SNP","BP","beta","SE","pval_nominal","varbeta")]
rm(Depression2,Depression3,Depression_qc1,Depression)
EAT<-Depression_qc3
rm(Depression_qc3)
