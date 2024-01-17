rm(list = ls())
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(ClumpMR)

#ALS
setwd("/Volumes/ZHAOYUJIE/COMMO/ALS/summary_statistics/")
ALS_filelist <- list.files()
ALS_all = list()
for(i in c(1:length(ALS_filelist))){
  print(i)
  ALS_current<-fread(ALS_filelist[i])
  ALS_all<-rbind(ALS_all,ALS_current)
}
ALS_all<-as.data.frame(ALS_all)
write.csv(ALS_all,file = "ALS_all.csv")
ALS_all<-fread("/Volumes/ZHAOYUJIE/COMMO/ALS/ALS_all.csv")
ALS_all<-ALS_all[,c(2:16)]
names(ALS_all)[1]<-"CHR"
names(ALS_all)[2]<-"SNP"
names(ALS_all)[3]<-"BP"
names(ALS_all)[4]<-"A1"
names(ALS_all)[5]<-"A2"
names(ALS_all)[7]<-"beta"
names(ALS_all)[8]<-"SE"
names(ALS_all)[9]<-"pval_nominal"
ALS<-ALS_all[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal")]
rm(ALS_all,ALS_current)

#chronicpain
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/chronic_pain-bgen.stats",header=T)
names(Depression)[5]<-'A1'
names(Depression)[6]<-'A2'
names(Depression)[7]<-'MAF'
names(Depression)[11]<-"beta"
names(Depression)[16]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
chronicpain<-Depression_qc3
rm(Depression_qc3)

#Asthma
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/Asthma/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv",header=T)
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"SNP"
names(Depression)[3]<-"BP"
names(Depression)[4]<-'A1'
names(Depression)[5]<-'A2'
names(Depression)[15]<-'beta'
names(Depression)[16]<-"SE"
names(Depression)[17]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
Asthma<-Depression_qc3
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
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
T2D<-Depression_qc3
T2D$pval_nominal<-as.numeric(T2D$pval_nominal)
rm(Depression_qc3)


#AD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/AD/AD_GCST90027158_buildGRCh38.tsv",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[4]<-"BP"
names(Depression)[3]<-"CHR"
names(Depression)[5]<-'A1'
names(Depression)[6]<-'A2'
names(Depression)[7]<-'MAF'
names(Depression)[11]<-"beta"
names(Depression)[12]<-"SE"
names(Depression)[2]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
AD<-Depression_qc3
AD$pval_nominal<-as.numeric(AD$pval_nominal)
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
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
CHD<-Depression_qc3
rm(Depression_qc3)


##lungfunction
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/lungfunction/Shrine_30804560_FEV1_to_FVC_RATIO_meta-analysis.txt",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"CHR"
names(Depression)[3]<-"BP"
names(Depression)[4]<-'A1'
names(Depression)[5]<-'A2'
names(Depression)[7]<-'MAF'
names(Depression)[10]<-"pval_nominal"
names(Depression)[11]<-"sdd"
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
lungfunction<-Depression_qc3
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
Depression_qc3<-Depression[,c("SNP","A1","A2","MAF","beta","SE","pval_nominal")]
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
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2)
arthritis<-Depression_qc3
rm(Depression_qc3)



#Stroke
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/Stroke/MEGASTROKE.1.AS.EUR.out",header=T)
names(Depression)[1]<-"SNP"
names(Depression)[2]<-"A1"
names(Depression)[3]<-"A2"
names(Depression)[4]<-"MAF"
names(Depression)[5]<-"beta"
names(Depression)[6]<-"SE"
names(Depression)[7]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
Stroke<-Depression_qc3
rm(Depression_qc3)

###MS
rm(MS)
Depression<-fread('/Volumes/ZHAOYUJIE/COMMO/MS/MS2019lnOR.csv',header=T)
names(Depression)[3]<-"beta"
names(Depression)[4]<-"SE"
names(Depression)[9]<-"pval_nominal"
names(Depression)[1]<-"SNP"
names(Depression)[5]<-"A1"
names(Depression)[6]<-"A2"
names(Depression)[7]<-"MAF"
Depression_qc3<-Depression[,c("SNP","A1","A2","MAF","beta","SE","pval_nominal")]
rm(Depression,Depression_qc2,Depression2,Depression3)
MS<-Depression_qc3
rm(Depression_qc3)

#ANX
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ANX/Anxietydisorder.meta.full.cc.txt",header=T)
names(Depression)[7]<-"beta"
names(Depression)[9]<-"pval_nominal"
names(Depression)[1]<-"SNP"
names(Depression)[8]<-"SE"
names(Depression)[4]<-"A1"
names(Depression)[5]<-"A2"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
ANX<-Depression_qc3
rm(Depression_qc3)

#OCD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/OCD/OCD_aug2017.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression)
OCD<-Depression_qc3
rm(Depression_qc3)

#ASD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ASD/ASD_iPSYCH-PGC_ASD_Nov2017.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression)
ASD<-Depression_qc3
rm(Depression_qc3)


#SCZ
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/SCZ/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv",header=T)
names(Depression)[9]<-"beta"
names(Depression)[11]<-"pval_nominal"
names(Depression)[1]<-"CHR"
names(Depression)[3]<-"BP"
names(Depression)[2]<-"SNP"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
SCZ<-Depression_qc3
rm(Depression_qc3)


#BP
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/BP/BP_pgc-bip2021-all.txt",header=T)
names(Depression)[6]<-"beta"
names(Depression)[8]<-"pval_nominal"
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
BP<-Depression_qc3
rm(Depression_qc3)


#ADHD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/ADHD/ADHD_daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression,Depression_qc2)
ADHD<-Depression_qc3
rm(Depression_qc3)


#MDD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/MDD/MDD_daner_pgc_mdd_meta_w2_no23andMe_rmUKBB_N.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","BP","A1","A2","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
MDD<-Depression_qc3
rm(Depression_qc3)


#PTSD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PTSD/pts_eur_freeze2_overall.results",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[11]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","BP","beta","A1","A2","CHR","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PTSD<-Depression_qc3
rm(Depression_qc3)


#Eat
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/EAT/pgcAN2.2019-07_EAT.txt",header=T)
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
names(Depression)[4]<-"A1"
names(Depression)[5]<-"A2"
names(Depression)[6]<-"beta"
names(Depression)[8]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","BP","beta","A1","A2","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression)
EAT<-Depression_qc3
rm(Depression_qc3)


#BP
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/BP/BP_pgc-bip2021-all.txt",header=T)
names(Depression)[6]<-"beta"
names(Depression)[8]<-"pval_nominal"
names(Depression)[1]<-"CHR"
names(Depression)[2]<-"BP"
names(Depression)[3]<-"SNP"
names(Depression)[10]<-"MAF"
Depression_qc3<-Depression[,c("SNP","BP","MAF","A1","A2","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
BP<-Depression_qc3
rm(Depression_qc3)

#Alcohol
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/Alcohol_dependency_pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","A1","A2","BP","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression)
Alcohol<-Depression_qc3
rm(Depression_qc3)


#TS
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/TS/Tourette_syndrome_Oct2018.txt",header=T)
Depression$beta<-log(Depression$OR)
names(Depression)[9]<-"pval_nominal"
Depression_qc3<-Depression[,c("SNP","CHR","A1","A2","BP","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
TS<-Depression_qc3
rm(Depression_qc3)

#HTN
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/HTN/BP-ICE_EUR_HTN_15-04-2020.txt",header=T)
names(Depression)[11]<-"pval_nominal"
names(Depression)[1]<-"SNP"
names(Depression)[3]<-"A1"
names(Depression)[4]<-"A2"
Depression_qc3<-Depression[,c("SNP","A1","A2","MaxFreq","Zscore","pval_nominal","N")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
HTN<-Depression_qc3
HTN$beta<-HTN$Zscore/sqrt(2*HTN$MaxFreq*(1-HTN$MaxFreq)*(HTN$N+HTN$Zscore^2))
HTN$SE<-1/sqrt(2*HTN$MaxFreq*(1-HTN$MaxFreq)*(HTN$N+HTN$Zscore^2))
rm(Depression_qc3)

#PD
Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PD/PD_2019Nalls_includeUKB/PD_2019Nalls.txt",header=T)
names(Depression)[3]<-"beta"
Depression_qc3<-Depression[,c("SNP","A1","A2","beta","SE","P")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PD<-Depression_qc3
rm(Depression_qc3)
names(PD)[6]<-"pval_nominal"

Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PD/PD_meta_analysis_PDGENE_PDWBS_top10k_exUKB/PD_meta_analysis_PDGENE_PDWBS.txt",header=T)
Depression$beta= log(Depression$OR.META)
names(Depression)[12]<-"pval_nominal"
names(Depression)[14]<-"SE"
Depression_qc3<-Depression[,c("SNP","CHR","BP","A1","A2","beta","SE","pval_nominal")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PD<-Depression_qc3
rm(Depression_qc3)

#ALS

Depression<-fread("/Volumes/ZHAOYUJIE/COMMO/PD/PD_2019Nalls_includeUKB/PD_2019Nalls.txt",header=T)
names(Depression)[3]<-"beta"
Depression_qc3<-Depression[,c("SNP","A1","A2","beta","SE","P")]
rm(Depression2,Depression3,Depression_qc1,Depression,Depression_qc2)
PD<-Depression_qc3
rm(Depression_qc3)

#pathway<-"/Volumes/ZHAOYUJIE/COMMO/"

#filename_cc<-"DIAMANTE-EUR.sumstat.txt"

#filename_bd<-"ADHD_daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"
  
#Chroniccondition<- fread(paste(pathway,"/T2D/",filename_cc,sep = ""),,header=T)
                         
#Braindisorder<-fread(paste(pathway,"/ADHD/",filename_bd,sep = ""),header=T)


Chroniccondition_name<-"TS"
Braindisorder_name<-"CHD"
Chroniccondition<-TS
Braindisorder<-CHD


Chroniccondition$phenotype <- 'exposure' 
Braindisorder$phenotype <- 'outcome' 

Chroniccondition<-Chroniccondition[Chroniccondition$pval_nominal<5e-6,]
Chroniccondition <- format_data(Chroniccondition,
                             type ='exposure',
                             snp_col = "SNP",
                             #z_col = "Zscore",
                             beta_col = "beta",
                             se_col = "SE",
                             effect_allele_col = "A1",
                             other_allele_col = "A2",
                             pval_col = "pval_nominal",
                             chr_col = "CHR",
                             phenotype_col = "phenotype")

Braindisorder<-format_data(Braindisorder,
                           type='outcome',
                           snp_col = "SNP",
                           beta_col = "beta",
                           se_col = "SE",
                           effect_allele_col = "A1",
                           other_allele_col = "A2",
                           pval_col = "pval_nominal",
                           chr_col = "CHR",
                           phenotype_col = "phenotype")
#library(remotes)
#remotes::install_github("huangyebao/ClumpMR")
library(ClumpMR)
rm(Depression_qc3)
Chroniccondition0 <- clump_data(Chroniccondition,clump_r2=0.01,clump_kb=1000)


MR_function <- function (exposure.data,outcome.data){
  final.data = harmonise_data(exposure.data, outcome.data,action = 3) #harmonise就是一个对齐和翻转正负链的过程,action=3, We droped all the palindromic SNPs so the results were more conservative
  final.data = final.data[which(!duplicated(final.data)),]
  #fwrite(final.data,paste0(levels(factor(final.data$exposure.data)),".csv"))
  results = mr(final.data,method_list= c("mr_ivw_mre","mr_weighted_median","mr_weighted_mode","mr_ivw_fe","mr_simple_median"))
  results$SNPs = dim(final.data)[1]
  results$OR = exp(results$b)
  results$LCI = exp(results$b - 1.96*results$se)
  results$UCI = exp(results$b + 1.96*results$se)
  #result.more <- generate_odds_ratios(results) #funnel plot
  heterogeneity = mr_heterogeneity(final.data) #heterogeneity test
  pleiotropy = mr_pleiotropy_test(final.data)
  results$IVW.Qpval <- heterogeneity$Q_pval[2]
  results$IVW.Q <- heterogeneity$Q[2]
  results$pleiotropy.pval<- pleiotropy$pval#pleiotropy test 1
  results$pleiotropy.intercept <- pleiotropy$egger_intercept
  results$pleiotropy.se <- pleiotropy$se
  return(results)
  #return(final.data)
  #return(result.more)
}

results<-MR_function(Chroniccondition0,Braindisorder)
results
final = harmonise_data(Chroniccondition0, Braindisorder,action = 3) 

radial_data<-format_radial(final$beta.exposure,final$beta.outcome,
                           final$se.exposure,final$se.outcome,
                           final$SNP)

ivw.model<-ivw_radial(radial_data,0.05,3,0.0001)
ivw.model$outliers

Chroniccondition1<-setdiff(Chroniccondition0$SNP,ivw.model$outliers$SNP)
Chroniccondition11<-Chroniccondition[which(is.element(Chroniccondition$SNP,Chroniccondition1)),]

final.data = harmonise_data(Chroniccondition11, Braindisorder,action = 3)

results0<-MR_function(Chroniccondition11,Braindisorder)

results0

results0$id.exposure<-Chroniccondition_name
results0$id.outcome<-Braindisorder_name

results$id.exposure<-Chroniccondition_name
results$id.outcome<-Braindisorder_name
setwd("/Users/yujiezhao/Desktop/COMMO/MR/")
write.csv(results0,file = paste(Chroniccondition_name,"_",Braindisorder_name,"_mrresults.csv",sep = ""))

save(results0,Chroniccondition11,final.data,file = paste(Chroniccondition_name,"_",Braindisorder_name,"_mrresults.RData",sep = ""))

save(Chroniccondition0,results,final,file = paste(Chroniccondition_name,"_",Braindisorder_name,"_mrresults.RData",sep = ""))

rm(Chroniccondition,final.data,Chroniccondition11,heterogeneity, pleiotropy,Chroniccondition1,Chroniccondition0,final,ivw.model,radial_data,results,results0)
rm(Braindisorder,final.data,Chroniccondition11,heterogeneity, pleiotropy,Chroniccondition1,final,ivw.model,radial_data,results,results0)



#reverse
rm(list = ls())
library(TwoSampleMR)
library(data.table)

#install.packages("remotes")
#remotes::install_github("WSpiller/RadialMR")
library(RadialMR)

Chroniccondition<- fread('GWAS_chr_full_glm_linear_Chroniccondition20230408.txt',header=T)
Braindisorder<-fread('Braindisorder_GWAS_daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.txt',header=T)
Braindisorder$BETA <- log(Braindisorder$OR) #

#reverse causality
Chroniccondition$phenotype <- 'outcome' 
Braindisorder$phenotype <- 'exposure' 

Braindisorder<-Braindisorder[Braindisorder$P<1e-7,]

Chroniccondition <- format_data(Chroniccondition,
                         type ='outcome',
                         snp_col = "ID",
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "ALT",
                         other_allele_col = "REF",
                         pval_col = "P",
                         chr_col = "CHROM",
                         phenotype_col = "phenotype")

Braindisorder<-format_data(Braindisorder,type='exposure',
                        snp_col = "SNP",beta_col = "BETA",
                        se_col = "SE",effect_allele_col = "A1",
                        other_allele_col = "A2",
                        pval_col = "P",
                        chr_col = "CHR",
                        phenotype_col = "phenotype")


Braindisorder0 <- clump_data(Braindisorder,clump_r2=0.01,clump_kb=1000)

results<-MR_function(Braindisorder0,Chroniccondition)


final = harmonise_data(Braindisorder0, Chroniccondition,action = 3) 
radial_data<-format_radial(final$beta.exposure,final$beta.outcome,
                           final$se.exposure,final$se.outcome,
                           final$SNP)

ivw.model<-ivw_radial(radial_data,0.05,3,0.0001)
ivw.model$outliers

Braindisorder1<-setdiff(Braindisorder0$SNP,ivw.model$outliers$SNP)
Braindisorder11<-Braindisorder[which(is.element(Braindisorder$SNP,Braindisorder1)),]

final.data = harmonise_data(Braindisorder11, Chroniccondition,action = 3)

results0<-MR_function(Braindisorder11, Chroniccondition)


rm(list = ls())
library(ggplot2)
res_single <- mr_singlesnp(final.data,all_method = c("mr_simple_median","mr_weighted_mode","mr_weighted_median","mr_ivw_fe"))
pp<-mr_forest_plot(res_single)

pp$NksTEY.ZoyTHV + xlab('MR effect size for Chroniccondition on Braindisorder,β')  +
  theme_classic(base_size = 17)


#mr_funnel_plot(res_single)


results0<-read.csv("/Users/yujiezhao/Desktop/Chroniccondition_Braindisorder/GWAS/MR20230408/Inverse_MR_20230414.csv")
final.data<-read.csv("/Users/yujiezhao/Desktop/Chroniccondition_Braindisorder/GWAS/MR20230408/Harmonised_20230414.csv")
save(results0,final.data,file = "/Users/yujiezhao/Desktop/Chroniccondition_Braindisorder/GWAS/MR20230408/mr_results_inversed_20230414.RData")

library(ggplot2)
p<-mr_scatter_plot(results0[c(2:5),],final.data)
p$nn4Uf6.S2Frtf
p$nn4Uf6.S2Frtf +
  #scale_x_continuous(limits= c(0,0.075), breaks= seq(0, 0.075, 0.025)) +
  xlab('SNP effect on Braindisorder') + ylab('SNP effect on Chroniccondition ') +
  theme_classic(base_size = 17)

