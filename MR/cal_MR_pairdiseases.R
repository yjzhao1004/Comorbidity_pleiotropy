rm(list = ls())
library(TwoSampleMR)
library(data.table)
library(RadialMR)
library(ClumpMR)


Somaticcondition_name<-"T2D"
Braindisorder_name<-"ALD"
Somaticcondition<-T2D
Braindisorder<-ALD


Somaticcondition$phenotype <- 'exposure' 
Braindisorder$phenotype <- 'outcome' 

Somaticcondition<-Somaticcondition[Somaticcondition$pval_nominal<5e-8,]
Somaticcondition <- format_data(Somaticcondition,
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
Somaticcondition0 <- clump_data(Somaticcondition,clump_r2=0.01,clump_kb=1000)


MR_function <- function (exposure.data,outcome.data){
  final.data = harmonise_data(exposure.data, outcome.data,action = 3) 
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

results<-MR_function(Somaticcondition0,Braindisorder)
results
final = harmonise_data(Somaticcondition0, Braindisorder,action = 3) 

radial_data<-format_radial(final$beta.exposure,final$beta.outcome,
                           final$se.exposure,final$se.outcome,
                           final$SNP)

ivw.model<-ivw_radial(radial_data,0.05,3,0.0001)
ivw.model$outliers

Somaticcondition1<-setdiff(Somaticcondition0$SNP,ivw.model$outliers$SNP)
Somaticcondition11<-Somaticcondition[which(is.element(Somaticcondition$SNP,Somaticcondition1)),]

final.data = harmonise_data(Somaticcondition11, Braindisorder,action = 3)

results0<-MR_function(Somaticcondition11,Braindisorder)

results0

results0$id.exposure<-Somaticcondition_name
results0$id.outcome<-Braindisorder_name

results$id.exposure<-Somaticcondition_name
results$id.outcome<-Braindisorder_name
setwd("/Users/yujiezhao/Desktop/COMMO/MR/")
write.csv(results0,file = paste(Somaticcondition_name,"_",Braindisorder_name,"_mrresults.csv",sep = ""))

save(results0,Somaticcondition11,final.data,file = paste(Somaticcondition_name,"_",Braindisorder_name,"_mrresults.RData",sep = ""))

save(Somaticcondition0,results,final,file = paste(Somaticcondition_name,"_",Braindisorder_name,"_mrresults.RData",sep = ""))

rm(Somaticcondition,final.data,Somaticcondition11,heterogeneity, pleiotropy,Somaticcondition1,Somaticcondition0,final,ivw.model,radial_data,results,results0)
rm(Braindisorder,final.data,Somaticcondition11,heterogeneity, pleiotropy,Somaticcondition1,final,ivw.model,radial_data,results,results0)



#### reversed causality ####
rm(list = ls())
library(TwoSampleMR)
library(data.table)

#install.packages("remotes")
#remotes::install_github("WSpiller/RadialMR")
library(RadialMR)

Somaticcondition_name<-"T2D"
Braindisorder_name<-"ALD"
Somaticcondition<-T2D
Braindisorder<-ALD

Somaticcondition$phenotype <- 'outcome' 
Braindisorder$phenotype <- 'exposure' 

Braindisorder<-Braindisorder[Braindisorder$P<5e-8,]

Somaticcondition <- format_data(Somaticcondition,
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

results<-MR_function(Braindisorder0,Somaticcondition)


final = harmonise_data(Braindisorder0, Somaticcondition,action = 3) 
radial_data<-format_radial(final$beta.exposure,final$beta.outcome,
                           final$se.exposure,final$se.outcome,
                           final$SNP)

ivw.model<-ivw_radial(radial_data,0.05,3,0.0001)
ivw.model$outliers

Braindisorder1<-setdiff(Braindisorder0$SNP,ivw.model$outliers$SNP)
Braindisorder11<-Braindisorder[which(is.element(Braindisorder$SNP,Braindisorder1)),]

final.data = harmonise_data(Braindisorder11, Somaticcondition,action = 3)

results0<-MR_function(Braindisorder11, Somaticcondition)

write.csv(results0,file = paste(Braindisorder_name,"_",Somaticcondition_name,"_mrresults.csv",sep = ""))



