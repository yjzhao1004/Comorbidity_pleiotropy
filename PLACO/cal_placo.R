rm(list = ls())

require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")

#Before this step we have obtained disease pairs with significant genetic correlation (sig-rg)
#fileorder_disease1: one disease in the sig-rg disease pairs
#fileorder_disease2: the other disease in the sig-rg disease pairs
fileorder_disease1<-c("MCP",	"MCP",	"T2D",	"MCP",	"PUD",
                      "PUD",	"T2D",	"PUD",	"MCP",	"CHD",	"MCP",
                      "T2D",	"PUD",	"MCP",	"T2D",	"MCP",	"CHD",
                      "MCP",	"T2D",	"CHD",	"T2D",	"stroke",	"CHD",	"CLD",
                      "MCP",	"T2D",	"T2D",	"T2D",	"CHD",	"MCP",	"PUD",
                      "stroke",	"T2D",	"RA",	"RA")
fileorder_disease2<-c("ADHD",	"MDD",	"ADHD",	"PTSD",	"ADHD",	"MDD",	"Eat",
                      "ALD",	"ALD",	"ADHD",	"TS",	"MDD",	"PTSD",
                      "BP",	"PTSD",	"SCZ",	"MDD",	"ANX",	"MS",	"AD",
                      "OCD",	"ADHD",	"PTSD",	"SCZ",	"ASD",	"ALD",
                      "ANX",	"AD",	"ALD",	"MS",	"SCZ",	"BP",	"SCZ",	"BP",	"MS")

  
library(data.table)

for (j in c(1:length(fileorder_disease2))){
  
  print(j)
  filename_disease1<-paste('/Users/yujiezhao/Desktop/COMMO/LDSC/data/sum_gwas/',fileorder_disease1[j],'.sumstats',sep = '')
  disease1<-fread(filename_disease1,header=T)
  
  disease1<-na.omit(disease1)
  disease1<-subset(disease1,Z^2<=80)
  
  filename_disease1_p<-paste('/Users/yujiezhao/Desktop/COMMO/LDSC/data/QC_GWAS_new/',fileorder_disease1[j],'_gwas_qc.txt',sep = '')
  disease1_p<-fread(filename_disease1_p,header=T)
  
  filename_disease2<-paste('/Users/yujiezhao/Desktop/COMMO/LDSC/data/sum_gwas/',fileorder_disease2[j],'.sumstats',sep = '')
  disease2<-fread(filename_disease2,header=T)
  
  disease2<-na.omit(disease2)
  disease2<-subset(disease2,Z^2<=80)
  
  filename_disease2_p<-paste('/Users/yujiezhao/Desktop/COMMO/LDSC/data/QC_GWAS_new/',fileorder_disease2[j],'_gwas_qc.txt',sep = '')
  disease2_p<-fread(filename_disease2_p,header=T)
  
  names(disease1)[4]<-'Z_disease1'
  names(disease2)[4]<-'Z_disease2'
  
  Z.matrix<-merge(disease1[,c(1,4)],disease2[,c(1,4)],by = "SNP")
  
  
  names(disease1_p)[which(colnames(disease1_p) == 'P')]<-'p_value_disease1'
  names(disease2_p)[which(colnames(disease2_p) == 'P')]<-'p_value_disease2'
  
  
  P.matrix<-merge(disease1_p[,c("SNP","p_value_disease1")],disease2_p[,c("SNP","p_value_disease2")],by = "SNP")
  
  
  zp_matrix<-merge(Z.matrix,P.matrix,by = 'SNP')
  

  var<-var.placo(as.matrix(zp_matrix[,c(2,3)]), as.matrix(zp_matrix[,c(4,5)]), p.threshold=1e-4)
  
  
  placolist = list()
  
  
  for (i in c(1:dim(zp_matrix)[1])){
    print(i)
    a<-placo(as.matrix(zp_matrix[i,c(2,3)]),var)
    placolist$SNP[i]<-zp_matrix$SNP[i]
    placolist$T.placo[i]<-a$T.placo
    placolist$p.placo[i]<-a$p.placo
  }
  
  
  placolist<-as.data.frame(placolist)
  placolist_sig<-subset(placolist,placolist$p.placo<5e-8)
  
  filesigplaconame<-paste( "/Users/yujiezhao/Desktop/COMMO/PLACO/sig_placoresults/",fileorder_disease1[j],"_",fileorder_disease2[j],".txt",sep = '')
  write.table(placolist_sig,file = filesigplaconame,row.names = F,col.names = T,sep = " ",quote=F)
  
  fileplaconame<-paste( "/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults/",fileorder_disease1[j],"_",fileorder_disease2[j],".txt",sep = '')
  write.table(placolist,file = fileplaconame,row.names = F,col.names = T,sep = " ",quote=F)
  
  rm(placolist,placolist_sig,var,a,Z.matrix,P.matrix,zp_matrix)
}



  # Here Z.matrix is the pxk matrix of Z-scores where p is total no. of variants in the dataset, and k is the no. of traits 
  # Similarly, P.matrix is the corresponding pxk matrix of p-values where p is total no. of variants in the dataset, and k is the no. of traits
  # p.threshold determines which variants are deemed to have no association marginally (default: 1e-4)
  # checks
  
#placo
# Z: vector of Z-scores of size k=2 (i.e., collection of Z-scores of a particular SNP for k=2 traits)
# VarZ: vector of variances of Z-scores (covariance assumed 0; so need to be independent traits)



#collect all placo results txt in all pairwise diseases
j = 1
rm(qq_disease_file1,qq_disease_file)
qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/PLACO/sig_placoresults/',fileorder_disease1[j],"_",fileorder_disease2[j],'.txt',sep = '')
qq_disease_file1<-fread(qq_disease,header=T)
qq_disease_file1$name<-paste(fileorder_disease1[j],"_",fileorder_disease2[j],sep = "")

for( j in c(2:length(fileorder_disease2))){
  qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/PLACO/sig_placoresults/',fileorder_disease1[j],"_",fileorder_disease2[j],'.txt',sep = '')
  qq_disease_file<-fread(qq_disease,header=T)
  qq_disease_file$name<-paste(fileorder_disease1[j],"_",fileorder_disease2[j],sep = "")
  
  qq_disease_file1<-rbind(qq_disease_file1,qq_disease_file)
  
}
write.csv(qq_disease_file1,file = '/Users/yujiezhao/Desktop/COMMO/PLACO/all_sig_placo_bon_20230607.csv')


#### draw qq plot ####
##qq plot:https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
#Quantile-quantile plots (qq-plots) can be useful for verifying that a set of values 
#come from a certain distribution. For example in a genome-wide association study, 
#we expect that most of the SNPs we are testing not to be associated with the disease. 
#Under the null, this means that the p-values we get from tests where no true association 
#exists should follow a uniform(0,1) distribution. 
#Since we're usually most interested in really small p-values, 
#we generally transform the p-values by -log10 so that the smallest values 
#near zero become the larger values and are thus easier to see.

##draw qq plot in each disease pair
for( j in c(1:length(fileorder_disease2))){
  print(j)
  qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults/',fileorder_disease1[j],"_",fileorder_disease2[j],'.txt',sep = '')
  qq_disease_file<-fread(qq_disease,header=T)
  
  qq_disease_file$sigif<-NA
  qq_disease_file$sigif[qq_disease_file$p.placo<5e-8]<-"Sig."
  qq_disease_file$sigif[qq_disease_file$p.placo>=5e-8]<-"Non-Sig."
  
  
  qq_disease_file2<-qq_disease_file[order(qq_disease_file$p.placo,decreasing=FALSE),]
  
  qq_disease_file2$obs<--log10(qq_disease_file2$p.placo)
  qq_disease_file2$exp<--log10(ppoints(length(qq_disease_file2$p.placo)))
  
  qq_disease_file2$sigif <- factor(qq_disease_file2$sigif,levels=c("Non-Sig.","Sig."),labels=c("Non-Sig.","Sig."))
  
  library(qqman)
  jpeg(file = paste('/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults_plot/',fileorder_disease1[j],"_",fileorder_disease2[j],'.png',sep = ""))
  
  qq(qq_disease_file2$p.placo, main = paste("Q-Q plot of pleiotropic p-values in ", fileorder_disease1[j],"_",fileorder_disease2[j],sep = ""),xlim = c(0, 7), ylim = c(0,12), pch = 18, col = qq_disease_file2$sigif, cex = 1, las = 1)
  
  dev.off();
  
}

