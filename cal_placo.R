rm(list = ls())

require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.1.1.R?raw=TRUE")

fileorder_disease1<-c("CHD","CHD","CHD","CHD","CHD","CHD","CHD","Asthma","Asthma","Asthma",
                      "Asthma","chronicpain","chronicpain","chronicpain",
                      "chronicpain","chronicpain","chronicpain","chronicpain",
                      "chronicpain","chronicpain","chronicpain","chronicpain",
                      "chronicpain","stroke","stroke","stroke","stroke","stroke",
                      "PUD","PUD","PUD","PUD","PUD","PUD","PUD","PUD","T2D","T2D",
                      "T2D","T2D","T2D","T2D","T2D","T2D","T2D","T2D","HTN","arthritis",
                      "arthritis","arthritis","cancerlung","cancerlung","cancerlung","cancerlung",
                      "cancerlung","lungfunction","lungfunction")

fileorder_disease2<-c("ADHD","Alcohol","ANX","MDD","Insomnia","OCD","PTSD","BP","MDD","Insomnia","AD",
                      "ADHD","Alcohol","ANX","ASD","BP","MDD","Insomnia","OCD",
                      "PTSD","SCZ","TS","AD","ADHD","BP","Insomnia","OCD","TS",
                      "ADHD","Alcohol","ANX","MDD","Insomnia","PTSD","SCZ","AD",
                      "ADHD","Alcohol","Eat","ANX","MDD","Insomnia","OCD","PTSD",
                      "SCZ","AD","Insomnia","ADHD","ASD","BP","ADHD","Eat","ANX",
                      "MDD","Insomnia","BP","SCZ")

  
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
  
  #disease2_p<-disease2_p[,c("SNP","CHR","BP","A1","A2","OR","P")]
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


##qq plot:https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R
#Quantile-quantile plots (qq-plots) can be useful for verifying that a set of values 
#come from a certain distribution. For example in a genome-wide association study, 
#we expect that most of the SNPs we are testing not to be associated with the disease. 
#Under the null, this means that the p-values we get from tests where no true association 
#exists should follow a uniform(0,1) distribution. 
#Since we're usually most interested in really small p-values, 
#we generally transform the p-values by -log10 so that the smallest values 
#near zero become the larger values and are thus easier to see.

rm(list = ls())
fileorder_disease1<-c("chronicpain","chronicpain","T2D","chronicpain","PUD","PUD","CHD","T2D",
                      "PUD","chronicpain","chronicpain","CHD","T2D","PUD","chronicpain","T2D",
                      "chronicpain","chronicpain","CHD","CHD","T2D","stroke",
                      "T2D","lungfunction","chronicpain",	"T2D","T2D","T2D","chronicpain","CHD","PUD","stroke",	"T2D",	"arthritis",	"arthritis",	"CHD",	"Asthma")
                      
fileorder_disease2<-c("ADHD","MDD","ADHD","PTSD","ADHD","MDD","ADHD","Eat","Alcohol","Alcohol",
                      "TS","MDD","MDD","PTSD","BP","PTSD","SCZ","ANX","PTSD","AD2","MS2","ADHD","OCD",
                      "SCZ","ASD","Alcohol","ANX",	"AD",	"MS2",	"Alcohol",	"SCZ",	"BP",	"SCZ",	"BP",	"MS2",	"ANX",	"MDD")
  

fileorder_disease1<-c("chronicpain",	"chronicpain",	"T2D",	"chronicpain",	"PUD",
                      "PUD",	"T2D",	"PUD",	"chronicpain",	"CHD2",	"chronicpain",
                      "T2D",	"PUD",	"chronicpain",	"T2D",	"chronicpain",	"CHD2",
                      "chronicpain",	"T2D",	"CHD2",	"T2D",	"stroke",	"CHD2",	"lungfunction",
                      "chronicpain",	"T2D",	"T2D",	"T2D",	"CHD2",	"chronicpain",	"PUD",
                      "stroke",	"T2D",	"arthritis",	"arthritis")
fileorder_disease2<-c("ADHD",	"MDD",	"ADHD",	"PTSD",	"ADHD",	"MDD",	"Eat",
                      "Alcohol",	"Alcohol",	"ADHD",	"TS",	"MDD",	"PTSD",
                      "BP",	"PTSD",	"SCZ",	"MDD",	"ANX",	"MS2",	"AD2",
                      "OCD",	"ADHD",	"PTSD",	"SCZ",	"ASD",	"Alcohol",
                      "ANX",	"AD2",	"Alcohol",	"MS2",	"SCZ",	"BP",	"SCZ",	"BP",	"MS2")



library(data.table)
j = 1
rm(qq_disease_file1,qq_disease_file)
qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/PLACO/sig_placoresults/',fileorder_disease1[j],"_",fileorder_disease2[j],'.txt',sep = '')
qq_disease_file1<-fread(qq_disease,header=T)
qq_disease_file1$name<-paste(fileorder_disease1[j],"_",fileorder_disease2[j],sep = "")

for( j in c(2:35)){
  qq_disease<-paste('/Users/yujiezhao/Desktop/COMMO/PLACO/sig_placoresults/',fileorder_disease1[j],"_",fileorder_disease2[j],'.txt',sep = '')
  qq_disease_file<-fread(qq_disease,header=T)
  qq_disease_file$name<-paste(fileorder_disease1[j],"_",fileorder_disease2[j],sep = "")
  
  qq_disease_file1<-rbind(qq_disease_file1,qq_disease_file)
  
}

write.csv(qq_disease_file1,file = '/Users/yujiezhao/Desktop/COMMO/PLACO/all_sig_placo_bon_20230607.csv')


for( j in c(4:57)){
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



rm(list = ls())
pleioloci<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci20230619.csv")

trait_unique_pair<-unique(pleioloci$Trait_pair)

setwd("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/")
pleioloci$beta_chronic<-0
pleioloci$beta_brain<-0
for (i in c(1:length(trait_unique_pair))){
  
  filename<-paste(trait_unique_pair[i],"_pleiosnp.csv",sep = "")
  pleiosnp<-read.csv(filename)
  
  pleioloci$beta_chronic[which(pleioloci$Trait_pair==trait_unique_pair[i])]<-pleiosnp[match(pleioloci[which(pleioloci$Trait_pair==trait_unique_pair[i]),]$TopSNP,pleiosnp$rsID),17]
  pleioloci$beta_brain[which(pleioloci$Trait_pair==trait_unique_pair[i])]<-pleiosnp[match(pleioloci[which(pleioloci$Trait_pair==trait_unique_pair[i]),]$TopSNP,pleiosnp$rsID),14]
  rm(pleiosnp,filename)
  
}
pleioloci$beta2<-pleioloci$beta_chronic*pleioloci$beta_brain
pleioloci$direction<-0
pleioloci$direction[pleioloci$beta2<0]="-"
pleioloci$direction[pleioloci$beta2>0]="+"
write.csv(pleioloci,file = "all_fuma_pleioloci_direction_20230717.csv")

multi_pleioloci<-pleioloci[duplicated(pleioloci$Locus_boundary),]
all_multipleioloci = list()
for (i in c(1:length(multi_pleioloci$X))){
  current<-subset(pleioloci,Locus_boundary==multi_pleioloci[i,14])
  all_multipleioloci<-rbind(all_multipleioloci,current)
}
write.csv(all_multipleioloci,file = "all_multi_fuma_pleioloci_direction_20230717.csv")

rm(list = ls())
all_multipleioloci<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_multi_fuma_pleioloci_direction_20230717.csv")


#ANNOVAR
intronic<-subset(pleioloci,Functional_annotation=="ncRNA_exonic")
unique<-unique(intronic$TopSNP)




