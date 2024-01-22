library(data.table)
rm(list = ls())
#take T2D and ALD for example 
#input and quality control for the gwas data of these two diseases

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
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
T2D2<-T2D[T2D$CHR==6,]
T2D3<-T2D2[T2D2$BP>25000000,]
T2D3<-T2D3[T2D3$BP<35000000,]
T2D_qc1<-subset(T2D,SNP %in% setdiff(T2D$SNP,T2D3$SNP))
#keeping biallelic SNPs with minor allele frequency (MAF) > 0.01
T2D_qc2<-subset(T2D_qc1,MAF>0.01)
rm(T2D2,T2D3,T2D_qc1)
T2D_qc3<-T2D_qc2[,c("SNP","CHR","BP","beta","SE","pval_nominal")]
rm(T2D,T2D_qc2)
T2D<-T2D_qc3
T2D$pval_nominal<-as.numeric(T2D$pval_nominal)
rm(T2D_qc3)

#ALD
ALD<-fread("/Volumes/ZHAOYUJIE/COMMO/Alcohol/Alcohol_dependency_pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt",header=T)
ALD$beta<-log(ALD$OR)
names(ALD)[9]<-"pval_nominal"
#excluding SNPs in major histocompatibility complex region (MHC, chr 6: 25–35 Mb) 
#due to its complex LD structure
ALD2<-ALD[ALD$CHR==6,]
ALD3<-ALD2[ALD2$BP>25000000,]
ALD3<-ALD3[ALD3$BP<35000000,]
ALD_qc1<-subset(ALD,SNP %in% setdiff(ALD$SNP,ALD3$SNP))
ALD_qc3<-ALD_qc1[,c("SNP","CHR","BP","beta","SE","pval_nominal")]
rm(ALD2,ALD3,ALD_qc1,ALD)
ALD<-ALD_qc3
rm(ALD_qc3)



##
placo_snp<-fread("/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults/T2D_ALD.txt")

fuma_placo_loci_T2D_ALD<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_ALD_pleioloci.csv")

a<-strsplit(fuma_placo_loci_T2D_ALD$Locus_boundary, " - ")
a<-as.matrix(a)

input <- merge(T2D, ALD, by="SNP", all=FALSE, suffixes=c("_T2D","_ALD"))
head(input)

placo_snp_inloci<-subset(placo_snp,SNP%in%input$SNP)

dir.create("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/T2D_ALD")
setwd("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/T2D_ALD")

for(i in c(1:length(fuma_placo_loci_T2D_BP$Locus_boundary))){
  inputcurrent<-subset(input,input$BP_T2D<as.numeric(a[[i]][2])&input$BP_T2D>as.numeric(a[[i]][1]))
  inputcurrent<-subset(inputcurrent,inputcurrent$BP_ALD<as.numeric(a[[i]][2])&inputcurrent$BP_ALD>as.numeric(a[[i]][1]))
  
  placo_snp_inloci<-subset(placo_snp,SNP%in%inputcurrent$SNP)
  write.table(placo_snp_inloci[,c("SNP","p.placo")],file =paste("placo_gwas_snpp_loci",as.character(i),".txt",sep =""),quote = FALSE,row.names = FALSE,col.names = FALSE)
}


#use magma: For each pleiotropic loci of disease pairs, we identified the genes significantly associated with the pleiotropic SNPs based on PLACO outputs.
"/Users/yujiezhao/Desktop/COMMO/MAGMA/magma_py/magma_v1/magma --bfile /Users/yujiezhao/Desktop/COMMO/MAGMA/magma_py/g1000_eur/g1000_eur --pval /Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/T2D_ALD/placo_gwas_snpp_loci1.txt  N=492191 --gene-annot /Users/yujiezhao/Desktop/COMMO/MAGMA/magma_py/g1000_eur.genes.annot --out /Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/T2D_ALD/genebased_loci1"


##

rm(list = ls())
setwd("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/T2D_ALD")
fuma_placo_loci_T2D_ALD<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/T2D_ALD_pleioloci.csv")

T2D_magma_gene<-read.table("/Volumes/ZHAOYUJIE/COMMO/T2D/FUMA_jobT2D/magma.genes.out",header = TRUE)
ALD_magma_gene<-read.table("/Volumes/ZHAOYUJIE/COMMO/BP/FUMA_jobALD/magma.genes.out",header = TRUE)

magma_T2D_ALD <- merge(T2D_magma_gene, ALD_magma_gene, by="GENE", all=FALSE, suffixes=c("_T2D","_ALD"))
head(magma_T2D_ALD)


i = 1
loci_placo_gene<-read.table(paste("genebased_loci",as.character(i),".genes.out",sep = ""),header = TRUE)
loci_placo_gene<-subset(loci_placo_gene,START%in%magma_T2D_ALD$START_T2D)
names(magma_T2D_ALD)[3]<-"START"
names(loci_placo_gene)[9]<-"P_placo"
loci_placo_gene_twod<-merge(loci_placo_gene,magma_T2D_ALD[,c(3,9,10,18)],by = "START")
loci_placo_gene_twod<-loci_placo_gene_twod[,c(2,3,1,4,8,9,10,12,11)]
loci_placo_gene_twod$trait_pairs<-"T2D_ALD"
loci_placo_gene_twod$No_loci<-i
loci_placo_gene_twod$P_placo_bonf<-loci_placo_gene_twod$P_placo*length(loci_placo_gene_twod$GENE)
loci_placo_gene_twod$STOP<-magma_T2D_ALD$STOP_T2D[match(loci_placo_gene_twod$SYMBOL_T2D,magma_T2D_ALD$SYMBOL_T2D)]
loci_placo_gene_twod_collection<-loci_placo_gene_twod


for(i in c(2:length(fuma_placo_loci_T2D_ALD$Locus_boundary))){
  loci_placo_gene<-read.table(paste("genebased_loci",as.character(i),".genes.out",sep = ""),header = TRUE)
  loci_placo_gene<-subset(loci_placo_gene,START%in%magma_T2D_ALD$START)
  names(loci_placo_gene)[9]<-"P_placo"
  names(magma_T2D_ALD)[3]<-"START"
  loci_placo_gene_twod<-merge(loci_placo_gene,magma_T2D_ALD[,c(3,9,10,18)],by = "START")
  loci_placo_gene_twod<-loci_placo_gene_twod[,c(2,3,1,4,8,9,10,12,11)]
  if (length(loci_placo_gene_twod$GENE)==0){
    next
  }
  else{
    loci_placo_gene_twod$trait_pairs<-"T2D_ALD"
    loci_placo_gene_twod$No_loci<-i
    loci_placo_gene_twod$P_placo_bonf<-loci_placo_gene_twod$P_placo*length(loci_placo_gene_twod$GENE)
    loci_placo_gene_twod$STOP<-magma_T2D_ALD$STOP_T2D[match(loci_placo_gene_twod$SYMBOL_T2D,magma_T2D_ALD$SYMBOL_T2D)]
    loci_placo_gene_twod_collection<-rbind(loci_placo_gene_twod_collection,loci_placo_gene_twod)
  }
  }
write.csv(loci_placo_gene_twod_collection,file = "loci_placo_gene_twod_collection.csv")

loci_placo_gene_twod_collection_sig<-subset(loci_placo_gene_twod_collection,P_placo_bonf<0.05)
loci_placo_gene_twod_collection_sig<-subset(loci_placo_gene_twod_collection_sig,P_T2D<0.05)
loci_placo_gene_twod_collection_sig<-subset(loci_placo_gene_twod_collection_sig,P_ALD<0.05)

write.csv(loci_placo_gene_twod_collection_sig,file = "loci_placo_gene_twod_collection_sig.csv")


### gene_collection
rm(list = ls())
fileorder_disease1<-c("MCP","MCP","MCP","MCP","MCP","MCP",
                      "MCP","MCP","MCP","CAD","CAD","CAD","T2D","T2D","T2D","T2D","T2D",
                      "T2D","T2D","T2D","T2D","RA","RA","CLD","PUD","PUD","PUD","Stroke","Stroke")

fileorder_disease2<-c("ADHD","ANX","ASD","BP","MDD","MS","PTSD","SCZ","TS","ADHD","ALD","MDD","ADHD","ANX","ALD",
                      "EAT","MDD","MS","OCD","PTSD","SCZ","BP","MS","SCZ","ADHD","ALD","SCZ","ADHD","BP")

j=1
gene_sig<-read.csv(paste("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/",fileorder_disease1[j],"_",fileorder_disease2[j],"/loci_placo_gene_twod_collection_sig.csv",sep = ""))

gene_sig$Gene_position<-paste(gene_sig$START,"-",gene_sig$STOP,sep = " ")

fuma_placo_loci_pairs<-read.csv(paste("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/",fileorder_disease1[j],"_",fileorder_disease2[j],"_pleioloci.csv",sep = ""))

gene_sig$Loci_position<-fuma_placo_loci_pairs$Locus_boundary[match(gene_sig$No_loci,fuma_placo_loci_pairs$X)]
variable.names(gene_sig)
names(gene_sig)[which(variable.names(gene_sig)==paste("P_",fileorder_disease1[j],sep = ""))]<-"P_chroniccondition"
names(gene_sig)[which(variable.names(gene_sig)==paste("P_",fileorder_disease2[j],sep = ""))]<-"P_braindisorder"
names(gene_sig)[which(variable.names(gene_sig)==paste("SYMBOL_",fileorder_disease1[j],sep = ""))]<-"Gene_SYMBOL"
current<-gene_sig[,c("trait_pairs","No_loci","Loci_position","Gene_SYMBOL","GENE","Gene_position","ZSTAT","P_placo","P_placo_bonf","P_chroniccondition","P_braindisorder")]
collection<-current


for (j in c(2:length(fileorder_disease1))){
  gene_sig<-read.csv(paste("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/",fileorder_disease1[j],"_",fileorder_disease2[j],"/loci_placo_gene_twod_collection_sig.csv",sep = ""))
  
  if (length(gene_sig$GENE)==0){
    next
  }
  else{
    gene_sig$Gene_position<-paste(gene_sig$START,"-",gene_sig$STOP,sep = " ")
    fuma_placo_loci_pairs<-read.csv(paste("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/",fileorder_disease1[j],"_",fileorder_disease2[j],"_pleioloci.csv",sep = ""))
    gene_sig$Loci_position<-fuma_placo_loci_pairs$Locus_boundary[match(gene_sig$No_loci,fuma_placo_loci_pairs$X)]
    variable.names(gene_sig)
    names(gene_sig)[which(variable.names(gene_sig)==paste("P_",fileorder_disease1[j],sep = ""))]<-"P_chroniccondition"
    names(gene_sig)[which(variable.names(gene_sig)==paste("P_",fileorder_disease2[j],sep = ""))]<-"P_braindisorder"
    names(gene_sig)[which(variable.names(gene_sig)==paste("SYMBOL_",fileorder_disease1[j],sep = ""))]<-"Gene_SYMBOL"
    current<-gene_sig[,c("trait_pairs","No_loci","Loci_position","Gene_SYMBOL","GENE","Gene_position","ZSTAT","P_placo","P_placo_bonf","P_chroniccondition","P_braindisorder")]
    collection<-rbind(collection,current)
  }
}

write.csv(collection,file = "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_20230709.csv")

#collection_descendZ<-collection[order(collection$ZSTAT),]
#write.csv(collection_descendZ,file = "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_descendZ_20230709.csv")


##pleiotropy gene in multiple disease pairs
genepleio<-read.csv("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_20230709.csv")
genepleio<-genepleio[order(genepleio$P_placo_bonf,decreasing = FALSE),]
uniquegene<-unique(genepleio$GENE)
gene_multi=list()
genepleio_order=list()
for (i in c(1:length(uniquegene))){
  currentgenepleio<-subset(genepleio,GENE==uniquegene[i])
  gene_multi$geneid[i]<-uniquegene[i]
  gene_multi$genesymbol[i]<-paste(unique(currentgenepleio$Gene_SYMBOL),collapse = ",")
  gene_multi$trait_pairs[i]<-paste(unique(currentgenepleio$trait_pairs),collapse = ",")
  gene_multi$trait_pairsnum[i]<-length(unique(currentgenepleio$trait_pairs))
  genepleio_order<-rbind(genepleio_order,currentgenepleio)
  
}
gene_multi<-as.data.frame(gene_multi)
setwd("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/")
write.csv(genepleio_order,file = "magma_gene_ordergene_20230718.csv")
write.csv(gene_multi,file = "gene_multi_collection_20230718.csv")


###
##match genes in AHBA
setwd( "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/")
magma_gene<-read.csv("magma_gene_20230623.csv")
unique_gene<-unique(magma_gene$Gene_SYMBOL)
geneinfo<-read.csv("~/AHBA_gem/probeinfo.csv")
gene_select<-subset(magma_gene,magma_gene$EntrezID%in%geneinfo$X1)
unique_gene_select<-unique(gene_select$EntrezID)
write.csv(unique_gene_select,file = "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/unique_gene_ahba.csv",quote = FALSE,row.names = FALSE)





