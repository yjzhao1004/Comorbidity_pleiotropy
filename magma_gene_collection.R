###gene_collection
rm(list = ls())
fileorder_disease1<-c("chronicpain","chronicpain","chronicpain","chronicpain","chronicpain","chronicpain",
                      "chronicpain","chronicpain","chronicpain","CHD","CHD","CHD","T2D","T2D","T2D","T2D","T2D",
                      "T2D","T2D","T2D","T2D","arthritis","arthritis","lungfunction","PUD","PUD","PUD","Stroke","Stroke")

fileorder_disease2<-c("ADHD","ANX","ASD","BP","MDD","MS","PTSD","SCZ","TS","ADHD","Alcohol","MDD","ADHD","ANX","Alcohol",
                      "EAT","MDD","MS","OCD","PTSD","SCZ","BP","MS","SCZ","ADHD","Alcohol","SCZ","ADHD","BP")

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

collection_descendZ<-collection[order(collection$ZSTAT),]
write.csv(collection_descendZ,file = "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_descendZ_20230709.csv")

###
##match genes
setwd( "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/")

magma_gene<-read.csv("magma_gene_20230623.csv")

unique_gene<-unique(magma_gene$Gene_SYMBOL)

geneinfo<-read.csv("/Volumes/ZHAOYUJIE/gem/gem_new20230624/probeinfo.csv")

gene_select<-subset(magma_gene,magma_gene$EntrezID%in%geneinfo$X1)


unique_gene_select<-unique(gene_select$EntrezID)

write.csv(unique_gene_select,file = "/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/unique_gene_ahba.csv",quote = FALSE,row.names = FALSE)

getwd()


##
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


