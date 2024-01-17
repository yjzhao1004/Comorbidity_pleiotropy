rm(list = ls())
library(data.table)
setwd("/Users/yujiezhao/Desktop/COMMO/Enrichment/phenotype/")
MGI<-fread("informatics.jax.org_downloads_reports_HMD_HumanPhenotype.rpt.txt")
MGI<-MGI[,c(1,2,5)]
names(MGI) <- c("Human_marker_symbol", "EntrezID","PhenotypeID" )

MGI<-subset(MGI, PhenotypeID!="")

pleiogene<-fread("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_20230623.csv")
pleiogene<-fread("/Users/yujiezhao/Desktop/COMMO/Enrichment/phenotype/kjj.csv")


pleiogene_MGI<-subset(MGI,EntrezID%in%pleiogene$EntrezID)
pleiogene_MGI<-pleiogene_MGI[!duplicated(pleiogene_MGI$Human_marker_symbol),]

All_phenotype<-c("MP:0005381","MP:0010768","MP:0005384",
                 "MP:0005387","MP:0005386","MP:0003631",
                 "MP:0005380","MP:0005378","MP:0005385",
                 "MP:0005397","MP:0005376","MP:0005377",
                 "MP:0005382","MP:0005388","MP:0010771",
                 "MP:0005391","MP:0005389","MP:0005369",
                 "MP:0005370","MP:0005367","MP:0001186",
                 "MP:0002006","MP:0005379","MP:0005375",
                 "MP:0005390","MP:0005371","MP:0005394",
                 "MP:0005381","MP:0005386","MP:0005381",
                 "MP:0005386")
All_phenotype<-unique(All_phenotype)

#pleiogene
Phenotypes <- vector(mode = "list",length = 27)
for (i in 1:27) {
  Phenotypes[[i]]<-grep(All_phenotype[i],pleiogene_MGI$PhenotypeID)
  print(i)
}

ANTI0 <- pleiogene_MGI[Phenotypes[[1]],]
ANTI0$PHETYPE<-All_phenotype[1]
for (i in 2:27) {
  ANTI<-pleiogene_MGI[Phenotypes[[i]],]
  ANTI$PHETYPE<-All_phenotype[i]
  ANTI0 <- rbind(ANTI,ANTI0)
  print(i)
}
pleiogene_MGIphe<-ANTI0

##non-pleiotropic genes
MGI_nonpleiogenes<-subset(MGI,EntrezID %in% setdiff(MGI$EntrezID,pleiogene_MGI$EntrezID))
MGI_nonpleiogenes<-subset(MGI,EntrezID %in% setdiff(MGI$EntrezID,pleiogene_MGI$EntrezID))

MGI_nonpleiogenes<-MGI_nonpleiogenes[!duplicated(MGI_nonpleiogenes$EntrezID),]

Phenotypes_2 <- vector(mode = "list",length = 27)
for (i in 1:27) {
  Phenotypes_2[[i]]<-grep(All_phenotype[i],MGI_nonpleiogenes$PhenotypeID)
  #print(i)
}
ANTI0 <- MGI_nonpleiogenes[Phenotypes_2[[1]],]
ANTI0$PHETYPE<-All_phenotype[1]
for (i in 2:27) {
  ANTI<-MGI_nonpleiogenes[Phenotypes_2[[i]],]
  ANTI$PHETYPE<-All_phenotype[i]
  ANTI0 <- rbind(ANTI,ANTI0)
  #print(i)
}
  
nonpleiogene_MGIphe<-ANTI0

rm(ANTI,ANTI0)  
  

#collection
collection <- vector(mode = "list",length = 7)

for(i in c(1:27)){
  print(All_phenotype[i])
  collection[[1]][i]<-All_phenotype[i]
  collection[[4]][i]<-length(which(nonpleiogene_MGIphe$PHETYPE==All_phenotype[i]))
  collection[[5]][i]<-length(MGI_nonpleiogenes$EntrezID)-length(which(nonpleiogene_MGIphe$PHETYPE==All_phenotype[i]))
  collection[[2]][i]<-length(which(pleiogene_MGIphe$PHETYPE==All_phenotype[i]))
  collection[[3]][i]<-length(pleiogene_MGI$EntrezID)-length(which(pleiogene_MGIphe$PHETYPE==All_phenotype[i]))
  fisher_matrix<-matrix(c(as.numeric(collection[[2]][i]),as.numeric(collection[[3]][i]),as.numeric(collection[[4]][i]),as.numeric(collection[[5]][i])),nrow=2)
  results<-fisher.test(x = fisher_matrix, alternative = "two.sided")
  collection[[6]][i]<-results$p.value
  collection[[7]][i]<-paste(pleiogene_MGIphe$Human_marker_symbol[which(pleiogene_MGIphe$PHETYPE==All_phenotype[i])],collapse = ",")
  #Disease_pair<-subset(pleiogene[,c(2,5)],Gene_SYMBOL%in%pleiogene_MGIphe$Human_marker_symbol[which(pleiogene_MGIphe$PHETYPE==All_phenotype[i])])
  #collection[[8]][i]<-paste(Disease_pair$trait_pairs,collapse = ",")
  }
names(collection)<-c("PhenotypeID", "Pleiotropic_gene_group_1","Pleiotropic_gene_group_0","NonPleiotropic_gene_group_1","NonPleiotropic_gene_group_0","P_fisher","Genes")

write.csv(collection,file = "phenotype_enrichment_results_kjj.csv")



##plot
rm(list = ls())
Phenotypes_enrichment<-read.csv("phenotype_enrichment_results_20230704.csv")
Phenotypes_enrichment$term<-paste(Phenotypes_enrichment$PhenotypeID,Phenotypes_enrichment$Phenotype,sep = " ")
Phenotypes_enrichment_sig<-subset(Phenotypes_enrichment,P_fisher<0.05)

Phenotypes_enrichment_sig<-Phenotypes_enrichment_sig[order(Phenotypes_enrichment_sig$P_fisher,decreasing = FALSE),]
names(Phenotypes_enrichment_sig)[3]<-"No.of_genes"
Phenotypes_enrichment_sig$term <- factor(Phenotypes_enrichment_sig$term, levels = Phenotypes_enrichment_sig$term)
Phenotypes_enrichment_sig$LogP <- -log10(Phenotypes_enrichment_sig$P_fisher)


ggplot(Phenotypes_enrichment_sig,aes(No.of_genes,term))+
  geom_bar(aes(y=reorder(term,-P_fisher),x=No.of_genes,fill=LogP)
           ,stat='identity')+
  scale_fill_gradient(low="#FAF0E6",high="#CD0000")+
  #facet_grid(Category~., scale = 'free_y', space = 'free_y') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
)




