##
rm(list=ls())
setwd("/Users/yujiezhao/Desktop/COMMO/Enrichment/tissue_fuma/FUMA_gene2func115815")
tissue_specificity<-read.csv("gtex_v8_ts_DEG.csv",header = TRUE)
tissue_specificity_DEG_2side<-subset(tissue_specificity,Category=="DEG.twoside")
tissue_specificity_DEG_2side_sig<-subset(tissue_specificity_DEG_2side,adjP<0.05)


library(ggplot2)
tissue_specificity_DEG_2side_sig$LogP <- -log10(tissue_specificity_DEG_2side_sig$adjP)
tissue_specificity_DEG_2side_sig$GeneSet <- factor(tissue_specificity_DEG_2side_sig$GeneSet, levels = tissue_specificity_DEG_2side_sig$GeneSet)

ggplot(tissue_specificity_DEG_2side_sig,aes(N_overlap,GeneSet))+
  geom_bar(aes(y=reorder(GeneSet,-adjP),x=N_overlap,fill=LogP)
           ,stat='identity')+
  scale_fill_gradient(low="#FAF0E6",high="#CD0000")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))

gene<-read.table("geneIDs.txt",header = TRUE)
for(i in c(1:length(tissue_specificity_DEG_2side$Category))){
  gene_current<-tissue_specificity_DEG_2side$genes[i]
  gene_current_list <- strsplit(gene_current, ":")
  gene_current_select<-subset(gene,ensg%in%gene_current_list[[1]])
  tissue_specificity_DEG_2side$Gene_symbol[i]<-paste(gene_current_select$symbol,collapse = ",")
  rm(gene_current,gene_current_list,gene_current_select)
  }

write.csv(tissue_specificity_DEG_2side[,c(2,3,4,5,6,8)],file = "tissue_specificity_gtex_v8_results_20230711
          .csv")

###GS
rm(list = ls())
GS_enrichment<-read.csv("GS.csv")
GS_enrichment_BP<-subset(GS_enrichment,Category == "GO_bp")
GS_enrichment_BP<-GS_enrichment_BP[order(GS_enrichment_BP$adjP,decreasing = FALSE),]
#GS_enrichment_BP<-GS_enrichment_BP[c(1:20),]

GS_enrichment_BP$LogP <- -log10(GS_enrichment_BP$adjP)

GS_enrichment_BP<-subset(GS_enrichment_BP[c(1:5),],LogP>2)
ggplot(GS_enrichment_BP,aes(N_overlap,GeneSet))+
  geom_bar(aes(y=reorder(GeneSet,N_overlap),x=N_overlap,fill=LogP)
           ,stat='identity')+
  scale_fill_gradient(low="#FAF0E6",high="#8B0000")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))


GS_enrichment<-read.csv("GS.csv")
GS_enrichment_MF<-subset(GS_enrichment,Category == "GO_mf")
GS_enrichment_MF<-GS_enrichment_MF[order(GS_enrichment_MF$adjP,decreasing = FALSE),]
GS_enrichment_MF$LogP <- -log10(GS_enrichment_MF$adjP)


GS_enrichment<-read.csv("GS.csv")
GS_enrichment_CC<-subset(GS_enrichment,Category == "GO_cc")
GS_enrichment_CC$LogP <- -log10(GS_enrichment_CC$adjP)

GS_enrichment<-read.csv("GS.csv")
GS_enrichment_kegg<-subset(GS_enrichment,Category == "KEGG")
GS_enrichment_kegg$LogP <- -log10(GS_enrichment_kegg$adjP)


GS_enrichment4<-rbind(GS_enrichment_BP,GS_enrichment_MF,GS_enrichment_CC,GS_enrichment_kegg)
write.csv(GS_enrichment4,file = "GS_enrichment4.csv")

GS_enrichment4<-read.csv("/Users/yujiezhao/Desktop/COMMO/Enrichment/tissue_fuma/FUMA_gene2func115815/GS_enrichment4.csv")
GS_enrichment4<-GS_enrichment4[order(GS_enrichment4$adjP,decreasing = FALSE),]
GS_enrichment4<-GS_enrichment4[c(1:20),]
#GS_enrichment4<-subset(GS_enrichment4,LogP>2)
GS_enrichment4$term<-paste(GS_enrichment4$Category,GS_enrichment4$GeneSet,sep = ": ")

GS_enrichment4$term <- factor(GS_enrichment4$term, levels = GS_enrichment4$term)


ggplot(GS_enrichment4,aes(N_overlap,term))+
  geom_bar(aes(y=reorder(term,-adjP),x=N_overlap,fill=LogP)
           ,stat='identity')+
  scale_fill_gradient(low="#FAF0E6",high="#CD0000")+
  #geom_text(aes(label = adjP))+
  #facet_grid(Category~., scale = 'free_y', space = 'free_y') +
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
)

