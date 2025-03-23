BiocManager::install("gep2pep")
library(gep2pep)
install.packages("WriteXLS")
library("WriteXLS")
rm(list = ls())
setwd("/Volumes/yjzhao/COMMO_new/enrichmentnew/")
# download.file("http://dsea.tigem.it/data/Cmap_MSigDB_v6.1_PEPs.tar.gz",
# "Cmap_MSigDB_v6.1_PEPs.tar.gz")
# untar("/Volumes/yjzhao/COMMO_new/enrichmentnew/Cmap_MSigDB_v6.1_PEPs.tar")
rpBig <- openRepository("/Volumes/yjzhao/COMMO_new/enrichmentnew/Cmap_MSigDB_v6.1_PEPs")
gene<-fread("/Volumes/yjzhao/COMMO_new/magma/magma_loci/unique_gene_symbolSep9.txt")

pws <- gene2pathways(rpBig, "ZNRF3")
psea <- PathSEA(rpBig, pws)
collect<-list()

collect$C4_CGN<-getResults(psea, "C4_CGN")
collect$C4_CGN$condition<-colnames(getDetails(psea, "C4_CGN"))
collect$C4_CGN<-collect$C4_CGN%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C4_CM<-getResults(psea, "C4_CM")
collect$C4_CM$condition<-colnames(getDetails(psea, "C4_CM"))
collect$C4_CM<-collect$C4_CM%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C3_TFT<-getResults(psea, "C3_TFT")
collect$C3_TFT$condition<-colnames(getDetails(psea, "C3_TFT"))
collect$C3_TFT<-collect$C3_TFT%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C3_MIR<-getResults(psea, "C3_MIR")
collect$C3_MIR$condition<-colnames(getDetails(psea, "C3_MIR"))
collect$C3_MIR<-collect$C3_MIR%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$`C2_CP:KEGG`<-getResults(psea, "C2_CP:KEGG")
collect$`C2_CP:KEGG`$condition<-colnames(getDetails(psea, "C2_CP:KEGG"))
collect$`C2_CP:KEGG`<-collect$`C2_CP:KEGG`%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$`C2_CP`<-getResults(psea, "C2_CP")
collect$`C2_CP`$condition<-colnames(getDetails(psea, "C2_CP"))
collect$`C2_CP`<-collect$`C2_CP`%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$`C2_CP:REACTOME`<-getResults(psea, "C2_CP:REACTOME")
collect$`C2_CP:REACTOME`$condition<-colnames(getDetails(psea, "C2_CP:REACTOME"))
collect$`C2_CP:REACTOME`<-collect$`C2_CP:REACTOME`%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$`C2_CGP`<-getResults(psea, "C2_CGP")
collect$`C2_CGP`$condition<-colnames(getDetails(psea, "C2_CGP"))
collect$`C2_CGP`<-collect$`C2_CGP`%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$`C7_C7`<-getResults(psea, "C7_C7")
collect$`C7_C7`$condition<-colnames(getDetails(psea, "C7_C7"))
collect$`C7_C7`<-collect$`C7_C7`%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C5_BP<-getResults(psea, "C5_BP")
collect$C5_BP$condition<-colnames(getDetails(psea, "C5_BP"))
collect$C5_BP<-collect$C5_BP%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C5_CC<-getResults(psea, "C5_CC")
collect$C5_CC$condition<-colnames(getDetails(psea, "C5_CC"))
collect$C5_CC<-collect$C5_CC%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$C5_MF<-getResults(psea, "C5_MF")
collect$C5_MF$condition<-colnames(getDetails(psea, "C5_MF"))
collect$C5_MF<-collect$C5_MF%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$ARCHIVED_C5_CC<-getResults(psea, "ARCHIVED_C5_CC")
collect$ARCHIVED_C5_CC$condition<-colnames(getDetails(psea, "ARCHIVED_C5_CC"))
collect$ARCHIVED_C5_CC<-collect$ARCHIVED_C5_CC%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$ARCHIVED_C5_BP<-getResults(psea, "ARCHIVED_C5_BP")
collect$ARCHIVED_C5_BP$condition<-colnames(getDetails(psea, "ARCHIVED_C5_BP"))
collect$ARCHIVED_C5_BP<-collect$ARCHIVED_C5_BP%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

collect$ARCHIVED_C5_MF<-getResults(psea, "ARCHIVED_C5_MF")
collect$ARCHIVED_C5_MF$condition<-colnames(getDetails(psea, "ARCHIVED_C5_MF"))
collect$ARCHIVED_C5_MF<-collect$ARCHIVED_C5_MF%>%
  filter(PV<1*10^-6)%>%
  filter(ES>0.5)%>%
  arrange(-ES)

save(collect,file = "./drug/ZNRF3.RData")



###########################
##collect
###########################
rm(list = ls())
setwd("/Volumes/yjzhao/COMMO_new/enrichmentnew/drug/")

files <- list.files(pattern = "\\.RData$", full.names = TRUE)
filenames <- sub("^.*/(.*)\\.RData$", "\\1", files)

results <- data.frame(condition = character(), ES = numeric(), PV = numeric(), stringsAsFactors = FALSE)


for (i in c(1:length(filenames))) {
  load(files[i])  
  for (sublist in collect) {
    if ("condition" %in% names(sublist) && "ES" %in% names(sublist) && "PV" %in% names(sublist)) {
      results <- rbind(results, data.frame(condition = sublist$condition, ES = sublist$ES, PV = sublist$PV))
    }
  }
  results<-results%>%
    arrange(-ES)
  unique_results <- results[!duplicated(results$condition), ]
  
  if (nrow(unique_results)>=10){
    unique_results<-unique_results[c(1:10),]
    unique_results$gene<-filenames[i]
  }else{
    unique_results$gene<-filenames[i]
  }
  
  write.csv(unique_results,file = paste0("/Volumes/yjzhao/COMMO_new/enrichmentnew/drug_gene_max10/",filenames[i],"_drug.csv"),row.names = F)
}



###combine
subfolder <- "/Volumes/yjzhao/COMMO_new/enrichmentnew/drug_gene_max10" 
file_list <- list.files(path = subfolder, pattern = "*.csv", full.names = TRUE)
combined_data <- data.frame()
for (file in file_list) {
  data <- read.csv(file)        
  combined_data <- rbind(combined_data, data) 
}

write.csv(combined_data, "/Volumes/yjzhao/COMMO_new/enrichmentnew/combined_drug_gene_max10.csv", row.names = FALSE)  # 替换为你想要的输出文件名

###unique drug
combined_data<-fread("/Volumes/yjzhao/COMMO_new/enrichmentnew/combined_drug_gene_max10.csv")
uniquedrug_combined_data <- combined_data %>%
  group_by(condition) %>%
  arrange(desc(ES), PV) %>%
  slice(1) %>%
  ungroup()

gene_disease<-fread("/Volumes/yjzhao/COMMO_new/magma/magma_loci/magma_gene_loci.csv")
gene_disease_drug<-subset(gene_disease,SYMBOL%in%B)
gene_disease_drug<-gene_disease_drug%>%
  select(SYMBOL,diseasepair)%>%
  distinct()%>%
  rename(gene = SYMBOL)

combined_data<-combined_data%>%
  left_join(gene_disease_drug,by = "gene")
write.csv(combined_data, "/Volumes/yjzhao/COMMO_new/enrichmentnew/combined_drug_gene_max10.csv", row.names = FALSE)  


combined_data<-fread("/Volumes/yjzhao/COMMO_new/enrichmentnew/combined_drug_gene_max10.csv")

new_data <- combined_data %>%
  group_by(condition) %>%  
  summarise(
    gene = paste(gene, collapse = "/"),  
    count = n() 
  )%>%
  left_join(uniquedrug_combined_data[,c(1,2,3)],by = "condition")%>%
  rename(Description = condition,GeneRatio = ES, pvalue = PV,geneID = gene, Count = count)%>%
  select(Description,GeneRatio,pvalue,geneID,Count)

write.csv(new_data,file = "/Volumes/yjzhao/COMMO_new/enrichmentnew/drug_draw_max10.csv",row.names = F)

