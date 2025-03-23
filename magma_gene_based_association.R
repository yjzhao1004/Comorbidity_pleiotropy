rm(list = ls())
### Select significant genes from FUMA magma results and Bonferroni adjusted for gene number
library(dplyr)
library(data.table)
diseasename <- basename(list.dirs("/Volumes/yjzhao/COMMO_new/FUMA_results/"))[-1]

## Traits with pairwise pleiotropy loci
for (i in c(1:length(diseasename))) {
  path.file <- paste0("/Volumes/yjzhao/COMMO_new/FUMA_results/", diseasename[i], "/", "magma.genes.out")
  current_magmafile <- fread(path.file)    
  current_magmafile <- current_magmafile %>%
    filter(P < 0.05 / nrow(current_magmafile))
  current_magmafile$trait <- diseasename[i]
  write.csv(current_magmafile, file = paste0("/Volumes/yjzhao/COMMO_new/magma/gene_sig/", diseasename[i], "_sigbon_gene.csv"), row.names = F)
}

# Combine all traits results
setwd("/Volumes/yjzhao/COMMO_new/magma/gene_sig")
folder_path <- "/Volumes/yjzhao/COMMO_new/magma/gene_sig/"
file_suffix <- "_sigbon_gene.csv" 
file_list <- list.files(folder_path, pattern = file_suffix, full.names = TRUE)
combined_data <- data.frame()
for (file in file_list) {
  data <- read.csv(file, header = TRUE)  
  combined_data <- rbind(combined_data, data)
}
write.csv(combined_data, "combined_sigbon_gene.csv", row.names = FALSE)


## Gene count
gene_collection <- matrix(0, ncol = length(diseasename), nrow = 1)
colnames(gene_collection) <- diseasename
for (i in c(1:length(diseasename))) {
  path.file <- paste0("/Volumes/yjzhao/COMMO_new/FUMA_results/", diseasename[i], "/", "magma.genes.out")
  current_magmafile <- fread(path.file)    
  gene_collection[, diseasename[i]] <- nrow(current_magmafile)
}
gene_collection <- as.data.frame(gene_collection)
write.csv(gene_collection, "gene_num.csv", row.names = FALSE)

### Separate into different loci
rm(list = ls())
magma_gene <- fread("/Volumes/yjzhao/COMMO_new/magma/combined_sigbon_gene.csv")
loci_fuma <- fread("/Volumes/yjzhao/COMMO_new/loci_fuma/combined_conjfdr_snp_FUMA.csv")

process_row <- function(magma_row, loci_fuma) {
  # 1) Extract the CHR value from the first row of magma_gene
  chr_value <- magma_row$CHR
  
  # 2) Filter loci_fuma rows where CHR value matches
  loci_chr_match <- loci_fuma[loci_fuma$chr == chr_value, ]
  
  # 3) Further filter loci_fuma based on pos column conditions
  loci_filtered <- subset(loci_chr_match, pos > magma_row$START & pos < magma_row$STOP)
  
  # Return NULL if no matching rows found
  if (nrow(loci_filtered) == 0) return(NULL)
  
  # 4) Merge information for multiple matching rows
  loci_filtered$GenomicLocusstart <- loci_filtered$start
  loci_filtered$GenomicLocusend <- loci_filtered$end
  loci_filtered$topsnp <- loci_filtered$SNP
  loci_filtered$topsnpnearest <- loci_filtered$NearestGene
  loci_filtered$diseasepair <- loci_filtered$diseasepair
  loci_filtered$loci_num <- loci_filtered$locusnum
  
  # Combine results into a new row
  merged_row <- magma_row
  merged_row$GenomicLocusstart <- paste(loci_filtered$GenomicLocusstart, collapse = ",")
  merged_row$GenomicLocusend <- paste(loci_filtered$GenomicLocusend, collapse = ",")
  merged_row$topsnp <- paste(loci_filtered$topsnp, collapse = ",")
  merged_row$topsnpnearest <- paste(loci_filtered$topsnpnearest, collapse = ",")
  merged_row$diseasepair <- paste(loci_filtered$diseasepair, collapse = ",")
  merged_row$loci_num <- paste(loci_filtered$loci_num, collapse = ",")
  
  return(merged_row)
}

magma_gene_loci <- data.frame()

# 5) Iterate through each row of magma_gene and apply function
for (i in 1:nrow(magma_gene)) {
  result_row <- process_row(magma_gene[i, ], loci_fuma)
  
  # Store only non-NULL results
  if (!is.null(result_row)) {
    magma_gene_loci <- rbind(magma_gene_loci, result_row)
  }
}

magma_gene_loci <- magma_gene_loci %>%
  arrange(SYMBOL)

library(tidyr)
combined_data_snp <- magma_gene_loci %>%
  separate_rows(topsnp, GenomicLocusstart, GenomicLocusend, diseasepair, loci_num, sep = ",") %>%
  distinct()
unique_locustopsnp_genemap <- unique(combined_data_snp$topsnp)
unique_pairs_genemap <- unique(combined_data_snp$diseasepair)
unique_locus_genemap <- unique(combined_data_snp$loci_num)

for (i in c(1:nrow(combined_data_snp))) {
  if (combined_data_snp$START[i] >= combined_data_snp$GenomicLocusstart[i]) {
    if (combined_data_snp$STOP[i] <= combined_data_snp$GenomicLocusend[i]) {
      combined_data_snp$location[i] <- "within"
    } else {
      combined_data_snp$location[i] <- "end_outside"
    }
  } else {
    if (combined_data_snp$STOP[i] <= combined_data_snp$GenomicLocusend[i]) {
      combined_data_snp$location[i] <- "start_outside"
    } else {
      combined_data_snp$location[i] <- "both_outside"
    }
  }
}

setwd("/Volumes/yjzhao/COMMO_new/magma/magma_loci/")
write.csv(magma_gene_loci, "magma_gene_loci.csv", row.names = FALSE)

magma_gene_loci<-fread("magma_gene_loci.csv")
magma_gene_loci$gene_position<-paste0(magma_gene_loci$START,"-",magma_gene_loci$STOP)
magma_gene_loci$locus_boundary<-paste0(magma_gene_loci$GenomicLocusstart,"-",magma_gene_loci$GenomicLocusend)

magma_gene_loci$GenomicLocusstart_split <- strsplit(as.character(magma_gene_loci$GenomicLocusstart), ",")
magma_gene_loci$GenomicLocusend_split <- strsplit(as.character(magma_gene_loci$GenomicLocusend), ",")

# create new column 'locus_boundary'，connect first and second value

magma_gene_loci$locus_boundary <- mapply(function(start, end) {
  paste(start, end, sep = "-", collapse = ",")
}, magma_gene_loci$GenomicLocusstart_split, magma_gene_loci$GenomicLocusend_split)

magma_gene_loci<-as.data.frame(magma_gene_loci)
magma_gene_loci<-magma_gene_loci[,c(-20,-21)]

write.csv(magma_gene_loci, "magma_gene_loci.csv", row.names = FALSE)

unique_gene<-magma_gene_loci%>%
  select(SYMBOL,GENE, CHR, diseasepair,gene_position,loci_num,locus_boundary)%>%
  distinct()
unique_gene<-as.data.frame(unique_gene)
write.csv(unique_gene,file = "/Volumes/yjzhao/COMMO_new/magma/magma_loci/unique_gene_symbolNov15.csv",row.names = F)


unique_gene$Disease_pair_split <- strsplit(as.character(unique_gene$diseasepair), ",")

# create two dataframe
multi <- unique_gene[sapply(unique_gene$Disease_pair_split, length) >= 2, ]  # > 2 items
single <- unique_gene[sapply(unique_gene$Disease_pair_split, length) == 1, ]  # only 1 item

write.csv(multi[,c(-8)],file = "/Volumes/yjzhao/COMMO_new/magma/magma_loci/multi_diseasepair_gene_Nov16.csv",row.names = F)
write.csv(single[,c(-8)],file = "/Volumes/yjzhao/COMMO_new/magma/magma_loci/single_diseasepair_gene_Nov16.csv",row.names = F)

multi$Disease_pair_split <- NULL
single$Disease_pair_split <- NULL

combined_data_snp$pos<-loci_fuma$pos[match(combined_data_snp$topsnp,loci_fuma$SNP)]
combined_data_snp<-combined_data_snp%>%
  select(SYMBOL,CHR,GenomicLocusstart,GenomicLocusend,topsnp,diseasepair,pos)%>%
  distinct()

combined_data_snp_plot<-combined_data_snp%>%
  select(topsnp,CHR,pos,SYMBOL)%>%
  rename(snp=topsnp,chr = CHR,phenotype = SYMBOL)
write.table(combined_data_snp_plot, "combined_data_snp_plot.txt", quote=F,row.names = F,sep = "\t")


library(tidyr)
combined_data_snp_plot2 <- combined_data_snp %>%
  separate(diseasepair, into = c("brainD", "SomaticC"), sep = "-")
combined_data_snp_plot3 <- combined_data_snp_plot2 %>%
  pivot_longer(cols = c(brainD, SomaticC), 
               names_to = "Type", 
               values_to = "diseasepair") %>%
  select(-Type)%>%
  select(topsnp,CHR,pos,diseasepair)

combined_data_snp_plot3$ethnicity<-rep(c("brainD","SomaticC"),129)

combined_data_snp_plot3<-combined_data_snp_plot3%>%
  rename(position = pos,phenotype = diseasepair,snp = topsnp,chr = CHR)

write.table(combined_data_snp_plot3, "combined_data_snp_plot3.txt", quote=F,row.names = F,sep = "\t")

