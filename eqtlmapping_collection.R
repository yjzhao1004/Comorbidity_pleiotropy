#eqtl mapping
rm(list = ls())
pleioloci<-read.csv("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci20230619.csv")

trait_unique_pair<-unique(pleioloci$Trait_pair)


setwd("/Users/yujiezhao/Desktop/COMMO/eQTLmapping/")

library(data.table)
eqtl_allpairs_collection= list()
for (j in c(1:length(trait_unique_pair))){
  print(paste(as.character(j),trait_unique_pair[j],sep = ":"))
  currentpleio<-pleioloci[pleioloci$Trait_pair==trait_unique_pair[j],]
  traits<-strsplit(trait_unique_pair[j],"_")
  eqtlfilename1<-paste("eqtl_",traits[[1]][1],".txt",sep = "")
  eqtlfile1<-fread(eqtlfilename1)
  eqtlfile1<-subset(eqtlfile1,db %in%c("GTEx/v8","BRAINEAC"))
  
  eqtlfilename2<-paste("eqtl_",traits[[1]][2],".txt",sep = "")
  eqtlfile2<-fread(eqtlfilename2)
  eqtlfile2<-subset(eqtlfile2, db%in%c("GTEx/v8","BRAINEAC"))
  
  eqtl_collection = list()
  for (i in c(1:length(currentpleio$X))){
    
    eqtlfile1_current<-subset(eqtlfile1,pos==currentpleio$BP[i])
    eqtlfile2_current<-subset(eqtlfile2,pos==currentpleio$BP[i])
    
    eqtlfile1_current<-subset(eqtlfile1_current,symbol==currentpleio$Nearest_Gene[i])
    eqtlfile2_current<-subset(eqtlfile2_current,symbol==currentpleio$Nearest_Gene[i])
    
    eqtlfile1_current$orgnolocus<-currentpleio$V1[i]
    eqtlfile1_current$snp<-currentpleio$TopSNP[i]
    eqtlfile1_current$trait_pairs<-currentpleio$Trait_pair[i]
    
    eqtl_collection<-rbind(eqtl_collection,eqtlfile1_current)
    rm(eqtlfile1_current,eqtlfile2_current)
  }
  
  eqtl_allpairs_collection<-rbind(eqtl_allpairs_collection,eqtl_collection)
  
  rm(eqtl_collection)
}

eqtl_allpairs_collection<-subset(eqtl_allpairs_collection,db == "GTEx/v8")
eqtl_allpairs_collection<-subset(eqtl_allpairs_collection,FDR < 0.05)

write.csv(eqtl_allpairs_collection,file = "eqtl_allpairs_collection_20230717.csv")

uniquegene_pleio<-unique(pleioloci$Nearest_Gene)

uniquegene_eqtl<-unique(eqtl_allpairs_collection$symbol)


uniquesnp_eqtl<-unique(eqtl_allpairs_collection$snp)



