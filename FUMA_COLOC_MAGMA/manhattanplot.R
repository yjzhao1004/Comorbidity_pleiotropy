#### 1. product Pplaco min for each SNP txt ####
setwd("/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults_recollect_new/")
rm(list = ls())
files<-dir(path = "/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults_recollect_new",
           full.names = T,
           pattern = ".txt")

df<-map(files,fread)

for (i in c(1:35)){
  df[[i]]<-df[[i]][,c("SNP","p.placo")]
}

df1 = list()
for (i in c(1:35)){
  df1[[i]]<-df[[i]][,c("SNP")]
}

df1<-reduce(df1,union)
df2<-reduce(df,union)

df2<-df2[order(df2$SNP,df2$p.placo,decreasing = F),]
df2<-df2[!duplicated(df2$SNP),]


write.csv(df3,"union_placosnp.csv",row.names = FALSE)

library(qqman)
library(RColorBrewer)

cp<-fread("/Volumes/ZHAOYUJIE/COMMO/chronicpain/chronicpain_gwas.txt")
df3<-merge(df2,cp[,c("SNP","CHR","BP")],by = "SNP")

rm(cp)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df1$SNP,df3$SNP)
otherdflist = list()
otherdflist$b<-otherdf
otherdflist<-as.data.frame(otherdflist)
otherdflist$a<-"dbsnp"
write.table(otherdflist[,c(2,1)],"setdiff_placosnp.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

otherdflist2<-fread("/Users/yujiezhao/Desktop/COMMO/PLACO/placoresults_recollect_new/gen_coords_557cf4cf.txt")
names(otherdflist2)[1]<-"SNP"
names(otherdflist2)[3]<-"CHR"
names(otherdflist2)[4]<-"BP"
otherdflist2<-otherdflist2[,c(1,3,4)]
otherdflist2<-merge(otherdflist2,df2,by = "SNP")
df3<-rbind(df3,otherdflist2)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df2$SNP,df3$SNP)
otherdflist = list()
otherdflist$b<-otherdf
otherdflist<-as.data.frame(otherdflist)
otherdflist$a<-"dbsnp"
write.table(otherdflist[,c(2,1)],"setdiff_placosnp2.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)

names(otherdflist)[1]<-"SNP"
cp<-fread("/Volumes/ZHAOYUJIE/COMMO/ADHD/ADHD_daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta")
otherdf3<-merge(otherdflist,cp[,c("SNP","CHR","BP")],by = "SNP")
otherdf3<-otherdf3[,c("SNP","CHR","BP")]
otherdflist3<-merge(otherdf3,df2,by = "SNP")
df3<-rbind(df3,otherdflist3)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df2$SNP,df3$SNP)
otherdflist = list()
otherdflist$SNP<-otherdf
otherdflist<-as.data.frame(otherdflist)

cp<-fread("/Volumes/ZHAOYUJIE/COMMO/PTSD/pts_eur_freeze2_overall.results")
otherdf3<-merge(otherdflist,cp[,c("SNP","CHR","BP")],by = "SNP")
otherdflist3<-merge(otherdf3,df2,by = "SNP")
df3<-rbind(df3,otherdflist3)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df2$SNP,df3$SNP)
otherdflist = list()
otherdflist$SNP<-otherdf
otherdflist<-as.data.frame(otherdflist)

cp<-fread("/Volumes/ZHAOYUJIE/COMMO/MDD/MDD_daner_pgc_mdd_meta_w2_no23andMe_rmUKBB_N.txt")
otherdf3<-merge(otherdflist,cp[,c("SNP","CHR","BP")],by = "SNP")
otherdflist3<-merge(otherdf3,df2,by = "SNP")
df3<-rbind(df3,otherdflist3)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df2$SNP,df3$SNP)
otherdflist = list()
otherdflist$SNP<-otherdf
otherdflist<-as.data.frame(otherdflist)

cp<-fread("/Volumes/ZHAOYUJIE/COMMO/T2D/t2d_gwas.csv")
otherdf3<-merge(otherdflist,cp[,c("SNP","CHR","BP")],by = "SNP")
otherdflist3<-merge(otherdf3,df2,by = "SNP")
df3<-rbind(df3,otherdflist3)
df3<-df3[order(df3$SNP),]
df3<-df3[!duplicated(df3$SNP),]
otherdf<-setdiff(df2$SNP,df3$SNP)
otherdflist = list()
otherdflist$SNP<-otherdf
otherdflist<-as.data.frame(otherdflist)


write.table(df3[,c(1,3,4,2)],"/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/plotplacosnp.txt",quote = FALSE,row.names = FALSE,col.names = TRUE)


#### 2. qqman package plot ####
library(data.table)
library(tidyverse)
library(dplyr)

df3<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/plotplacosnp.txt")
setwd("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/")

manhattan(df3, chr = "CHR",bp = "BP",snp = "SNP",p = "p.placo")
#annotation
loci_anno<-fread("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci_direction_20230717.csv")
manhattan(df3, chr = "CHR",bp = "BP",snp = "SNP",p = "p.placo",highlight = loci_anno$TopSNP)

#### 3. design by ggplot2 ####
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
df3<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/plotplacosnp.txt")
setwd("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/")
#annotation
loci_anno<-fread("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci_direction_20230717.csv")

# 1)计算chr长度
chr_len <- df3 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP))
# 2） 计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len)
#3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(df3, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)

#查看转化后的数据
head(Snp_pos[which(Snp_pos$CHR==2)],8)

p <- ggplot(Snp_pos, aes(x=BPcum, y=-log10(p.placo)))

X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)

p+
  #设置点的大小，透明度
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.2) +
  #设置颜色
  scale_color_manual(values = rep(c("#99CCFF", "#003399"), 22 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  geom_hline(yintercept = c(-log10(5*10^-8)), color = c('red'),linewidth = 0.5, linetype = c("dashed")) + 
  #设置主题
  theme_bw() +
  theme(legend.position="none",
    panel.border = element_blank(),
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())


#11 multitraits loci
library(ggrepel)

multi_traitsloci<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/11_multitraits_loci.csv")

#准备数据
rm(data)
data <- Snp_pos %>%
  # 添加高亮和注释信息：snpsOfInterest中的rs编号
  mutate( is_highlight=ifelse(SNP %in% loci_anno$TopSNP, "yes", "no"))

data$is_annotate<-"no"
uniquesnp<-multi_traitsloci[!duplicated(multi_traitsloci$TopSNP),c(1:2)]
data$is_annotate[which(data$SNP%in%uniquesnp$TopSNP)]<-uniquesnp$uniquenum
#unique(data$is_annotate)
highlightdata=subset(data, is_annotate!="no")


p <- ggplot(data, aes(x=BPcum, y=-log10(p.placo)))

X_axis <-  data %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)

p+
  #设置点的大小，透明度
  geom_point(data=subset(data, is_highlight=="no"), aes(color=as.factor(CHR)), alpha=0.8, size=1.2) +
  #设置颜色
  scale_color_manual(values = rep(c("#99CCFF", "#003399"), 22 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  geom_hline(yintercept = c(-log10(5*10^-8)), color = c('red'),linewidth = 0.5, linetype = c("dashed")) + 
  # 添加高亮点
  geom_point(data=subset(data, is_highlight=="yes"), color="red", size=2,shape = 17) +
  # 添加高亮label，且防止重叠
  geom_label_repel( data=subset(data, is_annotate!="no"), aes(label=SNP), size=2)+
  geom_point( data=highlightdata, color="orange", size=2.5,shape = 17)+
  
  #设置主题
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#### 4. prepare gene data ####
rm(list = ls())
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)
genedata<-fread("/Users/yujiezhao/Desktop/COMMO/MAGMA/prepare/magma_gene_20230709.csv")
locidata<-fread("/Users/yujiezhao/Desktop/COMMO/FUMA/pleioloci/pleio_fuma_results/all_fuma_pleioloci_direction_20230717.csv")

gene_topsnplocidata<-subset(locidata[,c("TopSNP","CHR","V1","BP","P.placo","Trait_pair","Locus_boundary")],locidata$Locus_boundary%in%genedata$Loci_position)

#write.csv(gene_topsnplocidata,file = "gene_topsnplocidata.csv")

gene_topsnplocidata<-read.csv("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/gene/gene_topsnplocidata.csv")

a<-match(genedata$Loci_position,gene_topsnplocidata$Locus_boundary)

genedata$TopSNP<-gene_topsnplocidata[a,1]
genedata$TopSNPBP<-gene_topsnplocidata[a,4]
genedata$TopSNPCHR<-gene_topsnplocidata[a,2]
genedata$TopSNPPplaco<-gene_topsnplocidata[a,5]

write.csv(genedata,file = "/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/gene/gene_topsnplocidata241gene.csv")

#plot
rm(list = ls())
df3<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/plotplacosnp.txt")
genedata<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/gene/gene_topsnplocidata241gene.csv")
# 1)计算chr长度
chr_len <- df3 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP))
# 2） 计算每条chr的初始位置
chr_pos <- chr_len  %>% 
  mutate(total = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len)
#3)计算累计SNP的位置
Snp_pos <- chr_pos %>%
  left_join(df3, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + total)
rm(chr_len,chr_pos)

unique(genedata$TopSNP)
unique(genedata$Loci_position)

#准备数据
rm(data)
data <- Snp_pos %>%
  # 添加高亮和注释信息：snpsOfInterest中的rs编号
  mutate( is_highlight=ifelse(SNP %in% genedata$TopSNP, "yes", "no"))

#multi-traits_genes
genes<-c("APOPT1","RP11-73M18.2","DCC","BAG5","PITPNM2","SBNO1","FUT2","FAM57B","TBX6","5-Mar","SLC9B1","UBE2D3","UQCR10","ALDOA","CDHR4","MST1R","RAI1","DNAJC11","RILPL2","AMIGO3","BSN","CDKAL1","NCOA5","LRRC48","ZZEF1","KLHL21","BTAF1","TCF4","RBM5","SGK223","MAPK3","MYO15A","MFHAS1","EMB","SSH2","WBP1L","FMNL3")

selectmultigenes<-subset(genedata,genedata$Gene_SYMBOL%in%genes)

#selectSNP<-unique(genedata$TopSNP[genedata$Gene_SYMBOL%in%genes])

data$is_annotate<-"no"
uniquesnp<-multi_traitsloci[!duplicated(multi_traitsloci$TopSNP),c(1:2)]

data$is_annotate[which(data$SNP%in%selectSNP)]<-selectSNP

#unique(data$is_annotate)
highlightdata=subset(data, is_annotate!="no")


p <- ggplot(data, aes(x=BPcum, y=-log10(p.placo)))

X_axis <-  data %>% group_by(CHR) %>% summarize(center=(max(BPcum)+min(BPcum))/2)

p+
  #设置点的大小，透明度
  geom_point(data=subset(data, is_highlight=="no"), aes(color=as.factor(CHR)), alpha=0.8, size=1.2) +
  #设置颜色
  scale_color_manual(values = rep(c("#99CCFF", "#003399"), 22 )) +
  #设定X轴
  scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
  #去除绘图区和X轴之间的gap
  scale_y_continuous(expand = c(0, 0) ) +  
  #添加阈值线
  geom_hline(yintercept = c(-log10(5*10^-8)), color = c('red'),linewidth = 0.5, linetype = c("dashed")) + 
  # 添加高亮点
  geom_point(data=subset(data, is_highlight=="yes"), color="red", size=2,shape = 17) +
  # 添加高亮label，且防止重叠
  geom_label_repel( data=subset(data, is_annotate!="no"), aes(label=SNP), size=2)+
  geom_point( data=highlightdata, color="orange", size=2.5,shape = 17)+
  
  #设置主题
  theme_bw() +
  theme(legend.position="none",
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#### 5. another plot for gene ####
rm(list = ls())
genedata<-fread("/Users/yujiezhao/Desktop/COMMO/Plot/loci&gene/gene/gene_topsnplocidata241gene.csv")
genedata<-genedata[,-1]
X_axis <-  genedata%>% group_by(trait_pairs)%>%summarize(center=(max(V1)+min(V1))/2)
X_axis<-X_axis[order(X_axis$center),]
X_axis$color<-factor(1:22)


XX_axis <-  genedata%>% group_by(trait_pairs)%>%summarize(center=max(V1))
XX_axis<-XX_axis[order(XX_axis$center),]
XX_axis$color<-factor(1:22)


genes<-c("APOPT1","RP11-73M18.2","DCC","BAG5","PITPNM2","SBNO1","FUT2","FAM57B","TBX6","5-Mar","SLC9B1","UBE2D3","UQCR10","ALDOA","CDHR4","MST1R","RAI1","DNAJC11","RILPL2","AMIGO3","BSN","CDKAL1","NCOA5","LRRC48","ZZEF1","KLHL21","BTAF1","TCF4","RBM5","SGK223","MAPK3","MYO15A","MFHAS1","EMB","SSH2","WBP1L","FMNL3")
selectmultigenes<-subset(genedata,genedata$Gene_SYMBOL%in%genes)

genedata$is_annotate<-"no"
genedata$is_annotate[which(genedata$Gene_SYMBOL%in%genes)]<-"yes"


#3.6,1
max(-log10(genedata$P_placo))
min(-log10(genedata$P_placo))

mydata <- data.frame(
  xstart = c(0,XX_axis$center[1:21]),
  xend   = XX_axis$center[1:22], 
  ystart = replicate(22,-0.1),
  yend   = replicate(22,-0.5),
  size   = replicate(22,3),
  classfy  = XX_axis$trait_pairs)

labeldata<-subset(genedata, is_annotate!="no")

ggplot(genedata, aes(x=V1, y=-log10(genedata$P_placo),size=-log10(genedata$P_placo))) +
  coord_cartesian(xlim = c(0,250),ylim = c(-1,15))+
  geom_hline(yintercept = -log10(0.05), color = 'black',size = 0.5,linetype="dashed") +
  geom_point(aes(size = -log10(genedata$P_placo), fill = factor(is_annotate)), alpha = 0.9,shape=21) +
  scale_fill_manual(values=c("yes" = "#CD3333", "no" = "#003399"))+
  scale_x_continuous(limits=c(0,241),label =X_axis$trait_pairs, breaks= X_axis$center) +
  scale_y_continuous(limits=c(-1,15),breaks = c(0,3,6,9,12,15)) +
  #geom_label_repel(data=labeldata, aes(label=labeldata$Gene_SYMBOL), size=2)+
  theme_bw() +
  theme(text=element_text(family="Arial"),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        #axis.line.x = element_line(),
        axis.text.x = element_text(color="black", size=10,angle=90,hjust=1),
        axis.text.y = element_text(color="black", size=10),
        axis.title.y=element_text(color="black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #plot.margin = unit(c(top_margin,right_margin,bottom_margin,left_margin),'lines')
  )+
  geom_rect(aes(xmin = mydata$xstart[1],xmax = mydata$xend[1] , ymin = mydata$ystart[1] , ymax = mydata$yend[1]) ,fill='#8ed4c8',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[2],xmax = mydata$xend[2] , ymin = mydata$ystart[2] , ymax = mydata$yend[2]) ,fill='#ffffb5',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[3],xmax = mydata$xend[3] , ymin = mydata$ystart[3] , ymax = mydata$yend[3]) ,fill='#bfbbdb',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[4],xmax = mydata$xend[4] , ymin = mydata$ystart[4] , ymax = mydata$yend[4]) ,fill='#f98577',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[5],xmax = mydata$xend[5] , ymin = mydata$ystart[5] , ymax = mydata$yend[5]) ,fill='#81b3d4',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[6],xmax = mydata$xend[6] , ymin = mydata$ystart[6] , ymax = mydata$yend[6]) ,fill='#fdb562',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[7],xmax = mydata$xend[7] , ymin = mydata$ystart[7] , ymax = mydata$yend[7]) ,fill='#fccee6',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[8],xmax = mydata$xend[8] , ymin = mydata$ystart[8] , ymax = mydata$yend[8]) ,fill='#dadada',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[9],xmax = mydata$xend[9] , ymin = mydata$ystart[9] , ymax = mydata$yend[9]) ,fill='#bd81be',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[10],xmax = mydata$xend[10] , ymin = mydata$ystart[10] , ymax = mydata$yend[10]) ,fill='#cdecc6',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[11],xmax = mydata$xend[11] , ymin = mydata$ystart[11] , ymax = mydata$yend[11]) ,fill='#8ed4c8',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[12],xmax = mydata$xend[12] , ymin = mydata$ystart[12] , ymax = mydata$yend[12]) ,fill='#ffffb5',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[13],xmax = mydata$xend[13] , ymin = mydata$ystart[13] , ymax = mydata$yend[13]) ,fill='#bfbbdb',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[14],xmax = mydata$xend[14] , ymin = mydata$ystart[14] , ymax = mydata$yend[14]) ,fill='#f98577',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[15],xmax = mydata$xend[15] , ymin = mydata$ystart[15] , ymax = mydata$yend[15]) ,fill='#81b3d4',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[16],xmax = mydata$xend[16] , ymin = mydata$ystart[16] , ymax = mydata$yend[16]) ,fill='#fdb562',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[17],xmax = mydata$xend[17] , ymin = mydata$ystart[17] , ymax = mydata$yend[17]) ,fill='#fccee6',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[18],xmax = mydata$xend[18] , ymin = mydata$ystart[18] , ymax = mydata$yend[18]) ,fill='#dadada',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[19],xmax = mydata$xend[19] , ymin = mydata$ystart[19] , ymax = mydata$yend[19]) ,fill='#bd81be',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[20],xmax = mydata$xend[20] , ymin = mydata$ystart[20] , ymax = mydata$yend[20]) ,fill='#cdecc6',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[21],xmax = mydata$xend[21] , ymin = mydata$ystart[21] , ymax = mydata$yend[21]) ,fill='#8ed4c8',alpha = 0.8)+
  geom_rect(aes(xmin = mydata$xstart[22],xmax = mydata$xend[22] , ymin = mydata$ystart[22] , ymax = mydata$yend[22]) ,fill='#ffffb5',alpha = 0.8)
  





