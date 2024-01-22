rm(list = ls())
ldsc<-read.csv("/Users/yujiezhao/Desktop/COMMO/LDSC/ldsc_results_20231031.csv")
library("tidyverse")
library(ggplot2)

names(ldsc)[3]<-"r"
names(ldsc)[1]<-"CP"
names(ldsc)[2]<-"BD"
names(ldsc)[4]<-"pvalue"


ldsc$BD<-factor(ldsc$BD,levels=c("MS","ALS","PD","AD","TS","SCZ","PTSD","OCD","MDD","EAT","BP","ASD","ANX","ALD","ADHD"))



ldsc$sig_label2<-NA
ldsc$sig_label2[ldsc$FDR<0.05] = 1
ldsc$sig_label2[ldsc$FDR>=0.05] = 2
#ldsc$sig_label2[ldsc$FDR>=0.05&ldsc$pvalue<0.05] = 3

ldsc$sig_label2 <- factor(ldsc$sig_label2,levels=c(1,2),labels=c("Sig_fdr","Non-sig_fdr"))

ldsc$logfdr<--1*log(ldsc$FDR)
ldsc$logfdrs<-scale(ldsc$logfdr)
ldsc$logfdrs<-scale(ldsc$logfdr)+3

ldsc$logfdrs[1]<-9
ldsc$logfdrs[2]<-7

top_margin<-1
right_margin<-1
bottom_margin<-1
left_margin<-1
ggplot(ldsc,aes(x=CP, y=BD))+
  geom_point(aes(size = logfdrs,fill=r,colour = sig_label2,shape=sig_label2)) +
  scale_size_area(max_size = 9) +
  scale_shape_manual(values = c("Sig_fdr" = 21,"Non-sig_fdr" = 21))+
  scale_colour_manual(values = c("Sig_fdr" = "black","Non-sig_fdr" = "white"))+
  scale_fill_gradient2(low='#0000FF',mid = "white",high='#8B0000',midpoint = 0,limits=c(min(ldsc$r),max(ldsc$r)))+
  theme_bw() +
  theme(axis.text.x = element_text(color="black", size=10,angle=0,hjust=0.5),
        axis.text.y = element_text(color="black", size=10),
        axis.title.y=element_text(color="black", size=10),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.margin = unit(c(top_margin,right_margin,bottom_margin,left_margin),'lines')
  )+
  labs(y = 'Brain Disorders', x = 'Chronic Physical Conditions')

  
  
