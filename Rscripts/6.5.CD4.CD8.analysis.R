library(ggplot2)
library(reshape2)
library(gridExtra)
library(DataCombine)
library(colorspace)
library(dplyr)
library(plotrix)
source("Rscripts/Pcorrection.R")
source("Rscripts/label_scientific.R")


colors<-qualitative_hcl(6, palette="Dark3")

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]
Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")


#Remove SIV only 
samples<-samples[samples$Cohort!="SIV only",]


### CD4/CD8 and diversity in Lung gran vs. non-gran
rna<-Sum21[,c("Cohort","Monkey","Tissue3","Week","mean","SIV.RNA.per.tissue",
               "CD4.percent","CD8.percent", "Granuloma")]

#CD4
cd4<-rna[!is.na(rna$CD4.percent),] #26 samples
colnames(cd4)[colnames(cd4)=="Tissue3"]<-"Tissue"

#Lung/MtbR only
cd4<-cd4[cd4$Cohort=="Mtb R"&cd4$Tissue=="Lung",]
cor.test(cd4$CD4.percent, cd4$mean, method = "spearman") 
#p-value = 0.04129 rho=0.5506621

Plots<-list()
Plots[[1]]<-ggplot(cd4, aes(x=CD4.percent*100, y=mean*100))+
    geom_point(size=2.5, alpha=0.8, color=colors[c(1)])+ylab("% Average diversity")+
    theme_bw()+xlab("% CD4")+
    annotate("text", x=42, y=0.33, label="rho=0.55*")


#####
#CD8
cor.test(cd4$CD8.percent, cd4$mean, method = "spearman") 
#p-value = 0.004201 rho=-0.7130312

Plots[[2]]<-ggplot(cd4, aes(x=CD8.percent*100, y=mean*100))+
    geom_point(size=2.5, alpha=0.8, color=colors[c(1)])+ylab("_")+
    theme_bw()+xlab("% CD8")+
    annotate("text", x=60, y=0.33, label="rho=-0.71**", hjust=0)+
    theme(axis.title.y = element_text(color="white"))


pdf("Output/Figures/CD4.CD8.MtbR.Lung.only.pdf", width = 8, height = 3.8)
do.call(grid.arrange, c(Plots, ncol=2))
dev.off()

cd4$ratio<-cd4$CD4.percent/cd4$CD8.percent
cor.test(cd4$ratio, cd4$mean, method="spearman")
#p-value = 0.001811
#rho= 0.7546759 

ggplot(cd4, aes(x=ratio, y=mean*100))+
    geom_point(size=3.5, alpha=0.8, color=colors[c(1)])+ylab("% Average diversity")+
    theme_bw()+xlab("CD4/CD8 ratio")+
    annotate("text", x=1.4, y=0.16, label="rho=-0.75**", hjust=0)
ggsave("Output/Figures/Cd4.Cd8.ratio.pdf",width = 4, height=3.8)