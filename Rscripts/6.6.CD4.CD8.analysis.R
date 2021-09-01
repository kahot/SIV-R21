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
#col_light<-qualitative_hcl(6, palette="Set3")
#hcl_palettes(plot = TRUE)#
#MFcolors<-c("#EB4E61","#FF9300","#9437FF")

SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]
samples$Tissue2[samples$Tissue2=="plasma"]<-"Plasma"
samples$Tissue2[samples$Tissue2=="Plasma"&samples$Week>5]<-"Plasma PM"
samples$Tissue2[samples$Tissue2=="Thoracic LN"]<-"tLN"

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]
Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")


#Remove SIV only 
samples<-samples[samples$Cohort!="SIV only",]


### CD4/8 and diversity copy numbers in tLN gran vs. non-gran
rna<-Sum21[,c("Cohort","Monkey","Tissue2","Week","mean","SIV.RNA.per.tissue",
               "CD4.percent","CD8.percent", "Granuloma")]
#Eliminate SIV only 
rna<-rna[rna$Cohort!="SIV only",]

#CD4
cd4<-rna[!is.na(rna$CD4.percent),] #26 samples
table(cd4$Cohort,cd4$Granuloma)
#        N  Y
#Mtb NR  2  4
#Mtb R   3 17

#select granuloma only
cd4<-cd4[cd4$Granuloma=="Y",] #21 samples
table(cd4$Cohort,cd4$Tissue2)
#        Lung Thoracic LN
#Mtb NR    0           4
#Mtb R    14           3


cor.test(cd4$CD4.percent, cd4$mean, method = "spearman") #P=0.03464
#rho=0.462788

colnames(cd4)[which(colnames(cd4)=="Tissue2")]<-"Tissue"
cd4$Cohort<-factor(cd4$Cohort,levels=c("Mtb NR", "Mtb R"))
ggplot(cd4, aes(x=CD4.percent*100, y=mean*100, color=Cohort, shape=Tissue))+
    geom_point(size=2.5, alpha=0.7)+ylab("% Average diversity")+
    theme_bw()+xlab("% CD4")+
    scale_color_manual(values=colors[c(3,1)])+
    annotate("text", x=45, y=0.33, label="rho=0.46*")
ggsave("Output/Figures/CD4.pdf", height = 4,width = 5)

#Lung/MtbR only
cd4.2<-cd4[cd4$Cohort=="Mtb R",]
cor.test(cd4.2$CD4.percent, cd4.2$mean, method = "spearman") #p-value = 0.06231 rho=0.46135

cd4.3<-cd4.2[cd4.2$Tissue=="Lung",]
cor.test(cd4.3$CD4.percent, cd4.3$mean, method = "spearman") #p-value = 0.04129 rho=0.5506621

Plots<-list()
Plots[[1]]<-ggplot(cd4.3, aes(x=CD4.percent*100, y=mean*100))+
    geom_point(size=2.5, alpha=0.8, color=colors[c(1)])+ylab("% Average diversity")+
    theme_bw()+xlab("% CD4")+
    annotate("text", x=42, y=0.33, label="rho=0.55*")
#ggsave("Output/Figures/CD4_MtbR_Lung.pdf", height = 3.8,width = 4)


#####
#CD8
cor.test(cd4$CD8.percent, cd4$mean, method = "spearman") #p-value = 0.00563
#rho=-0.5821146 

#cd4$Cohort<-factor(cd8$Cohort,levels=c("Mtb NR", "Mtb R"))
ggplot(cd4, aes(x=CD8.percent*100, y=mean*100, color=Cohort, shape=Tissue))+
    geom_point(size=2.5, alpha=0.7)+ylab("% Average diversity")+
    theme_bw()+xlab("% CD8")+
    scale_color_manual(values=colors[c(3,1)])+
    annotate("text", x=60, y=0.33, label="rho=-0.58***", hjust=0)
ggsave("Output/Figures/CD8.pdf", height = 4,width = 5)

#Lung /Mtb Ronly
cor.test(cd4.3$CD8.percent, cd4.3$mean, method = "spearman") #p-value = 0.004201 rho=-0.7130312

Plots[[2]]<-ggplot(cd4.3, aes(x=CD8.percent*100, y=mean*100))+
    geom_point(size=2.5, alpha=0.8, color=colors[c(1)])+ylab("_")+
    theme_bw()+xlab("% CD8")+
    annotate("text", x=60, y=0.33, label="rho=-0.71**", hjust=0)+
    theme(axis.title.y = element_text(color="white"))
#ggsave("Output/Figures/CD8_MtbR_Lung.pdf", height = 3.8,width = 4)


pdf("Output/Figures/CD4.CD8.MtbR.Lung.only.pdf", width = 8, height = 3.8)
do.call(grid.arrange, c(Plots, ncol=2))
dev.off()

cd4.3$ratio<-cd4.3$CD4.percent/cd4.3$CD8.percent
cor.test(cd4.3$ratio, cd4.3$mean, method="spearman")
#p-value = 0.001811
#rho= 0.7546759 
ggplot(cd4.3, aes(x=ratio, y=mean*100))+
    geom_point(size=3.5, alpha=0.8, color=colors[c(1)])+ylab("% Average diversity")+
    theme_bw()+xlab("CD4/CD8 ratio")+
    annotate("text", x=1.4, y=0.16, label="rho=-0.75**", hjust=0)
    
ggsave("Output/Figures/Cd4.Cd8.ratio.pdf",width = 4, height=3.8)
