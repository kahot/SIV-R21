library(ggplot2)
library(reshape2)
library(gridExtra)
library(DataCombine)
library(colorspace)
library(dplyr)
library(plotrix)
library(cowplot)
source("Rscripts/Pcorrection.R")
source("Rscripts/label_scientific.R")


colors<-qualitative_hcl(6, palette="Dark3")

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum21<-summary[1:69,]
Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")


### CD4/CD8 and diversity in Lung gran vs. non-gran
rna<-Sum21[,c("filename","Cohort","Monkey","Tissue3","Week","mean","SIV.RNA.per.tissue",
               "CD4.percent","CD8.percent", "Granuloma")]

#CD4
cd4<-rna[!is.na(rna$CD4.percent),] 
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




cd4$ratio<-cd4$CD4.percent/cd4$CD8.percent
cor.test(cd4$ratio, cd4$mean, method="spearman")
#p-value = 0.001811
#rho= 0.7546759 

ggplot(cd4, aes(x=ratio, y=mean*100))+
    geom_point(size=3.5, alpha=0.8, color=colors[c(1)])+ylab("% Average diversity")+
    theme_bw()+xlab("CD4/CD8 ratio")+
    annotate("text", x=1.4, y=0.16, label="rho=-0.75**", hjust=0)
ggsave("Output/Figures/Cd4.Cd8.ratio.pdf",width = 4, height=3.8)





### granzyme B analysis  ###
cd8<-read.csv("Data/granzyme_cd8.csv")
cd8<-merge(cd8, cd4[,c("filename","mean")], by="filename")

cor.test(cd8$mean, cd8$GranzymeB, method="spearman")
#data:  cd8$mean and cd8$GranzymeB
#S = 464.63, p-value = 0.02992
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#    rho 
#-0.6245652 

Plots[[3]]<-ggplot(cd8, aes(x=GranzymeB*100, y=mean*100))+
    geom_point(size=2.5, alpha=0.8, color=colors[c(1)])+ylab("_")+
    theme_bw()+xlab("% CD8 expressing grandzyme B")+
    annotate("text", x=57, y=0.33, label="rho=-0.62*", hjust=0)+
    theme(axis.title.y = element_text(color="white"))



png("Output/Figures/CD4.CD8.Lung.only.png", width = 12, height = 3.8, res=300, unit="in")
ggdraw()+
    draw_plot(Plots[[1]],x=0,y=0,width=0.33,height=1)+
    draw_plot(Plots[[2]],x=0.33,y=0,width=0.33,height=1)+
    draw_plot(Plots[[3]],x=0.66,y=0,width=0.33,height=1)+
    draw_plot_label(c("A", "B","C"), c(0, 0.33,0.66), c(1, 1,1), size = 15)
dev.off()





