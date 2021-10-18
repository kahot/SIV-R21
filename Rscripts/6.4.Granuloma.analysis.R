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
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]
summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)

Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]
Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")

### Granuloma vs. non-granuloma diversity comparison ####
### 1. Gran vs. No-gran  
gran<-Sum21[,c("Cohort","Monkey","Tissue2","Tissue3","Week","mean","Granuloma","SIV.RNA.per.granuloma", "SIV.RNA.per.tissue")]
gran<-gran[!is.na(gran$Granuloma),] #50 samples

gran$Granuloma[gran$Granuloma=="Y"]<-"Granuloma"
gran$Granuloma[gran$Granuloma=="N"]<-"Non-granuloma"


by.gran<- aggregate(gran["mean"],by=list(gran$Granuloma),mean,na.rm=T )

p1<-ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100), color="lightslategray",size=2, position=position_jitter(width=0.03), alpha=0.5)+
    geom_point(data=by.gran, aes(x=Group.1, y=mean*100), color="#023e8a",size=2)+
    xlab('')+ylab('% Average diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())

wilcox.test(gran$mean[gran$Granuloma=="Granuloma"], gran$mean[gran$Granuloma=="Non-granuloma"], alternative="two.sided")
#W = 332, p-value = 0.6489

# 2. By cohort
by.granC<- aggregate(gran["mean"],by=list(gran$Granuloma, gran$Cohort),mean,na.rm=T )

gran$Cohort<-factor(gran$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.granC$Group.2<-factor(by.granC$Group.2, levels=c("SIV only", "Mtb NR", "Mtb R"))

wilcox.test(gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
#W = 10, p-value = 0.0257

#non-granuloma comparison
gra.res<-data.frame(test=c("nongran_R.vs.nongran_NR","nongran_R.vs.nongran_SIV","nongran_SIV.vs.nongran_NR") )
re1<-wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
re2<-wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="SIV only"], alternative="two.sided")
re3<-wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="SIV only"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")

for (i in 1:3){
    wilc<-get(paste0("re",i))
    gra.res$rawP[i]<-wilc[[3]]
}
gra.res<-Pcorrection(gra.res)
#non are significantly different

p2<-ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100, color=Cohort),size=2, position=position_dodge(width=0.6))+
    geom_point(data=by.granC, aes(x=Group.1, y=mean*100, fill=Group.2), color="#023e8a",size=2, position=position_dodge(width=0.6))+
    xlab('')+ylab('% Average diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    scale_color_manual(values = paste0(colors[c(5,3,1)],99))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    annotate("segment", x = 0.85, xend = 0.85, y = 0.35, yend = 0.4, colour = "gray60", size=0.3) +
    annotate("segment", x = 1.15, xend = 1.15, y = 0.31, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.85, xend = 1.15, y = 0.4, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.41, colour = "gray60", label="*", size=5)

#Test gran vs. non-gran within each cohort
cohorts<-c("Mtb NR", "Mtb R")
gran.state<-c("Granuloma","Non-granuloma")
tissues<-c("Peripheral LN","Thoracic LN","Lung")

test2<-data.frame(test=c("R","NR"))
for (i in 1: length(cohorts)){
    r1<-wilcox.test(gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort==cohorts[i]], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort==cohorts[i]])
    test2$rawP[i]<-r1[[3]]
}
test2<-Pcorrection(test2)
#none are significant


#########
# 2. Diversity by tissue
G<-gran[gran$Granuloma=="Granuloma",]
N<-gran[gran$Granuloma=="Non-granuloma",]

#test in non-granuloma
comb<-combn(tissues,2)
comb<-t(comb)
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

test3<-data.frame(test=co_pairs)
for (i in 1:length(co_pairs)){
    re<-wilcox.test(N$mean[N$Tissue3==comb[i,1]],N$mean[N$Tissue3==comb[i,2]])
    test3$rawP[i]<-re[[3]]
}
test3<-Pcorrection(test3)
#non are significantly different

#Test in granulomas
wilcox.test(G$mean[G$Tissue3=="Lung"],G$mean[G$Tissue3=="Thoracic LN"], alternative="two.sided")
#W = 30, p-value = 0.0817



# Test thoracic LN/Lung granuloma between R and NR
gran.tissue<-data.frame(test=c("NR-tLN.vs.R-tLN","NR-tLN.vs.R-Lung","R-tLN.vs.R-Lung"))
r1<-wilcox.test(G$mean[G$Tissue3=="Thoracic LN"&G$Cohort=="Mtb NR"],G$mean[G$Tissue3=="Thoracic LN"&G$Cohort=="Mtb R"], alternative="two.sided")
r2<-wilcox.test(G$mean[G$Tissue3=="Thoracic LN"&G$Cohort=="Mtb NR"],G$mean[G$Tissue3=="Lung"&G$Cohort=="Mtb R"], alternative="two.sided")
r3<-wilcox.test(G$mean[G$Tissue3=="Thoracic LN"&G$Cohort=="Mtb R"], G$mean[G$Tissue3=="Lung"&G$Cohort=="Mtb R"], alternative="two.sided")
gran.tissue$rawP[1]<-r1[[3]]
gran.tissue$rawP[2]<-r2[[3]]
gran.tissue$rawP[3]<-r3[[3]]

gran.tissue<-Pcorrection(gran.tissue)
gran.tissue[,c(1:3,5)]
#              test      rawP Bonferroni      Holm
#1  NR-tLN.vs.R-tLN 0.2000000  0.6000000 0.4000000
#2 NR-tLN.vs.R-Lung 0.0248366  0.0745098 0.0745098
#3  R-tLN.vs.R-Lung 0.6450980  1.0000000 0.6450980

#plot
gran2<-gran[gran$Cohort!="SIV only",]
by.granA<- aggregate(gran2["mean"],by=list(gran2$Granuloma, gran2$Tissue3, gran2$Cohort),mean,na.rm=T )
by.granA$Group.2<-factor(by.granA$Group.2, levels=c("Peripheral LN","Thoracic LN","Lung"))
gran2$Tissue3<-factor(gran2$Tissue3, levels=c("Peripheral LN","Thoracic LN","Lung"))
gran2$Cohort<-factor(gran2$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.granA$Group.3<-factor(by.granA$Group.3, levels=c("SIV only", "Mtb NR", "Mtb R"))

p3<-ggplot()+
    geom_point(data=gran2, aes(x=Granuloma, y=mean*100, color=Cohort, shape=Tissue3),size=2.8, 
               position=position_dodge(width=0.6))+
    geom_point(data=by.granA, aes(x=Group.1, y=mean*100, color=Group.3, fill=Group.3, shape=Group.2),
               position=position_dodge(width=0.6), color="#023e8a",size=1.8)+
    xlab('')+ylab('% Average diversity')+
    scale_color_manual(values = paste0(colors[c(3,1)],99))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    scale_shape_manual(values=c(16,17,15),labels=c("Peripheral LN", "Thoracic LN", "Lung"))+
    guides(color = guide_legend(title = NULL))+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())


pdf("Output/Figures/FigureS5_Granuloma.pdf",width=9, height=6)
png("Output/Figures/FigureS5_Granuloma.png",width=9, height=6, units="in",res = 300)

ggdraw()+
    draw_plot(p1,0,0.5,0.36,0.5)+
    draw_plot(p2,0.36,0.5,0.58,0.5)+
    draw_plot(p3,0,0,0.75,0.5)+
    draw_plot_label(c("A", "B", "C"), c(0, 0.36, 0), c(1, 1, 0.5), size = 15)
dev.off()    


######################################################
## Do graunulomas within an animal differ in diversity?

#remove SIV only animals since no granuloma 
gran3<-Sum21[!is.na(Sum21$Granuloma)&Sum21$Granuloma=="Y",]
monkeys=unique(gran3$Monkey)
#Remove Monkeys with only 1 granuloma samples
monkeys<-monkeys[!(monkeys=="3116"|monkeys=="3216")]

gran3<-merge(gran3,samples[,c("filename","Sample")], by="filename")
Mean<-gran3[,c("Monkey","Cohort","Tissue3","mean","se","Sample","filename")]
colnames(Mean)[colnames(Mean)=="Sample"]<-"Location"
Mean$Location<-gsub("\\*","",Mean$Location)
Mean$Location<-gsub(" gran","",Mean$Location)
Mean$Location<-gsub("w/","",Mean$Location)
Mean$Location<-trimws(Mean$Location)

#Mean$Label<-paste(Mean$Tissue,Mean$Location)

library(tidyverse)
max_length<- max(str_length(Mean$Location))
Mean <- mutate(Mean, Location = str_pad(Location, max_length, side = "left"))


P<-list()
for (i in 1:length(monkeys)){
    df<-Mean[Mean$Monkey==monkeys[i],]
    df<-df[order(df$Tissue3),]
    df$Location<-factor(df$Location,levels=paste(df$Location))
    df$Tissue3<-factor(df$Tissue3, levels=c("Thoracic LN","Lung"))
    p<-ggplot(df, aes(x=Location, y=mean*100, shape=Tissue3))+
            geom_errorbar(mapping=aes(ymin=(mean-se)*100, ymax=(mean+se)*100), width=.18, color="gray50")+
            geom_point(size=3)+ylim(0.1,0.45)+
            ggtitle(df$Monkey[1])+ylab('')+
            scale_shape_manual(values=c(17,16))+
            theme(axis.text.x = element_text(angle=90, hjust=1,family = "Monaco"),axis.title.x = element_blank())+
            theme(legend.position = "none")
    if (i==1) P[[i]]<-p+ylab("% Average diversity")
    else P[[i]]<-p
    
 }

get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
l<-ggplot(df, aes(x=Location, y=mean, shape=Tissue3))+
    geom_errorbar(mapping=aes(ymin=mean-se, ymax=mean+se), width=.18, color="gray50")+
    geom_point(size=3)+ylim(0.001,0.005)+
    ggtitle(df$Monkey[1])+
    scale_shape_manual(values=c(17,16))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90, hjust=1, family="Menlo"),axis.title.x = element_blank())+
    theme(legend.title =element_blank())
    
legend1<-get_legend(l)


P[[6]]<-legend1


png("Output/Figures/Granuloma_diversity_withinAnimal1.png", width=13, height =4.5, units="in",res=300)
do.call(grid.arrange, c(P, ncol=6))
dev.off()

#Granuloma box plot
AllFreq<-read.csv("Output/Diversity/Allfile_frequency.csv")
AllFreq<-AllFreq[,-1]
gsamples<-unique(Mean$filename)
AllFreq<-AllFreq[AllFreq$filename %in% gsamples,]
AllFreq<-merge(AllFreq,Mean[,c("filename","Location")], by="filename")

P2<-list()
for (i in 1:length(monkeys)){
    fnames<-Mean$filename[Mean$Monkey==monkeys[i]]
    df<-AllFreq[AllFreq$filename %in% fnames,]
    df<-df[order(df$Tissue3),]
    df$Location<-factor(df$Location,levels=paste(unique(df$Location)))
    df$Tissue3<-factor(df$Tissue3, levels=c("Thoracic LN","Lung"))
    p<-ggplot()+
        geom_boxplot(data=df, aes(x=Location, y=freq.mutations*100, color=Tissue3),
                    alpha=0.2)+
        scale_y_continuous(trans = 'log10',breaks = c(0.01,1,50,100), labels=c(0.01,1,50,100), limits = c(0.0002,140),expand = c(0,0))+
        xlab('')+ylab('')+
        theme_bw()+theme(legend.position = "none")+
        theme(panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle=90, hjust=1,family = "Monaco"))
    if (i==1) P2[[i]]<-p+ylab("% Diversity")
    else P2[[i]]<-p
    
}

l<-ggplot()+
    geom_boxplot(data=df, aes(x=Location, y=freq.mutations*100, color=Tissue3),
                 alpha=0.2)+
    scale_y_continuous(trans = 'log10',breaks = c(0.01,1,50,100), labels=c(0.01,1,50,100), limits = c(0.0002,140),expand = c(0,0))+
    xlab('')+ylab('')+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank(),axis.text.x = element_text(angle=90, hjust=1,family = "Monaco"))

    
    
legend1<-get_legend(l)


P2[[6]]<-legend1

png("Output/Figures/Granuloma_diversity_withinAnimal_boxplot.png", width=13, height =4.5, units="in",res=300)
do.call(grid.arrange, c(P2, ncol=6))
dev.off()




# Age of granuloma related to diversity?
age<-read.csv("Data/SIVR21_granulomaAge..csv")
Ag<-Sum21[Sum21$filename %in% age$filename,]
grandiv<-age[,c("filename","gran.age.weeks")]
for (i in 1:nrow(age)){
    grandiv$diversity[i]<-Ag$mean[Ag$filename==age$filename[i]]
    grandiv$Tissue[i]<-Ag$Tissue2[Ag$filename==age$filename[i]]
    grandiv$Cohort[i]<- Ag$Cohort[Ag$filename==age$filename[i]]
}
grandiv$gran.age.weeks<-as.numeric(grandiv$gran.age.weeks)
colnames(grandiv)[2]<-"age"

cor.test(grandiv$age,grandiv$diversity, method="spearman")
#p-value = 0.3708 rho =-0.2592543 

ggplot(grandiv, aes(x=age, y=diversity*100))+
    geom_point(color="steelblue", size=2)+xlab("Age of granuloma (weeks)")+ylab("Average diversity %")+
    annotate("text", x=15,y=0.28, label= "P=0.37", color="darkblue")+
    theme_bw()
ggsave("Output/Figures/Granuloma.vs.age.pdf", width = 4.5, height = 4)

