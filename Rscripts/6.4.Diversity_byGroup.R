library(ggplot2)
library(reshape2)
library(gridExtra)
library(DataCombine)
library(colorspace)
library(dplyr)
library(FSA)
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

### 1. By Cohort (SIV only vs. Latent NR vs. Latent R)
#Wilcoxon test
coh<-c("SIV only","Mtb NR","Mtb R")

Co<-data.frame(cohort=coh)
for (j in 1:3){
    DF<-Sum21[Sum21$Cohort==coh[j],]
    Co$mean[j]<-mean(DF$mean)
    Co$se[j]<-mean(DF$se)
}

comb<-combn(coh,2)
comb<-t(comb)
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

wil.res<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    v1<-Sum21$mean[Sum21$Cohort==comb[i,1]]
    v2<-Sum21$mean[Sum21$Cohort==comb[i,2]]
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res$rawP[i]<-re1[[3]]
    wil.res$mean1[i]<-mean(v1, na.rm=T)
    wil.res$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
}

wil.res<-Pcorrection(wil.res)
write.csv(wil.res, paste0("Output/Diversity/Tests/Cohorts_div_wilcox_results.csv"))


#2. Diverstiy by Tissue types by Cohort
infec<-Sum21[,c("Cohort","Tissue3","Week","mean", "se")]
colnames(infec)[2]<-"Tissue"

infec$Cohort<-factor(infec$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.Coh<- aggregate(infec["mean"],by=list(infec$Cohort),mean,na.rm=T )
by.Coh$Group.1<-factor(by.Coh$Group.1, levels=c("SIV only", "Mtb NR", "Mtb R"))

#Plot cohort results
p1<-ggplot()+
    geom_point(data=infec, aes(x=Cohort, y=mean*100, color=Cohort),size=2.3, position=position_jitter(width=0.05))+
    geom_point(data=by.Coh, aes(x=Group.1, y=mean*100), color="black",size=2.3, shape=4)+
    xlab('')+ylab('% Average diversity')+
    scale_color_manual(values=paste0(colors[c(5,3,1)],"66"), guide='none')+
    theme(legend.title = element_blank())+
    scale_y_continuous(limits = c(0.1, 0.55))+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/Cohort_diversity.png", width = 3,height = 3,units="in", dpi=300)




#### Wilcoxon test on diversity between cohort by tissues ###

#comparison pairs 
comb<-combn(coh,2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

tissues<-unique(Sum21$Tissue2)

Results<-data.frame()
for (j in 1:length(tissues)){
    df<-Sum21[Sum21$Tissue2==tissues[j],]

    #wicoxon test
    wil.res<-data.frame(Test=comb_pairs)
    
    for (i in 1:nrow(comb)){
        v1<-df$mean[df$Cohort==comb[i,1]]
        v2<-df$mean[df$Cohort==comb[i,2]]
        if (length(v1)==0|length(v2)==0){
            wil.res$rawP[i]<-NA
            wil.res$mean1[i]<-NA
            wil.res$mean2[i]<-NA
            wil.res$higher[i]<-NA
        }
        else{
            re1<-wilcox.test(v1, v2, alternative ="two.sided")
            wil.res$rawP[i]<-re1[[3]]
            wil.res$mean1[i]<-mean(v1, na.rm=T)
            wil.res$mean2[i]<-mean(v2, na.rm=T)
            #which is higher in diversity
            wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
        }
    }
    wil.res$Tissue<-tissues[j]
    wil.res<-wil.res[!is.na(wil.res$rawP),]
    Results<-rbind(Results, wil.res)
}

Re<-Pcorrection(Results)
write.csv(Re, "Output/Diversity/Tests/Cohort_comparison_Wilcox_results_all.csv")



#calcualte the averages
sivMean<-aggregate(infec[,"mean"],by=list(infec$Cohort, infec$Tissue),mean,na.rm=T )
colnames(sivMean)<-c("Cohort", "Tissue","Mean")

#Order the factors
infec$Tissue<-factor(infec$Tissue, levels = c("Plasma","Plasma nex","Peripheral LN","Thoracic LN","Lung"))
sivMean$Tissue<-factor(sivMean$Tissue, levels = c("Plasma","Plasma nex","Peripheral LN","Thoracic LN","Lung"))
infec$Cohort<-factor(infec$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
sivMean$Cohort<-factor(sivMean$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))

#add jitter to one of the lung 
which(infec$Tissue=="Lung"& infec$Cohort=="Mtb NR")
#64, 67
infec$mean[67]<-0.00195
infec$mean[64]<-0.00205

p2<-ggplot()+
    geom_point(data=infec, aes(x=Tissue, y=mean*100, color=Cohort),alpha=0.6,
               position=position_jitterdodge(jitter.width = 0.1,jitter.height = 0,
                                             dodge.width = 0.5), size=2)+
    geom_point(data=sivMean, aes(x=Tissue, y=Mean*100, fill=Cohort), size=1.7,position = position_dodge(width = 0.5), shape=4)+
    xlab('')+ylab('% Average diversity')+
    theme(legend.title = element_blank())+
    scale_y_continuous(limits = c(0.16,0.49))+
    scale_color_manual(values=colors[c(5,3,1)])+
    scale_fill_manual(values=c(1,1,1),guide="none")+
    theme(theme(legend.title = element_blank()))+theme_bw()+
    theme(legend.title = element_blank(), panel.grid.major.x = element_blank())

ggsave("Output/Figures/Cohort_noletters.png", width = 5.5,height = 3.5, dpi=300)



#######
### Test if diversity in one tissue type is different from another tissue within each monkey

list.animal<-split(samples, samples$Monkey)
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
    if (nrow(list.animal[[i]])>1){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
    }
}
monkeys<-names(monkeyList)

#comparison pairs 
tis<-unique(Sum21$Tissue2)
tis[tis=="Plasma nex"]<-"Plasma_nex"
comb<-combn(tis,2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

Results<-data.frame()
for (j in 1:length(monkeys)){
    sample<-monkeyList[[j]]
    monkey<-names(monkeyList)[j]
    
    df<-Sum21[Sum21$Monkey==monkey,]

    #wicoxon test
    wil.res<-data.frame(Test=comb_pairs)
    
    for (i in 1:nrow(comb)){
        v1<-df$mean[df$Tissue2==comb[i,1]]
        v2<-df$mean[df$Tissue2==comb[i,2]]
        if (length(v1)<=1|length(v2)<=1){
            wil.res$rawP[i]<-NA
            wil.res$mean1[i]<-NA
            wil.res$mean2[i]<-NA
            wil.res$higher[i]<-NA
        }
        else{
            re1<-wilcox.test(v1, v2, alternative ="two.sided")
            wil.res$rawP[i]<-re1[[3]]
            wil.res$mean1[i]<-mean(v1, na.rm=T)
            wil.res$mean2[i]<-mean(v2, na.rm=T)
            #which is higher in diversity
            wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
        }
    }
    wil.res$Monkey<-monkey
    
    
    wil.res<-Pcorrection(wil.res)
    wil.res<-wil.res[!is.na(wil.res$rawP),]
    #write.csv(wil.res, paste0("Output/Diversity/Tests/",monkey,"_wilcox_results_tissues.csv"))
    Results<-rbind(Results, wil.res)
}
write.csv(Results, "Output/Diversity/Tests/Wilcox_results_all.csv")

################ Fig 5C ##############
## Plot all diversity (Figure 5C)


## with mean
all_mean<-Sum21[,c("Cohort","Monkey","Tissue2","Tissue3","Granuloma","Week","mean","median")]
all_mean<-InsertRow(all_mean,c(0,0,0,0,0,0, stks[1,c("mean","median")]),1)
all_mean[1,1:4]<-c("Stock","Stock","Stock","Stock")
#change plasma_nex to plasma for this plot
all_mean$Tissue3[all_mean$Tissue3=="Plasma nex"]<-"Plasma"
all_mean$label<-paste0(all_mean$Monkey,'_', all_mean$Tissue2,'_', all_mean$Week)
all_mean$Tissue3<-factor(all_mean$Tissue3, levels=c("Stock", "Plasma", "Peripheral LN","Thoracic LN","Lung"))

#order by ->stock, siv only, NR, N
all_mean$Cohort<-factor(all_mean$Cohort, levels=c("Stock", "SIV only", "Mtb NR", "Mtb R"))
all_mean<-all_mean[order(all_mean$Cohort, all_mean$Monkey,all_mean$Week, all_mean$Tissue3),]
labels<-unique(all_mean$label)
all_mean$label<-factor(all_mean$label, levels=paste(labels))

#plot all freq
AllFreq<-read.csv("Output/Diversity/Allfile_frequency.csv")
AllFreq<-AllFreq[,-1]
AllFreq$label<-factor(AllFreq$label, levels=paste(labels))
AllFreq$Tissue3[AllFreq$Tissue3=="Plasma nex"]<-"Plasma"
AllFreq$Tissue3<-factor(AllFreq$Tissue3, levels=c("Stock", "Plasma","Peripheral LN","Thoracic LN","Lung"))

#data for cohort and monkey info for the plot
dat<-distinct(all_mean, label, .keep_all=T) #50 rows
numb<-data.frame(table(dat$Monkey))
colnames(numb)<-c("Monkey","Freq")
morder<-unique(all_mean$Monkey)
#write.csv(morder,"Output/monkey_order.csv")
numb$Monkey<-factor(numb$Monkey, levels=paste(morder))
numb<-numb[order(numb$Monkey),]

numb$col<-numb$Freq
for (i in 2:nrow(numb)){
    numb$col[i]<-numb$Freq[i]+numb$col[i-1]
}
numb$lab_pos<-numb$Freq/2    
for (i in 2:nrow(numb)){
    numb$lab_pos[i]<-(numb$col[i]-numb$col[i-1])/2+numb$col[i-1]
}

pos1<-(numb$col[5]-1)/2+0.5
pos2<-(numb$col[8]-numb$col[5])/2+numb$col[5]+0.5
pos3<-(50-numb$col[8])/2+numb$col[8]+0.5

xlabels<-dat$Week
xlabels[1]<-"Stk"

#samples with granuloma
all_mean$mean2<-all_mean$mean
all_mean$mean2[all_mean$Granuloma!="Y"]<-NA
all_mean$mean2[is.na(all_mean$Granuloma)]<-NA
all_mean<-transform(all_mean,id=as.numeric(factor(label)))

#as.numeric(factor(all_mean$label))

#Points only
xlabels2<-xlabels
xlabels2[1]<-""
p4<-ggplot()+
    geom_point(data=all_mean, aes(x=label, y=mean*100, color=Tissue3),size=2,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Avearge diversity')+
    geom_point(data=all_mean, aes(x=label, y=mean2*100,color=Tissue3),size=2.2,shape=21, color="gray10",
              position = position_dodge(width = .5))+
    scale_y_continuous(limits = c(0.1,0.5))+
    scale_x_discrete(labels=xlabels2)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.1, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=0.5, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=0.5, label="Latent NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=0.5, label="Latent R", hjust=0.5, size=3.5, color='gray20')+
    #annotate("text",x=1, y=0.1, label="stock" , angle=90, hjust=1, size=3, color="gray30")+
    geom_text(aes(x=1, y=-Inf, label="stock") , angle=90, hjust=1.1, size=3, color="gray30")+
    #annotate("point", x=1:0, y=all_mean$mean2*100, shape=21, color="gray10", size=2.2)+
    coord_cartesian(clip = "off")
#ggsave("Output/Figures/Diversity_allSamples.png", width = 8,height = 3,units="in", dpi=300)
ggsave("Output/Figures/Diversity_allSamples.pdf", width = 8,height = 3)


#Violin plot

ggplot()+
    geom_violin(data=AllFreq, aes(x=label, y=freq.mutations*100, color=Tissue3), size=.2,trim=F, show.legend = F)+
    scale_y_continuous(trans = 'log10',breaks = c(0.01,1,100), labels=c(0.01,1,100), limits = c(0.0002,140),expand = c(0,0))+
    geom_point(data=all_mean, aes(x=label, y=mean*100, color=Tissue3),size=1,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Diversity')+
    geom_point(data=all_mean, aes(x=label, y=mean2*100,color=Tissue3),size=1.2,shape=21, color="gray10",
               position = position_dodge(width = .5))+
    scale_x_discrete(labels=xlabels)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.0003, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=0.0002, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=0.0002, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=0.0002, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=100, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=100, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=100, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    #annotate("point", x=all_mean$id, y=all_mean$mean2*100, shape=21, color="gray10", size=1.2)+
    geom_text(aes(x=1, y=-Inf, label="stock"), angle=90, hjust=1.1, size=3, color="gray30")+
    coord_cartesian(clip = "off")
ggsave("Output/Figures/Diversity_allsamples_Violinplot_mean.pdf",width = 10,height = 4.5)
ggsave("Output/Figures/Diversity_allSamples_Violinplot.pdf", width =  10,height = 4.5)




ggplot()+
    geom_boxplot(data=AllFreq, aes(x=label, y=freq.mutations*100, color=Tissue3),
                 position=position_dodge(width = 0.5), alpha=0.2)+
    scale_y_continuous(trans = 'log10',breaks = c(0.01,1,50,100), labels=c(0.01,1,50,100), limits = c(0.0002,140),expand = c(0,0))+
    geom_point(data=all_mean, aes(x=label, y=mean*100, fill=Tissue3),size=1.6,
               position = position_dodge(width = .5), shape=21, color="white")+
    xlab('')+ylab('% Diversity')+
    scale_x_discrete(labels=xlabels)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.0003, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=0.0002, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=0.0002, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=0.0002, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=100, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=100, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=100, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    coord_cartesian(clip = "off")+
    annotate("point", x=all_mean$id, y=all_mean$mean2*100, shape=21, color="gray10", size=1.6, stroke = 0.5)

ggsave("Output/Figures/Diversity_Allsamples_Boxplot_.pdf",width = 14,height = 4.5)




###  Use median  ####

#samples with granuloma
all_mean$median2<-all_mean$median
all_mean$median2[all_mean$Granuloma!="Y"]<-NA
all_mean$median2[is.na(all_mean$Granuloma)]<-NA

p3<-ggplot()+
    geom_point(data=all_mean, aes(x=label, y=median*100, color=Tissue3),size=2,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Median diversity')+
    #geom_point(data=all_mean, aes(x=label, y=mean2*100,color=Tissue3),size=2.2,shape=21, color="gray10",
    #           position = position_dodge(width = .5))+
    scale_y_continuous(limits = c(0.1,0.24))+
    scale_x_discrete(labels=xlabels2)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.1, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=0.24, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=0.24, label="Latent NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=0.24, label="Latent R", hjust=0.5, size=3.5, color='gray20')+
    #annotate("text",x=1, y=0.1, label="stock" , angle=90, hjust=1, size=3, color="gray30")+
    geom_text(aes(x=1, y=-Inf, label="stock") , angle=90, hjust=1.1, size=3, color="gray30")+
    annotate("point", x=all_mean$id, y=all_mean$median2*100, shape=21, color="gray10", size=2.2)+
    coord_cartesian(clip = "off")
#ggsave("Output/Figures/Diversity_allSamples.png", width = 8,height = 3,units="in", dpi=300)
ggsave("Output/Figures/Diversity_allSamples_median.pdf", width = 8,height = 3)


library(cowplot)
grid.arrange(arrangeGrob(p1,p2, ncol=2),p3, nrow=2 )

#draw_plot(plot, x, y, width, height)
pdf("Output/Figures/Figure5.pdf",width=9, height=6)
png("Output/Figures/Figure5.png",width=9, height=6, units="in",res=300)

ggdraw()+
    draw_plot(p1,0,0.5,0.36,0.5)+
    draw_plot(p2,0.36,0.5,0.6,0.5)+
    draw_plot(p3,0,0,1,0.5)+
    draw_plot_label(c("A", "B", "C"), c(0, 0.36, 0), c(1, 1, 0.5), size = 15)
dev.off()    

#mean
pdf("Output/Figures/Figure5_Mean.pdf",width=9, height=6)
#png("Output/Figures/Figure5.png",width=9, height=6, units="in",res=300)

ggdraw()+
    draw_plot(p1,0,0.5,0.36,0.5)+
    draw_plot(p2,0.36,0.5,0.6,0.5)+
    draw_plot(p4,0,0,1,0.5)+
    draw_plot_label(c("A", "B", "C"), c(0, 0.36, 0), c(1, 1, 0.5), size = 15)
dev.off()    


###########
#Diversity by monkeys

#Wilcoxon test
ids<-c("SIV only","Mtb NR","Mtb R")

Co<-data.frame(cohort=coh)
for (j in 1:3){
    DF<-Sum21[Sum21$Cohort==coh[j],]
    Co$mean[j]<-mean(DF$mean)
    Co$se[j]<-mean(DF$se)
}

comb<-combn(coh,2)
comb<-t(comb)
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

wil.res<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    v1<-Sum21$mean[Sum21$Cohort==comb[i,1]]
    v2<-Sum21$mean[Sum21$Cohort==comb[i,2]]
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res$rawP[i]<-re1[[3]]
    wil.res$mean1[i]<-mean(v1, na.rm=T)
    wil.res$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
}

wil.res<-Pcorrection(wil.res)
write.csv(wil.res, paste0("Output/Diversity/Tests/Cohorts_div_wilcox_results.csv"))


#2. Diverstiy by Tissue types by Cohort
infec<-Sum21[,c("Cohort","Tissue3","Week","mean", "se")]
colnames(infec)[2]<-"Tissue"

infec$Cohort<-factor(infec$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.Coh<- aggregate(infec["mean"],by=list(infec$Cohort),mean,na.rm=T )
by.Coh$Group.1<-factor(by.Coh$Group.1, levels=c("SIV only", "Mtb NR", "Mtb R"))

#Plot cohort results
p1<-ggplot()+
    geom_point(data=infec, aes(x=Cohort, y=mean*100, color=Cohort),size=2.3, position=position_jitter(width=0.05))+
    geom_point(data=by.Coh, aes(x=Group.1, y=mean*100), color="black",size=2.3, shape=4)+
    xlab('')+ylab('% Average diversity')+
    scale_color_manual(values=paste0(colors[c(5,3,1)],"66"), guide='none')+
    theme(legend.title = element_blank())+
    scale_y_continuous(limits = c(0.1, 0.55))+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/Cohort_diversity.png", width = 3,height = 3,units="in", dpi=300)





##### Other analyses: RNA copy numbers vs. diversity
rna<-Sum21[,c("Cohort","Monkey","Tissue3","Week","mean","SIV.RNA.per.granuloma","SIV.RNA.per.tissue",
              "CD4.percent","CD8.percent")]

rna1<-rna[!is.na(rna$SIV.RNA.per.granuloma),]
rna2<-rna[!is.na(rna$SIV.RNA.per.tissue),]

cor.test(rna1$SIV.RNA.per.granuloma, rna1$mean, method = "spearman") #P=0.2473
cor.test(rna2$SIV.RNA.per.tissue, rna2$mean, method = "spearman") #P= 0.3812

ggplot(rna1, aes(x=SIV.RNA.per.granuloma, y=mean*100))+
    geom_point(color="steelblue",alpha=0.7, size=2.5)+ylab("% Average diversity")+
    xlab("RNA copy number in granulooma")+
    theme_bw()
ggsave("Output/Figures/Correlation_RNAgranuloma_Div.pdf", height = 3.6, width = 4)
ggplot(rna2, aes(x=SIV.RNA.per.tissue  , y=mean))+
    geom_point(color="steelblue",alpha=0.7, size=2.5)+ylab("% Average diversity")+
    xlab("RNA copy number")+
    theme_bw()
ggsave("Output/Figures/Correlation_RNAtissue_Div.pdf", height = 3.6, width = 4)


#####
#Look at tLN and Lung
rna2<-rna2[rna2$Tissue3=="Lung"|rna2$Tissue3=="Thoracic LN",]

rna.sum<-aggregate(rna2[,c("SIV.RNA.per.tissue","mean")], by=list(rna2$Monkey,rna2$Tissue3), mean)
se<-aggregate(rna2[,c("SIV.RNA.per.tissue","mean")], by=list(rna2$Monkey,rna2$Tissue3), std.error)

#remove<-names(which(table(rna.sum$Group.1) != 2))
#rna.sum<-rna.sum[!(rna.sum$Group.1 %in% remove),]
ggplot(rna.sum, aes(x=SIV.RNA.per.tissue, y=mean, color=Group.2))+
    geom_point(size=3)+ggtitle("RNA vs. diversity in Lungs and thoracic LNs")
ggsave("Output/Figures/RNA.vs.diversity.in.tLN-Lungs.pdf", width = 6, height = 5)

cor.test(rna.sum$SIV.RNA.per.tissue, rna.sum$mean,method="spearman")
#p-value = 0.05958 rho=0.4421053  
cor.test(rna.sum$SIV.RNA.per.tissue[rna.sum$Group.2=="Lung"], rna.sum$mean[rna.sum$Group.2=="Lung"],method="spearman")
#p-value = 0.115
cor.test(rna.sum$SIV.RNA.per.tissue[rna.sum$Group.2=="Thoracic LN"], rna.sum$mean[rna.sum$Group.2=="Thoracic LN"],method="spearman")
#p-value = 0.4182

## RNA concentrations differ between tLN and Lung?
wilcox.test(rna2$SIV.RNA.per.tissue[rna2$Tissue3=="Lung"], rna2$SIV.RNA.per.tissue[rna2$Tissue3=="Thoracic LN"], alternative="two.sided")
# p-value = 9.991e-07 
#Yes, but Lung vs. tLN diversity did not as aggregate
wilcox.test(rna2$mean[rna2$Tissue3=="Lung"], rna2$mean[rna2$Tissue3=="Thoracic LN"], alternative="two.sided")
#p-value = 0.1267

#Look at each monkey separately to test if RNA concentrations differ between Lung and tLN
ms<-unique(rna.sum$Group.1)

sumtable<-data.frame(Monkey=ms)
for (i in 1:length(ms)){
    df<-rna2[rna2$Monkey==ms[i],]
    if (length(df$SIV.RNA.per.tissue[df$Tissue3=="Lung"])==0|length(df$SIV.RNA.per.tissue[df$Tissue3=="Thoracic LN"])==0){
        sumtable$RNA.pvalue[i]<-NA
        sumtable$Diversity.pvalue[i]<-NA
    }
    else {
        r1<-wilcox.test(df$SIV.RNA.per.tissue[df$Tissue3=="Lung"],df$SIV.RNA.per.tissue[df$Tissue3=="Thoracic LN"], alternative="two.sided")
        sumtable$RNA.pvalue[i]<-r1[[3]]
        r2<-wilcox.test(df$mean[df$Tissue3=="Lung"],df$mean[df$Tissue3=="Thoracic LN"], alternative="two.sided")
        sumtable$Diversity.pvalue[i]<-r2[[3]]
    }
}

#non are significant
write.csv(sumtable,file="Output/Diversity/Tests/RNA.conc_test.within.eachmonkey.csv")

ggplot(rna2, aes(x=Monkey, y=mean, color=Tissue3))+
    geom_point(position=position_dodge(width = 0.5))
ggplot(rna2, aes(x=Monkey, y=SIV.RNA.per.tissue, color=Tissue3))+
    geom_point(position=position_dodge(width = 0.5))

#Rna copy number summary for each chort and tissue
rna<-rna[!(rna$Tissue3=="Plasma"|rna$Tissue3=="Plasma nex"),]

RNAsum<-data.frame(aggregate(rna[,"SIV.RNA.per.tissue"], by=list(rna$Cohort,rna$Tissue3), mean, na.rm=T))
colnames(RNAsum)<-c("Cohort", "Tissue", "Mean")
se<-data.frame(aggregate(rna[,"SIV.RNA.per.tissue"], by=list(rna$Cohort,rna$Tissue3), std.error, na.rm=T))
RNAsum$SE<-se$x

write.csv(RNAsum,"Output/RNA.summary.perCohortperTissue.csv")
t(RNAsum)
