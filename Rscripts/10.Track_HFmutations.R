#Look at 'genetic drift' at moving into tissues from plasma
#mutation freq differences between tissues comparison

library(ggplot2)
library(reshape)
library(gridExtra)
source("Rscripts/label_scientific.R")
library(colorspace)
colors<-qualitative_hcl(6, palette="Dark3")

# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}
#sample info
SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
#stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

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


###############
samples_org<-samples
samples$tis<-samples$Tissue2
### Look for hsigh mut freq sites in laster samples and compare 
samples$tis[samples$tis=="plasma" & samples$Week>5]<-"plasma_nex"
samples<-samples[samples$tis!="plasma",]

Results<-list()
for (i in 1:length(monkeys)){
    sample<-samples[samples$Monkey==monkeys[i],]
    monkey<-monkeys[i]
    Ov<-OvDF[sample$filename]
    
    high.pos<-c()
    pos.type<-c()
    for (j in 1: length(Ov)){
        df<-OvDF[[sample$filename[j]]]
        df[which(df$TotalReads<100),17:26]<-NA
        df<-df[df$ref251.pos!=395,]
        df<-df[df$ref251.pos!=400,]
        df<-df[df$freq.mutations.ref>0.008,]
        df<-df[!is.na(df$ref251.pos),] 
        df$MFtype<-apply(df[,c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref")],1,function(x) c("freq.Ts.ref","freq.transv1.ref","freq.transv2.ref")[which.max(x)])
        high.pos<-c(high.pos, df$ref251.pos)
        pos.type<-c(pos.type, df$MFtype)
    }
    
    high.po<-data.frame(cbind(high.pos,pos.type))
    high.po<-unique(high.po[,1:2] )
    
    res<-sample[,c("Monkey","tis","Week","Cohort","filename")]
    for (f in 1:length(Ov)){
        dat<-Ov[[res$filename[f]]]
        #dat<-dat[dat$ref251.pos %in% high.po$high.pos,]
        for (j in 1:nrow(high.po)){
            res[f,(j+5)]<-dat[dat$ref251.pos==high.po$high.pos[j],paste0(high.po$pos.type[j])]
        }
    }
    types<-gsub('.ref', '',high.po$pos.type)
    types<-gsub('freq.', '', types)
    colnames(res)[6:ncol(res)]<-paste0("pos.",high.po$high.pos,".",types)
    write.csv(res,paste0("Output/MF/High.mf.", monkey,".csv"))
    Results[[i]]<-res
    names(Results)[i]<-monkey
}


Diff<-list()
for (i in 1:length(Results)){
    freq<-Results[[i]]
    freq$tis<-factor(freq$tis, c("plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    monkey<-names(Results)[i]
    
    diff<-data.frame(comp=c("plasma-pLN","plasma-tLN",'plasma-Lung'))
    for (c in 1:(ncol(freq)-5)){
        df<-freq[,c(1:5,(c+5))]
        means<-data.frame(aggregate(df[6],by=list(df$tis), mean, na.rm=T))
        
        diff[1,(c+1)]<-means[means$Group.1=="plasma_nex",2]-means[means$Group.1=="pLN",2]
        diff[2,(c+1)]<-means[means$Group.1=="plasma_nex",2]-means[means$Group.1=="Thoracic LN",2]
        if (nrow(means[means$Group.1=='Lung',])!=0){
            diff[3,(c+1)]<-means[means$Group.1=="plasma_nex",2]-means[means$Group.1=="Lung",2]
        }
        colnames(diff)[c+1]<-colnames(df)[6]
    }    
    
    cols<-paste(diff$comp)
    diff<-diff[,-1]
    diffT<-data.frame(t(diff), stringsAsFactors = F)
    colnames(diffT)<-cols
    Diff[[i]]<-diffT
    names(Diff)[i]<-monkey
    
    
}

#export the file for Monkey 3116 since the difference is small but significat
write.csv(Diff[[4]], "Output/MF/Difference.in.MF.3116.csv")

summary<-data.frame(monkey=monkeys)
for (i in 1:length(Diff)){
    diffT<-Diff[[i]]
    summary[i,2:4]<-apply(diffT, 2, function(x) mean(abs(x), na.rm=T))
    r1<-wilcox.test(abs(diffT$`plasma-pLN`),abs(diffT$`plasma-tLN`), alternative="two.sided")
    summary[i,5]<-r1[[3]]
    lung<-diffT$`plasma-Lung`
    if (length(lung[!is.na(lung)])!=0){
        r2<-wilcox.test(abs(diffT$`plasma-pLN`),abs(diffT$`plasma-Lung`), alternative="two.sided")
        r3<-wilcox.test(abs(diffT$`plasma-tLN`),abs(diffT$`plasma-Lung`), alternative="two.sided")
        summary[i,6]<-r2[[3]]
        summary[i,7]<-r3[[3]]
    }

}

colnames(summary)[2:7]<-c("Mean_plasma-pLN","Mean_plasma-tLN","Mean_plasma-Lung",
                          "Pvalue_pLNvs.tLN","Pvalue_pLNvs.lung" , "Pvalue_tLNvs.Lung")
write.csv(summary, "Output/MF/Difference_MF_betwTissues_mean_p-val.csv")




### 3116 monkey -are they really significantly different?
diffT<-read.csv( "Output/MF/Difference.in.MF.3116.csv", stringsAsFactors = F, row.names = 1)
plot(abs(  diffT$plasma.tLN),   pch=16,  col="red") #Thoracic
points(abs(diffT$plasma.Lung),pch=16,  col="green") #Lung
points(abs(diffT$plasma.pLN), pch=16,  col="blue")
abline(h=mean(abs(diffT$plasma.tLN)), col="red")
abline(h=mean(abs(diffT$plasma.Lung)), col="green")
mean(abs(diffT$plasma.tLN))

wilcox.test(abs(diffT$plasma.Lung),abs(diffT$plasma.tLN), alternative="two.sided")
#
#Wilcoxon rank sum test with continuity correction
#
#data:  abs(diffT$plasma.Lung) and abs(diffT$plasma.tLN)
#W = 806, p-value = 0.005119
#alternative hypothesis: true location shift is not equal to 0
#
#Warning message:
#    In wilcox.test.default(abs(diffT$plasma.Lung), abs(diffT$plasma.tLN),  :
#                               cannot compute exact p-value with ties



#### Plot the difference
diff<-read.csv("Output/MF/Difference_MF_betwTissues_mean_p-val.csv", stringsAsFactors = F, row.names = 1)
coh<-unique(samples[,c('Monkey','Cohort')])
coh$monkey<-as.integer(coh$Monkey)
diff<-merge(diff, coh[,2:3],by="monkey")
m<-read.csv("Output/monkey_order.csv", stringsAsFactors = F)
morder<-m$x

diff$Cohort<-factor(diff$Cohort,levels=c("SIV only","Mtb NR","Mtb R"))
diff$monkey<-factor(diff$monkey, levels=morder[2:12])
diff<-diff[order(diff$monkey),]

colnames(diff)[2:4]<-c("plas.pLN","plas.tLN","plas.Lung")

diffM<-melt(diff[,c(1:4,8)], id.vars=c("Cohort","monkey"))
diffM$monkey<-factor(diffM$monkey, levels=morder)
cols<-c("#023FA5","#C4A12C", "#8E063B")

clables<-c("Plasma vs. Perif LN", "Plasma vs. Thoracic LN", "Plasma vs. Lung")
colnames(diffM)[3]<-"Comparison"

d2<-diff[,c(1,5:8)]

for (i in 1:nrow(d2)){
    p<-d2[i,2:4]
    sig<-which(p<0.05)
    sig[sig==1]<-'a'
    sig[sig==2]<-'b'
    sig[sig==3]<-'c'
    
    sig<-paste(sig, collapse = ",")
    d2$sign[i]<-sig
    
}


ggplot(diffM,aes(x=factor(monkey), y=value*100, color=Comparison))+
    geom_point(position=position_dodge(width=0.7))+ylab('Average difference in frequency (%)')+
    xlab('')+geom_vline(xintercept =c(1:10)+0.5,  
                        color = "gray80", size=.4)+
    scale_color_manual(values = cols, labels=clables)+
    
    annotate("rect", xmin=0.5, xmax=4.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=4.5, xmax=7.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=7.5, xmax=11.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=2.5, y=9.7, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=6,   y=9.7, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=9.5, y=9.7, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    scale_y_continuous(limits=c(0,9.9),expand = c(0, 0.1))+
    theme_bw()+
    annotate("text", x=1:11, y=0.2, label=d2$sig, hjust=0.5, size=2.4, color='gray50')+
    
    theme(panel.grid.major.x = element_blank())+
    theme(axis.ticks.x = element_blank())
ggsave("Output/Figures/Difference.in.highMF.drift.png", width =7.5 , height =3 , dpi=300, units = "in")



aggregate(diff[,c(2:4)], by=list(diff$Cohort), mean, na.rm=T)
#   Group.1   plas.pLN   plas.tLN  plas.Lung
#1 SIV only 0.02058346 0.04151433 0.05201125
#2   Mtb NR 0.02695502 0.03587039 0.04228236
#3    Mtb R 0.02261113 0.03437687 0.02148160

hcl_palettes(plot = TRUE)

mean(diff$plas.pLN, na.rm = T) #0.02305849
mean(diff$plas.tLN, na.rm = T) # 0.03737963
mean(diff$plas.Lung, na.rm = T) # 0.0343142
mean(c(diff$plas.pLN,diff$plas.tLN,diff$plas.Lung), na.rm = T)
#[1] 0.0313111

wilcox.test(diff$plas.pLN,diff$plas.tLN)
#W = 43, p-value = 0.2703
wilcox.test(diff$plas.pLN,diff$plas.Lung)
#W = 39, p-value = 0.7168
wilcox.test(diff$plas.Lung,diff$plas.tLN)
#W = 41, p-value = 0.8404

