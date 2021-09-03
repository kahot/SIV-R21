#Look at 'genetic drift' at moving into tissues from plasma
#mutation freq differences between tissues comparison

library(ggplot2)
library(reshape2)
library(gridExtra)
library(colorspace)
source("Rscripts/label_scientific.R")
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
SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
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
#remove plasma samples from early time points
samples<-samples[samples$Tissue2!="Plasma",]

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
    
    res<-sample[,c("Monkey","Tissue3","Week","Cohort","filename","SIV.RNA.per.tissue")]
    for (f in 1:length(Ov)){
        dat<-Ov[[res$filename[f]]]
        #dat<-dat[dat$ref251.pos %in% high.po$high.pos,]
        for (j in 1:nrow(high.po)){
            res[f,(j+6)]<-dat[dat$ref251.pos==high.po$high.pos[j],paste0(high.po$pos.type[j])]
        }
    }
    types<-gsub('.ref', '',high.po$pos.type)
    types<-gsub('freq.', '', types)
    colnames(res)[7:ncol(res)]<-paste0("pos.",high.po$high.pos,".",types)
    #write.csv(res,paste0("Output/MF/High.mf.", monkey,".csv"))
    Results[[i]]<-res
    names(Results)[i]<-monkey
}


Diff<-list()
rna<-list()
for (i in 1:length(Results)){
    freq<-Results[[i]]
    freq$Tissue3<-factor(freq$Tissue3, c("Plasma nex", "Peripheral LN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$Tissue3),]
    monkey<-names(Results)[i]
    
    rna[[i]]<-data.frame(aggregate(freq[6],by=list(freq$Tissue3), mean, na.rm=T))
    names(rna)[i]<-monkey
    
    diff<-data.frame(comp=c("Plasma-pLN","Plasma-tLN",'Plasma-Lung'))
    for (c in 1:(ncol(freq)-6)){
        df<-freq[,c(1:6,(c+6))]
        means<-data.frame(aggregate(df[7],by=list(df$Tissue3), mean, na.rm=T))
        
        diff[1,(c+1)]<-means[means$Group.1=="Plasma nex",2]-means[means$Group.1=="Peripheral LN",2]
        diff[2,(c+1)]<-means[means$Group.1=="Plasma nex",2]-means[means$Group.1=="Thoracic LN",2]
        if (nrow(means[means$Group.1=='Lung',])!=0){
            diff[3,(c+1)]<-means[means$Group.1=="Plasma nex",2]-means[means$Group.1=="Lung",2]
        }
        colnames(diff)[c+1]<-colnames(df)[7]
    }    
    
    cols<-paste(diff$comp)
    diff<-diff[,-1]
    diffT<-data.frame(t(diff), stringsAsFactors = F)
    colnames(diffT)<-cols
    Diff[[i]]<-diffT
    names(Diff)[i]<-monkey
}


#Summarize of the frequency differences and RNA copy numbers

## mean (each monkey is an independent experiment)
summary<-data.frame(monkey=monkeys)
rna_sum<-data.frame(monkey=monkeys)
for (i in 1:length(Diff)){
    diffT<-Diff[[i]]
    monkey<-names(Diff)[i]
    summary[i,2:4]<-apply(diffT, 2, function(x) mean(abs(x), na.rm=T))
    lung<-samples$Tissue2[samples$Monkey==monkey&samples$Tissue2=="Lung"]
    stats<-data.frame(test=c("pLN.vs.tLN","pLN.vs.Lung","tLN.vs.Lung"))
    if (length(lung[!is.na(lung)])==0){
        r1<-wilcox.test(abs(diffT$`Plasma-pLN`),abs(diffT$`Plasma-tLN`), alternative="two.sided")
        summary[i,5]<-r1[[3]]
        summary[i,6]<-NA
        summary[i,7]<-NA
    }
    else{
        r1<-wilcox.test(abs(diffT$`Plasma-pLN`),abs(diffT$`Plasma-tLN`))
        r2<-wilcox.test(abs(diffT$`Plasma-pLN`),abs(diffT$`Plasma-Lung`))
        r3<-wilcox.test(abs(diffT$`Plasma-tLN`),abs(diffT$`Plasma-Lung`))
        stats$rawP[1]<-r1[[3]]        
        stats$rawP[2]<-r2[[3]]       
        stats$rawP[3]<-r3[[3]]
        stats<-Pcorrection(stats)
        summary[i,5]<-stats$Holm[1]
        summary[i,6]<-stats$Holm[2]
        summary[i,7]<-stats$Holm[3]
        
    }
    
    rn<-rna[[i]]
    rna_sum[i, 2:4]<-rn$SIV.RNA.per.tissue[2:4]
}

colnames(summary)[2:7]<-c("Mean_plasma-pLN","Mean_plasma-tLN","Mean_plasma-Lung",
                          "Pvalue_pLNvs.tLN","Pvalue_pLNvs.lung", "Pvalue_tLNvs.Lung")
colnames(rna_sum)[2:4]<-c("RNA.pLN","RNA.tLN","RNA.Lung")

write.csv(summary, "Output/MF/Difference_MF_betwTissues_mean_p-val.csv")
write.csv(rna_sum, "Output/MF/rna.copy.no_summary.csv")


#Test differences by tissues as aggregate
tissue<-c("pLN","tLN","Lung")
for (i in 1:length(tissue)){
    x<-lapply(Diff,'[', paste0("Plasma-",tissue[i]))
    x<-unlist(x)
    assign(paste(tissue[i]),abs(x))
}

comb<-combn(tissue,2)
comb<-t(comb)
co_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

wil.res<-data.frame(Test=co_pairs)
for (i in 1:nrow(comb)){
    x<-get(paste(comb[i,1]))
    y<-get(paste(comb[i,2]))
    r1<-wilcox.test(x,y)
    wil.res$rawP[i]<-r1[[3]]        
}

stats<-Pcorrection(wil.res)
#None are significant

# Assess RNA copy numbers vs. diversity relationship
rna.sum1<-summary[,c(1:4)]
rna.summ1<-melt(rna.sum1, id.var="monkey")
rna.summ1$tissue<-c(rep("pLN", time=11), rep("tLN", times=11), rep("Lung", times=11))
colnames(rna.summ1)[3]<-"Difference"

#RNA copy numbers
rna.sum2<-rna_sum[,c(1:4)]
rna.summ2<-melt(rna.sum2)
colnames(rna.summ2)[3]<-"RNA"
rna.sum<-cbind(rna.summ1, rna.summ2[,"RNA"])
colnames(rna.sum)[5 ]<-"RNA"

#plot(rna.sum$Difference, rna.sum$RNA, pch=16, col="blue")
cor.test(rna.sum$Difference, rna.sum$RNA, method="spearman")
#p-value =  0.6313   rho = 0.09098999 

#within tissue correlation 
cor.test(rna.sum$Difference[rna.sum$tissue=="pLN"], rna.sum$RNA[rna.sum$tissue=="pLN"], method="spearman")
#S = 166, p-value = 0.4682
cor.test(rna.sum$Difference[rna.sum$tissue=="tLN"], rna.sum$RNA[rna.sum$tissue=="tLN"], method="spearman")
#S = 204, p-value = 0.8388
cor.test(rna.sum$Difference[rna.sum$tissue=="Lung"], rna.sum$RNA[rna.sum$tissue=="Lung"], method="spearman")
#S = 72, p-value = 0.752
#Non are significant


#### Calculate the median
summary2<-data.frame(monkey=monkeys)
for (i in 1:length(Diff)){
    diffT<-Diff[[i]]
    summary2[i,2:4]<-apply(diffT, 2, function(x) median(abs(x), na.rm=T))
}
colnames(summary2)[2:4]<-c("Median_plasma.pLN","Median_plasma.tLN","Median_plasma.Lung")
summary<-cbind(summary, summary2[,2:4])
write.csv(summary, "Output/MF/Difference_MF_betwTissues_withMdian.csv")

#### Plot the difference
diff<-read.csv("Output/MF/Difference_MF_betwTissues_withMdian.csv", stringsAsFactors = F, row.names = 1)
coh<-unique(samples[,c('Monkey','Cohort')])
coh$monkey<-as.integer(coh$Monkey)
diff<-merge(diff, coh[,2:3],by="monkey")
m<-read.csv("Output/monkey_order.csv", stringsAsFactors = F)
morder<-m$x

diff$Cohort<-factor(diff$Cohort,levels=c("SIV only","Mtb NR","Mtb R"))
diff$monkey<-factor(diff$monkey, levels=morder[2:12])
diff<-diff[order(diff$monkey),]

#colnames(diff)[2:4]<-c("plas.pLN","plas.tLN","plas.Lung")

#### Plot median
diffM2<-melt(diff[,c(1,8:11)], id.vars=c("Cohort","monkey"))
diffM2$monkey<-factor(diffM2$monkey, levels=morder)
cols<-c("#023FA5","#C4A12C", "#8E063B")

clables<-c("Plasma vs. Peripheral LN", "Plasma vs. Thoracic LN", "Plasma vs. Lung")
colnames(diffM2)[3]<-"Comparison"


ggplot(diffM2,aes(x=factor(monkey), y=value*100, color=Comparison))+
    geom_point(position=position_dodge(width=0.7))+
    ylab(expression(paste('Absolute difference \n   in frequency (%)')))+
    xlab('')+geom_vline(xintercept =c(1:10)+0.5,  
                        color = "gray80", size=.4)+
    scale_color_manual(values = cols, labels=clables)+
    
    annotate("rect", xmin=0.5, xmax=4.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=4.5, xmax=7.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=7.5, xmax=11.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=2.5, y=2.4, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=6,   y=2.4, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=9.5, y=2.4, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    scale_y_continuous(limits=c(0,2.5),expand = c(0, 0.1))+
    theme_bw()+
    #annotate("text", x=1:11, y=0.2, label=d2$sig, hjust=0.5, size=2.4, color='gray50')+
    theme(axis.title.y = element_text(hjust = 0.5))+
    theme(panel.grid.major.x = element_blank())+
    theme(axis.ticks.x = element_blank())+
    theme(plot.margin=unit(c(0.3,0.2,0,1),"cm"))
ggsave("Output/Figures/Difference.in.highMF.drift_median.png", width =8 , height =3 , dpi=300, units = "in")



