#Diversity comparison
library(ggplot2)
library(reshape)
library(gridExtra)
source("Rscripts/label_scientific.R")

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
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]
Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")

#high mut freq in stock
stOv<-OvDF[[stock$filename]]
high<-stOv[stOv$freq.mutations>0.008,]
high<-high[!is.na(high$ref251.pos),]
high<-high[high$ref251.pos!=395,] 
high<-high[high$TotalReads>100,] #23 sites all are transition mutaitons

highpos<-high$ref251.pos

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

#look for the high mut positions in each monkey


results<-samples[,c("Monkey","Tissue2","Week","Cohort","filename")]
results$sample<-paste0("Week",results$Week,"_",results$Tissue2)

for (i in 1:nrow(results)){
    df<-OvDF[[results$filename[i]]]
     for (j in 1:length(highpos)){
        results[i,(j+6)]<-df$freq.Ts[df$ref251.pos==highpos[j]]
    }
}

colnames(results)[7:29]<-paste0("pos.",highpos)
write.csv(results, "Output/MF/BottleneckTrack.csv")


highPos<-paste0("pos.",highpos)

#create a list for each mutations
for (i in 1:length(highPos)){
    L<-list()
    assign(highPos[i],L)
}

for (j in 1:length(monkeyList)){
    monkey<-names(monkeyList)[j]
    
    freq<-results[results$Monkey==monkey,]
    freq<-freq[order(freq$Week),]
    freq$Tissue2<-factor(freq$Tissue2, c("plasma", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$Tissue2),]
    for (i in 1:length(highPos)){
        df<-freq[,c(1:6,(i+6))]
        colnames(df)[7]<-"Freq"
        df$sample<-factor(df$sample, levels=unique(df$sample))
        p1<-ggplot(df, aes(x=sample, y=Freq*100))+
            geom_point(color="blue")+
            #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100),  limits=c(0,1))+
            theme_bw()+
            ggtitle(monkey)+
            theme(plot.title = element_text(size=11))+
            xlab('')+ylab(paste0("Mut freq. at ",highPos[i]))+
            theme(axis.text.x = element_text(angle =45, hjust=1))   
        lis<-get(highPos[i])
        lis[[j]]<-p1
        names(lis)[j]<-monkey
        assign(highPos[i], lis)
    }
    
}

for (i in 1:length(highPos)){
    L<-get(highPos[i])
    png(paste0("Output/Figures/MF_",highPos[i],".png"), width = 6, height=15, units = "in",res=300)
    do.call(grid.arrange, c(L, ncol=2))
    dev.off()
    
}   


Positions<-data.frame(ref251.pos=highpos)
overview<-OvDF[[1]]
Positions$AA239pos<-overviews$AA239pos[overviews$ref251.pos %in% highpos]
Positions


######

results<-read.csv("Output/MF/BottleneckTrack.csv",stringsAsFactors = F, row.names = 1)
results2<-results
results2$Tissue2[results$Tissue2=="plasma"&results2$Week>5]<-"plasma_nex"
aggregate(results2[,c("pos.155")], by=list(results2$Tissue2),mean)
#      Group.1           x
#1        Lung 0.007766588
#2      plasma 0.002485597
#3  plasma_nex 0.003451250
#4         pLN 0.002073491
#5 Thoracic LN 0.001744538

#average-> Lung is higher than plasma, but pLN and LN are lower. 
re<-results2[,c(2,7:29)]
Sum<-aggregate(.~Tissue2, re, mean, na.action = na.omit)
mSum<-melt(Sum)

mSum$Tissue2<-factor(mSum$Tissue2,levels=c("plasma","plasma_nex","pLN","Thoracic LN", "Lung"))
ggplot(mSum, aes(x=Tissue2, y=value))+
    geom_point()+
    facet_wrap(~variable, scale='free',ncol=2)
ggsave("Output/MF/Track_high_mutations.pdf", width = 8,height = 20)    


all<-melt(re)
all$Tissue2<-factor(all$Tissue2,levels=c("plasma","plasma_nex","pLN","Thoracic LN", "Lung"))
ggplot()+
    geom_point(data=all, aes(x=Tissue2, y=value),position=position_jitter(width=0.1),color="dodgerblue4", alpha=0.5)+
    geom_point(data=mSum, aes(x=Tissue2, y=value), color="red")+
    facet_wrap(~variable, scale='free',ncol=2)
ggsave("Output/MF/Track_high_mutations_all2.pdf", width = 8,height = 20)    




###############
### Look for hsigh mut freq in Plasma_nex samples and compare with tissue
samples$tis<-samples$Tissue2
samples$tis[samples$tis=="plasma" & samples$Week>5]<-"plasma_nex"
samp<-samples[samples$tis=="plasma_nex",]
monk<-unique(samp$Monkey)

Results<-list()
for (i in 1:nrow(samp)){
    df<-OvDF[[samp$filename[i]]]
    monkey<-samp$Monkey[i]
    df[which(df$TotalReads<100),17:26]<-NA
    df<-df[df$ref251.pos!=395,]
    df<-df[df$freq.mutations>0.008&!is.na(df$ref251.pos),]
    df<-df[!is.na(df$ref251.pos),]
    
    df$MFtype<-apply(df[,c("freq.Ts","freq.transv1","freq.transv2")],1,function(x) c("freq.Ts","freq.transv1","freq.transv2")[which.max(x)])
    
    sample<-samples[samples$Monkey==monkey,]
    Ov<-OvDF[sample$filename]    
    
    res<-sample[,c("Monkey","tis","Week","Cohort","filename")]
    
    for (f in 1:nrow(res)){
        dat<-Ov[[res$filename[f]]]
        for (j in 1:nrow(df)){
            res[f,(j+5)]<-dat[dat$ref251.pos==df$ref251.pos[j],paste0(df$MFtype[j])]
        }
    }
    colnames(res)[6:ncol(res)]<-paste0("pos.",df$ref251.pos,".",df$MFtype)
    
    Results[[i]]<-res
    names(Results)[i]<-monkey
}



for (i in 1:length(Results)){
    freq<-Results[[i]]
    freq$tis<-factor(freq$tis, c("plasma","plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    monkey<-names(Results)[i]
    
    Plots<-list()
    for (c in 1:(ncol(freq)-5)){
        df<-freq[,c(1:5,(c+5))]
        colnames(df)[6]<-"Freq"
        Plots[[c]]<-ggplot(df, aes(x=tis, y=Freq*100))+
            geom_point(color="blue")+
            #scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100),  limits=c(0,1))+
            theme_bw()+
            ggtitle(monkey)+
            theme(plot.title = element_text(size=11))+
            xlab('')+ylab(colnames(freq)[c+5])+
            theme(axis.text.x = element_text(angle =45, hjust=1))   
        
    }
    
    pdf(paste0("Output/MF/TrackHighMF_",monkey,".pdf"), width = 8, height=ceiling((ncol(freq)-5)/2)*2)
    do.call(grid.arrange, c(Plots, ncol=2))
    dev.off()
    
}

Diff<-list()
for (i in 1:length(siv)){
    freq<-Results[[siv[i]]]
    freq<-freq[freq$tis!="plasma",]
    freq$tis<-factor(freq$tis, c("plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    
    diff<-c()
    for (j in 1:(ncol(freq)-5)){
        vec<-freq[,(j+5)]
        if (vec[1]>max(vec[2:length(vec)])){
            diff<-c(diff, vec[1]-min(vec[2:length(vec)]))
        }
    }
    Diff[[i]]<-diff
    names(Diff)[i]<-siv[i]
    
}

lapply(Diff, function(x) max(x))

difference<-unlist(Diff) #40 mutations out of 46 showed bottleneck 
median(difference)
hist(difference)

diff.large<-difference[difference<0.01] 
#19 mutations out of 40 < 1% diff
#27 are less than 1.5% difference


# Latent R
siv2<-unique(samples$Monkey[samples$Cohort=="Latent R"])
Diff2<-list()
for (i in 1:length(siv2)){
    freq<-Results[[siv2[i]]]
    freq<-freq[freq$tis!="plasma",]
    freq$tis<-factor(freq$tis, c("plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    
    diff<-c()
    for (j in 1:(ncol(freq)-5)){
        vec<-freq[,(j+5)]
        vec<-vec[!is.na(vec)]
        if (vec[1]>max(vec[2:length(vec)])){
            diff<-c(diff, vec[1]-min(vec[2:length(vec)]))
        }
    }
    Diff2[[i]]<-diff
    names(Diff2)[i]<-siv[i]
}

difference2<-unlist(Diff2) # total 25 sites bottleneck out of 43 sites
median(difference2) #25 sites  #0.009646302
hist(difference2)

diff.large<-difference2[difference2<0.01] 
#13 sites under 1%
diff.large<-difference2[difference2<0.015] 
#16 sites under 1.5%

mean(difference2) #0.03154329


## Latent R
siv3<-unique(samples$Monkey[samples$Cohort=="Latent NR"])
Diff3<-list()
for (i in 1:length(siv3)){
    freq<-Results[[siv3[i]]]
    freq<-freq[freq$tis!="plasma",]
    freq$tis<-factor(freq$tis, c("plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    
    diff<-c()
    for (j in 1:(ncol(freq)-5)){
        vec<-freq[,(j+5)]
        vec<-vec[!is.na(vec)]
        if (vec[1]>max(vec[2:length(vec)])){
            diff<-c(diff, vec[1]-min(vec[2:length(vec)]))
        }
    }
    Diff3[[i]]<-diff
    names(Diff3)[i]<-siv[i]
}

difference3<-unlist(Diff3) # total 5 sites bottleneck out of 13 sites
median(difference3) #25 sites  #0.009065534
hist(difference3)

length(difference3[difference3<0.01]) 
#4 sites under 1%
length(difference3[difference3<0.015])
#4 sites under 1.5%



#remove the si
hist(diff)
median(diff)

apply(freq[,6:ncol(freq)], 2, function(x){ if (max(x)==max(x[x$tis=="plasma_nex"]) })
    freq[,6:ncol(freq)]
    
    d<-apply(freq[,6:ncol(freq)], 2, function(x) {max(x)-min(x)})
    diff<-c(diff,d)


difference[[i]]
    freq$tis<-factor(freq$tis, c("plasma","plasma_nex", "pLN","Thoracic LN", "Lung"))
    freq<-freq[order(freq$tis),]
    monkey<-names(Results)[i]
    
    Plots<-list()
    