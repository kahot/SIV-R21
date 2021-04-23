#Look at the mutation frequency at AA112 in each animal

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


######################
## AA112
AA112<-list()
for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF[sample$filename]    
    
    #Eliminate sites with reads<100
    summ1<-sample[,c("Monkey","Tissue2","Week","Cohort","filename")]
    summ1$sample<-paste0("Week",summ1$Week,"_",summ1$Tissue2)
    
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        if (is.na(df$TotalReads[df$ref251.pos==335])|df$TotalReads[df$ref251.pos==335]<100){
            summ1$AA112[i]<-NA
        }
        
        else{
            summ1$AA112[i]<-df$a[df$ref251.pos==335]/df$TotalReads[df$ref251.pos==335]
        }
    }
    AA112[[j]]<-summ1
    
    summ1$sample<-factor(summ1$sample, levels=unique(summ1$sample))
    ggplot(summ1, aes(x=sample, y=AA112))+
        geom_point(color="blue")+
        scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100),  limits=c(0,1))+
        theme_bw()+
        ggtitle(monkey)+
        theme(plot.title = element_text(size=11))+
        xlab('')+ylab("Freq. of A at AA112 (nt335)")+
        theme(axis.text.x = element_text(angle =45, hjust=1))   
    ggsave(paste0("Output/MF/AA112/AA112_",monkey,".pdf"), height = 3,width =4.5 )
    }

aa112<-do.call(rbind,AA112)
write.csv(aa112,"Output/MF/AA112.csv")

aa112<-read.csv("Output/MF/AA112.csv", row.names = 1, stringsAsFactors = F)
aa112<-aa112[!is.na(aa112$AA112),]

aggregate(aa112[,c("AA112")], by=list(aa112$Cohort), mean, na.omit=T)
#    Group.1           x
#1 Latent NR 0.209514391
#2  Latent R 0.063485557
#3  SIV only 0.001152281

# Groups
gp1<-c(20615,30816)

aa<-aa112
aa$Group<-"No sweep"
aa$Group[aa$Monkey %in% gp1]<-"Sustained AA sub"

aa$Tissue2[aa$Tissue2=="plasma"]<-"Plasma"
aa$Tissue2[aa$Tissue2=="Plasma"&aa$Week>5]<-"Plasma PM"
aa$Tissue2[aa$Tissue2=="pLN"]<-"Perif LN"
aa$Group<-factor(aa$Group, levels=c("Sustained AA sub","No sweep"))

#aa130 in stock?
sk<-OvDF[[stock$filename]]
sk<-sk[sk$ref251.pos==335,]
f<-sk$a/sk$TotalReads

aa<-aa[,c("Monkey","Tissue2","Group","AA112")]
s<-data.frame(Monkey=rep("Stock", times=2), Tissue2=rep("Stock", times=2),
              Group=c("Sustained AA sub","No sweep"),
              AA112=f)
aa<-rbind(aa,s)
aa$Tissue2<-factor(aa$Tissue2, levels=c("Stock","Plasma","Plasma PM", "Perif LN","Thoracic LN", "Lung"))

m<-read.csv("Output/monkey_order.csv", stringsAsFactors = F)
morder<-m$x
aa$Monkey<-factor(aa$Monkey, levels=morder)
aa<-aa[order(aa$Monkey),]

ggplot(aa, aes(x=Tissue2, y=AA112, color=Monkey))+
    geom_point(position=position_dodge(width = 0.3))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
    theme_bw()+xlab('')+ylab("Freq. of A at AA112 (nt335)")+
    theme(axis.text.x = element_text(angle =45, hjust=1))+
    facet_wrap(~Group, ncol=2)+theme(legend.title = element_blank())
ggsave("Output/Figures/AA112_summary.png", width = 8, height = 3.5, dpi=300, unit="in")

