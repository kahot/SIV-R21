library(ggplot2)
library(reshape2)
library(gridExtra)

#library(plotly)
#library(tidyverse)
#library(htmlwidgets)

library(colorspace)
cols<-qualitative_hcl(6, palette="Dark3")

#source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}


SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stocks<-SampleSheet[SampleSheet$Monkey=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]
which(SampleSheet$Monkey=="stock_virus") #25,49
SampleSheet$filename[SampleSheet$Monkey=="stock_virus"]

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



#calculate the % indels from ref251 -> indels are rare
for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF[sample$filename]    
    
    Plots<-list()
    #where do the indels occur?
    Ins<-data.frame(ref251.pos=94:834)
    Del<-data.frame(ref251.pos=94:834)
    Deletion<-list()
    Indels<-list()
    for (i in 1:length(Ov)){
        DF<-Ov[[i]]
        DF$del.percent<-DF$deletion/DF$TotalReads
        DF$ins.percent<-DF$insertion/DF$TotalReads
        
        df2<-DF[,c("ref251.pos","del.percent","ins.percent")]
        colnames(df2)[2:3]<-c("deletion","insertion")
        dfm<-melt(df2, id.vars="ref251.pos")
        dfm[dfm==0]<-NA
        id<-paste0("Animal:", monkey," Week",sample$Week[i], ": ",sample$Sample[i])
        p<-ggplot(dfm, aes(x=ref251.pos, y=value, color=variable))+
            geom_point(size=.5)+
            scale_color_manual(values=cols[c(2,5)])+
            facet_grid(rows = vars(variable))+
            theme_bw()+ylab('Indel observed')+
            xlab('Env position')+ylim(0,1)+
            guides(color = guide_legend(title = NULL))+
            ggtitle(id)+theme(legend.position = "none") 
        Plots[[i]]<-p
        Indels[[i]]<-df2
        names(Indels)[[i]]<-id
            
        Ins[,id]<-df2[,"insertion"]
        Del[,id]<-df2[,"deletion"]
    }     
    
    pdf(paste0("Output/Indels/Indel_",monkey,".pdf"), width=8, height = length(Ov)*2)
    do.call(grid.arrange, c(Plots, ncol=2))
    dev.off()
    n<-ncol(Ins)
    Ins$mean<-rowMeans(Ins[,2:n], na.rm=T)
    Ins$count<-rowSums(Ins[,2:n]>0)
    Del$mean<-rowMeans(Del[,2:n], na.rm=T)
    Del$count<-rowSums(Del[,2:n]>0)
    
    write.csv(Ins, paste0("Output/Indels/Insertion_", monkey, ".csv"))  
    write.csv(Del, paste0("Output/Indels/Deletion_", monkey, ".csv"))  
}
        
#look at the insertion sites compared to ref239
#indel sites between codon 126 - 127 (6 bases)
# codon 141-142 (3 bases)
#Is there any files with different indel patterns? -> no
Deletion<-list()
Insertion<-list()
ins.files<-c()
del.files<-c()
k=1
g=1
pdf("Output/Indels/insertion_reads.pdf", width = 10, height = 40)
par(mfrow=c(15,5))

for (i in 1:length(OvDF)){
    print(names(OvDF[i]))
    DF<-OvDF[[i]]
#    Plots<-list()
    dfI<-DF[is.na(DF$ref239.pos)&!is.na(DF[,2]),]
    plot(dfI$TotalReads, main=names(OvDF[i]), pch=16, xlab='', ylab='# of reads', ylim=c(0,max(dfI$TotalReads)))

    if (length(is.na(dfI[,2]))!=9) ins.files<-c(ins.files, names(OvDF)[i])
    if (nrow(dfI)>0) {
        dfI2<-dfI[dfI$freq.mutations.ref>0.05,]
        dfI2<-dfI2[!is.na(dfI2$merged.pos),]
        if (nrow(dfI2)>0) {
            Insertion[[k]]<-dfI2
            names(Insertion)[k]<-names(OvDF)[i]
            k=k+1
        }
    }

    dfD<-DF[is.na(DF[,2]),]
    if (nrow(dfD)!=0) {
        del.files<-c(del.files, names(OvDF)[i])
        Deletion[[g]]<-dfD
        names(Deletion)[g]<-names(OvDF)[i]
        g=g+1
    }
}
dev.off()     
 
deletion<-data.frame(ID<-names(Deletion))
for (i in 1:length(Deletion)){
    d<-Deletion[[i]]
    deletion$nrow[i]<-nrow(d)
    
    deletion$positions[i]<-paste(d$ref251.pos, collapse = ',')
}
   


#24 files have less insertion than the 9 bases
#deletions are all at the same positions 430,431,432



##calculate indels for stock virus
DF<-OvDF[[stock$filename]]
    
    #where do the indels occur?
Ins<-data.frame(ref251.pos=94:834)
Del<-data.frame(ref251.pos=94:834)
Deletion<-list()
Indels<-list()

DF$del.percent<-DF$deletion/DF$TotalReads
DF$ins.percent<-DF$insertion/DF$TotalReads

df2<-DF[,c("ref251.pos","del.percent","ins.percent")]
colnames(df2)[2:3]<-c("deletion","insertion")
dfm<-melt(df2, id.vars="ref251.pos")
dfm[dfm==0]<-NA
id<-"Stock"
ggplot(dfm, aes(x=ref251.pos, y=value, color=variable))+
    geom_point(size=.5)+
    scale_color_manual(values=cols[c(2,5)])+
    facet_grid(rows = vars(variable))+
    theme_bw()+ylab('Indel observed')+
    xlab('Env position')+ylim(0,1)+
    guides(color = guide_legend(title = NULL))+
    ggtitle(id)+theme(legend.position = "none") 
ggsave(paste0("Output/Indels/Indel_",id,".pdf"), width=4, height = 2)

write.csv(df2,"Output/Indels/Indels_stockvirus.csv")

       