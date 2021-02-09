library(ggplot2)
library(reshape)
library(gridExtra)

#library(plotrix)
#library(grid)
#library(tidyverse)
#source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}


#Variable positions in V1 loop from Table 1
varAA<-c(120,132,135,136,138,139,198,201,202)


#nt positions
varNT<-c()
for (i in 1: length(varAA)){
    pos<-c(varAA[i]*3-2,varAA[i]*3-1,varAA[i]*3)
    varNT<-c(varNT,pos)
}

#List all available mutations 
df<-OvDF[[1]]
df<-df[df$ref239.pos %in% varNT,]
df<-df[!(df$Type=="syn"&df$Type.tv1=="syn"&df$Type.tv2=="syn"),]
df$AA.pos<-ceiling(df$ref239.pos/3)

mutations<-c()
for (j in 1:nrow(df)){
    if (df$Type[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$MUTAA[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv1[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$TVS1_AA[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv2[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$TVS2_AA[j])
        mutations<-c(mutations,mut)
    }        
}


#Table1<-list() 

muts<-data.frame(matrix(ncol=length(mutations)+1, nrow=length(OvDF)), stringsAsFactors = F) 
colnames(muts)<-c("Sample", mutations)
muts$Sample<-names(OvDF)

for (i in 1:length(OvDF)){
    df<-OvDF[[i]]
    df<-df[df$ref239.pos %in% varNT,]
    df<-df[!(df$Type=="syn"&df$Type.tv1=="syn"&df$Type.tv2=="syn"),]
    df$AA.pos<-ceiling(df$ref239.pos/3)
    
    
    mut<-c()
    for (j in 1:nrow(df)){
        if (df$Type[j]!="syn"){
            mut<-c(mut,df$freq.Ts.ref[j])
        }
        if (df$Type.tv1[j]!="syn"){
            mut<-c(mut, df$freq.transv1.ref[j])
        }
        if (df$Type.tv2[j]!="syn"){
            mut<-c(mut, df$freq.transv2.ref[j])
        }
    }
    muts[i,2:ncol(muts)]<-mut
}


SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stocks<-SampleSheet[SampleSheet$Monkey=="stock_virus",]
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

#save data frame with monkey info.
info<-SampleSheet[,c("filename","Monkey","Sample","Week")]
colnames(info)[1]<-"ID"
muts2<-muts
colnames(muts2)[1]<-"ID"
table1<-merge(info, muts2,by="ID", all.y = T)
write.csv(table1,"Output/Table1/AAchanges.summary.csv")




Plots<-list()    
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[i]
    
    freq<-muts[muts$Sample %in% sample$filename,]
    colnames(freq)[1]<-"ID"
    freq$Sample<-paste0("Animal_",monkey,"_wk",sample$Week,"_",sample$Sample)
    freq<-freq[,-1]
    freqM<-melt(freq, id.vars="Sample")
    colnames(freqM)[2:3]<-c("Mutation","Freq")
    
    highFs<-freqM$Mutation[freqM$Freq>=0.005]
    highF<-unique(as.character(highFs[!is.na(highFs)]))
    
    mutID<-levels(freqM$Mutation)
    colors<-sapply(mutID, function(x) if (x%in%highF) x=2 else 1)
    Plots[[i]]<-ggplot(freqM, aes(x=Mutation, y=Freq, color=Sample))+
        geom_point(position=position_dodge(width=0.8))+
        theme_bw()+ylim(0,0.4)+
        theme(axis.text.x=element_text(angle=90))+
        theme(axis.text.x=element_text(color=colors))+
        theme(
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)
        )+
        ggtitle(paste0("Animal ",monkey))+
        geom_hline(yintercept = 0.005, color="lightsteelblue3")
    #ggsave(paste0("Output/Table1/",monkey,".pdf"), width = 12, height = 3.5)
}    


pdf("Output/Table1/All_animals_2.pdf", width = 12, height=35)
do.call(grid.arrange, c(Plots, ncol=1))
dev.off()
