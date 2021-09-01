# Look for known AA substitutions reported in Ita et al (2018) Table 1.

library(ggplot2)
library(reshape2)
library(gridExtra)


# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}

OverviewFiles2<-list.files("Output/Overview2/",pattern=".csv")
OvDF2<-list()
for (i in 1:length(OverviewFiles2)){ 
    overviews<-read.csv(paste0("Output/Overview2/",OverviewFiles2[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF2[[i]]<-overviews
    names(OvDF2)[i]<-gsub(".csv",'',OverviewFiles2[i])
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
df<-df[!(df$Type.ref=="syn"&df$Type.tv1.ref=="syn"&df$Type.tv2.ref=="syn"),]

mutations<-c()
for (j in 1:nrow(df)){
    if (df$Type.ref[j]!="syn"){
        mut<-paste0(df$WTAA.ref[j],df$AA239pos[j],df$MUTAA.ref[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv1.ref[j]!="syn"){
        mut<-paste0(df$WTAA.ref[j],df$AA239pos[j],df$TVS1_AA.ref[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv2.ref[j]!="syn"){
        mut<-paste0(df$WTAA.ref[j],df$AA239pos[j],df$TVS2_AA.ref[j])
        mutations<-c(mutations,mut)
    }        
}


muts<-data.frame(matrix(ncol=length(mutations)+1, nrow=length(OvDF)), stringsAsFactors = F) 
colnames(muts)<-c("Sample", mutations)
muts$Sample<-names(OvDF)

for (i in 1:length(OvDF)){
    df<-OvDF[[i]]
    df<-df[df$ref239.pos %in% varNT,]
    df<-df[!(df$Type.ref=="syn"&df$Type.tv1.ref=="syn"&df$Type.tv2.ref=="syn"),]
    df$AA.pos<-ceiling(df$ref239.pos/3)

    mut<-c()
    for (j in 1:nrow(df)){
         if (df$Type.ref[j]!="syn"){
             if (is.na(df$TotalReads[j])|df$TotalReads[j]<100) mut<-c(mut, NA)
             else mut<-c(mut,df$freq.Ts.ref[j])
        }
        if (df$Type.tv1.ref[j]!="syn"){
            if (is.na(df$TotalReads[j])|df$TotalReads[j]<100) mut<-c(mut, NA)
            else mut<-c(mut, df$freq.transv1.ref[j])
        }
        if (df$Type.tv2.ref[j]!="syn"){
            if (is.na(df$TotalReads[j])|df$TotalReads[j]<100) mut<-c(mut, NA)
            else mut<-c(mut, df$freq.transv2.ref[j])
        }
    }
    muts[i,2:ncol(muts)]<-mut
}


SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
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
info<-SampleSheet[,c("filename","Monkey","Tissue3","Week", "Cohort")]
info$Monkey[info$filename=="Run_6_01_Animal_stock_virus"]<-"control"
info$Tissue2[info$filename=="Run_6_01_Animal_stock_virus"]<-"control"
info$Monkey [info$filename=="Run_5_01_Animal_stock_virus"]<-"stock"
info$Tissue2[info$filename=="Run_5_01_Animal_stock_virus"]<-"stock"
info$Cohort[70:71]<-c("Stock","Control")


muts2<-muts
colnames(muts2)[1]<-"filename"
table1<-merge(info, muts2,by="filename", all = T)
write.csv(table1,"Output/Table1/AAchanges.summary.csv")


Plots<-list()    
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[i]
    
    freq<-muts2[muts2$filename %in% sample$filename,]
    colnames(freq)[1]<-"ID"
    freq$Sample<-paste0("Week",sample$Week,"_",sample$Tissue2)
    freq<-freq[,-1]
    freqM<-melt(freq, id.vars="Sample")
    colnames(freqM)[2:3]<-c("Mutation","Freq")
    
    # find the frequncy higher than 0.006
    highFs<-freqM$Mutation[freqM$Freq>=0.006]
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
    ggsave(paste0("Output/Table1/",monkey,".pdf"), width = 12, height = 3.5)
}    


pdf("Output/Table1/All_Animals.pdf", width = 12, height=35)
do.call(grid.arrange, c(Plots, ncol=1))
dev.off()

Table1<-data.frame()
for (i in 1:length(mutations)){
    df<-table1[,c(1:5,(i+5))]s
    high<- which(df[,6]>0.005)
    if (length(high)==0) next
    else{
        df2<-df[which(df[,6]>0.005),]
        df2$AAsub<-colnames(df2)[6]
        colnames(df2)[6]<-"Freq"
        df2$Stock<-df[df$Monkey=="stock",6]
        Table1<-rbind(Table1, df2)
    }
}    
    
write.csv(Table1,"Output/Table1/Observed_frq_table1.csv")


#Create summary of Table 1

Table1<-Table1[Table1$Monkey!="control",]
Table1$Pos<-substring(Table1$AAsub,2,4)

aasubs<-unique(Table1$AAsub)

sum1<-data.frame(AA.sub=aasubs)
sum1$AA.pos<-substring(sum1$AA.sub,2,4)
sum1<-sum1[,c(2,1)]

for (i in 1:nrow(sum1)){
    df<-Table1[Table1$AAsub==aasubs[i],]
    #how many samples?
    sum1$count[i]<-nrow(df)
    
    #which tissues?
    tiss<-paste(unique(df$Tissue2), collapse = ", ")
    sum1$tissue[i]<-tiss
    #which cohort?
    sum1$Cohort[i]<-paste(unique(df$Cohort), collapse = ",")
    
    #Frequency
    sum1$min.freq[i]<-min(df$Freq)
    sum1$max.freq[i]<-max(df$Freq)
    sum1$mean.freq[i]<-mean(df$Freq)
    sum1$stock[i]<-Table1$Stock[Table1$AAsub==aasubs[i]]
    }


sum1$range<-paste0("(",round(sum1$min.freq, digits=4), " - ",round(sum1$max.freq, digits=4),")") 
sum1$freq_percent<-sum1$mean.freq*100
write.csv(sum1,"Output/Table1/Table1_Summary.csv")   


#how many % of samples had freq. higher than 0.5%?
df<-table1[table1$Tissue2!="stock",]
df<-df[df$Tissue2!="control",]

prop<-data.frame(AA.sub=aasubs)
for (i in 1:length(aasubs)){
    n<-length(!is.na(df[,aasubs[i]]))
    k<-sum1$count[sum1$AA.sub==aasubs[i]]
    prop$percent[i]<-k/n*100
}




