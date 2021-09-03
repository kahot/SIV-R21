#Diversity comparison
library(ggplot2)
library(reshape)
library(gridExtra)
source("Rscripts/label_scientific.R")

# read the files saved in Overview_output:
ItaFiles<-list.files("Output/Overview/Ita/",pattern="overview.csv")

ItaOv<-list()
for (i in 1:length(ItaFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/Ita/",ItaFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    ItaOv[[i]]<-overviews
    names(ItaOv)[i]<-gsub("_overview.csv",'',ItaFiles[i])
}


#sample info
SampleIta<-read.csv("Data/Ita.SampleInfo.csv", stringsAsFactors =F)
stock<-SampleIta[SampleIta$ID=="Stock",]
samples<-SampleIta[SampleIta$ID!="Stock",]

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


#calculate average mutation frequency (minor variant freq) for all files
Summary<-list()
cutoff<-0.009

for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-ItaOv[sample$SampleID]    
    
    summary<-sample[,c("Monkey","Week","SampleID")]
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #remove the sites with row read depth
        df[which(df$TotalReads<100),18:27]<-NA
        df<-df[df$AA251pos<276,]
        
        #calculate total ns freq at each position
        df$ns1<-as.numeric(apply(df[,c("Type","freq.Ts")],1, function(x) if (is.na(x["Type"])) NA else if (x["Type"]=="nonsyn") x["freq.Ts"] else 0))
        df$ns2<-as.numeric(apply(df[,c("Type.tv1","freq.transv1")],1, function(x) if (is.na(x["Type.tv1"])) NA else if (x["Type.tv1"]=="nonsyn") x=x["freq.transv1"] else 0))
        df$ns3<-as.numeric(apply(df[,c("Type.tv2","freq.transv2")],1, function(x) if (is.na(x["Type.tv2"])) NA else if (x["Type.tv2"]=="nonsyn") x=x["freq.transv2"] else 0))
        df$ns<-df$ns1+df$ns2+df$ns3
        #calculate total syn freq at each position
        df$syn1<-as.numeric(apply(df[,c("Type","freq.Ts")],1, function(x) if (is.na(x["Type"])) NA else if (x["Type"]=="syn") x["freq.Ts"] else 0))
        df$syn2<-as.numeric(apply(df[,c("Type.tv1","freq.transv1")],1, function(x) if (is.na(x["Type.tv1"])) NA else if (x["Type.tv1"]=="syn") x=x["freq.transv1"] else 0))
        df$syn3<-as.numeric(apply(df[,c("Type.tv2","freq.transv2")],1, function(x) if (is.na(x["Type.tv2"])) NA else if (x["Type.tv2"]=="syn") x=x["freq.transv2"] else 0))
        df$syn<-df$syn1+df$syn2+df$syn3
        #write.csv(df,paste0("Output/Overview2/Ita/", names(Ov[i]),".csv"))

        df2<-df[,c("AA239pos","ref251.pos","ns")]
        colnames(df2)[3]<-"freq"
        df2$Type<-"nonsyn"
        df2<-df2[!is.na(df2$freq),]
        
        df3<-df[,c("AA239pos","ref251.pos","syn")]
        colnames(df3)[3]<-"freq"
        df3$Type<-"syn"
        df3<-df3[!is.na(df3$freq),]
        
        mut<-rbind(df2,df3)

        missing<-df[is.na(df$freq.mutations),]
        missing<-missing[missing$ref251.pos>=95,]
        missing<-missing[missing$AA251pos>(min(missing$AA251pos)+1)&missing$AA251pos<max(missing$AA251pos),]
        

        #remove the center missing sections for average calculation
        mut<-mut[mut$ref251.pos<395|mut$ref251.pos>540,]
        alt<-df[df$MajNt!=df$Ita.stock,]
        alt<-alt[!is.na(alt$freq.mutations),]
        
        summary$no.sites[i]<-nrow(mut)
        summary$no.syn.sites[i]<-nrow(mut[mut$Type=="syn",])
        summary$no.ns.sites[i]<-nrow(mut[mut$Type=="nonsyn",])
        summary$mean[i]<-mean(df$freq.mutations[df$ref251.pos<395|df$ref251.pos>540], na.rm=T)
        summary$mean.syn[i]<-mean(mut$freq[mut$Type=="syn"], na.rm=T)
        summary$mean.ns[i]<-mean(mut$freq[mut$Type=="nonsyn"], na.rm=T)
        summary$mutated.sites.no[i]<-nrow(alt)
        if (nrow(alt)!=0) {summary$mutated.sites[i]<-paste(alt$ref251.pos, collapse = ',')}
        else summary$mutated.sites[i]<- NA
    }
    Summary[[j]]<-summary
    
}

Sum.divergence<-do.call(rbind,Summary)    
write.csv(Sum.divergence,"Output/Diversity/Diversity_summary_Ita.csv")
