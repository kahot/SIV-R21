#Diversity comparison of all samples, remove nt pos 395-540

library(ggplot2)
library(reshape)
library(gridExtra)
library(DataCombine)
library(dplyr)
source("Rscripts/label_scientific.R")
library(colorspace)

colors<-qualitative_hcl(6, palette="Dark3")
col_light<-qualitative_hcl(6, palette="Set3")

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


#calculate average MINOR VARIANT frequency (minor variant freq) for all files
summary<-SampleSheet[,c("Cohort", "Monkey","Tissue2","Week","filename")]
for (i in 1:length(OvDF)){
   
        df<-OvDF[[SampleSheet$filename[i]]]
        #remove the sites with row read depth
        df[which(df$TotalReads<100),17:26]<-NA
        df<-df[df$AA251pos<276,]
        
        #remove the center missing region
        df<-df[df$ref251.pos<395|df$ref251.pos>540,]
        
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
        #Save the overview with addtional info as Overview 2
        write.csv(df,paste0("Output/Overview2/", names(Ov[i]),".csv"))
                
        df2<-df[,c("AA251pos","ref251.pos","ns")]
        colnames(df2)[3]<-"freq"
        df2$Type<-"nonsyn"
        df2<-df2[!is.na(df2$freq),]
        
        df3<-df[,c("AA251pos","ref251.pos","syn")]
        colnames(df3)[3]<-"freq"
        df3$Type<-"syn"
        df3<-df3[!is.na(df3$freq),]
        
        mut<-rbind(df2,df3)
          
        #calculate how many sites with majority different from reference nt
        alt<-df[df$MajNt!=df$ref,]
        alt<-alt[!is.na(alt$freq.mutations),]
        
        summary$no.sites[i]<-nrow(mut)
        summary$no.syn.sites[i]<-nrow(mut[mut$Type=="syn",])
        summary$no.ns.sites[i]<-nrow(mut[mut$Type=="nonsyn",])
        summary$mean[i]<-mean(df$freq.mutations[df$ref251.pos<397|df$ref251.pos>540], na.rm=T)
        summary$mean.syn[i]<-mean(mut$freq[mut$Type=="syn"], na.rm=T)
        summary$mean.ns[i]<-mean(mut$freq[mut$Type=="nonsyn"], na.rm=T)
        summary$mutated.sites.no[i]<-nrow(alt)
        if (nrow(alt)!=0) {summary$mutated.sites[i]<-paste(alt$ref251.pos, collapse = ',')}
        else summary$mutated.sites[i]<- NA
        
          
}

write.csv(summary,"Output/Diversity_summary_R21.csv")



