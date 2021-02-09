#Script to analyse the frequency data and associate with features
library(dplyr)
source("Rscripts/baseRscript.R")

#get the file name
SIVFiles_SeqData<-list.files("Output/SeqData/Ita/",pattern="SeqData")

#create the type of mutations info for Sergio's stock virus 
RefDF<-read.csv(paste0("Output/SeqData/Ita/SeqData_mmStock.csv"),stringsAsFactors=FALSE, row.names = 1)

TypeOfSite<-c() 
TypeOfSite.tv1<-c()
TypeOfSite.tv2<-c()
for (codon in 1:(nrow(RefDF)/3)) { 
    positions <- c(codon*3-2,codon*3-1, codon*3)  
    WTcodon <- RefDF$Ita.stock[positions]  
    mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
    mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
    mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
    
    #transversion mutation to 'a' or 'c'
    mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
    mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
    mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
    #transversion mutation to 'g' or 't'
    mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
    mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
    mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
    
    #       }
    
    
    TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
    TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
    TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
    
    TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant1codon.tv1))
    TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant2codon.tv1))
    TypeOfSite.tv1<-c(TypeOfSite.tv1,typeofsitefunction(WTcodon,mutant3codon.tv1))
    
    TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant1codon.tv2))
    TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant2codon.tv2))
    TypeOfSite.tv2<-c(TypeOfSite.tv2,typeofsitefunction(WTcodon,mutant3codon.tv2))
    
    
} 

RefDF$Type<-TypeOfSite
RefDF$Type.tv1<-TypeOfSite.tv1
RefDF$Type.tv2<-TypeOfSite.tv2

for (k in 1:nrow(RefDF)){
    if (k%%3==1){
        RefDF$WTAA[k] = seqinr::translate(RefDF$Ita.stock[c(k,k+1,k+2)])
        RefDF$MUTAA[k] =seqinr::translate(c(transition(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
        RefDF$TVS1_AA[k] = seqinr::translate(c(transv1(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
        RefDF$TVS2_AA[k] = seqinr::translate(c(transv2(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
    } 
    if (k%%3==2){
        RefDF$WTAA[k]   = seqinr::translate(  RefDF$Ita.stock[c(k-1,k,k+1)])
        RefDF$MUTAA[k]  = seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transition(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
        RefDF$TVS1_AA[k]= seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transv1(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
        RefDF$TVS2_AA[k]= seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transv2(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
    }
    if (k%%3==0){
        RefDF$WTAA[k] =    seqinr::translate(  RefDF$Ita.stock[c(k-2,k-1,k)])
        RefDF$MUTAA[k] =   seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transition(RefDF$Ita.stock[k])))
        RefDF$TVS1_AA[k] = seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transv1(RefDF$Ita.stock[k])))
        RefDF$TVS2_AA[k] = seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transv2(RefDF$Ita.stock[k])))
        }
}
#Add whether AA change is drastic & makes CpG
RefDF$bigAAChange<-0
RefDF$bigAAChange.tv1<-0
RefDF$bigAAChange.tv2<-0
RefDF$makesCpG <- 0
RefDF$makesCpG.tv1 <- 0
RefDF$makesCpG.tv2 <- 0

for(j in 1:nrow(RefDF)){
    WT <- amCat(RefDF[j,'WTAA'])
    MUT<- amCat(RefDF[j,'MUTAA'])
    MUT1<-amCat(RefDF[j,'TVS1_AA'])
    MUT2<-amCat(RefDF[j,'TVS2_AA'])
    
    if (WT != MUT) RefDF$bigAAChange[j] <- 1
    if (WT != MUT1) RefDF$bigAAChange.tv1[j] <- 1
    if (WT != MUT2) RefDF$bigAAChange.tv2[j] <- 1
    
    trip <- RefDF$Ita.stock[c(j, j+1,j+2)]
    if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) next
    else {
        if (trip[1] == "c" & trip[2] == "a" ) RefDF$makesCpG[j] <- 1 
        if (trip[2] == "t" & trip[3] == "g")  RefDF$makesCpG[j] <- 1
        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) RefDF$makesCpG.tv2[j] <- 1                                
        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) RefDF$makesCpG.tv1[j] <- 1
        
    }
} 

RefDF<-RefDF[-c(1,7:21)]
RefDF$AA251pos<-ceiling(RefDF$ref251.pos/3)

write.csv(RefDF, "Output/Ita.stock_overview.csv")

#Attach the information to each SeqData file
#dir.create("Output/Overview/Ita")

for (i in 1:length(SIVFiles_SeqData)){   
        id<-gsub(".csv",'',paste(SIVFiles_SeqData[i]))
        id<-gsub("SeqData_",'',id)
        DF<-read.csv(paste0("Output/SeqData/Ita/",SIVFiles_SeqData[i]),stringsAsFactors=FALSE, row.names = 1)
                
        DF<-merge(DF,RefDF[c(1,3,6:19)], by="ref251.pos")
        DF<-DF[DF$ref251.pos>=94,]
        write.csv(DF,paste0("Output/Overview/Ita/",id,"_overview.csv"))
        print(id)
}        
#####################################################

### Read depths for all files ###
SIVFiles_SeqData<-list.files("Output/SeqData/Ita/",pattern="SeqData")

ReadsSummary<-data.frame(SampleID=matrix(nrow=length(SIVFiles_SeqData)))
ReadsSummary$MaxDepth<-""
ReadsSummary$AveDepth<-""

for (i in 1:length(SIVFiles_SeqData)){
        id<-gsub(".csv",'',paste(SIVFiles_SeqData[i]))
        id<-gsub("SeqData_",'',id)
        ReadsSummary$SampleID[i]<-id
        print(id)
        SeqData<-read.csv(paste("Output/SeqData/Ita/",SIVFiles_SeqData[i],sep=""))
        ReadsSummary$MaxDepth[i]<-max(SeqData$TotalReads,na.rm=T)
        ReadsSummary$AveDepth[i]<-mean(SeqData$TotalReads,na.rm=T)
        ReadsSummary$No.ofSites[i]<-nrow(SeqData[!is.na(SeqData$Ita.stock),])
        ReadsSummary$SE[i]<-std.error(SeqData$TotalReads, na.rm=T)
}

write.csv(ReadsSummary,"Output/ReadsSummary_Ita.csv")      


###
ReadsSummary<-read.csv("Output/ReadsSummary_Ita.csv", row.names = 1)
ggplot(ReadsSummary, aes(x=SampleID, y=AveDepth))+
    geom_point(size=0.8)+
    theme(axis.text.x=element_text(angle=90))+
    scale_y_continuous(labels = label_scientific())
ggsave("Output/average.readdepth_Ita.pdf", width = 8, height = 6)    

#plot together with R21 study data

ReadsSummary2<-read.csv("Output/ReadsSummary.csv", row.names = 1)

reads<-rbind(ReadsSummary2,ReadsSummary)
ggplot(reads, aes(x=SampleID, y=AveDepth))+
    geom_point(size=0.8)+
    theme(axis.text.x=element_text(angle=90))+
    scale_y_continuous(labels = label_scientific())
ggsave("Output/average.readdepth_R21andIta.pdf", width = 8, height = 6)    


#Create the average mf data

OverviewFiles<-list.files("Output/Overview/Ita/",pattern="overview.csv")
#remove the stock and control
#OverviewFiles<-OverviewFiles[-c(49,25)]

MF<-data.frame(sample=1:length(OverviewFiles))
for (i in 1:length(OverviewFiles)){ 
    df<-read.csv(paste0("Output/Overview/Ita/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    df<-df[df$TotalReads>=1000,]
    MF$SampleID[i]<-gsub("_overview.csv",'',OverviewFiles[i])
    MF$ave.mf[i]<-mean(df$freq.mutations.ref, na.rm=T)
}


Summary<-merge(MF, ReadsSummary[,c("SampleID","AveDepth")], by="SampleID")
write.csv(Summary,"Output/AveMF_Deapth_Ita.summary.csv")

ggplot(Summary,aes(x=ave.mf, y=AveDepth))+
    geom_point()+
    xlab("Average mutation freq")+ylab("Average read depth")
ggsave("Output/depth.vs.mf_Ita.pdf", width = 4, height = 4)

cor.test(Summary$ave.mf, Summary$AveDepth, method="spearman")
#p-value = 0.01804
#      rho 
#-0.5421053  


mean(ReadsSummary$AveDepth)
mean(ReadsSummary2$AveDepth)
