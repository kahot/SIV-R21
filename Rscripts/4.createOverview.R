#Script to analyse the frequency data and associate with features
library(dplyr)
source("Rscripts/baseRscript.R")

#get the file name
SIVFiles_SeqData<-list.files("Output/SeqData/",pattern="SeqData")

#create the type of mutations infor for ref251 
RefDF<-read.csv(paste0("Output/SeqData/",SIVFiles_SeqData[1]),stringsAsFactors=FALSE, row.names = 1)
TypeOfSite<-c() 
TypeOfSite.tv1<-c()
TypeOfSite.tv2<-c()

for (codon in 1:(nrow(DF)/3)) { 
    positions <- c(codon*3-2,codon*3-1, codon*3)  
    WTcodon <- DF$ref251[positions]  
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
        RefDF$WTAA[k] = seqinr::translate(RefDF$ref251[c(k,k+1,k+2)])
        RefDF$MUTAA[k] =seqinr::translate(c(transition(RefDF$ref251[k]), RefDF$ref251[c(k+1,k+2)]))
        RefDF$TVS1_AA[k] = seqinr::translate(c(transv1(RefDF$ref251[k]), RefDF$ref251[c(k+1,k+2)]))
        RefDF$TVS2_AA[k] = seqinr::translate(c(transv2(RefDF$ref251[k]), RefDF$ref251[c(k+1,k+2)]))
    } 
    if (k%%3==2){
        RefDF$WTAA[k]   = seqinr::translate(  RefDF$ref251[c(k-1,k,k+1)])
        RefDF$MUTAA[k]  = seqinr::translate(c(RefDF$ref251[c(k-1)],transition(RefDF$ref251[k]),RefDF$ref251[c(k+1)]))
        RefDF$TVS1_AA[k]= seqinr::translate(c(RefDF$ref251[c(k-1)],transv1(RefDF$ref251[k]),RefDF$ref251[c(k+1)]))
        RefDF$TVS2_AA[k]= seqinr::translate(c(RefDF$ref251[c(k-1)],transv2(RefDF$ref251[k]),RefDF$ref251[c(k+1)]))
    }
    if (k%%3==0){
        RefDF$WTAA[k] =    seqinr::translate(  RefDF$ref251[c(k-2,k-1,k)])
        RefDF$MUTAA[k] =   seqinr::translate(c(RefDF$ref251[c(k-2,k-1)],transition(RefDF$ref251[k])))
        RefDF$TVS1_AA[k] = seqinr::translate(c(RefDF$ref251[c(k-2,k-1)],transv1(RefDF$ref251[k])))
        RefDF$TVS2_AA[k] = seqinr::translate(c(RefDF$ref251[c(k-2,k-1)],transv2(RefDF$ref251[k])))
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
    
    trip <- RefDF$ref251[c(j, j+1,j+2)]
    if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) next
    else {
        if (trip[1] == "c" & trip[2] == "a" ) RefDF$makesCpG[j] <- 1 
        if (trip[2] == "t" & trip[3] == "g")  RefDF$makesCpG[j] <- 1
        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) RefDF$makesCpG.tv2[j] <- 1                                
        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) RefDF$makesCpG.tv1[j] <- 1
        
    }
} 

RefDF<-RefDF[-c(1, 6:14,17:26)]
write.csv(RefDF, "Output/Ref251_overview.csv")


for (i in 1:length(SIVFiles_SeqData)){   
        id<-gsub(".csv",'',paste(SIVFiles_SeqData[i]))
        id<-gsub("SeqData_",'',id)
        DF<-read.csv(paste0("Output/SeqData/",SIVFiles_SeqData[i]),stringsAsFactors=FALSE, row.names = 1)
                
        DF<-merge(DF,RefDF[c(1,7:19)], by="merged.pos")
        write.csv(DF,paste0("Output/Overview/",id,"_overview.csv"))
        print(id)
}        

#####################################################

### Read depths for all files ###
SIVFiles_SeqData<-list.files("Output/SeqData/",pattern="SeqData")

ReadsSummary<-data.frame(SampleID=matrix(nrow=length(SIVFiles_SeqData)))
ReadsSummary$MaxDepth<-""
ReadsSummary$AveDepth<-""

for (i in 1:length(SIVFiles_SeqData)){
        id<-gsub(".csv",'',paste(SIVFiles_SeqData[i]))
        id<-gsub("SeqData_",'',id)
        ReadsSummary$SampleID[i]<-id
        print(id)
        SeqData<-read.csv(paste("Output/SeqData/",SIVFiles_SeqData[i],sep=""))
        ReadsSummary$MaxDepth[i]<-max(SeqData$TotalReads,na.rm=T)
        ReadsSummary$AveDepth[i]<-mean(SeqData$TotalReads,na.rm=T)
        ReadsSummary$No.ofSites[i]<-nrow(SeqData[!is.na(SeqData$ref251),])
        ReadsSummary$SE[i]<-std.error(SeqData$TotalReads, na.rm=T)
}

write.csv(ReadsSummary,"Output/ReadsSummary.csv")      


###
ReadsSummary<-read.csv("Output/ReadsSummary.csv", row.names = 1)
ggplot(ReadsSummary, aes(x=SampleID, y=AveDepth))+
    geom_point()+
    theme(axis.text.x=element_text(angle=90))



