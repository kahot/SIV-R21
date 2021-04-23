#Script to analyse the frequency data and associate with features for Ita et al(2018) data
library(dplyr)
source("Rscripts/baseRscript.R")

#get the file name
ItaSeqData<-list.files("Output/SeqData/Ita/",pattern="SeqData")

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

RefDF$Type.ref<-TypeOfSite
RefDF$Type.tv1.ref<-TypeOfSite.tv1
RefDF$Type.tv2.ref<-TypeOfSite.tv2

for (k in 1:nrow(RefDF)){
    if (k%%3==1){
        RefDF$WTAA.ref[k] = seqinr::translate(RefDF$Ita.stock[c(k,k+1,k+2)])
        RefDF$MUTAA.ref[k] =seqinr::translate(c(transition(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
        RefDF$TVS1_AA.ref[k] = seqinr::translate(c(transv1(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
        RefDF$TVS2_AA.ref[k] = seqinr::translate(c(transv2(RefDF$Ita.stock[k]), RefDF$Ita.stock[c(k+1,k+2)]))
    } 
    if (k%%3==2){
        RefDF$WTAA.ref[k]   = seqinr::translate(  RefDF$Ita.stock[c(k-1,k,k+1)])
        RefDF$MUTAA.ref[k]  = seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transition(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
        RefDF$TVS1_AA.ref[k]= seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transv1(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
        RefDF$TVS2_AA.ref[k]= seqinr::translate(c(RefDF$Ita.stock[c(k-1)],transv2(RefDF$Ita.stock[k]),RefDF$Ita.stock[c(k+1)]))
    }
    if (k%%3==0){
        RefDF$WTAA.ref[k] =    seqinr::translate(  RefDF$Ita.stock[c(k-2,k-1,k)])
        RefDF$MUTAA.ref[k] =   seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transition(RefDF$Ita.stock[k])))
        RefDF$TVS1_AA.ref[k] = seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transv1(RefDF$Ita.stock[k])))
        RefDF$TVS2_AA.ref[k] = seqinr::translate(c(RefDF$Ita.stock[c(k-2,k-1)],transv2(RefDF$Ita.stock[k])))
        }
}
#Add whether AA change is drastic & makes CpG
RefDF$bigAAChange.ref<-0
RefDF$bigAAChange.tv1.ref<-0
RefDF$bigAAChange.tv2.ref<-0
RefDF$makesCpG.ref <- 0
RefDF$makesCpG.tv1.ref <- 0
RefDF$makesCpG.tv2.ref <- 0

for(j in 1:nrow(RefDF)){
    WT <- amCat(RefDF[j,'WTAA.ref'])
    MUT<- amCat(RefDF[j,'MUTAA.ref'])
    MUT1<-amCat(RefDF[j,'TVS1_AA.ref'])
    MUT2<-amCat(RefDF[j,'TVS2_AA.ref'])
    
    if (WT != MUT) RefDF$bigAAChange.ref[j] <- 1
    if (WT != MUT1) RefDF$bigAAChange.tv1.ref[j] <- 1
    if (WT != MUT2) RefDF$bigAAChange.tv2.ref[j] <- 1
    
    trip <- RefDF$Ita.stock[c(j, j+1,j+2)]
    if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) next
    else {
        if (trip[1] == "c" & trip[2] == "a" ) RefDF$makesCpG.ref[j] <- 1 
        if (trip[2] == "t" & trip[3] == "g")  RefDF$makesCpG.ref[j] <- 1
        if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) RefDF$makesCpG.tv2.ref[j] <- 1                                
        if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) RefDF$makesCpG.tv1.ref[j] <- 1
        
    }
} 

RefDF<-RefDF[-c(1,7:15,17:27)]
RefDF$AA251pos<-ceiling(RefDF$ref251.pos/3)
RefDF$AA239pos<-ceiling(RefDF$ref239.pos/3)
write.csv(RefDF, "Output/Ita.stock_overview.csv")

#Attach the information to each SeqData file
#dir.create("Output/Overview/Ita")
RefDF<-read.csv( "Output/Ita.stock_overview.csv", stringsAsFactors = F, row.names = 1)
for (i in 1:length(ItaSeqData)){   
        id<-gsub(".csv",'',paste(ItaSeqData[i]))
        id<-gsub("SeqData_",'',id)
        DF<-read.csv(paste0("Output/SeqData/Ita/",ItaSeqData[i]),stringsAsFactors=FALSE, row.names = 1)
        
        TypeOfSite<-c() 
        TypeOfSite.tv1<-c()
        TypeOfSite.tv2<-c()
        for (codon in 1:(nrow(DF)/3)) { 
            positions <- c(codon*3-2,codon*3-1, codon*3)  
            WTcodon <- DF$MajNt[positions]  
            mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])  
            mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
            mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
            if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
                WTcodon<-c('n','n','n')
                mutant1codon<-c('n','n','n')
                mutant2codon<-c('n','n','n')
                mutant3codon<-c('n','n','n')}
            else {
            #transversion mutation to 'a' or 'c'
                mutant1codon.tv1 <- c(transv1(WTcodon[1]), WTcodon[2:3]) 
                mutant2codon.tv1 <- c(WTcodon[1],transv1(WTcodon[2]), WTcodon[3])
                mutant3codon.tv1 <- c(WTcodon[1:2], transv1(WTcodon[3]))
                #transversion mutation to 'g' or 't'
                mutant1codon.tv2 <- c(transv2(WTcodon[1]), WTcodon[2:3])  
                mutant2codon.tv2 <- c(WTcodon[1],transv2(WTcodon[2]), WTcodon[3])
                mutant3codon.tv2 <- c(WTcodon[1:2], transv2(WTcodon[3]))
            }
            
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
        
        DF$Type<-TypeOfSite
        DF$Type.tv1<-TypeOfSite.tv1
        DF$Type.tv2<-TypeOfSite.tv2
        
        for (k in 1:nrow(DF)){
            if (k%%3==1){
                if (is.na(DF$MajNt[k])|is.na(DF$MajNt[k+1])|is.na(DF$MajNt[k+2])) { 
                    DF$WTAA[k]<-"NA"
                    DF$MUTAA[k]<-"NA"
                    DF$TVS1_AA[k]<-"NA"
                    DF$TVS2_AA[k]<-"NA"}
                else{
                    DF$WTAA[k] = seqinr::translate(DF$MajNt[c(k,k+1,k+2)])
                    DF$MUTAA[k] =seqinr::translate(c(transition(DF$MajNt[k]), DF$MajNt[c(k+1,k+2)]))
                    DF$TVS1_AA[k] = seqinr::translate(c(transv1(DF$MajNt[k]), DF$MajNt[c(k+1,k+2)]))
                    DF$TVS2_AA[k] = seqinr::translate(c(transv2(DF$MajNt[k]), DF$MajNt[c(k+1,k+2)]))}
            } 
            if (k%%3==2){
                if (is.na(DF$MajNt[k-1])|is.na(DF$MajNt[k])|is.na(DF$MajNt[k+1])) { 
                    DF$WTAA[k]<-"NA"
                    DF$MUTAA[k]<-"NA"
                    DF$TVS1_AA[k]<-"NA"
                    DF$TVS2_AA[k]<-"NA"}
                else{
                    DF$WTAA[k]   = seqinr::translate(  DF$MajNt[c(k-1,k,k+1)])
                    DF$MUTAA[k]  = seqinr::translate(c(DF$MajNt[c(k-1)],transition(DF$MajNt[k]),DF$MajNt[c(k+1)]))
                    DF$TVS1_AA[k]= seqinr::translate(c(DF$MajNt[c(k-1)],transv1(DF$MajNt[k]),DF$MajNt[c(k+1)]))
                    DF$TVS2_AA[k]= seqinr::translate(c(DF$MajNt[c(k-1)],transv2(DF$MajNt[k]),DF$MajNt[c(k+1)]))
                }}
            if (k%%3==0){
                if (is.na(DF$MajNt[k-2])|is.na(DF$MajNt[k-1])|is.na(DF$MajNt[k])) { 
                    DF$WTAA[k]<-"NA"
                    DF$MUTAA[k]<-"NA"
                    DF$TVS1_AA[k]<-"NA"
                    DF$TVS2_AA[k]<-"NA"}
                else{
                    DF$WTAA[k] =    seqinr::translate(  DF$MajNt[c(k-2,k-1,k)])
                    DF$MUTAA[k] =   seqinr::translate(c(DF$MajNt[c(k-2,k-1)],transition(DF$MajNt[k])))
                    DF$TVS1_AA[k] = seqinr::translate(c(DF$MajNt[c(k-2,k-1)],transv1(DF$MajNt[k])))
                    DF$TVS2_AA[k] = seqinr::translate(c(DF$MajNt[c(k-2,k-1)],transv2(DF$MajNt[k])))
                }}
        }
        #Add whether AA change is drastic & makes CpG
        DF$bigAAChange<-0
        DF$bigAAChange.tv1<-0
        DF$bigAAChange.tv2<-0
        DF$makesCpG <- 0
        DF$makesCpG.tv1 <- 0
        DF$makesCpG.tv2 <- 0
        
        for(j in 1:nrow(DF)){
            WT <- amCat(DF[j,'WTAA'])
            MUT<- amCat(DF[j,'MUTAA'])
            MUT1<-amCat(DF[j,'TVS1_AA'])
            MUT2<-amCat(DF[j,'TVS2_AA'])
            
            if (WT != MUT)  DF$bigAAChange[j] <- 1
            if (WT != MUT1) DF$bigAAChange.tv1[j] <- 1
            if (WT != MUT2) DF$bigAAChange.tv2[j] <- 1
            
            trip <- DF$MajNt[c(j, j+1,j+2)]
            if (is.na(trip[1])|is.na(trip[2])|is.na(trip[3])) next
            else {
                if (trip[1] == "c" & trip[2] == "a" ) DF$makesCpG[j] <- 1 
                if (trip[2] == "t" & trip[3] == "g")  DF$makesCpG[j] <- 1
                if (trip[1] == "c" & (trip[2]=="c"|trip[2]=='t')) DF$makesCpG.tv2[j] <- 1                                
                if (trip[3] == "g" & (trip[2]=="a"|trip[2]=="g")) DF$makesCpG.tv1[j] <- 1
                
            }
        } 
        
                
        DF<-merge(DF,RefDF[c(1,7:21)], by="ref251.pos")
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
    theme(axis.text.x=element_text(angle=90))+ylab("Average read depth")+
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
    df<-df[df$TotalReads>=100,]
    MF$SampleID[i]<-gsub("_overview.csv",'',OverviewFiles[i])
    MF$ave.mf[i]<-mean(df$freq.mutations[df$freq.mutations>0.005], na.rm=T)
}


Summary<-merge(MF, ReadsSummary[,c("SampleID","AveDepth")], by="SampleID")
write.csv(Summary,"Output/AveMF_Deapth_Ita.summary.csv")

ggplot(Summary,aes(x=ave.mf, y=AveDepth))+
    geom_point()+
    xlab("Average minor variant freq")+ylab("Average read depth")
ggsave("Output/depth.vs.mf_Ita.pdf", width = 4, height = 4)

cor.test(Summary$ave.mf, Summary$AveDepth, method="spearman")
#p-value = 0.005153
#      rho 
#-0.6245614  


mean(ReadsSummary$AveDepth)
mean(ReadsSummary2$AveDepth)
