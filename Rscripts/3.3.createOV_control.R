#Assess the control sample
library(dplyr)
source("Rscripts/baseRscript.R")

#get the file name

#create the type of mutations info for control (6_01)
RefDF<-read.csv("Output/SeqData/SeqData_Run_6_01_Animal_stock_virus.csv",stringsAsFactors=FALSE, row.names = 1)

TypeOfSite<-c() 
TypeOfSite.tv1<-c()
TypeOfSite.tv2<-c()
for (codon in 1:(nrow(RefDF)/3)) { 
    positions <- c(codon*3-2,codon*3-1, codon*3)  
    WTcodon <- RefDF$MajNt[positions] 
    
    if (is.na(WTcodon[1])|is.na(WTcodon[2])|is.na(WTcodon[3])){ 
        WTcodon<-c('n','n','n')
        mutant1codon<-c('n','n','n')
        mutant2codon<-c('n','n','n')
        mutant3codon<-c('n','n','n')}
    else{                        
    
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

RefDF$Type<-TypeOfSite
RefDF$Type.tv1<-TypeOfSite.tv1
RefDF$Type.tv2<-TypeOfSite.tv2

for (k in 1:nrow(RefDF)){
    if (k%%3==1){
        if (is.na(RefDF$MajNt[k])|is.na(RefDF$MajNt[k+1])|is.na(RefDF$MajNt[k+2])) {
            RefDF$WTAA[k]<-"NA"
            RefDF$MUTAA[k]<-"NA"
            RefDF$TVS1_AA[k]<-"NA"
            RefDF$TVS2_AA[k]<-"NA"}
        else{
            RefDF$WTAA[k] = seqinr::translate(RefDF$MajNt[c(k,k+1,k+2)])
            RefDF$MUTAA[k] =seqinr::translate(c(transition(RefDF$MajNt[k]), RefDF$MajNt[c(k+1,k+2)]))
            RefDF$TVS1_AA[k] = seqinr::translate(c(transv1(RefDF$MajNt[k]), RefDF$MajNt[c(k+1,k+2)]))
            RefDF$TVS2_AA[k] = seqinr::translate(c(transv2(RefDF$MajNt[k]), RefDF$MajNt[c(k+1,k+2)]))
        }
    } 
    if (k%%3==2){
        if (is.na(RefDF$MajNt[k-1])|is.na(RefDF$MajNt[k])|is.na(RefDF$MajNt[k+1])) {
            RefDF$WTAA[k]<-"NA"
            RefDF$MUTAA[k]<-"NA"
            RefDF$TVS1_AA[k]<-"NA"
            RefDF$TVS2_AA[k]<-"NA"}
        else{
            RefDF$WTAA[k]   = seqinr::translate(  RefDF$MajNt[c(k-1,k,k+1)])
            RefDF$MUTAA[k]  = seqinr::translate(c(RefDF$MajNt[c(k-1)],transition(RefDF$MajNt[k]),RefDF$MajNt[c(k+1)]))
            RefDF$TVS1_AA[k]= seqinr::translate(c(RefDF$MajNt[c(k-1)],transv1(   RefDF$MajNt[k]),RefDF$MajNt[c(k+1)]))
            RefDF$TVS2_AA[k]= seqinr::translate(c(RefDF$MajNt[c(k-1)],transv2(   RefDF$MajNt[k]),RefDF$MajNt[c(k+1)]))
        }
    }
    
    if (k%%3==0){
        if (is.na(RefDF$MajNt[k-2])|is.na(RefDF$MajNt[k-1])|is.na(RefDF$MajNt[k])) {
            RefDF$WTAA[k]<-"NA"
            RefDF$MUTAA[k]<-"NA"
            RefDF$TVS1_AA[k]<-"NA"
            RefDF$TVS2_AA[k]<-"NA"}
        else{
            RefDF$WTAA[k] =    seqinr::translate(  RefDF$MajNt[c(k-2,k-1,k)])
            RefDF$MUTAA[k] =   seqinr::translate(c(RefDF$MajNt[c(k-2,k-1)],transition(RefDF$MajNt[k])))
            RefDF$TVS1_AA[k] = seqinr::translate(c(RefDF$MajNt[c(k-2,k-1)],transv1(RefDF$MajNt[k])))
            RefDF$TVS2_AA[k] = seqinr::translate(c(RefDF$MajNt[c(k-2,k-1)],transv2(RefDF$MajNt[k])))
        }
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

colnames(RefDF)[1]<-"control.pos"
write.csv(RefDF, "Output/Control_overview.csv")

