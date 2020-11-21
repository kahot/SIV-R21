library(tidyverse)
source("Rscripts/baseRscript.R")


#######
SIVFiles<-list.files("Output/CSV/",pattern="csv")

start<-94
end<-834
no<-data.frame("merged.pos"=c(start:end))

csns<-read.csv("Consensus_merged.pos.csv", stringsAsFactors = F, row.names = 1)
csns<-csns[csns$merged.pos>=start & csns$merged.pos<=end,] #741 rows


for (i in 1:length(SIVFiles)){
    print(i)
    id<-gsub(".csv",'',paste(SIVFiles[i]))
    SeqData<-read.csv(paste0("Output/CSV/",SIVFiles[i]), row.names = 1, stringsAsFactors = F)
    SeqData<-SeqData[,-c(1,8)]
    colnames(SeqData)[1]<-"pos"
    colnames(SeqData)[7:8]<-c("deletion","insertion")
    colnames(SeqData)[2:5]<-c("a","c","g","t")
    
    #determine the majority nucleotide base at each site
    SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
    
    #Add the consensus ref and position info  
    cons<-csns[,c("merged.pos", "ref251.pos","ref239.pos", paste0(id,".pos"), "ref251")]
    
    SeqData<-merge(cons,SeqData,by.y="pos",by.x=paste0(id,".pos"), all.x=T)
    
    #nuceotides from transition mutations
    SeqData$transition.maj<-sapply(SeqData$MajNt, function(x) transition(x))        
    SeqData$transition.ref<-sapply(SeqData$ref251, function(x) transition(x))
    
    SeqData<-SeqData[order(SeqData$merged.pos),]
    #determine Transition mutation freq of every site.
    
    #i==49 Run6_01 StockVirus ->only one with long indel at merged.pos 412-423
    for (k in 1:nrow(SeqData)){
        if (is.na(SeqData$MajNt[k])) {
            SeqData$freq.Ts[k]<-NA #transition mutations
            SeqData$freq.Ts.ref[k]<-NA
            SeqData$freq.transv[k]<-NA #transversion mutations
            SeqData$freq.transv.ref[k]<-NA
            SeqData$freq.transv1[k]<-NA
            SeqData$freq.transv2[k]<-NA
            SeqData$freq.transv1.ref[k]<-NA
            SeqData$freq.transv2.ref[k]<-NA
            
            SeqData$freq.mutations.ref[k]<-NA #all mutations
            SeqData$freq.mutations[k]<-NA
            
        }
        else {
            MajNum <- SeqData[k,paste0(SeqData$MajNt[k])]
            MutNum1<- SeqData[k,paste0(SeqData$transition.maj[k])]
            WTNum <-SeqData[k,paste0(SeqData$ref251[k])]
            WTNum <-  SeqData[k,paste0(SeqData$ref251[k])]
            MutNum2<- SeqData[k,paste0(SeqData$transition.ref[k])]
            
            SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
            SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
            
            
            #mutation frequencies of all transversion mutataions
            if (SeqData$MajNt[k]=="a"|SeqData$MajNt[k]=='g'){
                TrvMutNum<-SeqData[k,"c"]+SeqData[k,"t"]}
            if (SeqData$MajNt[k]=="c"|SeqData$MajNt[k]=="t"){
                TrvMutNum<-SeqData[k,"a"]+SeqData[k,"g"]}
            SeqData$freq.transv[k]<-TrvMutNum/SeqData$TotalReads[k]
            if (SeqData$ref251[k]=="a"|SeqData$ref251[k]=='g'){
                TrvMutNum2<-SeqData[k,"c"]+SeqData[k,"t"]}
            if (SeqData$ref251[k]=="c"|SeqData$ref251[k]=="t"){
                TrvMutNum2<-SeqData[k,"a"]+SeqData[k,"g"]}
            SeqData$freq.transv.ref[k]<-TrvMutNum2/SeqData$TotalReads[k]
            
            #Frequenceis for specific transversion mutations (1 & 2)
            Tvs1Num<-SeqData[k,paste0(transv1(SeqData$MajNt[k]))]
            Tvs2Num<-SeqData[k,paste0(transv2(SeqData$MajNt[k]))]
            SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
            SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
            Tvs1rNum<-SeqData[k,paste0(transv1(SeqData$ref251[k]))]
            Tvs2rNum<-SeqData[k,paste0(transv2(SeqData$ref251[k]))]
            SeqData$freq.transv1.ref[k]<-Tvs1Num/SeqData$TotalReads[k]
            SeqData$freq.transv2.ref[k]<-Tvs2Num/SeqData$TotalReads[k]
            
            
            #Frequencies of all SNPs (no indels)
            AllMutNum<-SeqData$TotalReads[k]-MajNum
            AllMutNum2<-SeqData$TotalReads[k]-WTNum
            
            SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
            SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
            
        }
    }
    write.csv(SeqData,paste0("Output/SeqData/SeqData_",id,".csv"))
}



