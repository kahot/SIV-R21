#Read frequency tables created in Step 1 and calculate them to mutation freq (SeqData) for each site

library(tidyverse)
source("Rscripts/baseRscript.R")


#######
SIVFiles<-list.files("Output/CSV/",pattern="csv")

#Trim CSV files 


start<-94
end<-834
no<-data.frame("ref251.pos"=c(start:end))

csns<-read.csv("Output/Consensus_merged.pos.csv", stringsAsFactors = F, row.names = 1)
csns<-csns[csns$merged.pos>=start & csns$merged.pos<=end,] #741 rows


#StockVirus consensus vs. ref251 -how many site are different?
csns<-read.csv("Consensus_merged.pos.csv", stringsAsFactors = F, row.names = 1)
csns<-csns[csns$merged.pos>=94 ,] #741 rows
csns$con_count<-ifelse(csns$ref251==csns$Run_5_01_Animal_stock_virus,0,1)

sum(csns$con_count) #10
csns$ref251.pos[csns$con_count==1]
#400 825 826 828 829 830 831 832 833 834
#400 has only 2 reads, 826 has 3 reads, rest are ?
#only 1 site differnece @825

cs<-csns[,c("ref251.pos","Run_5_01_Animal_stock_virus.pos","ref251","ref239","Run_5_01_Animal_stock_virus")]

# Replace the pos 825 in of ref251 and call it ref 
csns$ref<-csns$ref251
csns$ref[csns$merged.pos==825]<-"c"
csns$ref[csns$merged.pos==826]<-"g"

#"merged.pos" and "ref251.pos" are same
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
    cons<-csns[,c("ref251.pos","ref239.pos", paste0(id,".pos"), "ref251","ref")]
    
    SeqData<-merge(cons,SeqData,by.y="pos",by.x=paste0(id,".pos"), all.x=T)
    
    #nuceotides from transition mutations
    SeqData$transition.maj<-sapply(SeqData$MajNt, function(x) transition(x))        
    SeqData$transition.ref<-sapply(SeqData$ref, function(x) transition(x))
    
    SeqData<-SeqData[order(SeqData$ref251.pos),]
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
            WTNum <-  SeqData[k,paste0(SeqData$ref[k])]
            MutNum2<- SeqData[k,paste0(SeqData$transition.ref[k])]
            
            SeqData$freq.Ts[k]<-MutNum1/SeqData$TotalReads[k]
            SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
            
            
            #mutation frequencies of all transversion mutataions
            if (SeqData$MajNt[k]=="a"|SeqData$MajNt[k]=='g'){
                TrvMutNum<-SeqData[k,"c"]+SeqData[k,"t"]}
            if (SeqData$MajNt[k]=="c"|SeqData$MajNt[k]=="t"){
                TrvMutNum<-SeqData[k,"a"]+SeqData[k,"g"]}
            SeqData$freq.transv[k]<-TrvMutNum/SeqData$TotalReads[k]
            if (SeqData$ref[k]=="a"|SeqData$ref[k]=='g'){
                TrvMutNum2<-SeqData[k,"c"]+SeqData[k,"t"]}
            if (SeqData$ref[k]=="c"|SeqData$ref[k]=="t"){
                TrvMutNum2<-SeqData[k,"a"]+SeqData[k,"g"]}
            SeqData$freq.transv.ref[k]<-TrvMutNum2/SeqData$TotalReads[k]
            
            #Frequenceis for specific transversion mutations (1 & 2)
            Tvs1Num<-SeqData[k,paste0(transv1(SeqData$MajNt[k]))]
            Tvs2Num<-SeqData[k,paste0(transv2(SeqData$MajNt[k]))]
            SeqData$freq.transv1[k]<-Tvs1Num/SeqData$TotalReads[k]
            SeqData$freq.transv2[k]<-Tvs2Num/SeqData$TotalReads[k]
            Tvs1rNum<-SeqData[k,paste0(transv1(SeqData$ref[k]))]
            Tvs2rNum<-SeqData[k,paste0(transv2(SeqData$ref[k]))]
            SeqData$freq.transv1.ref[k]<-Tvs1rNum/SeqData$TotalReads[k]
            SeqData$freq.transv2.ref[k]<-Tvs2rNum/SeqData$TotalReads[k]
            
            
            #Frequencies of all SNPs (no indels)
            AllMutNum<-SeqData$TotalReads[k]-MajNum
            AllMutNum2<-SeqData$TotalReads[k]-WTNum
            
            SeqData$freq.mutations[k]<-AllMutNum/SeqData$TotalReads[k]
            SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
            
        }
    }
    write.csv(SeqData,paste0("Output/SeqData/SeqData_",id,".csv"))
}



