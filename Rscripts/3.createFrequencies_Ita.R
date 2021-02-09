library(tidyverse)
source("Rscripts/baseRscript.R")


#######
SIVFiles<-list.files("Output/CSV/Ita/",pattern="csv")

#Trim CSV files 
Pos<-read.csv("Data/Ita.env.start.position.csv",stringsAsFactors = F)
Pos$ID<-gsub("-",".",Pos$ID)
Pos$ID<-paste0("mm",Pos$ID)

#start<-94
end<-834
no<-data.frame("merged.pos"=c(start:end))

csns<-read.csv("Consensus_merged.pos.ita.csv", stringsAsFactors = F, row.names = 1)
#csns<-csns[csns$merged.pos>=start & csns$merged.pos<=end,] #741 rows

for (i in 1:length(SIVFiles)){
    id<-gsub(".csv",'',paste(SIVFiles[i]))
    id<-gsub("-",".", id)
    id<-paste0("mm",id)
    cat(i," ",id)
    SeqData<-read.csv(paste0("Output/CSV/Ita/",SIVFiles[i]), row.names = 1, stringsAsFactors = F)
    SeqData<-SeqData[,-c(1,8)]
    
    colnames(SeqData)[7:8]<-c("deletion","insertion")
    colnames(SeqData)[2:5]<-c("a","c","g","t")
    
    #adjust the start to Env's start
    n<-Pos$env.start.pos[Pos$ID==id]
    SeqData<-SeqData[SeqData$start>=n,]
    
    #determine the majority nucleotide base at each site
    SeqData$MajNt<-apply(SeqData[,2:5],1,function(x) c("a","c","g","t")[which.max(x)])
    
    #check the start codon
    print(SeqData$MajNt[1:3]) #should bve atg
    
    #add the position number starting 1 and remove 'start column'
    SeqData$pos<-1:nrow(SeqData)
    SeqData<-SeqData[,-1]
    SeqData<-SeqData[1:end,]
    
    #Add the consensus ref and position info  
    #Sergio's sequences and ref251 have all same position numbers
    if (i==19){
        cons<-csns[,c("ref251.pos","ref239.pos", "Ita.stock","Ita.stock.pos", "codon","ref251")]
        SeqData<-merge(cons,SeqData,by.y="pos",by.x="Ita.stock.pos", all=T)
    } 
    else {
        cons<-csns[,c("ref251.pos","ref239.pos", paste0(id,".pos"), "Ita.stock", "codon","ref251")]
        SeqData<-merge(cons,SeqData,by.y="pos",by.x=paste0(id,".pos"), all=T)
    }
    #Sort based on merged.pos.Ita
    SeqData<-SeqData[order(SeqData$ref251.pos),]
    SeqData<-SeqData[1:end,]
    #nuceotides from transition mutations
    #SeqData$transition.maj<-sapply(SeqData$MajNt, function(x) transition(x))        
    SeqData$transition.ref<-sapply(SeqData$Ita.stock, function(x) transition(x))
    
    #determine Transition mutation freq of every site.
    
    for (k in 1:nrow(SeqData)){
        if (is.na(SeqData$MajNt[k])) {
            SeqData$freq.Ts.ref[k]<-NA
            SeqData$freq.transv.ref[k]<-NA
            SeqData$freq.transv1.ref[k]<-NA
            SeqData$freq.transv2.ref[k]<-NA
            SeqData$freq.mutations.ref[k]<-NA #all mutations
        }
        else {
            WTNum <-  SeqData[k,paste0(SeqData$Ita.stock[k])]
            MutNum2<- SeqData[k,paste0(SeqData$transition.ref[k])]
            
            SeqData$freq.Ts.ref[k]<-MutNum2/SeqData$TotalReads[k]
            
            #Frequenceis for specific transversion mutations (1 & 2)
            Tvs1Num<-SeqData[k,paste0(transv1(SeqData$Ita.stock[k]))]
            Tvs2Num<-SeqData[k,paste0(transv2(SeqData$Ita.stock[k]))]
            SeqData$freq.transv1.ref[k]<-Tvs1Num/SeqData$TotalReads[k]
            SeqData$freq.transv2.ref[k]<-Tvs2Num/SeqData$TotalReads[k]
            SeqData$freq.transv.ref[k]<-SeqData$freq.transv1.ref[k]+SeqData$freq.transv2.ref[k]
            
            #Frequencies of all SNPs (no indels)
            AllMutNum2<-SeqData$TotalReads[k]-WTNum
            SeqData$freq.mutations.ref[k]<-AllMutNum2/SeqData$TotalReads[k]
            
        }
    }
    write.csv(SeqData,paste0("Output/SeqData/Ita/SeqData_",id,".csv"))
}



#dir.create("Output/SeqData/Ita/")
