library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#dir.create("Output/CSV/")
#number of sampels to process
bamfiles<-list.files("Output/bam2/",pattern="bam$")

for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0("Output/bam2/",bam,'.bai')
        bf<-BamFile(paste0("Output/bam2/",bam), index=index)

        file.name<-paste(bam)
        file.name<-sub('_BWA_sort.bam','',file.name)
        p_param <- PileupParam(max_depth=500000,include_insertions=TRUE,include_deletions=TRUE)
        result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
        summary<-pileupFreq(result)

        summary$TotalReads<-rowSums(summary[3:6])
        maxr<-max(summary$TotalReads)
        print(file.name)
        cat("The maximum number of read depth is ", maxr)
        cat("\n")
        write.csv(summary, file=paste0("Output/CSV/",file.name,".csv",collapse=""))
        
  }



###############
###############

library(ape)
# Create a data frame with merged positions, ref (mac329) positions, and each consensus seq positions 
con<-read.dna("Data/R21_SIV251_Consensus_alignment.fasta", format = "fasta",as.character=T)
csns<-data.frame(t(con), stringsAsFactors = F)
csns$merged.pos<-1:nrow(csns)
colnames(csns)[3:73]<-gsub('_Consensus','',colnames(csns[3:73]))

#Name 251 Consensus as ref
colnames(csns)[1]<-"ref251"
colnames(csns)[2]<-"ref239"

#max position set to 834
csns<-csns[csns$merged.pos<=834,]


for (i in 1:73){
    nt<-csns[,i]
    k=1
    pos<-c()
    for (j in 1:length(nt)){
        x<-nt[j]
        if (x!="-") {
            pos<-c(pos,k)
            k=k+1}
        if (x=="-"){
            pos<-c(pos,NA)
        }
    }
    
    csns[paste0(colnames(csns)[i],".pos")]<-pos
}


#write.csv(csns,"Consensus_merged.pos.csv")
cs<-csns[,74:ncol(csns)]
write.csv(cs,"merged.pos.csv")

## codon positions
#Check insertions are with 3x nucetoides
for (i in 1:73){
    seq<-csns[,i]
    print(colnames(csns)[i])
    seq<-paste(seq,collapse = '')
    tri<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
    for (j in 1:length(tri)){
        cd<-tri[j]
        codon<-unlist(strsplit(cd,''))
        if ("-" %in% codon) {
            pos<-j*3+93
            print(pos)
            print(codon)
        }
    }
}

csns$codon<-rep(c(1,2,3),times=nrow(csns)/3)
write.csv(csns,"Consensus_merged.pos.csv")
