library(Rsamtools)
library(stringr)

source("Rscripts/pileupFreq.R")

#dir.create("Output/CSV/")
#number of sampels to process
bamfiles<-list.files("Output/bam2/Ita/",pattern="bam$")

for (i in 1:length(bamfiles)){
        bam<-bamfiles[i]
        index<-paste0("Output/bam2/Ita/",bam,'.bai')
        bf<-BamFile(paste0("Output/bam2/Ita/",bam), index=index)

        file.name<-paste(bam)
        file.name<-sub('_sort.bam','',file.name)
        p_param <- PileupParam(max_depth=500000,include_insertions=TRUE,include_deletions=TRUE)
        result<-pileup(bf, pileupParam = p_param, distinguish_strands=FALSE,ignore_query_Ns=FALSE)
        summary<-pileupFreq(result)

        summary$TotalReads<-rowSums(summary[3:6])
        maxr<-max(summary$TotalReads)
        print(file.name)
        cat("The maximum number of read depth is ", maxr)
        cat("\n")
        write.csv(summary, file=paste0("Output/CSV/Ita/",file.name,".csv",collapse=""))
        
  }


# Add majNT info to obrain consensus seq
csvs<-list.files("Output/CSV/Ita/")
for (i in 1:length(csvs)){
    df<-read.csv(paste0("Output/CSV/Ita/",csvs[i]), stringsAsFactors = F, row.names = 1)
    colnames(df)[3:6]<-c("a","c","g","t")
    df$MajNT<-apply(df[,3:6],1,function(x) c("a","c","g","t")[which.max(x)])
    consensus<-df$MajNT
    fname<-paste0(gsub(".csv",'',csvs[i]),"_consensus")
    write.fasta(consensus,fname, "Ita.consensus.fasta",open = "a",nbchar = 5000,as.string = FALSE)
}






###############
###############

library(ape)
# Create a data frame with merged positions, ref (mac329) positions, and each consensus seq positions 
con<-read.dna("Data/Sergio_consensus_aligned.fasta", format = "fasta", as.character=TRUE, , as.matrix =T)
csns<-data.frame(t(con), stringsAsFactors = F)
csns$merged.pos.Ita<-1:nrow(csns)
colnames(csns)<-gsub('_consensus','',colnames(csns))

#Name 251 Consensus as ref
colnames(csns)[3]<-"ref251"
colnames(csns)[4]<-"ref239"
colnames(csns)[2]<-"Ita.stock"
#max position set to 834
csns<-csns[csns$merged.pos<=834,]


for (i in 1:22){
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



## codon positions
#Check insertions are with 3x nucetoides
for (i in 1:22){
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

#cs<-csns[,23:ncol(csns)]
#rename the column names
#colnames(cs)<-gsub("X",'mm',colnames(cs))
#write.csv(cs,"merged.pos.ita.csv")

colnames(csns)<-gsub("X",'mm',colnames(csns))

write.csv(csns,"Consensus_merged.pos.ita.csv")
