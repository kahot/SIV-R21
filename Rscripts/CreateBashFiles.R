#create bash files to run bbmap and bwa
library(stringr)

#choose the fastq files to be prrocessed
SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stockrows<-which(SampleSheet$Monkey=="stock_virus")



#for (i in 1:nrow(SampleSheet)){
for (i in 1:nrow(SampleSheet)){

    AnimalNumber= SampleSheet$MiseqSample[i]
    if (AnimalNumber<10) AnimalNumber = paste(c("0",AnimalNumber),collapse="")
    SampleName = paste(c("Run",SampleSheet$SampleSheet[i], AnimalNumber, "Animal", as.character(SampleSheet$Monkey[i])), collapse = "_")
    SampleName<-gsub(pattern=" ", replace="",x=SampleName)
    SampleSheet$filename<-SampleName
    print(SampleName)
    #t1=Sys.time()
    
    FilePath = paste(c(SampleSheet$FilePathPart1[i], SampleSheet$FilePathPart2[i]),collapse="/")
    filelist<-list.files(FilePath, full.names = TRUE)
    FileIn1=filelist[1]
    FileIn2=filelist[2]
    
    BashLines<-readLines("Data/BashSampleSheetMac251.sh")
    BashLines = gsub(pattern="FASTQFile1", replace=FileIn1,x=BashLines)
    BashLines = gsub(pattern="FASTQFile2", replace=FileIn2,x=BashLines)
    BashLines = gsub(pattern="SAMPLE", replace=SampleName,x=BashLines)
    bashscriptName<-paste0("Data/BashScripts1/", SampleName,".sh")
    writeLines(BashLines,bashscriptName )
}

#write.csv(SampleSheet, "Data/SampleSheetMac251All.csv", row.names = F)

## Create second bash script files to map to sample's own consensus

SIVFiles<-list.files("Output/Consensus/",pattern=".fasta")

for (i in 1:length(SIVFiles)){
    BashLines<-readLines("Data/Bash_template_mac251_2.sh")
    SampleName<-gsub("_Consensus.fasta",'',SIVFiles[i])
    BashLines = gsub(pattern="SAMPLE", replace=SampleName,x=BashLines)
    bashscriptName<-paste0("Data/BashScripts2/", SampleName,".sh")
    writeLines(BashLines,bashscriptName )
    
}


### Sergio's bam files

Files<-list.files("Output/bam2/Ita_duprm_trim_bamFiles/",pattern=".bam")

for (i in 1:length(Files)){
    BashLines<-readLines("Data/Process_ita.sh")
    SampleName<-gsub(".duprm.trim.bam",'',Files[i])
    BashLines = gsub(pattern="SAMPLE.bam", replace=Files[i],x=BashLines)
    BashLines = gsub(pattern="SAMPLE1", replace=SampleName,x=BashLines)
    bashscriptName<-paste0("Data/Bash_ita/", SampleName,".sh")
    writeLines(BashLines,bashscriptName )
    
}
