library(ggplot2)
library(reshape)
library(gridExtra)

library(plotly)
library(tidyverse)
library(htmlwidgets)

#source("Rscripts/baseRscript.R")

# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}


SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stocks<-SampleSheet[SampleSheet$Monkey=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

list.animal<-split(samples, samples$Monkey)
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
    if (nrow(list.animal[[i]])>1){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
    }
}
monkeys<-names(monkeyList)





#calculate the % diversity from consensus per codon. Average over the 3 bases?
for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF[sample$filename]    
    
    Plots<-list()
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #df.org<-df
        #calculate total ns freq at each position
        df$ns1<-as.numeric(apply(df[,c("Type","freq.Ts.ref")],1, function(x) if (x["Type"]=="nonsyn") x["freq.Ts.ref"] else 0))
        df$ns2<-as.numeric(apply(df[,c("Type.tv1","freq.transv1.ref")],1, function(x) if (x["Type.tv1"]=="nonsyn") x=x["freq.transv1.ref"] else 0))
        df$ns3<-as.numeric(apply(df[,c("Type.tv2","freq.transv2.ref")],1, function(x) if (x["Type.tv2"]=="nonsyn") x=x["freq.transv2.ref"] else 0))
        df$ns<-df$ns1+df$ns2+df$ns3
        #calculate total syn freq at each position
        df$syn1<-as.numeric(apply(df[,c("Type","freq.Ts.ref")],1, function(x) if (x["Type"]=="syn") x["freq.Ts.ref"] else 0))
        df$syn2<-as.numeric(apply(df[,c("Type.tv1","freq.transv1.ref")],1, function(x) if (x["Type.tv1"]=="syn") x=x["freq.transv1.ref"] else 0))
        df$syn3<-as.numeric(apply(df[,c("Type.tv2","freq.transv2.ref")],1, function(x) if (x["Type.tv2"]=="syn") x=x["freq.transv2.ref"] else 0))
        df$syn<-df$syn1+df$syn2+df$syn3
    
        #df.org$ns<-df$ns
        #df.org$sy<-df$syn
        
        #add codon position based on ref239
        k=96/3
        df$codon.pos<-NA
        for (j in 1:nrow(df)){
            if (is.na(df$ref239.pos[j])) df$codon.pos[j]<-NA: next
            if (df$ref239.pos[j]%%3==1|df$ref239.pos[j]%%3==2) df$codon.pos[j]<-k
            if (df$ref239.pos[j]%%3==0) {
                df$codon.pos[j]<-k
                k=k+1
            }
        }
        
        #eliminate noise/seq errors
        df2<-df[df$ns>=0.05&df$TotalReads>=100,c("codon.pos","ref251.pos","ns")]
        colnames(df2)[3]<-"freq"
        df2$Type<-"nonsyn"
        df2<-df2[!is.na(df2$freq),]
        #if(nrow(df2)==0) df2[1,1:3]<-c(1,1,0): df2[1,4]<-"nonsyn"
        df3<-df[df$syn>=0.05,c("codon.pos","ref251.pos","syn")]
        colnames(df3)[3]<-"freq"
        df3$Type<-"syn"
        df3<-df3[!is.na(df3$freq),]
        #if(nrow(df3)==0) df3[1,1:3]<-c(1,1,0); df3[1,4]<-"syn"

        mut<-rbind(df2,df3)
        mut2<-melt(mut, id.vars=c("codon.pos","ref251.pos","Type"))
        mut2<-mut2[,-which(colnames(mut2)=="variable")]
        mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
        
        #nonsyn only : make all blue
        if (length(unique(mut2$Type))==1){
            if (unique(mut2$Type)=="nonsyn") col="blue"
            if (unique(mut2$Type)=="syn") col='red'
            if (i==1){
                p<-ggplot(data=mut2, aes(x=codon.pos, y=value))+
                    geom_bar(stat = "identity", width=0.1, color=col)+
                    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
                    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1.1))+
                    theme_bw()+
                    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
                    ylab("Diversity")+
                    xlab("Codon position")+
                    theme(legend.title = element_blank())+
                    geom_text(aes(label=codon.pos), hjust=0.3, vjust=-.8, size=2.5)+
                    geom_text(aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray")+
                    ggtitle(paste0("Animal:", monkey," Week",sample$Week[i] ))+
                    annotate(geom="text", x=36, y=1, hjust=0,label="nonsyn",color ='blue', size=2.5)+
                    annotate(geom="text", x=36,  y=.9, hjust=0, label="syn",color ='red', size=2.5)+
                    annotate("segment", x = 30, xend = 34, y = 1, yend = 1, colour = "blue") +
                    annotate("segment", x = 30, xend = 34, y = .9, yend = .9, colour = "red") 
                }
            else {
                p<-ggplot(data=mut2, aes(x=codon.pos, y=value))+
                    geom_bar(stat = "identity", width=0.1, color=col)+
                    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
                    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1.1))+
                    theme_bw()+
                    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
                    ylab("Diversity")+
                    xlab("Codon position")+
                    theme(legend.title = element_blank())+
                    geom_text(aes(label=codon.pos), hjust=0.3, vjust=-.8, size=2.5)+
                    geom_text(aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray")+
                    ggtitle(paste0("Animal:", monkey," Week",sample$Week[i] ))
            }
        }
        else{
            if (i==1){
                p<-ggplot(data=mut2, aes(x=codon.pos, y=value, color=Type, fill=Type))+
                    geom_bar(stat = "identity", width=0.1)+
                    scale_color_manual(values=c("red","blue"),guide = 'none')+
                    scale_fill_manual(values=c("red","blue"),guide = 'none')+
                    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
                    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1.1))+
                    theme_bw()+
                    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
                    ylab("Diversity")+
                    xlab("Codon position")+
                    theme(legend.title = element_blank())+
                    geom_text(aes(label=codon.pos), hjust=0.3, vjust=-.8, size=2.5)+
                    geom_text(aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray")+
                    ggtitle(paste0("Animal:", monkey," Week",sample$Week[i] ))+
                    annotate(geom="text", x=36, y=1, hjust=0,label="nonsyn",color ='blue', size=2.5)+
                    annotate(geom="text", x=36,  y=.9, hjust=0, label="syn",color ='red', size=2.5)+
                    annotate("segment", x = 30, xend = 34, y = 1, yend = 1, colour = "blue") +
                    annotate("segment", x = 30, xend = 34, y = .9, yend = .9, colour = "red") 
                }
            else {
                p<-ggplot(data=mut2, aes(x=codon.pos, y=value, color=Type, fill=Type))+
                    geom_bar(stat = "identity", width=0.1)+
                    scale_color_manual(values=c("red","blue"),guide = 'none')+
                    scale_fill_manual(values=c("red","blue"),guide = 'none')+
                    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
                    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1.1))+
                    theme_bw()+
                    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
                    ylab("Diversity")+
                    xlab("Codon position")+
                    theme(legend.title = element_blank())+
                    geom_text(aes(label=codon.pos), hjust=0.3, vjust=-.8, size=2.5)+
                    geom_text(aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray")+
                    ggtitle(paste0("Animal:", monkey," Week",sample$Week[i] ))
            }
        }
        #ggplotly(p)
        #htmlwidgets::saveWidget(as_widget(p), "OccupationWages.html")
        Plots[[i]]<-p
    }

    pdf(paste0("Output/MF/Animal",monkey,".pdf"), width = 4, height=nrow(sample)*2)
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()

}


#####################


#Variable positions in V1 loop from Table 1
varAA<-c(120,132,135,136,138,139,198,201,202)


#nt positions
varNT<-c()
for (i in 1: length(varAA)){
    pos<-c(varAA[i]*3-2,varAA[i]*3-1,varAA[i]*3)
    varNT<-c(varNT,pos)
}

#List all available mutations 
df<-OvDF[[1]]
df<-df[df$ref239.pos %in% varNT,]
df<-df[!(df$Type=="syn"&df$Type.tv1=="syn"&df$Type.tv2=="syn"),]
df$AA.pos<-ceiling(df$ref239.pos/3)

mutations<-c()
for (j in 1:nrow(df)){
    if (df$Type[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$MUTAA[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv1[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$TVS1_AA[j])
        mutations<-c(mutations,mut)
    }
    if (df$Type.tv2[j]!="syn"){
        mut<-paste0(df$WTAA[j],df$AA.pos[j],df$TVS2_AA[j])
        mutations<-c(mutations,mut)
    }        
}


#Table1<-list() 

muts<-data.frame(matrix(ncol=length(mutations)+1, nrow=length(OvDF)), stringsAsFactors = F) 
colnames(muts)<-c("Sample", mutations)
muts$Sample<-names(OvDF)

for (i in 1:length(OvDF)){
    df<-OvDF[[i]]
    df<-df[df$ref239.pos %in% varNT,]
    df<-df[!(df$Type=="syn"&df$Type.tv1=="syn"&df$Type.tv2=="syn"),]
    df$AA.pos<-ceiling(df$ref239.pos/3)
    
    
    mut<-c()
    for (j in 1:nrow(df)){
        if (df$Type[j]!="syn"){
            mut<-c(mut,df$freq.Ts.ref[j])
        }
        if (df$Type.tv1[j]!="syn"){
            mut<-c(mut, df$freq.transv1.ref[j])
        }
        if (df$Type.tv2[j]!="syn"){
            mut<-c(mut, df$freq.transv2.ref[j])
        }
    }
    muts[i,2:ncol(muts)]<-mut
}


SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stocks<-SampleSheet[SampleSheet$Monkey=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

list.animal<-split(samples, samples$Monkey)
monkeyList<-list()
k=1
for (i in 1:length(list.animal)){
    if (nrow(list.animal[[i]])>1){
        monkeyList[[k]]<-list.animal[[i]]
        names(monkeyList)[k]<-names(list.animal)[i]
        k=k+1
    }
}
monkeys<-names(monkeyList)


Plots<-list()    
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[i]
    
    freq<-muts[muts$Sample %in% sample$filename,]
    colnames(freq)[1]<-"ID"
    freq$Sample<-paste0("Animal_",monkey,"_wk",sample$Week,"_",sample$Sample)
    freq<-freq[,-1]
    freqM<-melt(freq, id.vars="Sample")
    colnames(freqM)[2:3]<-c("Mutation","Freq")
    
    highFs<-freqM$Mutation[freqM$Freq>=0.005]
    highF<-unique(as.character(highFs[!is.na(highFs)]))
    
    mutID<-levels(freqM$Mutation)
    colors<-sapply(mutID, function(x) if (x%in%highF) x=2 else 1)
    Plots[[i]]<-ggplot(freqM, aes(x=Mutation, y=Freq, color=Sample))+
        geom_point(position=position_dodge(width=0.8))+
        theme_bw()+
        theme(axis.text.x=element_text(angle=90))+
        theme(axis.text.x=element_text(color=colors))+
        theme(
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8)
        )+
        ggtitle(paste0("Animal ",monkey))+
        geom_hline(yintercept = 0.005, color="lightsteelblue3")
    #ggsave(paste0("Output/Table1/",monkey,".pdf"), width = 12, height = 3.5)
}    


pdf("Output/Table1/All_animals.pdf", width = 12, height=35)
do.call(grid.arrange, c(Plots, ncol=1))
dev.off()


