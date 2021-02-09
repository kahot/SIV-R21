library(ggplot2)
library(reshape)
library(gridExtra)

# read the files saved in Overview_output:
OverviewFiles<-list.files("Output/Overview/",pattern="overview.csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}

#sample info
SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
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

#amount of variation in control
control<-OvDF[["Run_6_01_Animal_stock_virus"]]
controlVar<-control[control$freq.mutations>0.005&control$TotalReads>100,] #236 positions
controlVar<-controlVar[!is.na(controlVar$ref251),] #5 positions
#pos 394 has high mf -> remove and calculate the background noise

control2<-control[control$freq.mutations<0.01&control$TotalReads>5000,] 
control2<-control2[!is.na(control2$freq.mutations),] #568 sites
#average diversity of control (% diversity/mf from the consensus)
mean(control2$freq.mutations, na.rm = T)
#0.002003343 (>5000 read depth)
#0.002015668 (>100 read depth)
mean(control2$TotalReads,na.rm = T)
#93884.1
control2<-control[control$freq.mutations<0.01&control$TotalReads>10000,] 
control2<-control2[!is.na(control2$freq.mutations),] #554 positions
mean(control2$freq.mutations, na.rm = T)
#0.001996905 (>10000 read depth)
mean(control2$TotalReads,na.rm = T)
#96077.11

#how many sites were above mf=0.05?
control3<-control2[control2$freq.mutations>0.005,] #382


# Look at the amount of variants in the stock
stockMf<-OvDF[[stock$filename]]
stockVar<-stockMf[stockMf$freq.mutations.ref>0.005&stockMf$TotalReads>1000,]
stockVar<-stockVar[!is.na(stockVar$freq.mutations.ref),] #31 positions


#% diversity of Stock
#replace the mut frea =NA for positios with total_reads<100
stockMf[which(stockMf$TotalReads<100),17:26]<-NA
ggplot(stockMf, aes(x=AA251pos, y=freq.mutations))+
    geom_point(position=position_dodge(width=0.8), size=0.5)+
    theme_bw()+ylim(0,0.5)+
    ylab("% Diversity from consensus")+
    xlab("Env codon position")+
    ggtitle(paste0("SIV251 stock diversity "))+
    geom_hline(yintercept = 0.005, color="lightsteelblue3")
#ggsave("Output/Diversity/Stock_diversity_R21.pdf", width = 5, height = 4)


## Sergio's stock virus:
stockMfS<-read.csv("Output/Overview/Ita/mmStock_overview.csv", stringsAsFactors = F, row.names = 1)
stockVarS<-stockMfS[stockMfS$freq.mutations.ref>0.005&stockMfS$TotalReads>1000,]
stockVarS<-stockVarS[!is.na(stockVarS$freq.mutations.ref),] #36 positions

stockMfS[which(stockMfS$TotalReads<100),17:26]<-NA


#separate syn vs. nonsyn and Plot the results & create a summary table
summary<-data.frame(study=c("R21","Ita"))

for (i in 1:2){
    if (i==1) df<-stockMf; title="R21"
    if (i==2) df<-stockMfS; title="Ita"
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
    
    df2<-df[,c("AA251pos","ref251.pos","ns")]
    colnames(df2)[3]<-"freq"
    df2$Type<-"nonsyn"
    df2<-df2[!is.na(df2$freq),]
    
    df3<-df[,c("AA251pos","ref251.pos","syn")]
    colnames(df3)[3]<-"freq"
    df3$Type<-"syn"
    df3<-df3[!is.na(df3$freq),]
    
    missing<-df[df$TotalReads<100|is.na(df$TotalReads),]
    missing<-missing[missing$ref251.pos>=95,]
    
    mut<-rbind(df2,df3)
    mut2<-melt(mut, id.vars=c("AA251pos","ref251.pos","Type"))
    mut2<-mut2[,-which(colnames(mut2)=="variable")]
    mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
    
    #freq<0.005=NA
    mut2$value[mut2$value<0.005]<-NA
    
    #Create a summary table
    #1 MF over 0.005    
    summary$no.sites.over.0.005[i]<-length(mut$freq[mut$freq>0.005])
    summary$mean.freq[i]<-mean(mut$freq,na.rm=T)
    summary$no.sites.syn[i]<-length(mut$freq[mut$Type=="syn" & mut$freq>0.005])
    summary$mean.syn.freq[i]<-mean(mut$freq[mut$Type=="syn"],na.rm=T)
    summary$no.sites.ns[i]<-length(mut$freq[mut$Type=="nonsyn" & mut$freq>0.005])
    summary$mean.ns.freq[i]<-mean(mut$freq[mut$Type=="nonsyn"],na.rm=T)
    
    if (i==1) {
        ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                      color="gray80",size=1,alpha=0.2)+
            geom_point(position=position_dodge(width=0.8), size=0.5)+
            scale_color_manual(values=c("red","blue"))+
            scale_fill_manual(values=c("red","blue"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(title,"SIV251 stock diversity"))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Diversity from consensus")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
    }
    if (i==2){
        ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            geom_point(position=position_dodge(width=0.8), size=0.5)+
            scale_color_manual(values=c("red","blue"))+
            scale_fill_manual(values=c("red","blue"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(title,"SIV251 stock diversity"))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Diversity from consensus")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
    }
    
    ggsave(paste0("Output/Diversity/",title,"_stock_diversity_syn.nonsyn.pdf"),width = 5.5, height = 4)
}

write.csv(summary,"Output/Diversity/Stock_virus_diversity_summary_R21_Ita.csv")




# Look at the amount of variants in the 1st sample of each animal (plasma)
#R21
Plots<-list()
for (i in 1:length(monkeyList)){
    sample<-monkeyList[[i]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[i]
    
    #select the first sample after infection
    df<-OvDF[[sample$filename[1]]]
    df[which(df$TotalReads<100),17:26]<-NA
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
    
    df2<-df[,c("AA251pos","ref251.pos","ns")]
    colnames(df2)[3]<-"freq"
    df2$Type<-"nonsyn"
    df2<-df2[!is.na(df2$freq),]
    
    df3<-df[,c("AA251pos","ref251.pos","syn")]
    colnames(df3)[3]<-"freq"
    df3$Type<-"syn"
    df3<-df3[!is.na(df3$freq),]
    
    missing<-df[df$TotalReads<100|is.na(df$TotalReads),]
    missing<-missing[missing$ref251.pos>=95,]
    
    mut<-rbind(df2,df3)
    mut2<-melt(mut, id.vars=c("AA251pos","ref251.pos","Type"))
    mut2<-mut2[,-which(colnames(mut2)=="variable")]
    mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
    
    #freq<0.005=NA
    mut2$value[mut2$value<0.005]<-NA
    
    p<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                      color="gray80",size=1,alpha=0.2)+
            geom_point(position=position_dodge(width=0.8), size=0.5)+
            scale_color_manual(values=c("red","blue"))+
            scale_fill_manual(values=c("red","blue"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
            theme_bw()+
            ggtitle(paste0(monkey," Week ",sample$Week[1]))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Divergence")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
    
    Plots[[i]]<-p
}

pdf("Output/Diversity/FounderDiversity.pdf", width = 11, height=20)
do.call(grid.arrange, c(Plots, ncol=2))
dev.off()

#sergio's data
ItaFiles<-list.files("Output/Overview/Ita/",pattern="overview.csv")

ItaOv<-list()
for (i in 1:length(ItaFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/Ita/",ItaFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    ItaOv[[i]]<-overviews
    names(ItaOv)[i]<-gsub("_overview.csv",'',ItaFiles[i])
}
#first sampling
first<-c("mm10.2","mm156.3","mm174.2","mm198.2" )

Plots2<-list()
for (i in 1:length(first)){
    df<-ItaOv[[first[i]]]
    if (i==1) monkey<-"TP1"; week=2
    if (i==2) monkey<-"RP1"; week=3
    if (i==3) monkey<-"SP1"; week=2
    if (i==4) monkey<-"TP2"; week=2
    
    
    df[which(df$TotalReads<100),17:26]<-NA
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
    
    df2<-df[,c("AA251pos","ref251.pos","ns")]
    colnames(df2)[3]<-"freq"
    df2$Type<-"nonsyn"
    df2<-df2[!is.na(df2$freq),]
    
    df3<-df[,c("AA251pos","ref251.pos","syn")]
    colnames(df3)[3]<-"freq"
    df3$Type<-"syn"
    df3<-df3[!is.na(df3$freq),]
    
    mut<-rbind(df2,df3)
    mut2<-melt(mut, id.vars=c("AA251pos","ref251.pos","Type"))
    mut2<-mut2[,-which(colnames(mut2)=="variable")]
    mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
    
    #freq<0.005=NA
    mut2$value[mut2$value<0.005]<-NA
    
    p<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
        geom_point(position=position_dodge(width=0.8), size=0.5)+
        scale_color_manual(values=c("red","blue"))+
        scale_fill_manual(values=c("red","blue"),guide = 'none')+
        scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
        scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
        theme_bw()+
        ggtitle(paste0(monkey," Week ",week))+
        theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
        ylab("% Divergence")+
        xlab("Env codon position")+
        theme(legend.title = element_blank())
    
    Plots2[[i]]<-p
}

pdf("Output/Diversity/FounderDiversity_Ita.pdf", width = 11, height=8)
do.call(grid.arrange, c(Plots2, ncol=2))
dev.off()

