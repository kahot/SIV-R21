#Diversity comparison

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

#sergio's data
ItaFiles<-list.files("Output/Overview/Ita/",pattern="overview.csv")

ItaOv<-list()
for (i in 1:length(ItaFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/Ita/",ItaFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    ItaOv[[i]]<-overviews
    names(ItaOv)[i]<-gsub("_overview.csv",'',ItaFiles[i])
}


#calculate average mutation frequency (divergence from stock consensus) for all files
Summary<-list()
for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF[sample$filename]    

    Plots<-list()
    Plots2<-list()
    summary<-sample[,c("Monkey","Tissue","Week","filename")]
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #remove the sites with row read depth
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
        
        summary$no.sites[i]<-nrow(mut[mut$freq>0.005,])
        summary$no.syn.sites[i]<-nrow(mut[mut$freq>0.005&mut$Type=="syn",])
        summary$no.ns.sites[i]<-nrow(mut[mut$freq>0.005&mut$Type=="nonsyn",])
        summary$mean[i]<-mean(mut$freq[mut$freq>0.005], na.rm=T)
        summary$mean.syn[i]<-mean(mut$freq[mut$freq>0.005&mut$Type=="syn"], na.rm=T)
        summary$mean.ns[i]<-mean(mut$freq[mut$freq>0.005&mut$Type=="nonsyn"], na.rm=T)
        
        p<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                      color="gray80",size=1,alpha=0.2)+
            geom_point(position=position_dodge(width=0.8), size=0.5)+
            scale_color_manual(values=c("red","blue"))+
            #scale_fill_manual(values=c("red","blue"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(monkey," Week ",sample$Week[i]," ",sample$Sample[i]))+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Divergence")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
        
        Plots[[i]]<-p
        
        #barplot
        p2<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                      color="gray80",size=1,alpha=0.2)+
            geom_bar(stat = "identity", width=0.1)+
            scale_color_manual(values=c("red","blue"),guide = 'none')+
            scale_fill_manual(values=c("red","blue"))+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(monkey," Week ",sample$Week[i]," ",sample$Sample[i]))+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Divergence")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
        Plots2[[i]]<-p2
        }    
    pdf(paste0("Output/Diversity/Divergence_",monkey,".pdf"), width = 4, height=nrow(sample)*2)
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()
    
    pdf(paste0("Output/Diversity/Divergence_barplot_",monkey,".pdf"), width = 4, height=nrow(sample)*2)
    do.call(grid.arrange, c(Plots2, ncol=1))
    dev.off()
    
    
    Summary[[j]]<-summary
    
}

Sum.divergence<-do.call(rbind,Summary)    
write.csv(Sum.divergence,"Output/Diversity/Divergence_summary_R21.csv")

#Plasma vs. lung vs. lymphnodes
Sum.divergence$Tissue[Sum.divergence$Tissue=="LML"]<-"Lung"
by.tissues<-aggregate(Sum.divergence,by=list(Sum.divergence$Tissue),mean )
by.tissues<-by.tissues[,c(1,6:11)]
#  Group.1 no.sites no.syn.sites no.ns.sites       mean   mean.syn    mean.ns
#1     HLN 21.00000     9.785714   11.214286 0.08233142 0.05290551 0.11221227
#2      LN 20.94118     9.823529   11.117647 0.06048309 0.01894168 0.09470716
#3    Lung 17.89474     8.000000    9.894737 0.06651460 0.02702838 0.11360183
#4  plasma 20.00000     7.578947   12.421053 0.10053080 0.07851937 0.14646764

#plasma's diversity is higher


#add Sergio's samples
# Too many sites between 0.005~0.01. Use the cut-off at mf>0.01 for Segio's file.
Ita.samples<-read.csv("Data/Ita.SampleInfo.csv", stringsAsFactors = F)
m<-unique(Ita.samples$Monkey)
m<-m[1:4]

Sum<-list()
for (j in 1:4){
    smp<-Ita.samples[Ita.samples$Monkey==m[j],]
    ov<-ItaOv[smp$SampleID]
    
    Plots<-list()
    Plots2<-list()
    summary<-smp[,c("Monkey","Week")]
    for (i in 1:length(ov)){
        df<-ov[[i]]
    
        #remove the sites with row read depth
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
        mut2$value[mut2$value<0.01]<-NA
        
        summary$no.sites[i]<-nrow(mut[mut$freq>0.005,])
        summary$no.syn.sites[i]<-nrow(mut[mut$freq>0.005&mut$Type=="syn",])
        summary$no.ns.sites[i]<-nrow(mut[mut$freq>0.005&mut$Type=="nonsyn",])
        summary$mean[i]<-mean(mut$freq[mut$freq>0.005], na.rm=T)
        summary$mean.syn[i]<-mean(mut$freq[mut$freq>0.005&mut$Type=="syn"], na.rm=T)
        summary$mean.ns[i]<-mean(mut$freq[mut$freq>0.005&mut$Type=="nonsyn"], na.rm=T)
        
        p<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            #geom_rect(data=missing, inherit.aes=FALSE,
            #          aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
            #          color="gray80",size=1,alpha=0.2)+
            geom_point(position=position_dodge(width=0.8), size=0.5)+
            scale_color_manual(values=c("red","blue"))+
            #scale_fill_manual(values=c("red","blue"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(m[j]," Week ",smp$Week[i]))+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Divergence")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
        if (nrow(missing)>0){
            p<-p+geom_rect(data=missing, inherit.aes=FALSE,
                aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                color="gray80",size=1,alpha=0.1)
        }
        
        Plots[[i]]<-p
        
        #barplot
        p2<-ggplot(data=mut2, aes(x=AA251pos, y=value, color=Type, fill=Type))+
            
            geom_bar(stat = "identity", width=0.1)+
            scale_color_manual(values=c("red","blue"),guide = 'none')+
            scale_fill_manual(values=c("red","blue"))+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(m[j]," Week ",smp$Week[i]))+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Divergence")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
        
        if (nrow(missing)>0){
            p2<-p2+geom_rect(data=missing, inherit.aes=FALSE,
                           aes(xmin=AA251pos,xmax=AA251pos,ymin=-Inf,ymax=Inf),
                           color="gray80",size=1,alpha=0.1)
        }
    
        Plots2[[i]]<-p2
    }    
    pdf(paste0("Output/Diversity/Divergence_Ita_",m[j],".pdf"), width = 4, height=nrow(smp)*2)
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()
    
    pdf(paste0("Output/Diversity/Divergence_Ita_barplot_",m[j],".pdf"), width = 4, height=nrow(smp)*2)
    do.call(grid.arrange, c(Plots2, ncol=1))
    dev.off()
    
    Sum[[j]]<-summary
}

sum.ita<-do.call(rbind,Sum)    
#add the mean to by.divergence summary
by.tissues[5,]<-c("Ita(plasma)",colMeans(sum.ita[,3:8]))

write.csv(Sum.divergence,"Output/Diversity/Divergence_summary_R21.csv")



colMeans(sum.ita[,3:8])


#Run statisitcal test by aggregating the numbers





#select the first sample after infection
df<-OvDF[[sample$filename[1]]]
df[which(df$TotalReads<100),17:26]<-NA


