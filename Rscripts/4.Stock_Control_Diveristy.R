# Assess the diversity of the stock and the control
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)

MFcolors<-c("#fb8072","#FF9300","#9437FF")

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



### 1. Amount of variation in control
control<-OvDF[["Run_6_01_Animal_stock_virus"]]
controlVar<-control[control$freq.mutations>0.005&control$TotalReads>100,] #236 positions
controlVar<-controlVar[!is.na(controlVar$ref251),] #5 positions
controlVar[,c("ref251.pos","freq.mutations", "TotalReads")]
# ref251.pos freq.mutations TotalReads
#        382    0.005953252      40986
#        394    0.012259326      13296
#        395    0.005123826       1171
#        524    0.005115090        391
#        819    0.005917160       9295

ctrl<-control[!is.na(control$freq.mutations),]
ctrl<-ctrl[ctrl$TotalReads>100,]
mean(ctrl$freq.mutations, na.rm = T)
#0.002032741
range(ctrl$freq.mutations, na.rm = T)
#0.00000000 0.01225933
median(ctrl$freq.mutations, na.rm = T)
#0.002012888

ggplot(ctrl, aes(x=AA239pos, y=freq.mutations*100))+
    geom_point(position=position_dodge(width=0.8), size=0.5, color="dodgerblue4", alpha=0.5)+
    theme_bw()+
    scale_y_continuous(limits=c(0,20))+
    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
    ylab("% Diversity")+
    xlab("Env codon position")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
    #ggtitle(paste0("SIV251 control "))
ggsave("Output/Diversity/Control.png",width = 4,height = 2.7, unit="in",dpi=300 )


ggplot(ctrl, aes(x=AA239pos, y=freq.mutations*100))+
    geom_point(position=position_dodge(width=0.8), size=0.5, color="dodgerblue4", alpha=0.5)+
    theme_bw()+
    scale_y_continuous(limits=c(0,1.23))+
    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
    ylab("% Diversity")+
    xlab("Env codon position")+
    geom_hline(yintercept = 0.5, color="red",size=0.2)+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave("Output/Diversity/Control2.png",width = 4,height = 2.7, unit="in",dpi=300 )




#################################
### 2) Diversity in stock 
# Look at the amount of variants in the stock
stockMf<-OvDF[[stock$filename]]

#replace the mut frea =NA for positios with total_reads<100
#% diversity of Stock
stockMf[which(stockMf$TotalReads<100),17:26]<-NA

#missing region aa251.pos 94, 397 to 540, 83, >826
stockMf[stockMf$ref251.pos>=395&stockMf$ref251.pos<=540,c(17:26)]<-NA


ggplot(stockMf, aes(x=AA239pos, y=freq.mutations))+
    geom_point(position=position_dodge(width=0.8), size=0.5)+
    theme_bw()+ylim(0,0.5)+
    ylab("% Diversity from consensus")+
    xlab("Env codon position")+
    ggtitle(paste0("SIV251 stock diversity "))+
    geom_hline(yintercept = 0.005, color="lightsteelblue3")
ggsave("Output/Diversity/Stock_diversity_R21.pdf", width = 5, height = 4)

mean(stockMf$freq.mutations, na.rm=T)
#0.002500958 
skpos<-stockMf[!is.na(stockMf$freq.mutations)&stockMf$freq.mutations>0.005,]

## Ita et al.(Sergio)'s stock virus:
stockMfS<-read.csv("Output/Overview/Ita/mmStock_overview.csv", stringsAsFactors = F, row.names = 1)
stockMfS[which(stockMfS$TotalReads<100),18:27]<-NA
#exclude the missing sites in R21 stock (pos397 to 540) & pos395
stockMfS2<-stockMfS[stockMfS$ref251.pos<395|stockMfS$ref251.pos>540,]
stockMfS2<-stockMfS2[stockMfS2$ref251.pos>94 & stockMfS2$ref251.pos<826,]
mean(stockMfS2$freq.mutations, na.rm=T)
#0.01543653


#Plot the stock diversity
df<-stockMf[stockMf$AA251pos<276,];title="SIV-Mtb"
cutoff=0.005
#missing sites
missing<-df[is.na(df$freq.mutations),]
missing<-missing[missing$ref251.pos>=95,]

#create buffers for plotting so that points don't appear in the gray rectangle
missing<-missing[missing$AA251pos>(min(missing$AA251pos)+1)&missing$AA251pos<max(missing$AA251pos),]

#calculate total ns freq at each position
df$ns1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="nonsyn") x["freq.Ts.ref"] else 0))
df$ns2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="nonsyn") x=x["freq.transv1.ref"] else 0))
df$ns3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="nonsyn") x=x["freq.transv2.ref"] else 0))
df$ns<-df$ns1+df$ns2+df$ns3
#calculate total syn freq at each position
df$syn1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="syn") x["freq.Ts.ref"] else 0))
df$syn2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="syn") x=x["freq.transv1.ref"] else 0))
df$syn3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="syn") x=x["freq.transv2.ref"] else 0))
df$syn<-df$syn1+df$syn2+df$syn3
#calculate total stop freq at each position
df$stop1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="stop") x["freq.Ts.ref"] else 0))
df$stop2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="stop") x=x["freq.transv1.ref"] else 0))
df$stop3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="stop") x=x["freq.transv2.ref"] else 0))
df$stop<-df$stop1+df$stop2+df$stop3

df2<-df[,c("AA239pos","ref251.pos","ns")]
colnames(df2)[3]<-"freq"
df2$Type<-"nonsyn"
df2<-df2[!is.na(df2$freq),]

df3<-df[,c("AA239pos","ref251.pos","syn")]
colnames(df3)[3]<-"freq"
df3$Type<-"syn"
df3<-df3[!is.na(df3$freq),]

#combine syn and nonsyn dataframes
mut<-rbind(df2,df3)
mut<-mut[mut$ref251.pos<395|mut$ref251.pos>540,]
mut2<-melt(mut, id.vars=c("AA239pos","ref251.pos","Type"))
mut2<-mut2[,-which(colnames(mut2)=="variable")]

#Plot only the sites with freq > cutoff to avoid too many points
mut2$value[mut2$value<cutoff]<-NA

ggplot(data=mut2, aes(x=AA239pos, y=value, color=Type, fill=Type))+
    geom_rect(data=missing, inherit.aes=FALSE,
              aes(xmin=AA239pos,xmax=AA239pos,ymin=-Inf,ymax=Inf),
              color="gray90",size=1,alpha=0.2)+
    geom_bar(stat = "identity", width=0.1)+
    scale_color_manual(values=MFcolors[2:3],guide = 'none')+
    scale_fill_manual(values=MFcolors[2:3], guide = 'none')+
    scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
    scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,10,20,30,40,50),  limits=c(0,0.5))+
    theme_bw()+
    #ggtitle(paste0(title," SIV251 stock diversity"))+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+
    xlab("Env codon position")+
    theme(legend.title = element_blank())+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos<250&mut2$AA239pos>62,], aes(label=AA239pos), hjust=0.3, vjust=-.8, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos<250&mut2$AA239pos>62,], aes(label=ref251.pos), hjust=-0.2, vjust=0, size=1.5, color="gray30",show.legend = FALSE)+
    #geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos>250,], aes(label=AA239pos), hjust=0.2, vjust=-1.5, size=2.5,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==254,], aes(label=AA239pos), hjust=0.5, vjust=-2.5, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==259,], aes(label=AA239pos), hjust=0.2, vjust=-1.7, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==261,], aes(label=AA239pos), hjust=0.2, vjust=-0.5, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==263,], aes(label=AA239pos), hjust=0.2, vjust=-2.5, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==272,], aes(label=AA239pos), hjust=-0.2, vjust=-1.5, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==57,], aes(label=AA239pos), hjust=0.2, vjust=-0.8, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==57,], aes(label=ref251.pos), hjust=-.7, vjust=0, size=1.5, color="gray30",show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==61,], aes(label=AA239pos), hjust=-0.2, vjust=-2.5, size=2,show.legend = FALSE)+
    geom_text(data=mut2[mut2$value>0.01&mut2$AA239pos==61,], aes(label=ref251.pos), hjust=-0.2, vjust=-1.5, size=1.5, color="gray30",show.legend = FALSE)+
    
    annotate(geom="text", x=36, y=.5, hjust=0,label="nonsyn",color =MFcolors[2], size=2.5)+
    annotate(geom="text", x=36,  y=.475, hjust=0, label="syn",color =MFcolors[3], size=2.5)+
    annotate("segment", x = 30, xend = 34, y = .5, yend = .5, colour = MFcolors[2]) +
    annotate("segment", x = 30, xend = 34, y = .475, yend = .475, colour = MFcolors[3]) 
ggsave(paste0("Output/Figures/R21_stock_diversity_barplot.png"),width =4, height = 3, units = "in", dpi=300)


#Create a table and plot the stock diversity (point plot) from both studies 
summary1<-data.frame(study=c("SIV-Mtb","Ita"))
cutoff<-0.005

for (i in 1:2){
    if (i==1) {df<-stockMf[stockMf$AA251pos<276&stockMf$ref251.pos>94,];title="SIV-Mtb"; cutoff=0.005}
    if (i==2) {df<-stockMfS[stockMfS$AA251pos<276&stockMfS$ref251.pos>94,];title="Ita";  cutoff=0.009}
    
    #missing sites
    missing<-df[is.na(df$freq.mutations),]
    missing<-missing[missing$ref251.pos>=95,]
    missing<-missing[missing$AA251pos>(min(missing$AA251pos)+1)&missing$AA251pos<max(missing$AA251pos),]
    
    #calculate total ns freq at each position
    df$ns1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="nonsyn") x["freq.Ts.ref"] else 0))
    df$ns2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="nonsyn") x=x["freq.transv1.ref"] else 0))
    df$ns3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="nonsyn") x=x["freq.transv2.ref"] else 0))
    df$ns<-df$ns1+df$ns2+df$ns3
    #calculate total syn freq at each position
    df$syn1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="syn") x["freq.Ts.ref"] else 0))
    df$syn2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="syn") x=x["freq.transv1.ref"] else 0))
    df$syn3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="syn") x=x["freq.transv2.ref"] else 0))
    df$syn<-df$syn1+df$syn2+df$syn3
    #calculate total stop freq at each position
    df$stop1<-as.numeric(apply(df[,c("Type.ref","freq.Ts.ref")],1, function(x) if (x["Type.ref"]=="stop") x["freq.Ts.ref"] else 0))
    df$stop2<-as.numeric(apply(df[,c("Type.tv1.ref","freq.transv1.ref")],1, function(x) if (x["Type.tv1.ref"]=="stop") x=x["freq.transv1.ref"] else 0))
    df$stop3<-as.numeric(apply(df[,c("Type.tv2.ref","freq.transv2.ref")],1, function(x) if (x["Type.tv2.ref"]=="stop") x=x["freq.transv2.ref"] else 0))
    df$stop<-df$stop1+df$stop2+df$stop3
    
    df2<-df[,c("AA239pos","ref251.pos","ns")]
    colnames(df2)[3]<-"freq"
    df2$Type<-"nonsyn"
    df2<-df2[!is.na(df2$freq),]
    
    df3<-df[,c("AA239pos","ref251.pos","syn")]
    colnames(df3)[3]<-"freq"
    df3$Type<-"syn"
    df3<-df3[!is.na(df3$freq),]
    
    sp<-df[,c("AA239pos","ref251.pos","stop")]
    colnames(sp)[3]<-"freq"
    sp$Type<-"stop"
    
    mut<-rbind(df2,df3)
    write.csv(mut,paste0("Output/Diversity/Div.mf.summary_",title,".csv"))
    mut2<-melt(mut, id.vars=c("AA239pos","ref251.pos","Type"))
    mut2<-mut2[,-which(colnames(mut2)=="variable")]
    #mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
    
    #Plot the sites with freq > cutoff to avoid too many points
    mut2$value[mut2$value<cutoff]<-NA
    
    #remove the center missing sections for average calculation
    mut<-mut[mut$ref251.pos<395|mut$ref251.pos>540,]
    
    
    #Create a summary table
    #1 all sites   
    summary1$no.sites[i]<-length(unique(mut$ref251.pos[mut$ref251.pos<395|mut$ref251.pos>540]))
    summary1$mean.freq.total[i]<-mean(df$freq.mutations.ref[df$ref251.pos<395|df$ref251.pos>540],na.rm=T)
    summary1$no.sites.syn[i]<-length(unique(mut$ref251.pos[mut$Type=="syn"&mut$freq!=0]))
    summary1$mean.syn.freq[i]<-mean(mut$freq[mut$Type=="syn"],na.rm=T)
    summary1$no.sites.ns[i]<-length(unique(mut$ref251.pos[mut$Type=="nonsyn"&mut$freq!=0]))
    summary1$mean.ns.freq[i]<-mean(mut$freq[mut$Type=="nonsyn"],na.rm=T)
    summary1$mean.freq.ns.plus.syn[i]<-summary1$mean.ns.freq[i]+summary1$mean.ns.freq[i]
    summary1$mean.stop.freq[i]<-mean(sp$freq, na.rm=T)
    summary1$no.sites.stop[i]<-length(unique(sp$ref251.pos[sp$freq!=0]))
    
    if (i==1) {
        p1<-ggplot(data=mut2, aes(x=AA239pos, y=value, color=factor(Type)))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA239pos,xmax=AA239pos-1,ymin=-Inf,ymax=Inf),
                      color="gray90",size=1,alpha=0.2, guide="none")+
            geom_point(position=position_dodge(width=0.8), size=0.8)+
            scale_color_manual(values=MFcolors[2:3], guide ='none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,10,20,30,40,50), limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(title))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Diversity")+
            xlab("Env codon position")
    }
    if (i==2){
        p2<-ggplot(data=mut2, aes(x=AA239pos, y=value, color=Type))+
            geom_point(position=position_dodge(width=0.8), size=0.6)+
            scale_color_manual(values=MFcolors[2:3])+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,10,20,30,40,50), limits=c(0,0.5))+
            theme_bw()+
            ggtitle("Ita et al. (2018)")+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())
    }
    
}


png("Output/Figures/Stock_diversity_bothStudies.png", unit="in",width = 8, height=3, res=300)
plot_grid(p1, p2, align = "v", ncol = 2, rel_widths = c(.75,1))
dev.off()

write.csv(summary1,"Output/Diversity/Stock_virus_diversity_summary_R21_Ita.csv")

#R21 numbers of site with diversity >0.5% 
mut1<-read.csv("Output/Diversity/Div.mf.summary_SIV-Mtb.csv", row.names = 1,stringsAsFactors = F)
mut1<-mut1[mut1$freq>0.005,]
table(mut1$Type)
#nonsyn    syn 
#12     13 
