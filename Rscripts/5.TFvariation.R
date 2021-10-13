#Assess Transmitted/Founder Variants in early plasma samples

library(ggplot2)
library(reshape)
library(gridExtra)
library(gtable)
library(grid)
library(colorspace)
colors<-qualitative_hcl(6, palette="Pastel1")
MFcolors<-c("#fb8072","#FF9300","#9437FF")

# read the overview files
OverviewFiles<-list.files("Output/Overview/",pattern=".csv")

OvDF<-list()
for (i in 1:length(OverviewFiles)){ 
    overviews<-read.csv(paste0("Output/Overview/",OverviewFiles[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF[[i]]<-overviews
    names(OvDF)[i]<-gsub("_overview.csv",'',OverviewFiles[i])
}


#sample info
SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]


# Plasma less than <4wk old files 
early<-samples[samples$Week<4,]
TF<-OvDF[early$filename]

#create a data frame to plot all together using facet_grid

freq<-data.frame()
miss<-data.frame()
for (i in 1:length(TF)){
    id<-early$filename[i]
    df<-TF[[id]]
    #compile syn and ns freq info
    df[which(df$TotalReads<100),17:26]<-NA
    df<-df[df$AA251pos<276,]
    
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
    
    df2<-df[,c("AA239pos","ref251.pos","ns")]
    colnames(df2)[3]<-"freq"
    df2$Type<-"nonsyn"
    df2<-df2[!is.na(df2$freq),]
    
    df3<-df[,c("AA239pos","ref251.pos","syn")]
    colnames(df3)[3]<-"freq"
    df3$Type<-"syn"
    df3<-df3[!is.na(df3$freq),]
    
    mut<-rbind(df2,df3)
    mut$Cohort<-early$Cohort[i]
    mut$Week<-early$Week[i]
    mut$Monkey<-early$Monkey[i]
    
    freq<-rbind(freq, mut)
    
    #missing sites
    missing<-df[is.na(df$freq.mutations),c(1,54)]
    missing<-missing[missing$ref251.pos>=95,]
    missing<-missing[missing$AA239pos>(min(missing$AA239pos, na.rm=T)+1)& missing$AA239pos<max(missing$AA239pos, na.rm=T),]
    #If no missing sites, create an empty data frame
    if (nrow(missing)==0){
        missing<-df[1,c(1,54)]
        missing$AA239pos[1]<-25
        }
    
    missing$Cohort<-early$Cohort[i]
    missing$Week<-early$Week[i]
    missing$Monkey<-early$Monkey[i]
    
    miss<-rbind(miss, missing)
}


#eliminate super low freq for a cleaner plot
freq$freq[freq$freq<=0.005]<-NA
freq<-freq[freq$AA239pos<272,]
freq<-freq[!is.na(freq$AA239pos),]

#remove nt395 position
freq<-freq[freq$ref251.pos!=395,]

#remove nt400 (not in stock and freq=1 in some)
freq<-freq[freq$ref251.pos!=400,]

freq$ID<-paste(freq$Monkey, 'week',freq$Week)

#order the plots 
freq$ID<-factor(freq$ID, levels=c("3616 week 2", "16314 week 3","3116 week 3","20615 week 3",  "3216 week 3" , "3516 week 3" ))
freq$Cohort<-factor(freq$Cohort,levels=c("SIV only", "Mtb NR", "Mtb R"))

miss$ID<-paste(miss$Monkey, 'week',miss$Week)
miss$ID<-factor(miss$ID, levels=c("3616 week 2", "16314 week 3","3116 week 3","20615 week 3",  "3216 week 3" , "3516 week 3" ))
miss$Cohort<-factor(miss$Cohort,levels=c("SIV only", "Mtb NR", "Mtb R"))

#color assignment
facetcols<-c(colors[5],colors[1],colors[3],colors[1],colors[3],colors[1])
facetcols<-paste0(facetcols, "66")
p1<-ggplot()+
    geom_rect(data=miss, inherit.aes=FALSE,
              aes(xmin=AA239pos,xmax=AA239pos,ymin=-Inf,ymax=Inf),
              color="gray90",size=1,alpha=0.2)+
    geom_bar(data=freq, aes(x=AA239pos, y=freq, color=Type, fill=Type), stat = "identity", width=0.1)+
    scale_color_manual(values=MFcolors[2:3],guide = 'none')+
    scale_fill_manual (values=MFcolors[2:3])+
    scale_x_continuous(breaks=c(50,100,150,200,250), limits=c(30,280))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),labels=c(0,25,50,75,100), limits=c(0,1))+
    theme_bw()+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Divergence")+
    xlab("Env codon position")+
    theme(legend.title = element_blank())+
    facet_wrap(~ID, scale='free',ncol=2)


gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
}

gtable_select <- function (x, ...) 
{
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
}


dummy <- ggplot(data=freq, aes(x=AA239pos, y=freq, color=Type, fill=Type))+
        facet_wrap(~ID, scale='free',ncol=2)+
        geom_rect(aes(fill=factor(Cohort)), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.3) +
        scale_fill_manual(values = colors[c(5,3,1)])+theme_minimal() 
    
    
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip-t", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips)
new_plot <- gtable_stack(g1, new_strips)


png("Output/Figures/T_F_divergence.png", width = 5.5, height=5, units="in", res=300) 
grid.newpage()
grid.draw(new_plot)
dev.off()


#create legend
ggplot(data=freq, aes(x=AA239pos, y=freq, color=Type, fill=Type))+
facet_wrap(~ID, scale='free',ncol=2)+
    geom_rect(aes(fill=factor(Cohort)), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    scale_fill_manual(values = colors[c(5,3,1)])+theme_minimal()+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Legend.png", width = 5.5, height=5, units="in", dpi=300)

ggplot(data=freq, aes(x=AA239pos, y=freq, color=Type, fill=Type))+
    facet_wrap(~ID, scale='free',ncol=2)+
    geom_line() +
    scale_color_manual(values = MFcolors[2:3])+theme_minimal()+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Legend2.png", width = 5.5, height=5, units="in", dpi=300)


###################################
### Initial Transmission bottleneck 
#using 6 samples at week 2/3

tf.div<-data.frame(file.name=early$filename)
for (i in 1:length(TF)){
    id<-early$filename[i]
    df<-TF[[id]]
    df[which(df$TotalReads<100),17:26]<-NA
    
    #remove the missing cetner region
    df<-df[df$AA251pos<276,]
    df<-df[df$ref251.pos<395|df$ref251.pos>540,]
    df<-df[!is.na(df$freq.mutations),]
    
    tf.div$Mean[i]<-mean(df$freq.mutations, na.rm=T)
    tf.div$SE[i]<-sqrt(tf.div$Mean[i]*(1-tf.div$Mean[i])/sum(df$TotalReads, na.rm=T))
    tf.div$read.mean[i]<-mean(df$TotalReads, na.rm=T)
}

#Mean of 6 samples
M<-mean(tf.div$Mean)
M*100
#0.2138419
SE<-mean(tf.div$SE)
SE*100
#0.0005083513


