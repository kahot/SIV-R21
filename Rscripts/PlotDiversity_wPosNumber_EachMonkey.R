#Diversity comparison
library(ggplot2)
library(reshape)
library(gridExtra)
source("Rscripts/label_scientific.R")

# read the files saved in Overview_output:
OverviewFiles2<-list.files("Output/Overview2/",pattern=".csv")

OvDF2<-list()
for (i in 1:length(OverviewFiles2)){ 
    overviews<-read.csv(paste0("Output/Overview2/",OverviewFiles2[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF2[[i]]<-overviews
    names(OvDF2)[i]<-gsub(".csv",'',OverviewFiles2[i])
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

#Plot Diversity with position info 
cutoff<-0.005
for (j in 1:length(monkeyList)){
    sample<-monkeyList[[j]]
    sample = sample[order(sample[,'Week']),]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF2[sample$filename]    

    Plots<-list()
    summary<-sample[,c("Monkey","Tissue","Week","filename")]
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #missing sites
        missing<-df[is.na(df$freq.mutations),]
        missing<-missing[missing$ref251.pos>=95,]
        missing<-missing[missing$AA251pos>(min(missing$AA251pos)+1)&missing$AA251pos<max(missing$AA251pos),]
        
        df2<-df[,c("AA239pos","ref251.pos","ns")]
        colnames(df2)[3]<-"freq"
        df2$Type<-"nonsyn"
        df2<-df2[!is.na(df2$freq),]
        
        df3<-df[,c("AA239pos","ref251.pos","syn")]
        colnames(df3)[3]<-"freq"
        df3$Type<-"syn"
        df3<-df3[!is.na(df3$freq),]
        
        mut<-rbind(df2,df3)
        mut2<-melt(mut, id.vars=c("AA239pos","ref251.pos","Type"))
        mut2<-mut2[,-which(colnames(mut2)=="variable")]
        #mut2$Type<-factor(mut2$Type, levels=c("syn", "nonsyn"))
        
        #freq<CUTOFF=NA
        mut2$value[mut2$value<cutoff]<-NA
        
      
        
        #barplot
        if (nrow(missing)==0){
            p<-ggplot(data=mut2, aes(x=AA239pos, y=value, color=Type, fill=Type))+
                geom_bar(stat = "identity", width=0.1)+
                scale_color_manual(values=c("blue","red"),guide = 'none')+
                scale_fill_manual(values=c("blue","red"),guide = 'none')+
                scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
                scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,10,20,30,40,50),  limits=c(0,0.5))+
                theme_bw()+
                ggtitle(paste0(monkey," Week ",sample$Week[i]," ",sample$Sample[i]))+
                theme(plot.title = element_text(size=11))+
                theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
                ylab("% Diversity")+
                xlab("Env codon position")+
                theme(legend.title = element_blank())+
                geom_text(data=mut2[mut2$value>0.04,], aes(label=AA239pos), hjust=0.3, vjust=-.8, size=2.5,show.legend = FALSE)+
                geom_text(data=mut2[mut2$value>0.04,], aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray30",show.legend = FALSE)+
                annotate(geom="text", x=36, y=.5, hjust=0,label="nonsyn",color ='blue', size=2.5)+
                annotate(geom="text", x=36,  y=.45, hjust=0, label="syn",color ='red', size=2.5)+
                annotate("segment", x = 30, xend = 34, y = .5, yend = .5, colour = "blue") +
                annotate("segment", x = 30, xend = 34, y = .45, yend = .45, colour = "red") 
            
        }
        else{
            p<-ggplot(data=mut2, aes(x=AA239pos, y=value, color=Type, fill=Type))+
            geom_rect(data=missing, inherit.aes=FALSE,
                      aes(xmin=AA239pos,xmax=AA239pos,ymin=-Inf,ymax=Inf),
                      color="gray90",size=1,alpha=0.2)+
            geom_bar(stat = "identity", width=0.1)+
            scale_color_manual(values=c("blue","red"),guide = 'none')+
            scale_fill_manual(values=c("blue","red"),guide = 'none')+
            scale_x_continuous(breaks=c(50,75,100,125,150,175,200,225,250,275), limits=c(30,280))+
            scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),labels=c(0,10,20,30,40,50),  limits=c(0,0.5))+
            theme_bw()+
            ggtitle(paste0(monkey," Week ",sample$Week[i]," ",sample$Sample[i]))+
            theme(plot.title = element_text(size=11))+
            theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
            ylab("% Diversity")+
            xlab("Env codon position")+
            theme(legend.title = element_blank())+
            geom_text(data=mut2[mut2$value>0.04,], aes(label=AA239pos), hjust=0.3, vjust=-.8, size=2.5,show.legend = FALSE)+
            geom_text(data=mut2[mut2$value>0.04,], aes(label=ref251.pos), hjust=-0.4, vjust=0, size=2.5, color="gray30",show.legend = FALSE)+
            annotate(geom="text", x=36, y=.5, hjust=0,label="nonsyn",color ='blue', size=2.5)+
            annotate(geom="text", x=36,  y=.45, hjust=0, label="syn",color ='red', size=2.5)+
            annotate("segment", x = 30, xend = 34, y = .5, yend = .5, colour = "blue") +
            annotate("segment", x = 30, xend = 34, y = .45, yend = .45, colour = "red") 
            
        }
        Plots[[i]]<-p
    }    
    png(paste0("Output/Diversity/Diversity_",monkey,".png"), width = 4, height=nrow(sample)*2, units = "in",res=300)
    do.call(grid.arrange, c(Plots, ncol=1))
    dev.off()
}

