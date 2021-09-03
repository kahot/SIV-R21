library(ggplot2)
library(reshape2)
library(DataCombine)
library(colorspace)

colors<-qualitative_hcl(6, palette="Dark3")
MFcolors<-c("#EB4E61","#FF9300","#9437FF")

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)

Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]

Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")
 

### 1) by week (all aggregate)
weeks<-Sum21[,c("Week", "mean", "mean.syn","mean.ns")]
weeks<-InsertRow(weeks,c(0, stks[1,"mean"],stks[1,"mean.syn"],stks[1,"mean.ns"]),1)
#weeks[1,c(1)]<-"Stock"
by.weeks<-aggregate(weeks[,2:4],by=list(weeks$Week),mean,na.rm=T )
colnames(by.weeks)[1]<-"Week"
weeksm<-melt(weeks, id.vars=c("Week"))
weekm<-melt(by.weeks, id.vars = "Week")

#cols<-c('#00BA38','#F8766D','#619CFF')

ggplot()+
    geom_point(data=weeksm, aes(x=Week, y=value, color=variable), position = position_dodge(width = 0.5), size=1.2,alpha = 0.4 )+
    geom_point(data=weekm, aes(x=Week, y=value, color=variable), position = position_dodge(width = 0.5),size=2)+
    geom_line(data=weekm, aes(x=Week, y=value, group=factor(variable),color=variable),position = position_dodge(width = 0.5),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0,0.001,0.002,0.003,0.004,0.005),labels=c(0,0.1,0.2,0.3,0.4,0.5),  limits=c(0,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8),labels=c("stock", 2, 4, 6, 8))+
    theme_bw()+
    #ggtitle("R21 Diversity")+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Average diversity")+xlab('Weeks post-infection')+
    scale_color_manual(values=MFcolors[c(1,3,2)], labels=c("Total", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byWeek_R21.png", units="in", width = 4.5, height = 2.5, dpi=300)


## By week with Sergio's data
#Sergio's data 
SumS<-read.csv("Output/Diversity/Diversity_summary_Ita.csv", stringsAsFactors = F, row.names = 1)
stks2<-read.csv("Output/Diversity/Stock_virus_diversity_summary_R21_Ita.csv", row.names = 1)

#by week
weeks2<-SumS[,c("Week", "mean", "mean.syn","mean.ns")]
weeks2<-InsertRow(weeks2,c(0, stks2[2,3],stks2[2,5],stks2[2,7]),1)
#weeks2[1,1]<-"Stock"
by.weeks2<-aggregate(weeks2[,2:4],by=list(weeks2$Week),mean,na.rm=T )
colnames(by.weeks2)[1]<-"Week"
weeksm2<-melt(weeks2, id.vars="Week")
weekm2<-melt(by.weeks2, id.vars = "Week")

#Combine R21 and Ita data
weeksm$Study<-"SIV-Mtb"
weeksm2$Study<-"Ita"
weekm$Study<-"SIV-Mtb"
weekm2$Study<-"Ita"

weeksm2<-weeksm2[weeksm2$Week<=10,]
weekm2<-weekm2[weekm2$Week<=10,]


Wksm<-rbind(weeksm, weeksm2)
Wkm<-rbind(weekm,weekm2)

ggplot(data=Wkm, aes(x=Week, y=value, color=variable,shape=Study,group=interaction(variable,Study)))+
    geom_point(data=Wksm, aes(x=Week, y=value, color=variable, shape=Study), position = position_dodge(width = 0.5), size=1.2,alpha = 0.4 )+
    geom_point(size=2.5, position = position_dodge(width = 0.5))+
    geom_line (position = position_dodge(width = 0.5),linetype = 1, size=0.2,stat = "identity")+
    scale_y_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c(0,0.5,1,1.5,2),  limits=c(0,0.02))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8,10))+
    theme_bw()+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Weeks post-infection')+
    scale_color_manual(values=MFcolors[c(1,3,2)], labels=c("Total", "Syn mean", "NS mean"))+
    scale_shape_manual(values=c(17,16))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byWeek_both.png", units="in", width = 6.5, height = 4, dpi=300)



## By week -Plasma only 
means<-Sum21[Sum21$Tissue2=="Plasma"|Sum21$Tissue2=="Plasma_nex",c("Monkey","Week","Cohort","mean")]

#add week 0 (stock)
monkeys<-unique(means$Monkey)
adds<-data.frame(Monkey=monkeys, Week=rep(0, times=length(monkeys)), mean=rep(stks[1,"mean"], times=length(monkeys)))
adds<-merge(adds, means[,c("Monkey","Cohort")], by="Monkey", all.x=F, all.y=F)
adds<-unique(adds)
means<-rbind(adds, means)

means$Cohort<-factor(means$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))

ggplot()+
    geom_point(data=means, aes(x=Week, y=mean, group=factor(Monkey), color=factor(Cohort)),  size=1.5 )+
    geom_line(data=means, aes(x=Week, y=mean, group=Monkey, color=factor(Cohort)),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0.002,0.003,0.004,0.005),labels=c(0.2,0.3,0.4,0.5),  limits=c(0.0015,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8), labels=c("stock", 2, 4, 6, 8))+
    geom_point(data=means, aes(x=0, y=adds[1,"mean"]), color="red",  size=1.5 )+
    scale_color_manual(values=colors[c(5,3,1)])+
    theme_bw()+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Average diversity")+xlab('Weeks post-infection')+
    #scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'), labels=c("Mean", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byMonkey_R21_plasma_only.png", units="in", width = 4.5, height = 2.5, dpi=300)


# By cohort

df<-means[means$Cohort=="SIV only",]
ggplot()+
    geom_point(data=df, aes(x=Week, y=mean,  color=factor(Monkey), shape=Cohort), size=1.5 )+
    geom_line(data=df, aes(x=Week, y=mean, color=factor(Monkey)),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0.002,0.003,0.004,0.005),labels=c(0.2,0.3,0.4,0.5),  limits=c(0.0015,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8), labels=c("stock", 2, 4, 6, 8))+
    geom_point(data=df, aes(x=0, y=adds[1,"mean"]), color="red",  size=1.5 )+
    theme(legend.title = element_blank())

df2<-means[means$Cohort=="Mtb NR",]
ggplot()+
    geom_point(data=df2, aes(x=Week, y=mean, group=factor(Monkey), color=factor(Monkey), shape=Cohort), size=1.5 )+
    geom_line(data=df2, aes(x=Week, y=mean, group=Monkey, color=factor(Monkey)),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0.002,0.003,0.004,0.005),labels=c(0.2,0.3,0.4,0.5),  limits=c(0.0015,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8), labels=c("stock", 2, 4, 6, 8))+
    geom_point(data=df2, aes(x=0, y=adds[1,"mean"]), color="red",  size=1.5 )+
    theme(legend.title = element_blank())

df3<-means[means$Cohort=="Mtb R",]
ggplot()+
    geom_point(data=df3, aes(x=Week, y=mean, group=factor(Monkey), color=factor(Monkey), shape=Cohort), size=1.5 )+
    geom_line(data=df3, aes(x=Week, y=mean, group=Monkey, color=factor(Monkey)),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0.002,0.003,0.004,0.005),labels=c(0.2,0.3,0.4,0.5),  limits=c(0.0015,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8), labels=c("stock", 2, 4, 6, 8))+
    geom_point(data=df3, aes(x=0, y=adds[1,"mean"]), color="red",  size=1.5 )+
    theme(legend.title = element_blank())
