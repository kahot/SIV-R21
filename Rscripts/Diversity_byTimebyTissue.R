library(ggplot2)
library(reshape)
library(gridExtra)
library(DataCombine)
source("Rscripts/label_scientific.R")

SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)

Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]

Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")




#by week
weeks<-Sum21[,c(2,4,8,9,10)]
weeks<-InsertRow(weeks,c(0,0, stks[1,3],stks[1,5],stks[1,7]),1)
weeks[1,1]<-"Stock"
by.weeks<-aggregate(weeks[,3:5],by=list(weeks$Week),mean,na.rm=T )
colnames(by.weeks)[1]<-"Week"
weeksm<-melt(weeks, id.vars=c("Week","Monkey"))
weekm<-melt(by.weeks, id.vars = "Week")

cols<-c('#00BA38','#F8766D','#619CFF')
    
ggplot()+
    geom_point(data=weeksm, aes(x=Week, y=value, color=variable), position = position_dodge(width = 0.5), size=1.2,alpha = 0.4 )+
    geom_point(data=weekm, aes(x=Week, y=value, color=variable), size=2.5)+
    geom_line(data=weekm, aes(x=Week, y=value, group=factor(variable),color=variable),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0,0.001,0.002,0.003,0.004,0.005),labels=c(0,0.1,0.2,0.3,0.4,0.5),  limits=c(0,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8))+
    theme_bw()+
    ggtitle("R21 Diversity")+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Weeks post-infection')+
    scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'), labels=c("Mean", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
#ggsave("Output/Figures/Diversity_overTime_byWeek_R21.pdf", width = 6, height = 4)
#ggsave("Output/Figures/Diversity_overTime_byWeek_R21.png", units="in", width = 6.5, height = 4, dpi=300)


#Sergio's data 
SumS<-read.csv("Output/Diversity/Diversity_summary_Ita.csv", stringsAsFactors = F, row.names = 1)

#by week
weeks2<-SumS[,c(1,2,7,8,9)]
weeks2<-InsertRow(weeks2,c(0,0, stks[2,3],stks[2,5],stks[2,7]),1)
weeks2[1,1]<-"Stock"
by.weeks2<-aggregate(weeks2[,3:5],by=list(weeks2$Week),mean,na.rm=T )
colnames(by.weeks2)[1]<-"Week"
weeksm2<-melt(weeks2, id.vars=c("Week","Monkey"))
weekm2<-melt(by.weeks2, id.vars = "Week")

#Combine R21 and Ita data
weeksm$Study<-"R21"
weeksm2$Study<-"Ita"
weekm$Study<-"R21"
weekm2$Study<-"Ita"

weeksm2<-weeksm2[weeksm2$Week<=10,]
weekm2<-weekm2[weekm2$Week<=10,]


Wksm<-rbind(weeksm, weeksm2)
Wkm<-rbind(weekm,weekm2)

ggplot()+
    geom_point(data=Wksm, aes(x=Week, y=value, color=variable, shape=Study), position = position_dodge(width = 0.5), size=1.2,alpha = 0.4 )+
    geom_point(data=Wkm, aes(x=Week, y=value, color=variable,shape=Study), size=2.5, position = position_dodge(width = 0.5))+
    geom_line(data=Wkm, aes(x=Week, y=value, group=interaction(variable,Study),color=variable),stat = "identity", position = position_dodge(width = 0.5),linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0,0.005,0.01,0.015,0.02),labels=c(0,0.5,1,1.5,2),  limits=c(0,0.02))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8,10))+
    theme_bw()+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Weeks post-infection')+
    scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'), labels=c("Mean", "Syn mean", "NS mean"))+
    scale_shape_manual(values=c(17,16))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byWeek_both.pdf", width = 6.5, height = 4)
ggsave("Output/Figures/Diversity_overTime_byWeek_both.png", units="in", width = 6.5, height = 4, dpi=300)


### All pointns, Mean only, R21 only 
means<-Sum21[,c("Monkey","Week","mean")]
#add week 0 (stock)
monkeys<-unique(means$Monkey)
adds<-data.frame(Monkey=monkeys, Week=rep(0, times=length(monkeys)), mean=rep(stks[1,3], times=length(monkeys)) )

means<-rbind(adds, means)
#fill in the missing number as 8 for now
means$Week[is.na(means$Week)]<-8

ggplot()+
    geom_point(data=means, aes(x=Week, y=mean, color=factor(Monkey)), position = position_dodge(width = 0.5), size=1.2 )+
    geom_line(data=means, aes(x=Week, y=mean, group=Monkey, color=factor(Monkey)),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0.002,0.003,0.004,0.005),labels=c(0.2,0.3,0.4,0.5),  limits=c(0.0015,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8))+
    theme_bw()+
    ggtitle("R21 Diversity")+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Weeks post-infection')+
    #scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'), labels=c("Mean", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byMonkey_R21.pdf", width = 6, height = 4)
ggsave("Output/Figures/Diversity_overTime_byMonkey_R21.png", units="in", width = 6, height = 4, dpi=300)


#########  By Tissue ###########
tissues<-Sum21[,c(1,2,3,8,9,10, )]
#tissues$Week[is.na(tissues$Week)]<-8
tissues$Tissue[tissues$Tissue=="plasma"]<-"Plasma"
tissues$Tissue[tissues$Tissue=="plasma"&tissues$Week>5]<-"Plasma_nex"
tissues<-InsertRow(tissues,c(0,0, 0,stks[1,3],stks[1,5],stks[1,7]),1)
tissues[1,1]<-"Stock"
tissues[1,2]<-"Stock"

by.tissues<- aggregate(tissues[c(4:6)],by=list(tissues$Tissue),mean,na.rm=T )
colnames(by.tissues)[1]<-"Tissue"
by.tissues$Tissue<-factor(by.tissues$Tissue, levels=c("Stock", "Plasma", "Plasma_nex", "LN",'HLN',"Lung"))
tissues$Tissue<-factor(tissues$Tissue, levels=c("Stock", "Plasma", "Plasma_nex", "LN",'HLN',"Lung"))

ggplot()+
    geom_point(data=tissues, aes(x=Tissue, y=mean*100), color="lightblue",size=2)+
    geom_point(data=by.tissues, aes(x=Tissue, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByTissues_diversity.pdf", width = 5.5,height = 4)
ggsave("Output/Figures/ByTissues_diversity.png", width = 5.5,height = 4,units="in", dpi=300)


### by tissue with throcic LN and Peripheral LN division
tissues<-Sum21[,c(13, 2,4,8,9,10)]
tissues$Tissue2[tissues$Tissue2=="plasma"]<-"Plasma"
tissues$Tissue2[tissues$Tissue2=="Plasma"&tissues$Week>5]<-"Plasma_nex"
tissues<-InsertRow(tissues,c(0,0, 0,stks[1,3],stks[1,5],stks[1,7]),1)
tissues[1,1]<-"Stock"
tissues[1,2]<-"Stock"

by.tissues<- aggregate(tissues[c(4:6)],by=list(tissues$Tissue),mean,na.rm=T )
colnames(by.tissues)[1]<-"Tissue"
by.tissues$Tissue<-factor(by.tissues$Tissue, levels=c("Stock", "Plasma", "Plasma_nex", "pLN","Thoracic LN","HLN","Lung"))
tissues$Tissue<-factor(tissues$Tissue, levels=c("Stock", "Plasma", "Plasma_nex",  "pLN","Thoracic LN",'HLN',"Lung"))



ggplot()+
    geom_point(data=tissues, aes(x=Tissue, y=mean*100), color="lightblue",size=2)+
    geom_point(data=by.tissues, aes(x=Tissue, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByTissues_diversity2.pdf", width = 5.5,height = 4)
ggsave("Output/Figures/ByTissues_diversity2.png", width = 5.5,height = 4,units="in", dpi=300)




###SIV only vs. co-infected
SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

samples$coinfection<-"Coinfected"
samples$coinfection[samples$Notes=="SIV only"]<-"SIV only"

infec<-Sum21[,c(1,2,3,4,8,9,10)]
infec$Tissue[infec$Tissue=="plasma"]<-"Plasma"
infec$Tissue[infec$Tissue=="Plasma"&infec$Week>5]<-"Plasma_nex"

Siv<-merge(infec,samples[,c("filename","coinfection")])

sivMean<-aggregate(Siv[,"mean"],by=list(Siv$coinfection, Siv$Tissue),mean,na.rm=T )
colnames(sivMean)<-c("infection", "Tissue","Mean")
Siv$Tissue<-factor(Siv$Tissue, levels = c("Plasma","Plasma_nex","LN","HLN","Lung"))
sivMean$Tissue<-factor(sivMean$Tissue, levels = c("Plasma","Plasma_nex","LN","HLN","Lung"))

ggplot()+
    geom_point(data=Siv, aes(x=Tissue, y=mean*100, color=coinfection),alpha=0.6,position = position_dodge(width = 0.5))+
    geom_point(data=sivMean, aes(x=Tissue, y=Mean*100, fill=infection), size=1.5,position = position_dodge(width = 0.5), shape=4)+
    xlab('')+ylab('% Diversity')+
    #scale_y_continuous(breaks=c(0.0024,0.0026,0.0028,0.003,0.0032),labels=c(0.24,0.26,0.28,0.3,0.32),  limits=c(0.0024,0.0032))+
    theme(legend.title = element_blank())+
    scale_fill_manual(values=c(1,1),guide="none")+
    theme(theme(legend.title = element_blank()))+theme_bw()+
    theme(legend.title = element_blank(), panel.grid.major.x = element_blank())
ggsave("Output/Figures/SIV.vs.coninfected.pdf", width = 5.5,height = 3.5)
ggsave("Output/Figures/SIV.vs.coninfected.png", width = 5.5,height = 3.5, dpi=300)

#######
### Test if one tissue is different from another

#1. Plasma_nex vs. LN
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

OverviewFiles2<-list.files("Output/Overview2/",pattern=".csv")
OvDF2<-list()
for (i in 1:length(OverviewFiles2)){ 
    overviews<-read.csv(paste0("Output/Overview2/",OverviewFiles2[i]),stringsAsFactors=FALSE, row.names = 1)
    OvDF2[[i]]<-overviews
    names(OvDF2)[i]<-gsub(".csv",'',OverviewFiles2[i])
}


Results<-data.frame()
for (j in 1:length(monkeys)){
    sample<-monkeyList[[j]]
    sample$Tissue[sample$Tissue=="plasma"]<-"Plasma"
    sample$Tissue[sample$Tissue=="Plasma"&sample$Week>5]<-"Plasma_nex"
    monkey<-names(monkeyList)[j]
    Ov<-OvDF2[sample$filename]    
    
    lis<-lapply(Ov, "[", "freq.mutations")
    for (i in 1:length(lis)){
        #colnames(lis[[i]])<-sample$Tissue[i]
        names(lis)[i]<-sample$Tissue[i]
    }
    
    LN<-c()
    HLN<-c()
    Lung<-c()
    Plasma_nex<-c()
    
    
    for (i in 1: length(lis)){
        df<-lis[[i]]
        if (names(lis)[i]=="LN"){
            LN<-c(LN, df$freq.mutation)
        }
        if (names(lis)[i]=="HLN"){
            HLN<-c(HLN, df$freq.mutation)
        }
        if (names(lis)[i]=="Lung"){
            Lung<-c(Lung, df$freq.mutation)
        }
        if (names(lis)[i]=="Plasma_nex"){
            Plasma_nex<-c(Plasma_nex, df$freq.mutation)
        }
  
    }
    
    #wicoxon test
    tis<-unique(sample$Tissue)
    tis<-tis[tis!="Plasma"]
    comb<-combn(tis,2)
    comb<-t(comb)
    
    wil.res<-data.frame(matrix(nrow=nrow(comb), ncol=0))
    for (i in 1:nrow(comb)){
        v1<-get(paste0(comb[i,1]))
        v2<-get(paste0(comb[i,2]))
        
        re1<-wilcox.test(v1, v2, alternative ="two.sided")
        
        wil.res$test[i]<-paste0(comb[i,1],".vs.", comb[i,2])
        wil.res$pvalue[i]<-re1[[3]]
        wil.res$mean1[i]<-mean(v1, na.rm=T)
        wil.res$mean2[i]<-mean(v2, na.rm=T)
        #which is higher in diversity
        wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
    }
    wil.res$Monkey<-monkey
    write.csv(wil.res, paste0("Output/Diversity/Tests/",monkey,"_wilcox_results_tissues.csv"))
    Results<-rbind(Results, wil.res)
}

write.csv(Results, "Output/Diversity/Tests/Wilcox_results_all.csv")

#summarize test results
tests<-unique(Results$test)
divsum<-data.frame(test=tests)
for (i in 1:length(tests)){
    res1<-Results[Results$test==tests[i],]
    res1<-res1[res1$pvalue<=0.05,]
    if (nrow(res1)==0) {divsum$Sig.sample.no[i]<-0
        divsum$higher1[i]<-0
        divsum$higher1_no[i]<-0
        divsum$higher2[i]<-0
        divsum$higher2_no[i]<0}
    else {
        divsum$Sig.sample.no[i]<-nrow(res1)
        counts<-data.frame(table(res1$higher))
        divsum$higher1[i]<-paste(counts[1,1])
        divsum$higher1_no[i]<-counts[1,2]
        divsum$higher2[i]<-paste(counts[2,1])
        divsum$higher2_no[i]<-counts[2,2]
    }
}

write.csv(divsum, "Output/Diversity/Tests/Wilcox.test.summary.byTest.csv")


#plasma vs. LN
tests<-unique(Results$test)
res1<-Results[Results$test==tests[1],]
res1[res1$pvalue<=0.05,]
#               test       pvalue       mean1       mean2     higher Monkey
#7  Plasma_nex.vs.LN 2.448235e-02 0.002806766 0.002674381 Plasma_nex  20615
#13 Plasma_nex.vs.LN 3.069376e-05 0.003642822 0.002893473 Plasma_nex  30816
#16 Plasma_nex.vs.LN 4.058204e-03 0.002389570 0.002448629         LN   3116
#25 Plasma_nex.vs.LN 4.586825e-02 0.002290543 0.002381386         LN  31316
#34 Plasma_nex.vs.LN 5.719643e-04 0.002640443 0.002672089         LN   3316
#43 Plasma_nex.vs.LN 1.070064e-02 0.002607587 0.002702677         LN   3616

res2<-Results[Results$test==tests[2],]
res2[res2$pvalue<=0.05,]
#                test       pvalue       mean1       mean2     higher Monkey
#18 Plasma_nex.vs.HLN 0.0002682888 0.002389570 0.001984572 Plasma_nex   3116
#35 Plasma_nex.vs.HLN 0.0004109602 0.002640443 0.002651330        HLN   3316
#38 Plasma_nex.vs.HLN 0.0009258810 0.003068068 0.002409447 Plasma_nex   3516

res2<-Results[Results$test==tests[2],]
res2[res2$pvalue<=0.05,]




### Granuloma or not ####
### by tissue with throcic LN and Peripheral LN division
gran<-Sum21[,c(13:15,2,4,8,9,10)]
gran<-gran[!is.na(gran$Granuloma),] #50 samples

by.gran<- aggregate(gran[c(4:6)],by=list(gran$Granuloma),mean,na.rm=T )

ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100), color="lightblue",size=2)+
    geom_point(data=by.gran, aes(x=Group.1, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByGranuloma_diversity.pdf", width = 3,height = 3)
ggsave("Output/Figures/ByGranuloma_diversity.png", width = 3,height = 3,units="in", dpi=300)

# Look at Latent R cohort
granR<-gran[gran$Cohort=="Latent R",]
ggplot()+
    geom_point(data=granR, aes(x=Granuloma, y=mean*100), color="darkturquoise",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/LatentR_Granuloma_div.pdf", width = 3,height = 3)

ggplot()+
    geom_point(data=granR, aes(x=Granuloma, y=mean*100, color=Tissue2), size=2,
               position=position_dodge(width = .7))+
    xlab('')+ylab('% Diversity')+
    ggtitle("Latent R Cohort: Granuloma or not")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank())
ggsave("Output/Figures/LatentR_Granuloma_byTissue_div.pdf", width = 3,height = 3)


granNR<-gran[gran$Cohort=="Latent NR",]
ggplot()+
    geom_point(data=granNR, aes(x=Granuloma, y=mean*100), color="darkturquoise",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/LatentNR_Granuloma_div.pdf", width = 3,height = 3)

ggplot()+
    geom_point(data=granNR, aes(x=Granuloma, y=mean*100, color=Tissue2), size=2,
               position=position_dodge(width = .7))+
    xlab('')+ylab('% Diversity')+
    ggtitle("Latent NR Cohort: Granuloma or not")+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank())
ggsave("Output/Figures/LatentNR_Granuloma_byTissue_div.pdf", width = 3,height = 3)


#Latent R vs. Latent NR
lat<-Sum21[,c(15,2,4,8,9,10)]
by.lat<- aggregate(lat[c(4:6)],by=list(lat$Cohort),mean,na.rm=T )

ggplot()+
    geom_point(data=lat, aes(x=Cohort, y=mean*100), color="lightblue",size=2)+
    geom_point(data=by.lat, aes(x=Group.1, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByLatentRNR_diversity.pdf", width = 3.5,height = 3)
ggsave("Output/Figures/ByLatentRNR_diversity.png", width = 3.5,height = 3,units="in", dpi=300)

##### Others
rna<-Sum21[,c(17:22,2,4,8,9,10)]
rna1<-rna[!is.na(rna$SIV.RNA.per.granuloma),]
rna2<-rna[!is.na(rna$SIV.RNA.per.tissue),]

plot(rna1$SIV.RNA.per.granuloma, rna1$mean, pch=16, col="blue")
plot(rna2$SIV.RNA.per.tissue, rna2$mean, pch=16, col="blue")
cor.test(rna1$SIV.RNA.per.granuloma, rna1$mean, method = "kendall") #P=1
cor.test(rna2$SIV.RNA.per.tissue, rna2$mean, method = "kendall") #P= 0.8869

ggplot(rna1, aes(x=SIV.RNA.per.granuloma, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/Correlation_RNA_Div.pdf", height = 3.6, width = 4)
ggplot(rna2, aes(x=SIV.RNA.per.tissue  , y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/Correlation_RNAtissue_Div.pdf", height = 3.6, width = 4)

cd4<-rna[!is.na(rna$CD4.percent),]
ggplot(cd4, aes(x=CD4.percent, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/CD4.pdf", height = 3.6,width = 4)
cor.test(cd4$CD4.percent, cd4$mean, method = "kendall") #P=0.7701

cd8<-rna[!is.na(rna$CD8.percent),]
ggplot(cd8, aes(x=CD8.percent, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/CD8.pdf", height = 3.6,width = 4)
cor.test(cd8$CD8.percent, cd8$mean, method = "kendall") # p=0.2985


# 31316, 3116 vs. 30816 vs. rest of ininfected

coinf<-Sum21[Sum21$Cohort!="SIV only",]
coinf$Type<-"typical"
coinf$Type[coinf$Monkey==31316|coinf$Monkey==3116]<-"high"
coinf$Type[coinf$Monkey==30816]<-"low"

aggregate(coinf$SIV.RNA.per.granuloma, by=list(coinf$Type), mean, na.rm=T )
#  Group.1         x
#1    high  726379.0
#2     low   70479.6
#3 typical 2393405.5
aggregate(coinf$SIV.RNA.per.tissue, by=list(coinf$Type), mean, na.rm=T )
#  Group.1        x
#1    high 80510182
#2     low 12845897
#3 typical 24964405

table(coinf[,c("Type","Granuloma")])
#         Granuloma
#Type       N  Y
#  high     5  4
#  low      2  2
#  typical 11 16

aggregate(coinf$CD4.percent, by=list(coinf$Type), mean, na.rm=T )
#  Group.1         x
#1    high 0.4100000
#2     low 0.3550000
#3 typical 0.2661111
aggregate(coinf$CD8.percent, by=list(coinf$Type), mean, na.rm=T )
#  Group.1         x
#1    high 0.4900000
#2     low 0.4700000
#3 typical 0.4766667

#CFU
cfu<-Sum21[,c(16,2,4,8,9,10)]
cfu<-cfu[!is.na(cfu$CFU),]
ggplot(cfu, aes(x=CFU, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/CFU_Div.pdf", height = 3.6, width = 4)

#######
##
all_mean<-Sum21[,c(13:15,2,4,8,9,10)]
all_mean<-InsertRow(all_mean,c(0,0,0,0,0, stks[1,3],stks[1,5],stks[1,7]),1)
all_mean[1,1:4]<-c("Stock","Stock","Stock","Stock")
all_mean$label<-paste0(all_mean$Monkey,'_', all_mean$Tissue2,'_', all_mean$Week)
all_mean$Tissue2[all_mean$Tissue2=="plasma"]<-"Plasma"
all_mean$Tissue2<-factor(all_mean$Tissue2, levels=c("Stock", "Plasma", "pLN","Thoracic LN","Lung"))

#order by 'stock,'siv only, R, NR
#library(Hmisc)
#library(gdata)

all_mean$Cohort<-factor(all_mean$Cohort, levels=c("Stock", "SIV only", "Latent R", "Latent NR"))

all_mean<-all_mean[order(all_mean$Cohort, all_mean$Monkey,all_mean$Week, all_mean$Tissue2),]

all_mean$label<-factor(all_mean$label, levels=unique(all_mean$label))


library(dplyr)
dat<-distinct(all_mean, label, .keep_all=T) #49 rows

numb<-data.frame(table(dat$Monkey))
colnames(numb)<-c("Monkey","Freq")
morder<-unique(all_mean$Monkey)
numb$Monkey<-factor(numb$Monkey, levels=morder)
numb<-numb[order(numb$Monkey),]

numb$col<-numb$Freq
for (i in 2:nrow(numb)){
    numb$col[i]<-numb$Freq[i]+numb$col[i-1]
}
numb$lab_pos<-numb$Freq/2    
for (i in 2:nrow(numb)){
    numb$lab_pos[i]<-(numb$col[i]-numb$col[i-1])/2+numb$col[i-1]
}

pos1<-(numb$col[5]-1)/2+0.5
pos2<-(numb$col[9]-numb$col[5])/2+numb$col[5]+0.5
pos3<-(50-numb$col[9])/2+numb$col[9]+0.5



library(colorspace)
hcl_palettes(plot = TRUE)
colors<-qualitative_hcl(5, palette="Pastel")
ggplot()+
    geom_point(data=all_mean, aes(x=label, y=mean*100, color=Tissue2),size=2,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Diversity')+
    scale_y_continuous(limits = c(0.1,0.5))+
    scale_x_discrete(labels=paste0(dat$Week))+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.1, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[9]+0.5, xmax=50.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=0.5, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=0.5, label="Latent R", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=0.5, label="Latent NR", hjust=0.5, size=3.5, color='gray20')
    #geom_text(aes(label="Week", x= 51, y=0.1), size=2.5, color="gray20",hjust=0)+
    #    coord_cartesian(clip = "off")+
    #theme(plot.margin = unit(c(2, 2, 1, 1), "lines"))
ggsave("Output/Figures/diversity_allSamples.pdf", width = 8,height = 3)
ggsave("Output/Figures/diversity_allSamples.png", width = 8,height = 3,units="in", dpi=300)
