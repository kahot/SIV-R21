library(ggplot2)
library(reshape)
library(gridExtra)
library(DataCombine)
source("Rscripts/label_scientific.R")
library(colorspace)
library(dplyr)

colors<-qualitative_hcl(6, palette="Dark3")
col_light<-qualitative_hcl(6, palette="Set3")
#hcl_palettes(plot = TRUE)
MFcolors<-c("#EB4E61","#FF9300","#9437FF")

SampleSheet<-read.csv("Data/SampleSheetMac251All.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)

Sum21<-summary[1:69,]
stks<-summary[summary$filename=="Run_5_01_Animal_stock_virus",]

Sum21<-merge(Sum21, samples[,c("filename","Granuloma","SIV.RNA.per.granuloma","SIV.RNA.per.tissue","CD4.percent","CD8.percent")], by="filename")


#by week (all aggregate)
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
    geom_point(data=weekm, aes(x=Week, y=value, color=variable), position = position_dodge(width = 0.5),size=2.5)+
    geom_line(data=weekm, aes(x=Week, y=value, group=factor(variable),color=variable),position = position_dodge(width = 0.5),stat = "identity",linetype = 1, size=0.2)+
    scale_y_continuous(breaks=c(0,0.001,0.002,0.003,0.004,0.005),labels=c(0,0.1,0.2,0.3,0.4,0.5),  limits=c(0,0.005))+
    scale_x_continuous(breaks=c(0,2, 4, 6,8),labels=c("stock", 2, 4, 6, 8))+
    theme_bw()+
    #ggtitle("R21 Diversity")+
    theme(plot.title = element_text(size=11))+
    theme(panel.grid.major.x=element_blank(),panel.grid.minor.x=element_blank())+
    ylab("% Diversity")+xlab('Weeks post-infection')+
    scale_color_manual(values=MFcolors[c(1,3,2)], labels=c("Total", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
ggsave("Output/Figures/Diversity_overTime_byWeek_R21.png", units="in", width = 6.5, height = 3.8, dpi=300)


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
ggsave("Output/Figures/Diversity_overTime_byWeek_both.pdf", width = 6.5, height = 4)
ggsave("Output/Figures/Diversity_overTime_byWeek_both.png", units="in", width = 6.5, height = 4, dpi=300)



#By week Plasma only 

means<-Sum21[Sum21$Tissue2=="plasma",c("Monkey","Week","Cohort","mean")]

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
    ylab("% Diversity")+xlab('Weeks post-infection')+
    #scale_color_manual(values=c('#00BA38','#F8766D','#619CFF'), labels=c("Mean", "Syn mean", "NS mean"))+
    theme(legend.title = element_blank())
#ggsave("Output/Figures/Diversity_overTime_byMonkey_R21_plasma_only.pdf", width = 6, height = 3.8)
ggsave("Output/Figures/Diversity_overTime_byMonkey_R21_plasma_only.png", units="in", width = 5, height = 3.2, dpi=300)



### By Cohort (SIV only vs. Latent NR vs. Latent R)
#### Wilcoxon test on diversity between cohort per tissues ###
samples2<-samples
samples2$Tissue2[samples2$Tissue2=="plasma"]<-"Plasma"
samples2$Tissue2[samples2$Tissue2=="Plasma"&samples2$Week>5]<-"Plasma PM"
samples2$Tissue2[samples2$Tissue2=="Thoracic LN"]<-"tLN"

OverviewFiles2<-list.files("Output/Overview2/",pattern=".csv")
OvDF2<-list()
for (i in 1:length(OverviewFiles2)){ 
    overviews<-read.csv(paste0("Output/Overview2/",OverviewFiles2[i]),stringsAsFactors=FALSE, row.names = 1)
    overviews<-overviews[overviews$ref251.pos!=395,]
    OvDF2[[i]]<-overviews
    names(OvDF2)[i]<-gsub(".csv",'',OverviewFiles2[i])
}
#comparison pairs 
cohorts<-c("SIV","NR","R")
comb<-combn(cohorts,2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))

tissues<-unique(samples2$Tissue2)

Results<-data.frame()
for (j in 1:length(tissues)){
    sample<-samples2[samples2$Tissue2==tissues[j],]
    Ov<-OvDF2[sample$filename] 
    
    SIV<-c()
    R<-c()
    NR<-c()
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #Remove the missing inner section from all files
        df<-df[df$ref251.pos<397|df$ref251.pos>540,]
        
        if (sample$Cohort[i]=="SIV only"){
            SIV<-c(SIV, df$freq.mutations)
        }
        if (sample$Cohort[i]=="Mtb NR"){
            NR<-c(NR, df$freq.mutations)
        }
        if (sample$Cohort[i]=="Mtb R"){
            R<-c(R, df$freq.mutations)
        }
    }
    
    #wicoxon test
    wil.res<-data.frame(Test=comb_pairs)
    
    for (i in 1:nrow(comb)){
        v1<-get(paste0(comb[i,1]))
        v2<-get(paste0(comb[i,2]))
        if (length(v1)==0|length(v2)==0){
            wil.res$pvalue[i]<-NA
            wil.res$mean1[i]<-NA
            wil.res$mean2[i]<-NA
            wil.res$higher[i]<-NA
        }
        else{
            re1<-wilcox.test(v1, v2, alternative ="two.sided")
            wil.res$pvalue[i]<-re1[[3]]
            wil.res$mean1[i]<-mean(v1, na.rm=T)
            wil.res$mean2[i]<-mean(v2, na.rm=T)
            #which is higher in diversity
            wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
        }
    }
    wil.res$Tissue<-tissues[j]
    wil.res<-wil.res[!is.na(wil.res$pvalue),]
    #write.csv(wil.res, paste0("Output/Diversity/Tests/",tissues[j],"_wilcox_results_tissues.csv"))
    Results<-rbind(Results, wil.res)
}
write.csv(Results, "Output/Diversity/Tests/Cohort_comparison_Wilcox_results_all.csv")

### Plot the significance symbols
Test.results<-data.frame(tissue=tissues)
for (i in 1:length(tissues)){
    
    df<-Results[Results$Tissue==tissues[i],]
    symb1<-NA
    symb2<-NA
    symb3<-NA
    if (df$pvalue[df$Test=='SIV.vs.NR']<0.05) symb1="a"
    if (df$pvalue[df$Test=='SIV.vs.R']<0.05) symb2="b"
    if (df$pvalue[df$Test=='NR.vs.R']<0.05) symb3="c"
    
    symb<-c(symb1, symb2, symb3) 
    symb<-symb[!is.na(symb)]
    symb<-paste(symb, collapse = ',' )
    if (length(symb)==0) symb<-''
    Test.results$Sig[i]<-symb
}


####### PLOT the RESULTS
###SIV only vs. Latent NR vs. Latent R
infec<-Sum21[,c("Cohort","Tissue2","Week","mean")]
colnames(infec)[2]<-"Tissue"
#Divide the plasam to earlier time points vs. at the end
infec$Tissue[infec$Tissue=="plasma"]<-"Plasma"
infec$Tissue[infec$Tissue=="Plasma"&infec$Week>5]<-"Plasma PM"
infec$Tissue[infec$Tissue=="pLN"]<-"Periferal LN"

#calcualte the averages
sivMean<-aggregate(infec[,"mean"],by=list(infec$Cohort, infec$Tissue),mean,na.rm=T )
colnames(sivMean)<-c("Cohort", "Tissue","Mean")

#Order the factors
infec$Tissue<-factor(infec$Tissue, levels = c("Plasma","Plasma PM","Periferal LN","Thoracic LN","Lung"))
sivMean$Tissue<-factor(sivMean$Tissue, levels = c("Plasma","Plasma PM","Periferal LN","Thoracic LN","Lung"))
infec$Cohort<-factor(infec$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
sivMean$Cohort<-factor(sivMean$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))

#add jitter to one of the lung 
which(infec$Tissue=="Lung"& infec$Cohort=="Mtb NR")
#64, 67
infec$mean[67]<-0.00195
infec$mean[64]<-0.00205

ggplot()+
    geom_point(data=infec, aes(x=Tissue, y=mean*100, color=Cohort),alpha=0.6,
               position=position_jitterdodge(jitter.width = 0.1,jitter.height = 0,
                dodge.width = 0.5), size=2)+
    geom_point(data=sivMean, aes(x=Tissue, y=Mean*100, fill=Cohort), size=1.7,position = position_dodge(width = 0.5), shape=4)+
    xlab('')+ylab('% Diversity')+
    #scale_y_continuous(breaks=c(0.0024,0.0026,0.0028,0.003,0.0032),labels=c(0.24,0.26,0.28,0.3,0.32),  limits=c(0.0024,0.0032))+
    theme(legend.title = element_blank())+
    scale_y_continuous(limits = c(0.16,0.49))+
    scale_color_manual(values=colors[c(5,3,1)])+
    scale_fill_manual(values=c(1,1,1),guide="none")+
    theme(theme(legend.title = element_blank()))+theme_bw()+
    theme(legend.title = element_blank(), panel.grid.major.x = element_blank())+
    annotate(geom="text", x=1:5, y=0.49,  label=Test.results$Sig, color="gray40", size=2.7 )

ggsave("Output/Figures/Cohort.png", width = 5.5,height = 3.3, dpi=300)


#######################################
### Diversity by cohort (not by tisseus) 
#wicoxon test
coh<-c("SIV only","Mtb NR","Mtb R")
for (j in 1:3){
    sample<-samples2[samples2$Cohort==coh[j],]
    Ov<-OvDF2[sample$filename] 
    vec<-c()
    #SIV<-c()
    #R<-c()
    #NR<-c()
    for (i in 1:length(Ov)){
        df<-Ov[[i]]
        #Remove the missing inner section from all files
        df<-df[df$ref251.pos<395|df$ref251.pos>540,]
        vec<-c(vec,df$freq.mutations)
    }
    assign(cohorts[j],vec)
}

wil.res<-data.frame(Test=comb_pairs)
for (i in 1:nrow(comb)){
    v1<-get(paste0(comb[i,1]))
    v2<-get(paste0(comb[i,2]))
    
    re1<-wilcox.test(v1, v2, alternative ="two.sided")
    wil.res$pvalue[i]<-re1[[3]]
    wil.res$mean1[i]<-mean(v1, na.rm=T)
    wil.res$mean2[i]<-mean(v2, na.rm=T)
    #which is higher in diversity
    wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
}
write.csv(wil.res, paste0("Output/Diversity/Tests/Cohorts_div_wilcox_results.csv"))



by.Coh<- aggregate(infec["mean"],by=list(infec$Cohort),mean,na.rm=T )
by.Coh$Group.1<-factor(by.Coh$Group.1, levels=c("SIV only", "Mtb NR", "Mtb R"))

ggplot()+
    geom_point(data=infec, aes(x=Cohort, y=mean*100, color=Cohort),size=2.3, position=position_jitter(width=0.05))+
    geom_point(data=by.Coh, aes(x=Group.1, y=mean*100), color="black",size=2.3, shape=4)+
    xlab('')+ylab('% Diversity')+
    scale_color_manual(values=paste0(colors[c(5,3,1)],"66"), guide='none')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    annotate(geom="text", x=3, y=0.49,  label="a,c", color="gray40", size=3 )
ggsave("Output/Figures/ByLatentRNR_diversity.png", width = 3,height = 3,units="in", dpi=300)





#######
### Test if one tissue is different from another in each monkey

#Change tissue names (Add plasam_nex and thorcic ->tLN)
samples2<-samples
samples2$Tissue2[samples2$Tissue2=="plasma"]<-"Plasma"
samples2$Tissue2[samples2$Tissue2=="Plasma"&samples2$Week>5]<-"Plasma_nex"
samples2$Tissue2[samples2$Tissue2=="Thoracic LN"]<-"tLN"

list.animal<-split(samples2, samples2$Monkey)
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
    overviews<-overviews[overviews$ref251.pos!=395,]
    OvDF2[[i]]<-overviews
    names(OvDF2)[i]<-gsub(".csv",'',OverviewFiles2[i])
}

#comparison pairs 
tis<-unique(samples2$Tissue2)
#tis<-tis[tis!="Plasma"]
comb<-combn(tis,2)
comb<-t(comb)
comb_pairs<-apply(comb, 1,function(x) paste0(x[1],".vs.",x[2]))
Results<-data.frame()
for (j in 1:length(monkeys)){
    sample<-monkeyList[[j]]
    monkey<-names(monkeyList)[j]
    Ov<-OvDF2[sample$filename]    
    
    lis<-lapply(Ov, "[", "freq.mutations")
    for (i in 1:length(lis)){
        #colnames(lis[[i]])<-sample$Tissue[i]
        names(lis)[i]<-sample$Tissue2[i]
    }
    
    pLN<-c()
    tLN<-c()
    Lung<-c()
    Plasma_nex<-c()
    Plasma<-c()
    
    for (i in 1: length(lis)){
        df<-lis[[i]]
        if (names(lis)[i]=="pLN"){
            pLN<-c(pLN, df$freq.mutation)
        }
        if (names(lis)[i]=="tLN"){
            tLN<-c(tLN, df$freq.mutation)
        }
        if (names(lis)[i]=="Lung"){
            Lung<-c(Lung, df$freq.mutation)
        }
        if (names(lis)[i]=="Plasma_nex"){
            Plasma_nex<-c(Plasma_nex, df$freq.mutation)
        }
        if (names(lis)[i]=="Plasma"){
            Plasma<-c(Plasma, df$freq.mutation)
        }
    }
    
    #wicoxon test

    
    wil.res<-data.frame(Test=comb_pairs)
    
    for (i in 1:nrow(comb)){
        v1<-get(paste0(comb[i,1]))
        v2<-get(paste0(comb[i,2]))
        if (length(v1)==0|length(v2)==0){
            wil.res$pvalue[i]<-NA
            wil.res$mean1[i]<-NA
            wil.res$mean2[i]<-NA
            wil.res$higher[i]<-NA
        }
        else{
            re1<-wilcox.test(v1, v2, alternative ="two.sided")
            wil.res$pvalue[i]<-re1[[3]]
            wil.res$mean1[i]<-mean(v1, na.rm=T)
            wil.res$mean2[i]<-mean(v2, na.rm=T)
            #which is higher in diversity
            wil.res$higher[i]<-ifelse((wil.res$mean1[i]-wil.res$mean2[i])>0, comb[i,1],comb[i,2])
        }
    }
    wil.res$Monkey<-monkey
    wil.res<-wil.res[!is.na(wil.res$pvalue),]
    write.csv(wil.res, paste0("Output/Diversity/Tests/",monkey,"_wilcox_results_tissues.csv"))
    Results<-rbind(Results, wil.res)
}

write.csv(Results, "Output/Diversity/Tests/Wilcox_results_all.csv")

#summarize test results
tests<-unique(Results$Test)
divsum<-data.frame(test=tests)
for (i in 1:length(tests)){
    res1<-Results[Results$Test==tests[i],]
    t<-nrow(res1)
    divsum$n[i]<-t
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

### 
#Attach metadata to Results for exporting
Re<-Results[Results$Test=="tLN.vs.Lung",]
sam<-unique(samples2[ , c("Monkey", "Cohort") ] )
Re<-merge(Re, sam, by="Monkey")

#sample info on the animals that had higher thoracic LN diversity than Lung
sa<-samples2[samples2$Monkey %in% Re$Monkey, ]
#sample size
size<-as.data.frame.matrix(table(sa$Monkey, sa$Tissue2))
size$Monkey<-rownames(size)
#attache the sample size
Re<-merge(Re, size, by="Monkey")
write.csv(Re, "Output/Diversity/Tests/tLN.Lung.Summary.csv")




### Granuloma or not ####
### by tissue 
gran<-Sum21[,c("Cohort","Tissue2","Week","mean","Granuloma","SIV.RNA.per.granuloma")]
gran<-gran[!is.na(gran$Granuloma),] #50 samples

gran$Granuloma[gran$Granuloma=="Y"]<-"Granuloma"
gran$Granuloma[gran$Granuloma=="N"]<-"Non-granuloma"


by.gran<- aggregate(gran["mean"],by=list(gran$Granuloma),mean,na.rm=T )

ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100), color="#00b4d899",size=2, position=position_jitter(width=0.03))+
    geom_point(data=by.gran, aes(x=Group.1, y=mean*100), color="#023e8a",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByGranuloma_diversity2.png", width = 3,height = 3,units="in", dpi=300)

ggplot()+
    geom_boxplot(data=gran, aes(x=Granuloma, y=mean*100), width=0.5,color="#00b4d899")+
    geom_point(data=by.gran, aes(x=Group.1, y=mean*100), color="#023e8a",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/ByGranuloma_diversity3.png", width = 3,height = 3,units="in", dpi=300)


wilcox.test(gran$mean[gran$Granuloma=="Granuloma"], gran$mean[gran$Granuloma=="Non-granuloma"], alternative="two.sided")
#W = 332, p-value = 0.6489

#1. Diversity of granuloma by cohort
by.granC<- aggregate(gran["mean"],by=list(gran$Granuloma, gran$Cohort),mean,na.rm=T )

gran$Cohort<-factor(gran$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.granC$Group.2<-factor(by.granC$Group.2, levels=c("SIV only", "Mtb NR", "Mtb R"))

ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100, color=Cohort),size=2, position=position_dodge(width=0.6))+
    geom_point(data=by.granC, aes(x=Group.1, y=mean*100, fill=Group.2), color="#023e8a",size=2, position=position_dodge(width=0.6))+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    scale_color_manual(values = paste0(colors[c(5,3,1)],99))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    annotate("segment", x = 0.85, xend = 0.85, y = 0.35, yend = 0.4, colour = "gray60", size=0.3) +
    annotate("segment", x = 1.15, xend = 1.15, y = 0.31, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.85, xend = 1.15, y = 0.4, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.41, colour = "gray60", label="*", size=5)
ggsave("Output/Figures/byGranuloma.byCohort_diversity.png", width = 4.5,height = 3,units="in", dpi=300)

wilcox.test(gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
#W = 10, p-value = 0.0257
wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
#W = 44, p-value = 0.761
wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="SIV only"], alternative="two.sided")
#W = 47, p-value = 0.8534
wilcox.test(gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="SIV only"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
#W = 45, p-value = 0.6965

#Test gran vs. non-gran within each cohort
wilcox.test(gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb R"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb R"], alternative="two.sided")
#W = 87, p-value = 0.9063
wilcox.test(gran$mean[gran$Granuloma=="Granuloma"&gran$Cohort=="Mtb NR"], gran$mean[gran$Granuloma=="Non-granuloma"&gran$Cohort=="Mtb NR"], alternative="two.sided")
#W = 24, p-value = 0.2141




#########
# 2. Diversity in Granuloma by tissue
by.gran2<- aggregate(gran["mean"],by=list(gran$Granuloma, gran$Tissue2),mean,na.rm=T )
by.gran2$Group.2<-factor(by.gran2$Group.2, levels=c("pLN","Thoracic LN","Lung"))
gran$Tissue2<-factor(gran$Tissue2, levels=c("pLN","Thoracic LN","Lung"))

ggplot()+
    geom_point(data=gran, aes(x=Granuloma, y=mean*100, color=Tissue2),size=2, 
               position=position_dodge(width=0.5))+
    geom_point(data=by.gran2, aes(x=Group.1, y=mean*100, color=Group.2, fill=Group.2),
                                  position=position_dodge(width=0.5), color="#023e8a",size=2)+
        xlab('')+ylab('% Diversity')+
    scale_color_manual(values = c("#00BF7D99", "#00B0F699", "#E76BF399"), labels=c("Peripheral LN", "Thoracic LN", "Lung"))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    guides(color = guide_legend(title = NULL))+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    #annotate("segment", x = 0.85, xend = 0.85, y = 0.35, yend = 0.4, colour = "gray60", size=0.3) +
    #annotate("segment", x = 1.15, xend = 1.15, y = 0.31, yend = 0.4, colour = "gray60", size=0.3)+
    #annotate("segment", x = 0.85, xend = 1.15, y = 0.4, yend = 0.4, colour = "gray60", size=0.3)+
    #annotate("text", x = 1, y = 0.41, colour = "gray60", label="*", size=5)

ggsave("Output/Figures/ByGranuloma.byTissue_diversity.png", width=4.5,height = 3,units="in", dpi=300)

G<-gran[gran$Granuloma=="Granuloma",]
N<-gran[gran$Granuloma=="Non-granuloma",]
wilcox.test(G$mean[G$Tissue2=="Lung"],G$mean[G$Tissue2=="Thoracic LN"], alternative="two.sided")
#W = 30, p-value = 0.0817

wilcox.test(N$mean[N$Tissue2=="Lung"],N$mean[N$Tissue2=="Thoracic LN"], alternative="two.sided")
#W = 26, p-value = 0.7214
wilcox.test(N$mean[N$Tissue2=="Lung"],N$mean[N$Tissue2=="pLN"], alternative="two.sided")
#W = 32, p-value = 0.6612
wilcox.test(N$mean[N$Tissue2=="pLN"],N$mean[N$Tissue2=="Thoracic LN"], alternative="two.sided")
#W = 48, p-value = 0.2875

# Non granuloma Latent R samples 
wilcox.test(N$mean[N$Tissue2=="Lung"&N$Cohort=="Mtb R"],N$mean[N$Tissue2=="Thoracic LN"&N$Cohort=="Mtb R"], alternative="two.sided")
#W = 5, p-value = 0.3333

# Test thoracic LN granuloma between R and NR
wilcox.test(G$mean[G$Tissue2=="Thoracic LN"&G$Cohort=="Mtb NR"],G$mean[G$Tissue2=="Thoracic LN"&G$Cohort=="Mtb R"], alternative="two.sided")
#W = 13, p-value = 0.2

wilcox.test(G$mean[G$Tissue2=="Thoracic LN"&G$Cohort=="Mtb NR"],G$mean[G$Tissue2=="Lung"&G$Cohort=="Mtb R"], alternative="two.sided")

gg<-G[G$Tissue2=="Thoracic LN",]

#####
#3. Diversity of granuloma by cohory and by tissue (SIV only has non-granuma only. Remove SIV only)
gran2<-gran[gran$Cohort!="SIV only",]

by.granA<- aggregate(gran2["mean"],by=list(gran2$Granuloma, gran2$Tissue2, gran2$Cohort),mean,na.rm=T )

by.granA$Group.2<-factor(by.granA$Group.2, levels=c("pLN","Thoracic LN","Lung"))
gran2$Tissue2<-factor(gran2$Tissue2, levels=c("pLN","Thoracic LN","Lung"))
gran2$Cohort<-factor(gran2$Cohort, levels=c("SIV only", "Mtb NR", "Mtb R"))
by.granA$Group.3<-factor(by.granA$Group.3, levels=c("SIV only", "Mtb NR", "Mtb R"))


ggplot()+
    geom_point(data=gran2, aes(x=Granuloma, y=mean*100, color=Cohort, shape=Tissue2, group=Cohort),size=2.5, 
               position=position_dodge(width=0.6))+
    geom_point(data=by.granA, aes(x=Group.1, y=mean*100, color=Group.3, fill=Group.3, shape=Group.2, group=Group.3),
               position=position_dodge(width=0.6), color="#023e8a",size=1.5)+
    xlab('')+ylab('% Diversity')+
    scale_color_manual(values = paste0(colors[c(3,1)],99))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    guides(color = guide_legend(title = NULL))+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())
    #annotate("segment", x = 0.85, xend = 0.85, y = 0.35, yend = 0.4, colour = "gray60", size=0.3) +
    #annotate("segment", x = 1.15, xend = 1.15, y = 0.31, yend = 0.4, colour = "gray60", size=0.3)+
    #annotate("segment", x = 0.85, xend = 1.15, y = 0.4, yend = 0.4, colour = "gray60", size=0.3)+
    #annotate("text", x = 1, y = 0.41, colour = "gray60", label="*", size=5)

ggsave("Output/Figures/ByGranuloma.byTissueByCohort_diversity1.png", width = 4.5,height = 3,units="in", dpi=300)

ggplot()+
    geom_point(data=gran2, aes(x=Granuloma, y=mean*100, color=Cohort, shape=Tissue2),size=2.5, 
               position=position_dodge(width=0.6))+
    geom_point(data=by.granA, aes(x=Group.1, y=mean*100, color=Group.3, fill=Group.3, shape=Group.2),
               position=position_dodge(width=0.6), color="#023e8a",size=1.5)+
    xlab('')+ylab('% Diversity')+
    scale_color_manual(values = paste0(colors[c(3,1)],99))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    scale_shape_manual(values=c(16,17,15),labels=c("Peripheral LN", "Thoracic LN", "Lung"))+
    guides(color = guide_legend(title = NULL))+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    annotate("segment", x = 0.8, xend = 0.85, y = 0.35, yend = 0.4, colour = "gray60", size=0.3) +
    annotate("segment", x = 1.2, xend = 1.15, y = 0.31, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("segment", x = 0.85, xend = 1.15, y = 0.4, yend = 0.4, colour = "gray60", size=0.3)+
    annotate("text", x = 1, y = 0.41, colour = "gray60", label="*", size=5)

ggsave("Output/Figures/ByGranuloma.byTissueByCohort_diversity2.png", width = 4.8,height = 3,units="in", dpi=300)


wilcox.test(N$mean[N$Tissue2=="Thoracic LN"&N$Cohort=="Mtb NR"],N$mean[N$Tissue2=="Lung"&N$Cohort=="Mtb NR"], alternative="two.sided")
#W = 2, p-value = 0.8
wilcox.test(N$mean[N$Tissue2=="Thoracic LN"&N$Cohort=="Mtb NR"],N$mean[N$Tissue2=="pLN"&N$Cohort=="Mtb NR"], alternative="two.sided")
#W = 3, p-value = 0.7
wilcox.test(N$mean[N$Tissue2=="Lung"&N$Cohort=="Mtb NR"],N$mean[N$Tissue2=="pLN"&N$Cohort=="Mtb NR"], alternative="two.sided")
#W = 0, p-value = 0.2

wilcox.test(N$mean[N$Tissue2=="Thoracic LN"&N$Cohort=="Mtb R"],N$mean[N$Tissue2=="Lung"&N$Cohort=="Mtb R"], alternative="two.sided")
#W = 0, p-value = 0.3333
wilcox.test(N$mean[N$Tissue2=="Thoracic LN"&N$Cohort=="Mtb R"],N$mean[N$Tissue2=="pLN"&N$Cohort=="Mtb R"], alternative="two.sided")
#W = 15, p-value = 0.2857
wilcox.test(N$mean[N$Tissue2=="Lung"&N$Cohort=="Mtb R"],N$mean[N$Tissue2=="pLN"&N$Cohort=="Mtb R"], alternative="two.sided")
#W = 4, p-value = 0.4


######
# Look at 1) Latent R cohort
granR<-gran[gran$Cohort=="Mtb R",]
by.granR<- aggregate(granR["mean"],by=list(granR$Granuloma),mean,na.rm=T )

ggplot()+
    geom_point(data=granR, aes(x=Granuloma, y=mean*100), color="darkturquoise",size=2)+
    geom_point(data=by.granR, aes(x=Group.1, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())
ggsave("Output/Figures/LatentR_Granuloma_div.png", width = 3,height = 3,units="in", dpi=300)

wilcox.test(granR$mean[granR$Granuloma=="Granuloma"],granR$mean[granR$Granuloma=="Non-granuloma"], alternative="two.sided")
#W = 87, p-value = 0.9063


granRm<-melt(granR[,c("Tissue2","Granuloma","mean")])
granRm$Tissue2<-factor(granRm$Tissue2, levels=c("pLN","Thoracic LN","Lung"))

ggplot()+
    geom_point(data=granRm, aes(x=Granuloma, y=value*100, color=Tissue2), position=position_dodge(width=0.5),size=2)+
    #geom_point(data=by.granR, aes(x=Group.1, y=mean*100), color="black",size=2)+
    xlab('')+ylab('% Diversity')+
    theme(legend.title = element_blank())+
    theme_bw()+
    ggtitle("Mtb R")+
    theme(panel.grid.major.x = element_blank())+    
    scale_color_manual(values = c("#00BF7DCC", "#00B0F6CC", "#E76BF3CC"), labels=c("Peripheral LN", "Thoracic LN", "Lung"))+
    scale_fill_manual(values=c(rep("#023e8a",times=3)), guide='none')+
    theme(legend.title = element_blank(), plot.title = element_text(size=12))
ggsave("Output/Figures/LatentR_Granuloma_byTissue_div.png", width = 3,height = 3,units="in", dpi=300)


### Latent NR cohort
granNR<-gran[gran$Cohort=="Mtb NR",]
ggplot()+
    geom_point(data=granNR, aes(x=Granuloma, y=mean*100), color="darkturquoise",size=2)+
    xlab('')+ylab('% Diversity')+
    ggtitle("Mtb NR")+
    theme(legend.title = element_blank())+
    theme_bw()+
    theme(panel.grid.major.x = element_blank(), plot.title = element_text(size=12))
ggsave("Output/Figures/LatentNR_Granuloma_div.png", width = 3,height = 3,units="in", dpi=300)

granRm<-melt(granR[,c("Tissue2","Granuloma","mean")])
granRm$Tissue2<-factor(granRm$Tissue2, levels=c("pLN","Thoracic LN","Lung"))

ggplot()+
    geom_point(data=granNR, aes(x=Granuloma, y=mean*100, color=Tissue2), size=2,
               position=position_dodge(width = .7))+
    xlab('')+ylab('% Diversity')+
    ggtitle("Mtb NR")+
    scale_color_manual(values = c("#00BF7DCC", "#00B0F6CC", "#E76BF3CC"), labels=c("Peripheral LN", "Thoracic LN", "Lung"))+
    theme_bw()+
    theme(panel.grid.major.x = element_blank())+
    theme(legend.title = element_blank(), plot.title = element_text(size=12))

ggsave("Output/Figures/LatentNR_Granuloma_byTissue_div.png", width = 3,height = 3,units="in", dpi=300)








##### Others
rna<-Sum21[,c("Cohort","Tissue2","Week","mean","SIV.RNA.per.granuloma","SIV.RNA.per.tissue",
              "CD4.percent","CD8.percent")]

rna1<-rna[!is.na(rna$SIV.RNA.per.granuloma),]
rna2<-rna[!is.na(rna$SIV.RNA.per.tissue),]

cor.test(rna1$SIV.RNA.per.granuloma, rna1$mean, method = "spearman") #P=0.2473
cor.test(rna2$SIV.RNA.per.tissue, rna2$mean, method = "spearman") #P= 0.3812

ggplot(rna1, aes(x=SIV.RNA.per.granuloma, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/Correlation_RNA_Div.pdf", height = 3.6, width = 4)
ggplot(rna2, aes(x=SIV.RNA.per.tissue  , y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/Correlation_RNAtissue_Div.pdf", height = 3.6, width = 4)

#### rna concetraiton by tissue by cohort
rna.conc<-data.frame(aggregate(rna2[,6], by=list(rna2$Cohort,rna2$Tissue2), mean, na.action=na.omit))
write.csv(rna.conc, "Output/Diversity/RNA_concentration_byTissues_byCohort.csv" )

plot(rna2$mean,rna2$SIV.RNA.per.tissue, pch=16)
cor.test(rna2$mean,rna2$SIV.RNA.per.tissue, method="spearman")
#S = 18196, p-value = 0.3812
#      rho 
#0.1262425


#####
cd4<-rna[!is.na(rna$CD4.percent),]
ggplot(cd4, aes(x=CD4.percent, y=mean))+
    geom_point()+
    theme_bw()+
    ggsave("Output/Figures/CD4.pdf", height = 3.6,width = 4)
cor.test(cd4$CD4.percent, cd4$mean, method = "spearman") #P=0.5474

cd8<-rna[!is.na(rna$CD8.percent),]
ggplot(cd8, aes(x=CD8.percent, y=mean*100))+
    geom_point()+
    theme_bw()+ylab("Diversity")
    ggsave("Output/Figures/CD8.pdf", height = 3.6,width = 4)
cor.test(cd8$CD8.percent, cd8$mean, method = "spearman") # p=0.03326, rhp=-0.3773173 



## Plot all diversity (Figure 4C)
all_mean<-Sum21[,c("Cohort","Monkey","Tissue2","Granuloma","Week","mean")]
all_mean<-InsertRow(all_mean,c(0,0,0,0,0, stks[1,"mean"]),1)
all_mean[1,1:3]<-c("Stock","Stock","Stock")
all_mean$label<-paste0(all_mean$Monkey,'_', all_mean$Tissue2,'_', all_mean$Week)
all_mean$Tissue2[all_mean$Tissue2=="plasma"]<-"Plasma"
all_mean$Tissue2[all_mean$Tissue2=="pLN"]<-"Periferal LN"
all_mean$Tissue2<-factor(all_mean$Tissue2, levels=c("Stock", "Plasma", "Periferal LN","Thoracic LN","Lung"))

#order by 'stock,'siv only, NR, N
all_mean$Cohort<-factor(all_mean$Cohort, levels=c("Stock", "SIV only", "Mtb NR", "Mtb R"))
all_mean<-all_mean[order(all_mean$Cohort, all_mean$Monkey,all_mean$Week, all_mean$Tissue2),]

all_mean$label<-factor(all_mean$label, levels=unique(all_mean$label))


dat<-distinct(all_mean, label, .keep_all=T) #50 rows

numb<-data.frame(table(dat$Monkey))
colnames(numb)<-c("Monkey","Freq")
morder<-unique(all_mean$Monkey)
write.csv(morder,"Output/monkey_order.csv")
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
pos2<-(numb$col[8]-numb$col[5])/2+numb$col[5]+0.5
pos3<-(50-numb$col[8])/2+numb$col[8]+0.5

xlabels<-dat$Week
xlabels[1]<-""

#samples with granuloma
all_mean$mean2<-all_mean$mean
all_mean$mean2[all_mean$Granuloma!="Y"]<-NA
all_mean$mean2[is.na(all_mean$Granuloma)]<-NA

all_mean<-transform(all_mean,id=as.numeric(factor(label)))

ggplot()+
    geom_point(data=all_mean, aes(x=label, y=mean*100, color=Tissue2),size=2.5,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Diversity')+
    scale_y_continuous(limits = c(0.1,0.5))+
    scale_x_discrete(labels=xlabels)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.1, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=0.5, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=0.5, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=0.5, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    #annotate("text",x=1, y=0.1, label="stock" , angle=90, hjust=1, size=3, color="gray30")+
    geom_text(aes(x=1, y=-Inf, label="stock") , angle=90, hjust=1.1, size=3, color="gray30")+
    coord_cartesian(clip = "off")+
    annotate("point", x=all_mean$id, y=all_mean$mean2*100, shape=21, color="gray10", size=2.7)
ggsave("Output/Figures/diversity_allSamples3.png", width = 8,height = 3,units="in", dpi=300)

ggplot()+
    geom_point(data=all_mean, aes(x=label, y=mean*100, color=Tissue2),size=2.5,
               position = position_dodge(width = .5))+
    xlab('')+ylab('% Diversity')+
    scale_y_continuous(limits = c(0.1,0.5))+
    scale_x_discrete(labels=xlabels)+
    theme_bw()+
    theme(legend.title = element_blank())+
    theme(panel.grid.major.x = element_blank())+
    geom_vline(xintercept = numb$col[1:11]+0.5,  
               color = "gray80", size=.4)+
    annotate("text", x=numb$lab_pos[2:12]+0.5, y=0.1, label= numb$Monkey[2:12], hjust=0.5 , size=3)+
    annotate("rect", xmin=1.5, xmax=numb$col[5]+0.5,ymax=Inf,ymin=-Inf, fill=colors[5], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[5]+0.5, xmax=numb$col[8]+0.5,ymax=Inf,ymin=-Inf, fill=colors[3], alpha=0.1 , color=NA)+
    annotate("rect", xmin=numb$col[8]+0.5, xmax=50.5,ymax=Inf,ymin=-Inf, fill=colors[1], alpha=0.1 , color=NA)+
    annotate("text", x=pos1, y=0.5, label="SIV only", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos2, y=0.5, label="Mtb NR", hjust=0.5, size=3.5, color='gray20')+
    annotate("text", x=pos3, y=0.5, label="Mtb R", hjust=0.5, size=3.5, color='gray20')+
    #annotate("text",x=1, y=0.1, label="stock" , angle=90, hjust=1, size=3, color="gray30")+
    geom_text(aes(x=1, y=-Inf, label="stock") , angle=90, hjust=1.1, size=3, color="gray30")+
    coord_cartesian(clip = "off")+
    annotate("point", x=all_mean$id, y=all_mean$mean2*100, shape=4, color="white", size=2.7)
ggsave("Output/Figures/diversity_allSamples2.png", width = 8,height = 3,units="in", dpi=300)
