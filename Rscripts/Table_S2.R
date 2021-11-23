ECreate Table S2

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
stock<-SampleSheet[SampleSheet$Tissue=="stock_virus",]
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]

summ<-samples[,c("Cohort","Monkey","Tissue2","Granuloma","SIV.RNA.per.tissue","CD4.count","CD4.percent", "Week")]

summ$Cohort<-factor (summ$Cohort, c("SIV only", "Mtb NR", "Mtb R"))
summ$Tissue2<-factor(summ$Tissue2, levels = c("Plasma","Plasma_nex","pLN","tLN","Lung"))

summ<-summ[order(summ$Cohort, summ$Monkey,summ$Tissue2),]
summ$Granuloma[is.na(summ$Granuloma)]<-""


vir<-read.csv("Data/viremia.csv")

pla<-summ[summ$Tissue2=="Plasma"|summ$Tissue2=="Plasma_nex",]
pla$Monkey<-as.integer(pla$Monkey)

pla<-merge(pla[,c(1:4,7:8)], vir, by=c("Monkey","Week"), all=T)

tis<-summ[!(summ$Tissue2=="Plasma"|summ$Tissue2=="Plasma_nex"),]
colnames(tis)[5:6]<-c("vRNA","CD4.counts")
sum<-rbind(pla, tis)
sum$Monkey<-as.integer(sum$Monkey)
sum<-sum[order(sum$Cohort, sum$Monkey,sum$Tissue2),]

write.csv(sum, "Output/Table_S2.csv")
