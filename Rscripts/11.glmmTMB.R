library(glmmTMB)
library(reshape2)
library(ggplot2)
library(bbmle) 
#library(DHARMa)
library(car)
library(effects)
#library(MuMIn)
#library(glmmLasso)
library(cowplot)
library(gridExtra)
library(multcomp)
library(emmeans)

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]
summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum<-summary[1:69,]
Sum<-merge(Sum[,c(1,2,3,5,6,10)], samples[,c("filename","Granuloma","SIV.RNA.per.tissue")], by="filename")

colnames(Sum)[colnames(Sum)=="Tissue2"]<-"Tissue"
#create ID
Sum$ID<-paste0(Sum$Monkey,"_", Sum$Week,"_",Sum$Tissue)

#replace NA to N for Granuloma
Sum$Granuloma[is.na(Sum$Granuloma)]<-"N"
Sum$Granuloma<- ifelse(Sum$Granuloma=="N", 0 ,1)

Sum$Cohort<-factor(Sum$Cohort, levels = c("SIV only", "Mtb NR", "Mtb R"))
Sum$Tissue<-factor(Sum$Tissue, levels = c("Plasma", "Plasma_nex", "pLN", "tLN","Lung"))


#glmmTMB

beta<-glmmTMB(mean~ Week+ Tissue+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta)
#Conditional model:
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      -6.026e+00  1.192e-01  -50.54   <2e-16 ***
#Week             -4.541e-02  2.970e-02   -1.53   0.1263    
#TissuePlasma_nex  2.340e-01  1.397e-01    1.68   0.0938 .  
#TissuepLN         1.950e-01  1.444e-01    1.35   0.1770    
#TissuetLN         2.911e-01  1.415e-01    2.06   0.0396 *  
#TissueLung        2.319e-01  1.440e-01    1.61   0.1072    
#CohortMtb NR      1.963e-01  9.661e-02    2.03   0.0422 *  
#CohortMtb R       9.237e-05  8.770e-02    0.00   0.9992  


#look at the residuals
beta1<- simulateResiduals(beta)
plot(beta1)
Anova(beta, type=3)
#Response: mean
#Chisq Df Pr(>Chisq)  
#Week    2.3377  1    0.12628  
#Tissue  8.4053  4    0.07781 .
#Cohort  5.4980  2    0.06399 .



#AIC comparison across models
#1. full model
beta<-glmmTMB(mean~ Week+ Tissue+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta)

#2. without random factor
beta2<-glmmTMB(mean~ Week+ Tissue+Cohort, data=Sum, family=beta_family())
summary(beta2)

#3. without week
beta3<-glmmTMB(mean~ Tissue+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta3)

#without week and the random factor
beta4<-glmmTMB(mean~ Tissue+Cohort, data=Sum, family=beta_family())
summary(beta4)

#5. with granulomas
gran<-glmmTMB(mean~ Week+ Tissue+Cohort+Granuloma +(1|Monkey), data=Sum, family=beta_family())
summary(gran)

#with interaction
betaI<-glmmTMB(mean~ Week+Tissue*Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(betaI)
#-887.3
#Don't worry about the interaction. not significant.


AICtab(beta, beta2,beta3,beta4, gran, betaI)
#      dAIC df
#beta   0.0 10
#beta3  0.3 9 
#gran   2.0 11
#betaI 10.7 18
#beta2 21.3 9 
#beta4 24.9 8 

#The full model -granuloma is the best


#dispersion parameter
sigma(beta)
#30572.82



#plot the estimated and actual
#EstimatedEffects
ae<-allEffects(beta)
fac<-c("Week","Tissue","Cohort")
Estimates<-list()
Plots<-list()
for (i in 1:length(fac)){
    f<-fac[i]
    ave<-data.frame(aggregate(Sum["mean"], by=list(Sum[,f]), mean))
    colnames(ave)<-c("Factor","Obs")
    est<-ae[[f]]
    Est<-data.frame(x=est$x, Est=exp(est$fit), SE_low=exp(est$lower), SE_high=exp(est$upper))
    colnames(Est)[1]<-"Factor"
    Est<-merge(Est, ave, by= "Factor", all=T)
    
    Estimates[[i]]<-Est

    if (f!="Week"){
        Est$Factor<-factor(Est$Factor, levels=paste(ave$Factor))
        
        Plots[[i]]<-ggplot(Est)+
            geom_point(aes(x=Factor, y=Est*100), shape=21, color="blue", size=3)+
            geom_errorbar(aes(x=Factor, y=Est*100,ymin=SE_low*100, ymax=SE_high*100), width=.2, size=.3, color="royalblue")+
            geom_point(aes(x=Factor, y=Obs*100), shape=23, size=3, color="maroon")+
            theme_bw()+
            #scale_y_continuous(limit=c(0.18,0.32))+
            xlab('')+ylab('% Average diversity')+
            ggtitle(paste(f))+
            theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
    }
    if (f=="Week"){
        Est1<-Est[!is.na(Est$Est),]
        Plots[[i]]<-ggplot()+
            geom_line(data=Est1, aes(x=Factor, y=Est*100), color="blue", size=1)+
            geom_point(data=Est,aes(x=Factor, y=Obs*100), shape=23, size=3, color="maroon")+
            geom_ribbon(data=Est1,aes(x=Factor, y=Est*100,ymin=SE_low*100,ymax=SE_high*100),alpha=0.2, fill="royalblue")+
            theme_bw()+
            #scale_y_continuous(limit=c(0.18,0.32))+
            xlab('')+ylab('% Average diversity')+
            ggtitle(paste(f))+
            theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
    }   
}

#add legend
get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}
df<-Estimates[[2]]
df<-df[,c("Factor","Est","Obs")]
colnames(df)[2:3]<-c("Estimated", "Observed")
dfm<-melt(df)
l<-ggplot()+
    geom_point(data=dfm,aes(x=Factor, y=value*100, fill=variable, color=variable, shape=variable), size=3, fill='white')+
    scale_color_manual(values = c("blue","maroon"))+
    scale_shape_manual(values=c(21,23))+
    theme_bw()+theme(legend.title = element_blank())
  
legend1<-get_legend(l)
Plots[[4]]<-legend1

#draw_plot(plot, x, y, width, height)
pdf("Output/Figures/GLMM_estimates.pdf",width=9, height=6)
ggdraw()+
    draw_plot(Plots[[1]],0,0.5,0.6,0.5)+
    draw_plot(Plots[[3]],0.6,0.5,0.4,0.5)+
    draw_plot(Plots[[2]],0,0,0.65,0.5)+
    draw_plot(Plots[[4]],0.65,0,0.2,0.5)+
    draw_plot_label(c("A", "B", "C"), c(0, 0.6, 0), c(1, 1, 0.5), size = 15)
dev.off()    



### Multicomparison  ###
emmeans(beta, ~ Cohort, at = list(F2 = levels(Sum$Tissue)), infer=T)

coh.emm<-emmeans(beta, "Cohort",infer = TRUE)
pairs(coh.emm)
# contrast           estimate     SE df t.ratio p.value
#SIV only - Mtb NR -1.96e-01 0.0966 59  -2.032  0.1135
#SIV only - Mtb R  -9.24e-05 0.0877 59  -0.001  1.0000
#Mtb NR - Mtb R     1.96e-01 0.0925 59   2.121  0.0943


tis.emm<-emmeans(beta, "Tissue")
pairs(tis.emm)
# contrast            estimate     SE df t.ratio p.value
#Plasma - Plasma_nex -0.23404 0.1397 59  -1.676  0.4565
#Plasma - pLN        -0.19500 0.1444 59  -1.350  0.6614
#Plasma - tLN        -0.29108 0.1415 59  -2.058  0.2522
#Plasma - Lung       -0.23192 0.1440 59  -1.611  0.4965
#Plasma_nex - pLN     0.03905 0.0517 59   0.756  0.9421
#Plasma_nex - tLN    -0.05704 0.0450 59  -1.266  0.7125
#Plasma_nex - Lung    0.00212 0.0484 59   0.044  1.0000
#pLN - tLN           -0.09608 0.0454 59  -2.117  0.2265
#pLN - Lung          -0.03693 0.0486 59  -0.760  0.9408
#tLN - Lung           0.05916 0.0404 59   1.466  0.5882


#Look at the effect sizes
source("Rscripts/coeff.matrix.R")
beta.model<-coeff.matrix(beta)
#save the effect size
write.csv(beta.model, "Output/GLMM/glmm.beta.csv")


#Plot the effect sizes
beta.model$Effect<-as.numeric(beta.model$Effect)
beta.model$Factor<-factor(beta.model$Factor, levels=paste(beta.model$Factor))

xlabels<-as.character(beta.model$Factor)
tissues<-c("Plasma nex","Peripheral LN", "Thoracic LN", "Lung")
xlabels[3:6]<-paste("Tissue: ", tissues)
xlabels[7:8]<-c("Cohort: Mtb NR", "Cohort: Mtb R")

beta.model<-beta.model[-1,]
#Horizontal
ggplot(beta.model, aes(Factor,Effect*100)) +
    geom_bar(stat="identity", color="#2C67AB", fill=paste0("#2C67AB","CC"))+
    theme_test() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")
ggsave("Output/Figures/GLMM_effectsizes.pdf", height = 5, width = 7)
