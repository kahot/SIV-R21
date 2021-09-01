library(glmmTMB)
library(reshape2)
library(ggplot2)
library("bbmle") 
library(DHARMa)
library(car)
library(effects)
library(MuMIn)
library(glmmLasso)
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

#Change the mean to standardized counts.
siv<-Sum
siv$Count<-round(siv$mean*1000000)

mean(siv$Count)
var(siv$Count)

#variance is greater than mean ==overdispersion
#need to sue negative bionomial distribution

po<-glmmTMB(Count~ Week+ Tissue+Cohort+(1|Monkey), data=siv, family=poisson)
summary(po)

nb1<-glmmTMB(Count~ Week+ Tissue+Cohort+(1|Monkey), data=siv, family=nbinom1)
summary(nb1)

nb2<-glmmTMB(Count~ Week+ Tissue+Cohort+(1|Monkey), data=siv, family=nbinom2)
summary(nb2)

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


#Nested random effects
#not converged
betaN<-glmmTMB(mean~ Week+ Tissue+Cohort+(1|Cohort/Tissue/Monkey), data=Sum, family=beta_family())
summary(betaN)


#look at the residuals
beta1<- simulateResiduals(beta)
plot(beta1)
Anova(beta, type=3)
#Response: mean
#Chisq Df Pr(>Chisq)  
#Week    2.3377  1    0.12628  
#Tissue  8.4053  4    0.07781 .
#Cohort  5.4980  2    0.06399 .


re<-nb2
sim<- simulateResiduals(re)
plot(sim)
Anova(re)

#nb1 and beta have very similar results
#Response: Count
#         Chisq Df Pr(>Chisq)  
#Week    2.3348  1    0.12651  
#Tissue  8.4275  4    0.07712 .
#Cohort  5.5031  2    0.06383 

#nb2
#Response: Count
#         Chisq Df Pr(>Chisq)  
#Week    2.5103  1    0.11310  
#Tissue  9.9616  4    0.04108 *
#Cohort  5.6449  2    0.05946 .


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

beta4<-glmmTMB(mean~ Tissue+Cohort, data=Sum, family=beta_family())
summary(beta4)

#5. with granulomas
gran<-glmmTMB(mean~ Week+ Tissue+Cohort+Granuloma +(1|Monkey), data=siv, family=beta_family())
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

##  using dredge to do model selection ## --somehow different from above
#Generate a model selection table of models with combinations (subsets) of fixed 
#effect terms in the global model, with optional rules for model inclusion.

#nb2.dredge<-MuMIn::dredge(nb2)
#print(nb2.dredge)

beta.dredge<-MuMIn::dredge(beta)
print(beta.dredge)

#can visualized the model selection
op <- par(mar=c(3,4,7,4))
plot(beta.dredge)

# another way to print out the details of models tested
get.models(beta.dredge, subset=T)


nb2.dredge<-MuMIn::dredge(nb2)
print(nb2.dredge)#
nb22<-glmmTMB(Count~ Week+ Tissue+Cohort, data=siv, family=nbinom2)
summary(nb22)
nb22.dredge<-MuMIn::dredge(nb22)
print(nb22.dredge)#




#### Effects ###
#Visualize the effects
ae<-allEffects(beta)
plot(ae)


#plot the estimated and actual
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
        ggsave(plot=Plots[[i]],filename=paste0("Output/GLMM/Estimated.effects.",f,".pdf"), width = 5,height = 3)
        
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
        ggsave(plot=Plots[[i]],filename=paste0("Output/GLMM/Estimated.effects.",f,".pdf"), width = 5,height = 3)
        
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


png("Output/GLMM/All_estimated_effects.png", width=, height =4.5, units="in",res=300)
do.call(grid.arrange, c(Plots, ncol=1))
dev.off()

#draw_plot(plot, x, y, width, height)
pdf("Output/Figures/GLMM_estimates.pdf",width=9, height=6)
#png("Output/Figures/Figure5.png",width=9, height=6, units="in",res=300)
ggdraw()+
    draw_plot(Plots[[1]],0,0.5,0.6,0.5)+
    draw_plot(Plots[[3]],0.6,0.5,0.4,0.5)+
    draw_plot(Plots[[2]],0,0,0.65,0.5)+
    draw_plot(Plots[[4]],0.65,0,0.2,0.5)+
    draw_plot_label(c("A", "B", "C"), c(0, 0.6, 0), c(1, 1, 0.5), size = 15)
dev.off()    


library(colorspace)
source("Rscripts/label_scientific.R")
colors<-qualitative_hcl(6, palette="Dark3")

## Add all poitns to all plots?
#1.Cohort
Est<-Estimates[[3]]

Est$Factor<-factor(Est$Factor, levels=paste(ave$Factor))
Estm<-melt(Est, id.vars="Factor")
#Plots[[i]]<-
    ggplot(Est)+
    geom_point(data=Sum, aes(x=Cohort, y=mean*100, color=Cohort),size=2.3, position=position_jitter(width=0.05))+
    scale_color_manual(values=paste0(colors[c(5,3,1)],"66"), guide='none')+
    
    geom_point(aes(x=Factor, y=Est*100), shape=21, color="blue", size=3, stroke =1.2)+
    geom_errorbar(aes(x=Factor, y=Est*100,ymin=SE_low*100, ymax=SE_high*100), width=.2, size=.3, color="royalblue")+
    #geom_point(aes(x=Factor, y=Obs*100), position=position_dodge(width=0.2),shape=23, size=3,  stroke = 1.2, color="maroon")+
    theme_bw()+
    #scale_y_continuous(limit=c(0.18,0.32))+
    xlab('')+ylab('% Average diversity')+
    ggtitle(paste(f))+
    theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+
    annotate("")



p<-ggplot+
    
p3+geom_point(data=Sum, aes(x=Cohort, y=mean*100, color=Cohort),size=2.3, position=position_jitter(width=0.05))+
    scale_color_manual(values=paste0(colors[c(5,3,1)],"66"), guide='none')

    
    
    
    
    
    
    
overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(beta)

### Multicomparison  ###
levels(Sum$Tissue )
emmeans(beta, ~ Cohort, at = list(F2 = levels(Sum$Tissue)), infer=T)

emmeans(beta,"mean", level = 0.90, infer = TRUE )
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





write.csv(df,"Example.csv")


coeff.matrix<-function(x, name){
    x2<-summary(x)
    coef<-x2$coefficients
    df<-coef[[1]]
    df<-data.frame(df)
    df$Factor<-rownames(df)
    df$Effect<-''
    for (g in 1:length(row.names(df)) ){
        if (g==1){
            df$Effect[1]<- exp(df[1,g])
        }
        else{
            df$Effect[g]<- (((exp(df[1,1] + df$Estimate[g]) - exp(df[1,1])) /exp(df[1,1])))#add estimate % column
        }
    }
    write.csv(df, paste0("Output/GLMM/",name,".csv"))
    return(data.frame(df))
    
}


beta.model<-coeff.matrix(beta, "glmm.beta")

beta.model$Effect<-as.numeric(beta.model$Effect)
beta.model$Factor<-factor(beta.model$Factor, levels=paste(beta.model$Factor))

#Horizontal
ggplot(beta.model, aes(Factor,Effect*100)) +
    geom_bar(stat="identity", color="#2C67AB", fill=paste0("#2C67AB","CC"))+
    theme_test() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")







ae1<-allEffects(glmm1)
plot(ae1)


#the emmeans package computes estimated marginal means (previously known as least-squares means) for the fixed effects of any component
library(emmeans)
library(glmmLasso)
emmeans(glmm1, ~Week |Tissue2)

#Generate a model selection table of models with combinations (subsets) of fixed 
#effect terms in the global model, with optional rules for model inclusion.

glmm1.dredge<-MuMIn::dredge(glmm1)
print(glmm1.dredge)

glmm1.dredge<-MuMIn::dredge(glmm1, fixed=Cohort)
print(glmm1.dredge)


op <- par(mar=c(3,4,7,4))
plot(glmm1.dredge)

get.models(glmm1.dredge, subset=T)

lasso1<-glmmLasso(mean~Week+as.factor(Cohort)+as.factor(Tissue2), rnd=list(Monkey=~1), data=Sum, 
          lambda=100, family=beta_family(),control = list(print.iter=TRUE,start=c(1,rep(0,29)),q.start=0.7))


## Tutorial
dater<-read.csv("~/Downloads/Example data.csv")
head(dater)
options(na.action = "na.fail")

mod1<-lm(density~distance+elev, data = dater) 
mod2<-lm(density~slope+pct.cover, data = dater) 
mod3<-lm(density~slope+distance, data = dater) 
mod4<-lm(density~slope+distance+elev, data = dater) 
out.put<-model.sel(mod1,mod2,mod3,mod4) 
out.put

#create a confidence set of models using the subset function 
# select models with delta AICc less than 5 
# IMPORTANT: Weights have been renormalized!! 
subset(out.put, delta <5) 

# select models using Royall's 1/8 rule for strength of evidence 
# IMPORTANT: Weights have been renormalized!! 
subset(out.put, 1/8 < weight/max(out.put$weight)) 

# select models 95% cumulative weight criteria 
# IMPORTANT: Weights have been renormalized!! 
subset(out.put, cumsum(out.put$weight) <= .95) 

