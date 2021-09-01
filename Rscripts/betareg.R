library(glmmTMB)
library(reshape2)
library(mltools)
library(data.table)
library(betareg)
library(ggplot2)

SampleSheet<-read.csv("Data/SampleSheet_Mac251.csv", stringsAsFactors =F)
samples<-SampleSheet[SampleSheet$Monkey!="stock_virus",]
summary<-read.csv("Output/Diversity_summary_R21.csv", stringsAsFactors = F, row.names = 1)
Sum<-summary[1:69,]
Sum<-merge(Sum[,c(1,2,3,5,6,10)], samples[,c("filename","Granuloma","SIV.RNA.per.tissue")], by="filename")
#create ID
Sum$ID<-paste0(Sum$Monkey,"_", Sum$Week,"_",Sum$Tissue2)

#replace NA to N for Granuloma
Sum$Granuloma[is.na(Sum$Granuloma)]<-"N"
Sum$Granuloma<- ifelse(Sum$Granuloma=="N", 0 ,1)

Sum$Cohort<-factor(Sum$Cohort, levels = c("SIV only", "Mtb NR", "Mtb R"))
Sum$Tissue2<-factor(Sum$Tissue2, levels = c("Plasma", "Plasma_nex", "pLN", "tLN","Lung"))

#Categorical variables for one hot encoding
cat<-Sum[,c("Monkey","Cohort","Tissue2")]
colnames(cat)[colnames(cat)=="Tissue2"]<-"Tissue"
cat<-data.frame(lapply(cat,factor))
new<-one_hot(as.data.table(cat))

#Combine all variables for beta regression
siv<-cbind(Sum[,c("ID","mean","Week","Granuloma")],new)
write.csv(siv,"Output/betareg.shpae.data.csv",row.names = F)


#categorical as serial numbers
cat2<-data.frame(lapply(cat,as.numeric))
siv2<-cbind(Sum[,c("filename","ID","mean","Week","Granuloma")],cat2)
write.csv(siv2,"Output/betareg.shpae.data2.csv",row.names = F)

#############

#### Read data #####
siv<-read.csv("Output/betareg.shpae.data2.csv", stringsAsFactors = F)

range(siv$mean)
# the number of sites used to calculate the frequencies
siteNo<-read.csv("Output/Diversity/NumberofUsed.sites.persample.csv",stringsAsFactors = F, row.names = 1) 
siv<-merge(siv, siteNo[,c("filename","n")], by="filename")

#transformation of zero inflated y values (recommended by betareg) 
#siv$y.trans<-unlist(apply(siv["mean"],1, function(x) y.trans.betareg(x,siv[,"n"]))
siv$y.trans<-(siv$mean*((siv$n)-1)+0.5)/siv$n


# y.trans.betareg <- function(y, n){
#    (y * (n - 1) + 0.5) /n
#}


siv<-read.csv("Output/betareg.shpae.data2.csv", stringsAsFactors = F)

range(siv$mean)
# the number of sites used to calculate the frequencies
siteNo<-read.csv("Output/Diversity/NumberofUsed.sites.persample.csv",stringsAsFactors = F, row.names = 1) 
siv<-merge(siv, siteNo[,c("filename","n")], by="filename")

#transformation of zero inflated y values (recommended by betareg) 
#siv$y.trans<-unlist(apply(siv["mean"],1, function(x) y.trans.betareg(x,siv[,"n"]))
siv$y.trans<-(siv$mean*((siv$n)-1)+0.5)/siv$n



vars<-colnames(siv)[3:ncol(siv)]
paste(vars,collapse='+ ') 
"Week+ Granuloma+ Monkey_16314+ Monkey_20615+ Monkey_30816+ Monkey_3116+ Monkey_31316+ Monkey_3216+ Monkey_3316+ Monkey_3516+ Monkey_3616+ Monkey_3816+ Monkey_4016+ Cohort_SIV only+ Cohort_Mtb NR+ Cohort_Mtb R+ Tissue_Plasma+ Tissue_Plasma_nex+ Tissue_pLN+ Tissue_tLN+ Tissue_Lung"



mod0<-betareg(y.trans~ Week+ Granuloma+ Monkey+ Cohort+ Tissue, data=siv)
summary(mod0)

#Coefficients (mean model with logit link):
#Estimate Std. Error  z value Pr(>|z|)    
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -5.6364109  0.0864869 -65.171   <2e-16 ***
#Week        -0.0102575  0.0105469  -0.973   0.3308    
#Granuloma    0.0298166  0.0381020   0.783   0.4339    
#Monkey       0.0002433  0.0063339   0.038   0.9694    
#Cohort      -0.0649059  0.0267105  -2.430   0.0151 *  
#Tissue       0.0221198  0.0106922   2.069   0.0386 * 
AIC(mod0)
#-868.5067

mod1<-betareg(mean~ Week+ Granuloma+ Monkey+ Cohort+ Tissue, data=siv)
summary(mod1)
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -5.912e+00  1.158e-01 -51.072   <2e-16 ***
#Week        -1.412e-02  1.414e-02  -0.999   0.3180    
#Granuloma    3.897e-02  5.096e-02   0.765   0.4444    
#Monkey      -1.803e-05  8.509e-03  -0.002   0.9983    
#Cohort      -8.547e-02  3.575e-02  -2.391   0.0168 *  
#Tissue       2.984e-02  1.435e-02   2.080   0.0376 * 


mod2 <- update(mod1, ~. -Monkey)
summary(mod2)
#             Estimate Std. Error  z value Pr(>|z|)    
#(Intercept) -5.91189    0.11427 -51.737  < 2e-16 ***
#Week        -0.01413    0.01336  -1.057  0.29037    
#Granuloma    0.03901    0.04851   0.804  0.42137    
#Cohort      -0.08551    0.03059  -2.795  0.00518 ** 
#Tissue       0.02984    0.01433   2.083  0.03725 *  

mod3 <- update(mod2, ~. -Granuloma)
summary(mod3)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -5.90711    0.11429 -51.686  < 2e-16 ***
#Week        -0.01084    0.01281  -0.846  0.39735    
#Cohort      -0.08811    0.03038  -2.900  0.00373 ** 
#Tissue       0.02626    0.01370   1.916  0.05533 .

mod4 <- update(mod3, ~. -Week)
summary(mod4)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -5.98052    0.07650 -78.173  < 2e-16 ***
#Cohort      -0.08646    0.03050  -2.835  0.00458 ** 
#Tissue       0.02348    0.01328   1.768  0.07707 .  

 
AIC(mod0,mod1,mod2,mod3,mod4)
#     df       AIC
#mod0  7 -868.5067
#mod1  7 -871.6073
#mod2  6 -873.6073
#mod3  5 -874.9732
#mod4  4 -876.2328

mod02 <- update(mod0, ~. -Monkey)
mod03 <- update(mod02, ~. -Granuloma)
mod04 <- update(mod03, ~. -Week)
AIC(mod0,mod02,mod03,mod04)


model1<-summary(mod1)
model3<-summary(mod3)
model4<-summary(mod4)


# full model
"Week+ Granuloma+ Monkey_16314+ Monkey_20615+ Monkey_30816+ Monkey_3116+ Monkey_31316+ Monkey_3216+ Monkey_3316+ Monkey_3516+ Monkey_3616+ Monkey_3816+ Monkey_4016+ Cohort_SIV only+ Cohort_Mtb NR+ Cohort_Mtb R+ Tissue_Plasma+ Tissue_Plasma_nex+ Tissue_pLN+ Tissue_tLN+ Tissue_Lung"
siv<-read.csv("Output/betareg.shpae.data.csv")
vars<-colnames(siv)[3:ncol(siv)]
paste(vars,collapse='+ ') 


md<-betareg(mean~ Cohort_SIV.only+ Cohort_Mtb.NR+ Cohort_Mtb.R, data=siv)



md<-betareg(mean~ Week+Tissue_tLN, data=siv)
summary(md)
AIC(md)

#(Intercept) -6.035854   0.097185 -62.107   <2e-16 ***
#Week        -0.009947   0.013384  -0.743   0.4573    
#Tissue_tLN   0.110491   0.048840   2.262   0.0237 *  
#-870.8768

md1<-betareg(mean~ Week+ Tissue_tLN+Cohort_Mtb.NR+Monkey_20615+Monkey_30816, data=siv)
summary(md1)
AIC(md1)

md1<-betareg(mean~ Week+ Tissue_tLN+Cohort_Mtb.NR+Tissue_Lung*Cohort_Mtb.NR, data=siv)
summary(md1)
AIC(md1)
#(Intercept)   -6.06068    0.08964 -67.608  < 2e-16 ***
#Week          -0.01188    0.01233  -0.963 0.335505    
#Tissue_tLN     0.09501    0.04556   2.085 0.037026 *  
#Cohort_Mtb.NR  0.15693    0.04496   3.490 0.000483 ***
#AIC= -879.6759

#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)     -6.06245    0.09100 -66.624  < 2e-16 ***
#Week            -0.01187    0.01233  -0.962 0.335871    
#Tissue_tLN       0.09500    0.04556   2.085 0.037055 *  
#Cohort_Mtb.NR    0.15863    0.04762   3.331 0.000865 ***
#Cohort_SIV.only  0.00578    0.05283   0.109 0.912880   
# -877.6879

#(Intercept)   -6.05667    0.09706 -62.400  < 2e-16 ***
#Week          -0.01187    0.01233  -0.962  0.33587    
#Tissue_tLN     0.09500    0.04556   2.085  0.03705 *  
#Cohort_Mtb.NR  0.15285    0.05831   2.622  0.00875 ** 
#Cohort_Mtb.R  -0.00578    0.05283  -0.109  0.91288    
#-877.6879

#(Intercept)   -6.070169   0.071253 -85.192  < 2e-16 ***
#Week          -0.010528   0.009808  -1.073   0.2831    
#Tissue_tLN     0.072772   0.035944   2.025   0.0429 *  
#Cohort_Mtb.NR  0.024662   0.043115   0.572   0.5673    
#Monkey_30816   0.403643   0.059131   6.826 8.71e-12 ***
#-911.5029


#(Intercept)   -6.201430   0.075978 -81.621  < 2e-16 ***
#Week           0.003102   0.009904   0.313 0.754111    
#Tissue_tLN     0.064897   0.032895   1.973 0.048510 *  
#Cohort_Mtb.NR  0.054717   0.040471   1.352 0.176381    
#Monkey_20615   0.171841   0.045443   3.781 0.000156 ***
#Monkey_30816   0.405897   0.054108   7.502 6.31e-14 ***
#-921.8835




#write coefficients

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
    write.csv(df, paste0("Output/Betareg/",name,".csv"))
    return(data.frame(df))
   
    }

model<-coeff.matrix(mod1, "Model1")


### Calculate the effects:

# Plot the effect size:
model<-read.csv("Output/Betareg/Model1.csv", row.names = 1, stringsAsFactors = F)

model$Factor<-factor(model$Factor, levels=paste(model$Factor))

#Horizontal
ggplot(model, aes(Factor,Effect*100)) +
    geom_bar(stat="identity", color="#2C67AB", fill=paste0("#2C67AB","CC"))+
    theme_test() +
    theme(axis.text=element_text(size=13), axis.title.y=element_text(size=13))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))+
    theme(panel.grid.major.y = element_line(color="gray80",linetype=5))+
    labs(x="", y="Estimated effects (%)")





#glmmTMB
library("bbmle") 
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  NCalls=SiblingNegotiation,
                  FT=FoodTreatment)


fit_zipoisson <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+
                             offset(log(BroodSize))+(1|Nest),
                         data=Owls,
                         ziformula=~1,
                         family=poisson)
summary(fit_zipoisson)


#Change the mean to standardized counts.
siv<-Sum
siv$Count<-round(siv$mean*1000000)

mean(siv$Count)
sd(siv$Count)^2
var(siv$Count)

#variance is greater than mean ==overdispersion
#need to sue negative bionomial distribution

po<-glmmTMB(Count~ Week+ Tissue2+Cohort+(1|Monkey), data=siv, family=poisson)
summary(po)

nb1<-glmmTMB(Count~ Week+ Tissue2+Cohort+(1|Monkey), data=siv, family=nbinom1)
summary(nb1)

nb2<-glmmTMB(Count~ Week+ Tissue2+Cohort+(1|Monkey), data=siv, family=nbinom2)
summary(nb2)

beta<-glmmTMB(mean~ Week+ Tissue2+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta)

library(DHARMa)
library(car)
library(effects)
library(MuMIn)
library(glmmLasso)
library(bbmle)

#look at the residuals
beta1<- simulateResiduals(beta)
plot(beta1)
Anova(beta)
#Response: mean
#Chisq Df Pr(>Chisq)  
#Week    2.3377  1    0.12628  
#Tissue2 8.4053  4    0.07781 .
#Cohort  5.4980  2    0.06399 .

#nb1 and beta have very similar results
#Response: Count
#         Chisq Df Pr(>Chisq)  
#Week    2.3348  1    0.12651  
#Tissue2 8.4275  4    0.07712 .
#Cohort  5.5031  2    0.06383 

#nb2
#Response: Count
#         Chisq Df Pr(>Chisq)  
#Week    2.5103  1    0.11310  
#Tissue2 9.9616  4    0.04108 *
#Cohort  5.6449  2    0.05946 .



re<-nb2
    
sim<- simulateResiduals(re)
plot(sim)
Anova(re)

#Generate a model selection table of models with combinations (subsets) of fixed 
#effect terms in the global model, with optional rules for model inclusion.

nb2.dredge<-MuMIn::dredge(nb2)
print(nb2.dredge)

beta.dredge<-MuMIn::dredge(beta)
print(beta.dredge)

#can visualized the model selection
op <- par(mar=c(3,4,7,4))
plot(beta.dredge)

# another way to print out the details of models tested
get.models(beta.dredge, subset=T)

#AIC comparison across models
beta<-glmmTMB(mean~ Week+ Tissue2+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta)

beta2<-glmmTMB(mean~ Week+ Tissue2+Cohort, data=Sum, family=beta_family())
summary(beta2)

beta3<-glmmTMB(mean~ Tissue2+Cohort, data=Sum, family=beta_family())
summary(beta3)

beta4<-glmmTMB(mean~ Tissue2+Cohort+(1|Monkey), data=Sum, family=beta_family())
summary(beta4)

AICtab(beta, beta2,beta3,beta4)
#      dAIC df
#beta   0.0 10
#beta4  0.3 9 
#beta2 21.3 9 
#beta3 24.9 8 

#Having the random factor (monkey)

formula(beta)[-2]
#dispersion parameter
sigma(beta)
#compare with no monkeys

nb1.2<-glmmTMB(Count~ Week+ Tissue2+Cohort, data=siv, family=nbinom1)
summary(nb1.2)

nb2.2<-glmmTMB(Count~ Week+ Tissue2+Cohort, data=siv, family=nbinom2)
summary(nb2.2)

beta2<-glmmTMB(mean~ Week+ Tissue2+Cohort, data=Sum, family=beta_family())
summary(beta2)




glmm1<-glmmTMB(mean~ Week+ Tissue2+Cohort+(1|Monkey), data=Sum)
summary(glmm1)

glmm2<-glmmTMB(mean~ Week+ Tissue2+Cohort, data=Sum, family=beta_family())
summary(glmm2)


glmm3<-glmmTMB(mean~ Week+Cohort+Tissue2+Granuloma, data=Sum, family=beta_family())
summary(glmm3)

library(DHARMa)
library(car)
library(effects)
library(MuMIn)

siv_beta1<- simulateResiduals(glmm2)
plot(siv_beta1)
Anova(glmm2)
#Response: mean
#          Chisq Df Pr(>Chisq)    
#Week     5.9175  1  0.0149914 *  
#Tissue2  9.8283  4  0.0434207 *  
#Cohort  16.9158  2  0.0002122 ***

Anova(glmm1, type=3)
#         Chisq Df Pr(>Chisq)  
#Week    2.3377  1    0.12628  
#Tissue2 8.4051  4    0.07782 .
#Cohort  5.4981  2    0.06399 .

#Visualize the effeects
ae<-allEffects(glmm2)
plot(ae)

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

