coeff.matrix<-function(x){
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
    return(data.frame(df))
    
}