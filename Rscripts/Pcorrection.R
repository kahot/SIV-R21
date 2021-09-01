Pcorrection<-function(x){
    df<-x
    df$Bonferroni <- p.adjust(df$rawP,method = "bonferroni")
    df$BH<-p.adjust(df$rawP,method = "BH")
    df$Holm<- p.adjust(df$rawP,method = "holm")
    df$Hochberg = p.adjust(df$rawP, method = "hochberg")
    df$Hommel = p.adjust(df$rawP, method = "hommel")
    df$BY = p.adjust(df$rawP, method = "BY")
    return(df)
}

