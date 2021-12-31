compare_and_overlay_TE <- function(m, TE, p.thresholds=c(RNA=0.01, Protein=0.001, TE=0.001), RNA.logFC.threshold=0.2, setTE0.pval=1, do.return=TRUE){ 
    TE <- homogenizeDEAresults(TE)
    m$logFC.TE <- TE[row.names(m),"logFC"]
    m$PValue.TE <- TE[row.names(m),"PValue"]
    m$FDR.TE <- TE[row.names(m),"FDR"]
    
    inComplex <- unique(unlist(complexes))
    m$inComplex <- factor(c("noComplex","inComplex")[as.numeric(row.names(m) %in% inComplex)+1],c("noComplex","inComplex"))

    m2 <- m[which( (m$PValue.RNA < p.thresholds["RNA"] & abs(m$logFC.RNA)>RNA.logFC.threshold) | m$PValue.Protein < p.thresholds["Protein"] | m$PValue.TE < p.thresholds["TE"]),]
    m2[which(is.na(m2$logFC.TE) | m2$PValue.TE>setTE0.pval ),"logFC.TE"] <- 0
    suppressPackageStartupMessages(library(plotly))
    cols <- colorRamp(c("red", "grey", "blue"))
    p <- plot_ly(m2,x=~logFC.RNA, y=~logFC.Protein, type="scatter", mode="markers", symbol=~inComplex, symbols=15:16, color=~logFC.TE, colors=cols, marker=list(size=10), hoverinfo="text", hovertext=paste(row.names(m2),round(m2$logFC.TE,2),sep="\n"))
    if(do.return)   return(list(m=m,p=p))
    p
}
