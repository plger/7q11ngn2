volcanoDEGs <- function(lrt, alpha = 0.05, lfc = log2(1.2), maxdot = 30, cols = c(viridis(25, option="B")[5:10], viridis(25, option="B")[15:20]), grid = T, labels = T, xlim=NULL, scalefactor = 1.5, sig.cex=0.6, ...){
  require(viridis)
  require(edgeR)
  require(plotrix)
  if(!is.data.frame(lrt)) lrt <- as.data.frame(topTags(lrt,Inf))
  if(!("FDR" %in% colnames(lrt)) & "adj.P.Val" %in% colnames(lrt)) lrt$FDR <- lrt$adj.P.Val
  lrt <- lrt[order(lrt$FDR, decreasing=F),]
  
  lrt <- lrt[which(lrt$FDR<0.9 | abs(lrt$logFC)>0.03),]
  w <- which(lrt$FDR==0)
  if(length(w)>0){
    lrt$FDR[w] <- min(lrt$FDR[-w],na.rm=T)/2
  }
  
  topsig <- rownames(lrt[which(lrt$FDR < alpha & abs(lrt$logFC) > lfc),])
  
  dotcol <- colorRampPalette(cols)(length(topsig) * scalefactor)
  
  dotcol_neg <- dotcol[1:length(which(lrt[topsig,"logFC"] < 0))]
  
  dotcol_pos <- dotcol[(length(dotcol) - length(which(lrt[topsig,"logFC"] > 0))) : length(dotcol)]
  
  if (length(topsig) < maxdot)
  { 
    md <-  length(topsig)
  } 
  else 
  {
    md = maxdot
  }
  
  if(all(grepl(".",row.names(lrt),fixed=T))){
    tmp <- sapply(strsplit(row.names(lrt),".",fixed=T),FUN=function(x){ x[2]})
    if(all(is.na(suppressWarnings(as.numeric(tmp))))) lrt$gene <- tmp
  }
  if ("gene" %in% colnames(lrt))
  {
    gene_label <- lrt[topsig,"gene"][1:md]
  }	
  else
  {
    gene_label <- rownames(lrt[topsig,])[1:md]
  }
  
  if(is.null(xlim)) xlim <- range(lrt$logFC,na.rm=T)
    plot(0,0,col="white",xlim=xlim,ylim=c(0,max(-log10(lrt$FDR),na.rm=T)),
        ylab = "-log10(FDR)",
        xlab = "log2(Fold Change)",
        ...
    )
    
  if (grid == T)
    
  {
    abline(
      h = seq(0, round(max(-log10(lrt$FDR), na.rm = T),0), by = 2), 
      lwd = 0.5,  
      col = "gray"
    )
    abline(
      v = seq(range(round(lrt$logFC, 0))[1], range(round(lrt$logFC, 0))[2], by = 0.5), 
      lwd = 0.5,  
      col = "gray"
    )
  }
  
  points(
    x = lrt$logFC, 
    y = -log10(lrt$FDR), 
    col = "#A9A9A964", 
    pch = 16, 
    cex = 0.5
  )
  points(
    x = lrt[topsig,"logFC"][order(lrt[topsig,"logFC"], decreasing=F)], 
    y = -log10(lrt[topsig,"FDR"][order(lrt[topsig,"logFC"], decreasing=F)]), 
    col = c(dotcol_neg, dotcol_pos), 
    pch = 16, 
    cex = 0.6
  )
  
  if (labels == T)
  {
    points(
      x = lrt[topsig,"logFC"][1:md], 	
      y = -log10(lrt[topsig,"FDR"])[1:md], 
      col = "black", 
      pch = 1, 
      cex = 1.2
    )
    if(length(topsig)>2){
        text(
        x = lrt[topsig,"logFC"][1:md], 
        y = -log10(lrt[topsig,"FDR"])[1:md],
        labels = gene_label,
        pos = thigmophobe(x = lrt[topsig,"logFC"][1:md],y = -log10(lrt[topsig,"FDR"])[1:md]) ,
        col = "black",
        cex = sig.cex,
        xpd=NA
        )
    }else{
        text(
        x = lrt[topsig,"logFC"][1:md], 
        y = -log10(lrt[topsig,"FDR"])[1:md],
        labels = gene_label,
        pos = 1,
        col = "black",
        cex = sig.cex,
        xpd=NA
        )
    }
  }
  
  abline(
    v = c(-lfc, lfc), 
    lwd = 1, 
    lty = 2, 
    col = "black"
  )
  
  abline(
    h = c(-log10(alpha)), 
    lwd = 1, 
    lty = 2, 
    col = "black"
  )
}
 
