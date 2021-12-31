crossLayer <- function(ll, minLfc=c(0.25,0.25,0.1), maxP=c(0.05,0.05,0.01), maxAgP=0.05, useFDR=c(T,T,F), full=FALSE){
  stopifnot(all(names(ll)==c("RNA","Protein","TE")))
  tt <- table(unlist(lapply(ll, FUN=row.names)))
  gu <- names(tt)[tt==max(tt)]
  degs <- lapply(setNames(1:3,names(ll)), FUN=function(i){
    if(is.null(ll[[i]])) return(NULL)
    d <- ll[[i]][gu,]
    sigf <- ifelse(useFDR[[i]], "FDR","PValue")
    row.names(d)[which(abs(d$logFC) >= minLfc[[i]] & d[[sigf]] <= maxP[[i]])]
  })
  m <- cbindDEAs(ll)
  ff <- paste(ifelse(useFDR, "FDR","PValue"),names(ll),sep=".")
  ap1 <- m[,ff[1:2]]
  ap1[which(sign(m$logFC.RNA)!=sign(m$logFC.Protein)),] <- c(1,1)
  m$agP <- apply(ap1,1,FUN=aggregation::fisher)
  degs$aggregated <- row.names(m)[which(p.adjust(m$agP) <= maxAgP)]
  degs$TE <- intersect(degs$TE, unlist(degs[-3]))
  degs$aggregated <- setdiff(degs$aggregated, unlist(degs[-4]))
  print(c(lengths(degs), "=total"=length(unique(unlist(degs)))))
  if(!full) m <- m[unique(unlist(degs)),]

  m$ProtRNA.residual <- m$logFC.Protein-m$logFC.RNA
  m$translational <- m[[ff[3]]] < maxP[3] & sign(m$ProtRNA.residual) == sign(m$logFC.TE) &
                abs(m$logFC.TE)>=minLfc[3] & abs(m$ProtRNA.residual)>0.2 &
                (abs(m$logFC.RNA) >= minLfc[1] | abs(m$logFC.Protein) >= minLfc[2]) &
                ( abs(m$ProtRNA.residual) >= 0.5 * apply(abs(as.matrix(m[,c("logFC.RNA","logFC.Protein")])),1,FUN=max)
                  | sign(m$logFC.RNA) != sign(m$logFC.Protein)) #& row.names(m) != "GTF2I"

  m$forwarded <- row.names(m) %in% unlist(degs[-3]) & sign(m$logFC.RNA) == sign(m$logFC.Protein) &
    abs(m$ProtRNA.residual) <= apply(abs(as.matrix(m[,c("logFC.RNA","logFC.Protein")])),1,FUN=min)

  m$buffered <- row.names(m) %in% degs$RNA & abs(m$logFC.Protein)<minLfc[2] & abs(m$ProtRNA.residual) >= minLfc[2] &
    abs(m$ProtRNA.residual) >= 0.5 * apply(abs(as.matrix(m[,c("logFC.RNA","logFC.Protein")])),1,FUN=max)
  m
}


plot3Layers <- function(m, genes=NULL, write=NULL){
  if(is.null(genes)) genes <- row.names(m)
  ggplot(m[unique(unlist(genes)),], aes(logFC.RNA,logFC.Protein,colour=logFC.TE,size=-log10(PValue.TE))) +
    geom_vline(xintercept=0, linetype="dashed", colour="grey") +
    geom_hline(yintercept=0, linetype="dashed", colour="grey") + geom_abline(slope=1) +
    geom_point() + scale_color_gradient2(low="blue", mid="grey", high="red")
}
