
WBSg1 <- c("GTF2IP1","NCF1B","GTF2IRD2P","POM121","NSUN5","TRIM50","FKBP6","FZD9","BAZ1B","BCL7B","TBL2","MLXIPL","VPS37D","DNAJC30","WBSCR22","STX1A","WBSCR26","ABHD11","CLDN3","CLDN4","WBSCR27","WBSCR28","ELN","LIMK1","EIF4H","LAT2","RFC2","CLIP2","GTF2IRD1","GTF2I","NCF1","GTF2IRD2","STAG3L2","PMS2P5","WBSCR16")
WBSg <- c("NSUN5","TRIM50","FKBP6","FZD9","BAZ1B","BCL7B","TBL2","MLXIPL","VPS37D","DNAJC30","WBSCR22","STX1A","WBSCR26","ABHD11","CLDN3","CLDN4","WBSCR27","WBSCR28","ELN","LIMK1","EIF4H","LAT2","RFC2","CLIP2","GTF2IRD1","GTF2I","NCF1")

#UPDATED SYMBOLS: USE THIS
WBSg2 <-c("GTF2IP1","NCF1B","GTF2IRD2P","POM121","NSUN5","TRIM50","FKBP6","FZD9","BAZ1B","BCL7B","TBL2","MLXIPL","VPS37D","DNAJC30","BUD23","STX1A","ABHD11-AS","ABHD11","CLDN3","CLDN4","METTL27","TMEM270","ELN","LIMK1","EIF4H","LAT2","RFC2","CLIP2","GTF2IRD1","GTF2I","NCF1","GTF2IRD2","STAG3L2","PMS2P5","RCC1L")


prepAssay <- function(x){
  if(all(grepl(".",head(row.names(x)),fixed=T))) x <- toSymbol(x);
  if(!("log2FC" %in% assayNames(x))){
    e <- SEtools:::.chooseAssay(x)
    assays(x)$log2FC <- e - rowMeans(e[,which(x$genotype=="CTRL")])
  }
  x
}


# max.agP is across layers and what determines the genes that get included in the analysis
# p.thres is the layer-specific p.value threshold
classifyGenes <- function(x, max.agP=0.01, p.thres=c(RNA=0.05,Protein=0.05,TE=0.1),abslogFC.thres=log2(1.2),agP.weights=c(RNA=3,Protein=2,TE=1),reportParams=TRUE,doPlot=TRUE,...){
    # we subset the matrix to genes that 'change in whatever'
    gclass <- rep(NA_character_,nrow(x))
    xp <- x[,paste0("PValue.",c("RNA","Protein","TE"))]
    library(aggregation)
    x$agP <- apply(cbind(apply(xp,1,na.rm=T,FUN=min), apply(xp,1,weights=c(2,2,1),FUN=lancaster)),1,na.rm=T,FUN=min)
    w <- which(x$agP<max.agP)
    if(length(w)==0){
        message("Nothing matches criteria...")
        return(NULL)
    }
    if(reportParams)    message(paste0("Selecting ",length(w)," genes with aggregated p<",max.agP,".
An alteration is considered if the absolute logFC is greater than ",abslogFC.thres," and if the p-value is respectively ",paste(paste(names(p.thres),"<",p.thres),collapse=", ")))
    x <- x[w,,drop=F]
    xp <- x[,paste0("PValue.",c("RNA","Protein","TE")),drop=F]
    x2 <- x[,paste0("logFC.",c("RNA","Protein","TE")),drop=F]
    # we set logFC to 0 when p is above the threshold
    for(i in 1:3) x2[which(xp[,i]>p.thres[i]),i] <- 0
    x2 <- x2[,c(1,3,2)]
    # we classify using logFC
    y <- apply(1+ (x2>abslogFC.thres) - (x2 < -abslogFC.thres), 1, collapse="", FUN=paste)
    # the digits represent respectively RNA, TE and Protein
    # 0 means down, 1 means unchanging, and 2 means up
    map <- c(   "222"="Amplified (+)",
                "212"="Forwarded (+)",
                "000"="Amplified (-)",
                "010"="Forwarded (-)",
                "201"="Buffered translationally",
                "021"="Buffered translationally",
                "202"="Buffered translationally (inefficient)",
                "020"="Buffered translationally (inefficient)",
                "200"="Over-corrected translationally",
                "022"="Over-corrected translationally",
                "211"="Buffered",
                "011"="Buffered",
                "001"="Buffered post-translationally",
                "221"="Buffered post-translationally",
                #"121"="Buffered",
                #"101"="Buffered",
                "112"="Post-translational (+)",
                "110"="Post-translational (-)",
                "122"="Translational (+)",
                "100"="Translational (-)",
                "220"="Complex",
                "210"="Complex",
                "012"="Complex",
                "002"="Complex" )
    gclass[w] <- map[y]
    if(!doPlot) return(gclass)

    # plotting...

    pchs= c(
        "Buffered"=4,
        "Buffered translationally"=15,
        "Buffered translationally (inefficient)"=7,
        "Forwarded"=16,
        "Translational"=17,
        "Post-translational"=3,
        "Complex"=8,
        "Amplified"=18,
        "Over-corrected translationally"=11 )

    x$class <- map[y]
    x$class2 <- gsub(" (+)","",x$class, fixed=T)
    x$class2 <- gsub(" (-)","",x$class2, fixed=T)
    x <- x[which(!is.na(x$class)),]
    ct <- table(x$class2)
    pchs <- pchs[names(ct)]
    names(pchs) <- paste0(names(pchs)," (",ct[names(pchs)],")")
    x$class3 <- paste0(x$class2," (",ct[x$class2],")")
    p <- suppressMessages(suppressWarnings(plot_ly(x,x=~logFC.RNA, y=~logFC.Protein, type="scatter", mode="markers", symbol=~class3, symbols=pchs, color=~class3, marker=list(size=10, apha=0.7), hoverinfo="text", hovertext=paste(row.names(x),round(x$logFC.TE,2),sep="\n"))))
    return(list(classes=gclass, p=p))
}

translateClasses <- function(x){
  x <- gsub(" \\(|\\)|\\+|\\-","",x)
  tr <- c("Amplified"="Forwarded",
          "Buffered translationally"="Buffered",
          "Buffered translationally (inefficient)"="Buffered",
          "Buffered posttranslationally"="Buffered",
          "Overcorrected translationally"="Post-transcriptionally regulated",
          "Overcorrected"="Post-transcriptionally regulated",
          "Translational"="Post-transcriptionally regulated",
          "Complex"="Post-transcriptionally regulated",
          "Posttranslational"="Post-transcriptionally regulated"
          )
  w <- which(x %in% names(tr))
  x[w] <- tr[x[w]]
  x

}


# m is the output of `classifyGenes`
reclassify <- function(m, minAbsLFC=log2(1.2), resth.fun=function(x){ 2*median(abs(x)) }, fw.coef=NULL ){
  library(MASS)
    mo <- m

    layout(matrix(1:2,nrow=1))

    # we first used the genes classified as forwarded in order to establish a
    # normal relationship between RNA and protein foldchanges
    m.f <- m[grep("Forwarded",m$class),,drop=F]
    if(nrow(m.f)<2){
        warning("Too few classified genes to bootstrap reclassification")
        return(m)
    }
    # robust lm forced through 0
    if(is.null(fw.coef)){
      mod <- rlm(logFC.Protein~0+logFC.RNA,data=m.f)
      if(coefficients(mod)[[1]]<0){
        warning("Negative model coefficient!")
        return(m)
      }
      fw.coef <- coefficients(mod)[[1]]
      resth <- resth.fun(mod$residuals)
    }else{
      resth <- resth.fun(m$logFC.Protein-(m$logFC.RNA*fw.coef))
    }

    # we then use the residuals to establish forwarded genes
    m$RNAprot.residual <- m$logFC.Protein-(m$logFC.RNA*fw.coef)
    w <- which( ( abs(m$RNAprot.residual)<=resth | abs(m$RNAprot.residual) < minAbsLFC) &
                abs(m$logFC.RNA) >= minAbsLFC &
                #(abs(m$logFC.Protein) >= minAbsLFC | ( (minAbsLFC-abs(m$logFC.Protein)) < abs(m$RNAprot.residual) )) &
                abs(m$logFC.Protein) >= minAbsLFC &
                sign(m$logFC.RNA)==sign(m$logFC.Protein) )
    mo[w,"class"] <- paste0("Forwarded (",c("-","+","+")[sign(m[w,"logFC.RNA"])+2],")")

    # buffered genes have residuals outside the threshold towards the y origin
    res2 <- m$RNAprot.residual*sign(m$logFC.RNA)
    wB <- which( res2 < -resth &
                 mo$class=="Buffered" &
                 abs(m$logFC.RNA) >= minAbsLFC  )
    # regions 'buffered' in the previous definition, which show an inversion in the
    # direction of the foldchange beyond that threshold are considered over-corrected
    wB <- which( res2 < -resth &
                 abs(m$logFC.RNA) >= minAbsLFC &
                 abs(m$logFC.Protein) >= minAbsLFC &
                 sign(m$logFC.RNA) != sign(m$logFC.Protein))
    mo[wB,"class"] <- "Over-corrected"

    # amplified genes have residuals outside the threshold and away from the y origin
    wA <- which(  res2 > resth &
                  abs(m$logFC.RNA) >= minAbsLFC &
                  abs(m$logFC.Protein) >= minAbsLFC  )
    mo[wA,"class"] <- paste0("Amplified (",c("-","+","+")[sign(m[wA,"logFC.RNA"])+2],")")

    plot(m.f$logFC.RNA,m.f$logFC.Protein,pch=20,xlab="RNA log2FC", ylab="Protein log2FC", main="Forwarded genes")
    ablines0()
    abline(a=0,b=fw.coef,lwd=2)
    abline(a=resth,b=fw.coef,lty="dashed",col="blue")
    abline(a=-resth,b=fw.coef,lty="dashed",col="blue")
    polygon(c(-minAbsLFC,minAbsLFC,minAbsLFC,-minAbsLFC,-minAbsLFC),c(minAbsLFC,minAbsLFC,-minAbsLFC,-minAbsLFC,minAbsLFC),lty="dashed",border="grey")
    legend("topleft",bty="n", legend=paste0("r^2=",round(cor(m.f$logFC.RNA,m.f$logFC.Protein),2)))

    m.t <- m[grep("^Translational|translationally",m$class),,drop=F]

    if(nrow(m.t)>2){

        plot(m.t$logFC.TE,m.t$RNAprot.residual,pch=16,xlab="TE log2FC", ylab="RNA-protein residuals", main="TE-buffered genes")
        ablines0()
        legend("topleft",bty="n", legend=paste0("r^2=",round(cor(m.t$logFC.TE,m.t$RNAprot.residual),2)))

        mod <- try(rlm(RNAprot.residual~0+logFC.TE,data=m.t),silent=T)

        if(is(mod,"try-error")){

            legend("topright",legend="no fit")

        }else{

            abline(mod,lwd=2)
            if(coefficients(mod)[[1]]<0){
                warning("Negative model coefficient!")
                return(mo)
            }

            m$TEprot.residual <- m$RNAprot.residual-(m$logFC.TE*coefficients(mod)[[1]])
            resth2 <- resth.fun(mod$residuals)
            w <- which( abs(m$TEprot.residual)<=resth2 &
                        abs(m$logFC.TE) >= minAbsLFC &
                        abs(m$logFC.RNA) >= minAbsLFC &
                        sign(m$logFC.RNA)!=sign(m$logFC.Protein) )
            w2 <- which(abs(m$logFC.Protein)>=minAbsLFC & sign(m$logFC.Protein)==sign(m$logFC.RNA))
            mo[intersect(w,w2),"class"] <- "Buffered translationally (inefficient)"
            w3 <- which(abs(m$logFC.Protein)>=minAbsLFC & sign(m$logFC.Protein)!=sign(m$logFC.RNA))
            mo[intersect(w,w3),"class"] <- "Over-corrected translationally"
            mo[setdiff(w,c(w2,w3)),"class"] <- "Buffered translationally"

            abline(a=resth2,b=coefficients(mod)[[1]],lty="dashed",col="blue")
            abline(a=-resth2,b=coefficients(mod)[[1]],lty="dashed",col="blue")
            polygon(c(-minAbsLFC,minAbsLFC,minAbsLFC,-minAbsLFC,-minAbsLFC),c(minAbsLFC,minAbsLFC,-minAbsLFC,-minAbsLFC,minAbsLFC),lty="dashed",border="grey")

        }
    }else{
      message("Insufficient number of high-confidence translationally-regulated genes for residuals model. Using distribution of TE logFC for forwarded genes as null hypothesis.")
      te.fw <- abs(mo$logFC.TE[grep("Forwarded",mo$class)])
      te.fw <- te.fw[!is.na(te.fw)]
      te.fw <- c(te.fw, -te.fw)
      d <- fitdistr(te.fw,"normal")
      p <- pnorm(mo$logFC.TE, d$estimate, d$sd)
      mo$RNAprot.residual <- res <- m[row.names(mo),"RNAprot.residual"]
      mo$TE.p <- p
      hist(te.fw)
      #LSD::heatscatter(res[p<0.01],mo$logFC.TE[p<0.01])
      w <- which( p<0.1 & sign(res)==sign(mo$logFC.TE) & !is.na(mo$class) &
                    !grepl("Forwarded",mo$class) &
                    abs(mo$logFC.TE) > minAbsLFC/4 & abs(res)>minAbsLFC/2 )
      if(length(w)>0) mo[w,"class"] <- paste(mo[w,"class"],"translationally")
    }
    mo
}

multileveltest <- function(prot, rna, te){
  i <- intersect(intersect(row.names(prot), row.names(te)), row.names(rna))
  prot <- assay(prot)[i,]
  rna <- assay(rna)[i,]
  te <- assay(te)[i,]
  e <- t(sapply(1:length(i), FUN=function(x){
    mod <- lm(p~1+r+te,data=data.frame(p=prot[x,],r=rna[x,],te=te[x,]))
    tryCatch(coef(summary(mod))["te",c(1,4)], error=function(e) c(0,1))
  }))
  row.names(e) <- i
  e
}

plotClassification <- function(m, classColumn="class", aggregateDirections=TRUE, dropUnclassified=TRUE, psize=8, use.plotly=TRUE){
    if(dropUnclassified)    m <- m[which(!is.na(m[[classColumn]])),]
    class2 <- m[[classColumn]]
    if(aggregateDirections){
        class2 <- gsub(" (+)","",class2, fixed=T)
        class2 <- gsub(" (-)","",class2, fixed=T)
    }
    ct <- table(class2)
    w <- which(!is.na(class2))
    class2[w] <- paste0(class2[w]," (",ct[class2[w]],")")
    if(use.plotly){
      library(plotly)
      return(suppressMessages(suppressWarnings({
        plot_ly( m, x=~logFC.RNA, y=~logFC.Protein, type="scatter", mode="markers",
                 symbol=class2, color=class2, marker=list(size=psize, apha=0.7),
                 hoverinfo="text", hovertext=paste(row.names(m),round(m$logFC.TE,2),sep="\n"))
      })))
    }else{
      library(ggplot2)
      return(ggablines(ggplot(m, aes(logFC.RNA, logFC.Protein, symbol=class2, colour=class2, size=-log10(PValue.TE))) + geom_point()))
    }
}

.pa <- function(x){
  x <- as.numeric(x)
  if(all(x<0.05)){
    return(fisher(x))
  }
  #return(2^mean(log2(x)))
  mean(x)
}
geom <- function(x){ 2^mean(log2(x)) }

cbindDEAs <- function(ll){
  ll <- lapply(ll, FUN=homogenizeDEAresults)
  ln <- names(ll)
  ln <- setdiff(ln, c("reg","cn"))
  for(n in ln){
    ll[[n]]$baseMean <- NULL
    colnames(ll[[n]]) <- paste0(colnames(ll[[n]]),".",n)
  }
  rn <- row.names(ll[[1]])
  ll <- lapply(ll,rn=rn,function(x,rn) x[rn,])
  names(ll) <- NULL
  res <- do.call(cbind, ll)
  if("PValue" %in% colnames(res)) res <- res[order(res$PValue),]
  res
}

mergeDEAs <- function(x,y,PValue=mean,FDR=geom,sameDirectionOnly=TRUE){
  library(aggregation)
  i <- intersect(row.names(x),row.names(y))
  x <- homogenizeDEAresults(x)[i,]
  y <- homogenizeDEAresults(y)[i,]
  d <- data.frame( row.names=i,
                   logFC=(x$logFC+y$logFC)/2,
                   PValue=apply(cbind(x$PValue,y$PValue),1,FUN=PValue),
                   FDR=NA_real_
                 )
  if(sameDirectionOnly){
    w <- which( sign(x$logFC)==sign(y$logFC) )
  }else{
    w <- 1:nrow(x)
  }
  if(is.character(FDR) && FDR=="recompute"){
    d$FDR[w] <- p.adjust(d$PValue[w],method="BY")
  }else{
    d$FDR[w] <- apply(cbind(x$FDR,y$FDR)[w,],1,FUN=FDR)
  }
  for(f in c("wbs","dup")){
    if( paste0("logFC.",f) %in% colnames(x) && paste0("logFC.",f) %in% colnames(x) ){
      m2 <- mergeDEAs(subset_mDEA(x,f),subset_mDEA(x,f))
      colnames(m2) <- paste0(colnames(m2),".",f)
      d <- cbind(d, m2[row.names(d),])
    }
  }
  d[order(d$FDR, d$PValue),]
}

mergeDEAs2 <- function(x, y, min.x.fdr=0.05, min.y.p=0.001, min.y.fdr=1, fc.con=1, min.lfc=log2(1.2)){
  i <- intersect(row.names(x),row.names(y))
  x <- homogenizeDEAresults(x)[i,]
  y <- homogenizeDEAresults(y)[i,]
  d <- data.frame( row.names=i,
                   logFC=(x$logFC+y$logFC)/2,
                   PValue=10^rowMeans(log10(cbind(x$PValue,y$PValue))),
                   FDR=10^rowMeans(log10(cbind(x$FDR,y$FDR))) )
  if("stat" %in% colnames(x) && "stat" %in% colnames(x)){
    d$stat <- (x$stat+y$stat)/2
  }
  d$sig <- sign(x$logFC)==sign(y$logFC) &
            x$FDR < min.x.fdr & y$FDR < min.y.fdr &
            y$PValue < min.y.p & abs(x$logFC)>min.lfc & abs(y$logFC)>min.lfc &
            abs(x$logFC-y$logFC)/(x$logFC+y$logFC) < fc.con
  d$FDR[which(!d$sig)] <- 1
  d$PValue[which(sign(x$logFC)!=sign(y$logFC) & (x$FDR > min.x.fdr | y$FDR > min.y.fdr | y$PValue > min.y.p))] <- 1
  for(f in c("WBS","DUP","DUPvsWBS")){
    if( paste0("logFC.",f) %in% colnames(x) && paste0("logFC.",f) %in% colnames(x) ){
      m2 <- mergeDEAs2(subset_mDEA(x,f),subset_mDEA(y,f), min.x.fdr=min.x.fdr, min.y.p=min.y.p, min.y.fdr=min.y.fdr, fc.con=fc.con, min.lfc=min.lfc)
      colnames(m2) <- paste0(colnames(m2),".",f)
      d <- cbind(d, m2[row.names(d),])
    }
  }
  d[order(d$FDR, d$PValue),]
}


# quick wrapper for aggregation
plag <- function(x, by, agFun=sum, ...){
    ag <- aggregate(x,by=list(by),...,FUN=agFun)
    row.names(ag) <- ag[,1]; ag[,1] <- NULL
    return(ag)
}

# wrapper around pheatmap, handles logicals and empty columns in the annotation
plhm <- function(x, scale="none", annotation_row=NA, color="by", annotation_col=NA, annotation_colors=NA, ...){
    library("viridis")
    library("pheatmap")
    scale <- match.arg(scale, c("none", "row", "column"))
    x <- x[which(apply(x, 1, FUN = function(y) {
        !all(is.na(y))
    })), which(apply(x, 2, FUN = function(y) {
        !all(is.na(y))
    }))]
    if(!is.list(annotation_colors)) annotation_colors <- list()
    if(is.list(annotation_row)){
        ar <- annotation_row
        for(i in colnames(ar)){
            if(is.logical(ar[[i]])) ar[[i]] <- as.numeric(ar[[i]])
            if(is.numeric(ar[[i]]) && all(ar[[i]] %in% c(0,1,NA))){
                if(all(ar[row.names(x),i]==0)){
                    ar[[i]] <- NULL
                }else{
                    ar[[i]] <- factor(c("No","Yes")[1+ar[[i]]])
                    annotation_colors[[i]] <- c(No="white",Yes="darkblue")
                }
            }
        }
        annotation_row <- ar
    }
    color <- switch(color,
        by=colorRampPalette(c("blue", "black", "yellow"))(29),
        viridis=viridis::viridis(50),
        magma=viridis::viridis(50,option="B"),
        color)
    pheatmap(x, scale = scale, color=color, drop_levels=T, border_color = NA, annotation_row=annotation_row, annotation_col=annotation_col, annotation_colors=annotation_colors, ...)
}

voomDupCor <- function(counts, mm, individuals, tcoefs=NULL){
    library(limma)
    library(edgeR)
    if(is.null(tcoefs)) tcoefs <- colnames(mm)[ncol(mm)]
    dds <- calcNormFactors(DGEList(counts))
    v <- voom(dds,mm)
    dc <- duplicateCorrelation(v,mm,block=individuals)
    v <- voom(dds,mm,block=individuals,correlation=dc$consensus.correlation)
    fit <- eBayes(lmFit(v,mm,block=individuals,correlation=dc$consensus.correlation))
    res <- as.data.frame(topTable(fit,tcoefs,50000))
    res$B <- NULL
    res$t <- NULL
    return(res)
}

limmaDupCor <- function(normIntensities, mm, individuals, tcoefs=NULL){
    library(limma)
    normIntensities <- log(normIntensities+1)
    if(is.null(tcoefs)) tcoefs <- colnames(mm)[ncol(mm)]
    dc <- duplicateCorrelation(normIntensities,mm,block=individuals)
    fit <- eBayes(lmFit(normIntensities,mm,block=individuals,correlation=dc$consensus.correlation))
    res <- as.data.frame(topTable(fit,tcoefs,50000))
    res$B <- NULL
    res$t <- NULL
    return(res)
}


subset_mDEA <- function(x, suffix=NULL){
  if(is.null(suffix) || suffix==""){
    x <- x[,intersect(c("logFC","PValue","FDR","stat"),colnames(x))]
  }else{
    i <- paste0(c("logFC.","PValue."),suffix)
    if(paste0("FDR.",suffix) %in% colnames(x)) i <- c(i,paste0("FDR.",suffix))
    if(paste0("stat.",suffix) %in% colnames(x)) i <- c(i,paste0("stat.",suffix))
    x <- x[,i]
    colnames(x) <- c("logFC","PValue","FDR","stat")[1:length(i)]
  }
  return(x)
}

plotFCcomp <- function(a, b, nameA="A", nameB="B", use.plotly=FALSE, units="log2FC", restrictTo=1, useFDR=FALSE,
                cols=c(A="#DDCC77", B="#CC6677", both="#4477AA", none="#00000020"),
                pchs=c(A=16,B=15,both=17,none=1),
                cexs=c(A=1.2,B=1.2,both=1.5,none=0.5),
                main=NULL, pthres=0.01, legPos="bottomright", bty="o",
                label.Alfc=0.5, label.Blfc=0.2, label.type="both", suppressPlot=FALSE,
                label.which.fun=function(x,label.Alfc,label.Blfc,label.type){ abs(x[["logFC.x"]])>label.Alfc & x[["type"]] %in% label.type & abs(x[["logFC.y"]])>label.Blfc },
                ...){
    if(length(pthres)==1) pthres <- c(pthres,pthres)
    if(!("PValue" %in% colnames(a)) && "P.Value" %in% colnames(a)) a$PValue <- a$P.Value
    if(!("PValue" %in% colnames(b)) && "P.Value" %in% colnames(b)) b$PValue <- b$P.Value
    if(length(pchs)==1) pchs <- c(A=pchs,B=pchs,both=pchs,none=pchs)
    m <- merge(a,b,by="row.names")
    m$type <- "none"
    if(length(useFDR)==1) useFDR <- rep(useFDR,2)
    if(length(pthres)==1) pthres <- rep(pthres,2)
    m <- m[which(m$PValue.x <= restrictTo | m$PValue.y <= restrictTo),]
    sigA <- which(m[[ifelse(useFDR[1],"FDR.x","PValue.x")]]<=pthres[1])
    sigB <- which(m[[ifelse(useFDR[2],"FDR.y","PValue.y")]]<=pthres[2])
    m$type[sigA] <- "A"
    m$type[sigB] <- "B"
    if(length(intersect(sigA,sigB))>0) m$type[intersect(sigA,sigB)] <- "both"
    tc <- table(m$type)
    labs <- paste(units, c(nameA,nameB))
    useFDR <- unique(useFDR)
    pthres <- unique(pthres)
    if(is.null(main)){
      t2 <- paste(ifelse(useFDR,"q","p"),"<=",pthres)
      if(length(t2)>1) t2 <- paste(c(nameA, nameB), t2)
      main <- paste0("Foldchange comparison (",paste(t2,collapse=", "),")")
    }
    if(use.plotly){
        suppressPackageStartupMessages(library(plotly))
        m$type <- factor(m$type, levels=c("A","B","both","none"))
        return( plot_ly(data = m, x = ~logFC.x, y = ~logFC.y, text=m[,1],
                        type="scatter", mode="markers", hoverinfo="text", alpha = 0.8,
                        color=~type, colors=as.character(cols), symbol=~type, symbols=as.numeric(pchs),
                        marker = list(size = 10), xaxis = list( title = labs[1] ), yaxis = list( title = labs[2] )) )
    }else{
      if(!is.null(label.which.fun) & is.function(label.which.fun)){
        w <- which(sapply(split(m,1:nrow(m),drop=F),label.Alfc=label.Alfc,label.Blfc=label.Blfc,label.type=label.type,FUN=label.which.fun))
        if(length(w)>0){
          m$toPlot <- FALSE
          m$toPlot[w] <- TRUE
        }
      }
      if(!suppressPlot){
        # eventually replace with fcPlot
          plot(m$logFC.x, m$logFC.y, col=cols[m$type], pch=pchs[m$type], cex=cexs[m$type], xlab=labs[1], ylab=labs[2], main=main, ...)
          if(length(w)>0){
            library(wordcloud)
            try(suppressWarnings(textplot(m$logFC.x[w],m$logFC.y[w],m[w,1],new=F,xpd=NA,pos=3)),silent=T)
          }
          ablines0()
          legend(legPos, text.col=cols[1:3], col=cols[1:3], pch=c(pchs[1:3]), legend=paste0(c(nameA,nameB,"both")," (",tc[c("A","B","both")],")"), bty=bty, text.font=2)
      }
    }
    row.names(m) <- m[,1]
    m[,1] <- NULL
    colnames(m) <- gsub("\\.x$",paste0(".",nameA),colnames(m))
    colnames(m) <- gsub("\\.y$",paste0(".",nameB),colnames(m))
    return(m)
}

fcPlot <- function(m, label.which=NA, nlabs=c(A="A", B="B"), cols=c(A="#DDCC77", B="#CC6677", both="#4477AA", none="#00000020"), pchs=c(A=16,B=15,both=17,none=20), cexs=c(A=1.6,B=1.6,both=2.2,none=0.5), write.cor=TRUE){
  n <- c(A="A", B="B", both="both", none="none")
  if(length(label.which)==1 && is.na(label.which)){
    label.which <- c()
    if(!is.null(m$toPlot)) label.which <- which(m$toPlot)
  }
  for(f in names(nlabs)) n[[f]] <- nlabs[[f]]
  m$type <- as.factor(m$type)
  nb <- paste0(n[levels(m$type)], " (", as.numeric(table(m$type)), ")")
  names(nb) <- levels(m$type)
  nb[["none"]] <- "none"
  nlabs <- paste0("logFC.",nlabs)
  p <- ggplot(m, aes_string(nlabs[[1]], nlabs[[2]], colour="type", shape="type", size="type"))
  if(suppressWarnings(require("ggraster", quietly=TRUE))){
    p <- p + rasterise(geom_point(), dpi=150)
  }else{
    p <- p + geom_point()
  }
  p <- p +scale_shape_manual(values=pchs, labels=nb) + scale_size_manual(values=cexs, labels=nb) +
    scale_color_manual(values=cols, labels=nb)
  if(length(label.which)>0){
    m$gene <- row.names(m)
    p <- p + geom_text_repel(data=m[label.which,,drop=FALSE], aes(logFC.RNA, logFC.Protein, label=gene), size=4, show.legend=FALSE)
  }
  if(write.cor){
    pc <- cor.test(m[which(m$type!="none"),nlabs[[1]]],m[which(m$type!="none"),nlabs[[2]]])
    lab <- paste0(" r=",round(pc$estimate,2),"\n p~",format(pc$p.value,digits = 1))
    p <- p + annotate("text", -Inf, Inf, label=lab, hjust=0, vjust=1)
  }
  ggablines(p)
}

fcPlots <- function(a,b,genes=NULL,ghighlight=WBSg,labels=c("RNA","Protein"),type=FALSE,hl_labs=FALSE,TE=NULL,write.cor=TRUE,...){
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  cntrst <- c("Regression on copynumber"="","WBS vs CTRL"="WBS", "DUP vs CTRL"="DUP")
  if(!is.null(genes)){
    if(!is.list(genes)){
      genes <- list(genes,genes,genes)
    }
    if(is.null(names(genes))) names(genes) <- names(cntrst)
    genes <- lapply(genes, y=intersect(row.names(a), row.names(b)), FUN=intersect)
  }
  if(!is.null(ghighlight)){
    if(!is.list(ghighlight)) ghighlight <- list(ghighlight, ghighlight, ghighlight)
    if(is.null(names(ghighlight))) names(ghighlight) <- names(cntrst)
  }
  lapply( names(cntrst), FUN=function(xn){
    if(!is.null(genes)){
      a <- a[genes[[xn]],]
      b <- b[genes[[xn]],]
    }
    if(nrow(a)<2) return(NULL)
    x <- cntrst[xn]
    if(x=="") x <- NULL
    m <- plotFCcomp(subset_mDEA(a, x), subset_mDEA(b, x), labels[1], labels[2], ..., suppressPlot=TRUE)
    labs <- paste0("logFC.",labels)
    m <- m[!apply(m[,labs],1,FUN=function(x) any(is.na(x))),]
    if(!is.null(TE)){
      te <- subset_mDEA(TE,x)
      colnames(te) <- paste0(colnames(te),".TE")
      m <- merge(m,te,by="row.names",all.x=TRUE)
      row.names(m) <- m[,1]
      colnames(m)[1] <- "gene"
      m <- m[order(m$PValue.TE, decreasing=TRUE),]
    }else{
      m$gene <- row.names(m)
    }
    names(labels) <- LETTERS[1:2]
    if(type){
      p <- fcPlot(m, c(), labels)
    }else{
      if(!is.null(TE)){
        tr <- scales::trans_new("signedsqrt", function(x) sign(x)*sqrt(abs(x)), function(x) sign(x)*x^2)
        p <- ggplot(m, aes_string(labs[[1]], labs[[2]], colour="logFC.TE", size="-log10(PValue.TE)")) +
          geom_point() + scale_size_continuous(range=c(0.3,4),name="TE\n-log10(p)")
        p <- tescale(p)
        p <- p + geom_point(data=m[which(is.na(m$logFC.TE)),], mapping=aes_string(labs[[1]], labs[[2]]), size=0.3, colour="lightgrey")
      }else{
        p <- ggplot(m, aes_string(labs[[1]], labs[[2]])) + geom_point()
      }
      if(write.cor){
        pc <- cor.test(m[,labs[[1]]],m[,labs[[2]]])
        lab <- paste0(" r^2=",round(pc$estimate,2),"\n p~",format(pc$p.value,digits = 1))
        p <- p + annotate("text", -Inf, Inf, label=lab, hjust=0, vjust=1)
      }
    }
    if(!is.null(ghighlight) & any(ghighlight[[xn]] %in% row.names(m))){
      labels <- paste0("logFC.",labels)
      m2 <- m[intersect(row.names(m),ghighlight[[xn]]),,drop=FALSE]
      if(hl_labs){
        p <- p + geom_label_repel(data=m2, mapping=aes_string(labels[[1]],labels[[2]], label="gene"))
      }else{
        p <- p + geom_point(data=m2, mapping=aes_string(labels[[1]],labels[[2]]), colour="blue")
      }
    }
    ggablines(p) + geom_abline(slope = 1, intercept = 0, colour="darkgrey") + ggtitle(xn)
  })
}

tescale <- function(p){
  tr <- scales::trans_new("signedsqrt", function(x) sign(x)*sqrt(abs(x)), function(x) sign(x)*x^2)
  p + scale_color_gradient2(trans=tr, low="darkred", mid="lightgrey", high="darkblue", na.value = "lightgrey", breaks=c(-0.5,0,0.5))
}

mergelayers <- function(rna, prot, te, keep.all=TRUE){
  colnames(rna) <- paste(colnames(rna),"RNA",sep=".")
  colnames(prot) <- paste(colnames(prot),"Protein",sep=".")
  m <- merge(rna, prot, by="row.names", all=keep.all)
  colnames(te) <- paste(colnames(te),"TE",sep=".")
  m <- merge(m, te, by.x="Row.names", by.y="row.names", all=keep.all)
  row.names(m) <- m[,1]
  m
}
layersplot <- function(m, genes, suffix=NULL, write.cor=FALSE){
  vars <- c("RNA", "Protein", "TE")
  if(!is.null(suffix)) vars <- paste0(suffix,".",vars)
  vars <- paste0("logFC.",vars)
  m <- m[intersect(genes,row.names(m)),]
  p <- ggplot(m, aes_string(vars[1], vars[2], colour=vars[3], size=paste0("-log10(",gsub("logFC","PValue", vars[3]),")"))) + geom_point() + geom_abline(slope = 1,intercept = 0)
  if(write.cor){
    pc <- cor.test(m[,vars[[1]]],m[,vars[[2]]])
    lab <- paste0(" r^2=",round(pc$estimate,2),"\n p~",format(pc$p.value,digits = 1))
    p <- p + annotate("text", -Inf, Inf, label=lab, hjust=0, vjust=1)
  }
  ggablines(tescale(p))
}


ggablines <- function(p){
  p + geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") + geom_vline(xintercept=0, linetype="dashed", colour="darkgrey")
}


whichGene4Integration <- function(x){
    which((  (m$PValue.RNA < 0.01 & abs(m$logFC.RNA)>0.2) +
            (m$PValue.Protein < 0.01) +
            (m$PValue.TE < 0.01) ) >1 )
}

ablines0 <- function(lty="dashed",col="grey",...){
    abline(v=0,lty=lty,col=col,...)
    abline(h=0,lty=lty,col=col,...)
}

mergeDEAs_old <- function(x,y,sig.ag.fun=mean){
    f <- function(p1,p2,lfc1,lfc2,sig.ag.fun=mean){
        lfc <- rowMeans(cbind(lfc1,lfc2))
        p <- apply(cbind(p1,p2),1,FUN=sig.ag.fun)
        p[which(sign(lfc1)!=sign(lfc2))] <- NA
        fdr <- p.adjust(p,method="fdr")
        cbind(logFC=lfc,PValue=p,FDR=fdr)
    }
    m <- merge(x,y,by="row.names")
    m2 <- as.data.frame(f(m$PValue.x,m$PValue.y,m$logFC.x,m$logFC.y,sig.ag.fun), row.names=m[,1])
    if( "logFC.wbs" %in% colnames(x) && "logFC.wbs" %in% colnames(y) ){
        wbs <- f(m$PValue.wbs.x,m$PValue.wbs.y,m$logFC.wbs.x,m$logFC.y,sig.ag.fun)[,-3]
        colnames(wbs) <- paste0(colnames(wbs),".wbs")
        dup <- f(m$PValue.dup.x,m$PValue.dup.y,m$logFC.dup.x,m$logFC.y,sig.ag.fun)[,-3]
        colnames(dup) <- paste0(colnames(dup),".dup")
        m2 <- cbind(m2,wbs,dup)
    }
    return(m2[order(m2$FDR,m2$PValue),])
}

filterGenesByCount <- function(x,conditions=NULL,minc=20,propSmallestGroup=0.75){
    if(is.null(conditions)){
        ns <- max(2,floor(0.2*ncol(x)))
    }else{
        tmp <- table(conditions)
        tmp <- tmp[which(tmp>1)]
    }
    ns <- floor(propSmallestGroup*min(tmp))
    x[which(apply(x,1,ns=ns,minc=minc,FUN=function(x,ns,minc){ sum(x>=minc)>=ns })),]
}

plot_nDEGs_by_threshold <- function(ll,field="FDR",logFC.thres=0,sig.thres=c(0.1,0.05,0.01,0.005,10^(-3:-10)),lwd=3,pch=20,type="b",ylim=NULL,...){
    if(!is(ll,"list")) ll <- list(ll)
    if(length(ll)>1){
        cols <- getQualitativePalette(length(ll))
        if(is.null(names(ll))) names(ll) <- 1:length(ll)
    }else{
        cols <- "black"
    }
    m <- as.data.frame(sapply(ll,field=field,sig.thres=sig.thres,logFC.thres=logFC.thres,FUN=function(x,field,sig.thres,logFC.thres){
        sapply(sig.thres,sig=x[[field]][which(abs(x$logFC)>logFC.thres)],FUN=function(y,sig){ sum(sig<=y,na.rm=T)})
    }))
    if(length(ll)>1){
        m$common <- sapply(sig.thres,field=field,logFC.thres=logFC.thres,ll=ll,FUN=function(x,ll,field,logFC.thres){
            ll2 <- lapply(ll,field=field,logFC.thres=logFC.thres,sig.thres=x,FUN=function(y,field,logFC.thres,sig.thres){
                row.names(y)[which(y[[field]]<=sig.thres & abs(y$logFC)>logFC.thres)]
            })
            g <- table(unlist(ll2))
            sum(g>=length(ll))
        })
        cols <- c(cols,"black")
    }

    sapply(ll,field=field,sig.thres=sig.thres,logFC.thres=logFC.thres,FUN=function(x,field,sig.thres,logFC.thres){
        sapply(sig.thres,sig=x[[field]][which(abs(x$logFC)>logFC.thres)],FUN=function(y,sig){ sum(sig<=y,na.rm=T)})
    })
    x <- 1:length(sig.thres)
    if(is.null(ylim)) ylim <- range(m)
    plot(x, m[,1], col=cols[1], xaxt="n", xlab=paste(field,"threshold"), ylab="Number of DEGs",ylim=ylim,type=type,lwd=lwd,pch=pch,...)
    axis(1,x,sig.thres)
    if(length(ll)>1){
        for(i in 2:ncol(m)) points(x, m[,i],col=cols[i],type=type,lwd=lwd,pch=pch)
        legend("topright",bty="n",legend=colnames(m), text.col=cols, text.font=2)
    }
}


homogenizeDEAresults <- function(x){
    if(is(x,"list")) return(lapply(x,homogenizeDEAresults))
    x <- as.data.frame(x)
    colnames(x)[which(colnames(x) %in% c("FDR","padj","adj.P.Val"))] <- "FDR"
    colnames(x)[which(colnames(x) %in% c("P.Value","pvalue","PValue"))] <- "PValue"
    colnames(x)[which(colnames(x) %in% c("t","stat"))] <- "stat"
    colnames(x)[which(colnames(x) %in% c("log2FoldChange","logFC"))] <- "logFC"
    for(v in c("logFC","PValue","FDR")) x[[v]] <- dround(x[[v]])
    if(is.null(x$stat)) x$stat <- sign(x$logFC)*-log10(x$PValue)
    return(x[order(x$FDR),])
}

DEP2dea <- function(reg, pw){
  ll <- lapply(c("WBS_vs_CTRL","DUP_vs_CTRL","DUP_vs_WBS"), pw=pw, FUN=function(x,pw){
    row.names(pw) <- pw$name
    pw[,paste0(x,c("_ratio","_p.val","_p.adj"))]
  })
  suffixes <- c(".WBS",".DUP",".DUPvsWBS")
  ll <- lapply(1:3, s=suffixes, ll=ll, FUN=function(x,s,ll){
    a <- ll[[x]]
    colnames(a) <- c("logFC","PValue","FDR")
    a <- homogenizeDEAresults(a)
    colnames(a) <- paste0(colnames(a),s[[x]])
    a
  })
  reg <- homogenizeDEAresults(reg)
  ll <- lapply(ll, rn=row.names(reg), FUN=function(x,rn){ x[rn,] })
  ll <- do.call(cbind, ll)
  cbind(reg,ll)
}

getDEGs <- function(x,logFC.thres=0,FDR.thres=0.05){
    row.names(x)[which(abs(x$logFC)>logFC.thres & x$FDR<=FDR.thres)]
}
getDEGs.any <- function(x,logFC.thres=0,FDR.thres=0.05, max=Inf){
  comps <- sapply(strsplit(grep("FDR",colnames(x),value=TRUE),".",fixed=T),FUN=function(x){
    if(length(x)==1) return("")
    x[[2]]
  })
  ll <- lapply(comps, FUN=function(y){
    if(y=="") return(subset_mDEA(x))
    subset_mDEA(x,y)
  })
  ll <- lapply(ll, FUN=function(x){
    head(row.names(x)[which(abs(x$logFC)>logFC.thres & x$FDR<=FDR.thres)],max)
  })
  unique(unlist(ll))
}



#' toSymbol
#'
#' @param x A DEA data.frame or a vector of ids
#' @param ambiguousUseFirst Logical; whether to use the first gene when a protein
#' group is ambiguous, and there is no non-ambiguous counterpart (default TRUE).
#' This works under the assumption that proteins with several shared peptides are
#' likely to be functionally related.
#'
#' @return A gene-level DEA data.frame.
toSymbol <- function(x, ambiguousUseFirst=TRUE){
  if(is.data.frame(x)){
    # assume dea
    x <- x[order(x$PValue),]
    s <- sapply(strsplit(row.names(x),".",fixed=T), FUN=function(x){ if(length(x)==1) return(x); x[[2]] })
    s <- sapply(strsplit(s,",",fixed=T),us=unique(s), f=ambiguousUseFirst, FUN=function(x, us, f){
      if(length(x)==1) return(x)
      if(f && !all(x %in% us)) return(setdiff(x,us)[1])
      return(NA)
    })
    w <- which(!duplicated(s) & !is.na(s))
    x <- x[w,]
    row.names(x) <- s[w]
    return(x)
  }
  if(is(x,"SummarizedExperiment")){
    library(SummarizedExperiment)
    s <- sapply(strsplit(row.names(x),".",fixed=T), FUN=function(x){
      if(length(x)>1) return(x[[2]]);
      x[[1]]
    })
    a <- lapply(assays(x[!is.na(s),]),s=s[!is.na(s)],FUN=function(x,s){
      ag <- aggregate(x,by=list(gene=s),FUN=sum)
      row.names(ag) <- ag[,1]
      ag[,1] <- NULL
      as.matrix(ag)
    })
    return(SummarizedExperiment(assays = a, colData=colData(x)))
  }
  if(is(x, "character")){
    s <- sapply(strsplit(x,".",fixed=T), FUN=function(x){
      if(length(x)>1) return(x[[2]]);
      NA
    })
    return(unique(s))
  }
  stop("Not implemented!")
}

plWriteDEA <- function(x, file){
    x$Feature <- row.names(x)
    x <- x[,c(ncol(x),1:(ncol(x)-1))]
    write.table(x,file,row.names=F,col.names=T,sep="\t",quote=F)
}

printDEA <- function(x){
    library(DT)
    x <- as.data.frame(cbind(row.names(x),x))
    row.names(x) <- NULL
    colnames(x)[1] <- "Feature"
    for(i in grep("logFC|logCPM",colnames(x))) x[,i] <- round(x[,i],2)
    for(i in grep("PValue|FDR|baseMean",colnames(x))) x[,i] <- format(x[,i],digits=2)
    if("lfcSE" %in% colnames(x)) x[,"lfcSE"] <- round(x[,"lfcSE"],3)
    if("lfcSE" %in% colnames(x)) x[,"lfcSE"] <- round(x[,"lfcSE"],3)
    datatable(x,filter="top",extensions=c("ColReorder","Buttons"), options=list(buttons=c('csv', 'excel', 'pdf', 'print')))
}

CDplot <- function(x, DEA=NULL, addBG=FALSE, cols=1:8, xlim=c(-2,2), xlab="log2(foldchange)", main="", sub=NULL, ...){
    if(!is.null(DEA)){
        x <- lapply(x,y=row.names(DEA),FUN=intersect)
        if(addBG){
            x$background <- unique(setdiff(row.names(DEA),unlist(x)))
            x <- x[c(length(x),1:(length(x)-1))]
        }
        x <- lapply(x,DEA=DEA,FUN=function(y,DEA){
            y <- DEA[unique(y),"logFC"]
            y[!is.na(y)]
        })
    }
    x2 <- lapply(x,FUN=ecdf)
    p <- sapply(x[-1],a=x[[1]],FUN=function(x,a){ suppressWarnings(ks.test(x,a)$p.value) })
    plot(x2[[1]],col=cols[1],xlab=xlab,ylab="Cumulative proportion",xlim=xlim,main=main,sub=sub,...)
    abline(v=0,lty="dashed",col="grey")
    for(i in 2:length(x)) lines(x2[[i]],col=cols[i])
    legend("bottomright",bty="n",lwd=3,col=cols,legend=paste(names(x), c("",paste0(" (p~",sapply(p,digits=2,FUN=format),")"))))
}



#' multintersect
#'
#' Performs pair-wise overlaps between the elements of one or two lists
#'
#' @param ll A list of vectors to be compared
#' @param ll2 An optional second list of vectors
#' @param universe An optional vector of the universe (all terms), used to calculate overlap probabilities and enrichments. If NULL, the union of all lists will be used. Only elements of `ll` and `ll2` that are in the `universe` will be considered.
#' @param keyCol The values to use for the heatmap's colors. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param keyWrite The values to be written in the cells. Possible values are: prob, enrichment, log2Enrichment, log10Prob, overlap, jaccard
#' @param addSetSize Logical; whether to add the set sizes to set names (default TRUE)
#' @param breakNames If not NULL (default), should be an number indicating the character length threshold above which set names should be split onto two lines.
#' @param margin Either a single number indicating the size of the left and bottom margins, or a list (l,r,b,t,pad) indicating the size of each margin.
#' @param returnTables Logical; whether to return tables in addition to the plot (default FALSE)
#'
#' @return Plots a heatmap and return a list of pair-wise overlap metrics
#'
#' @export
multintersect <- function(ll, ll2=NULL, universe=NULL, keyCol="log2Enrichment", keyWrite="overlap", addSetSize=TRUE, breakNames=NULL, margin=100, returnTables=FALSE){
  library(plotly)
  keyCol <- match.arg(keyCol, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  keyWrite <- match.arg(keyWrite, c("prob","enrichment","log2Enrichment","log10Prob","overlap","jaccard"))
  if(is.null(ll2)){
    symm <- TRUE
    ll2 <- ll
  }else{
    symm <- FALSE
  }
  if(is.null(universe)) universe <- unique(c(unlist(ll),unlist(ll2)))
  ll <- lapply(ll,y=universe,FUN=intersect)
  ll2 <- lapply(ll2,y=universe,FUN=intersect)
  ll <- ll[sapply(ll,FUN=length)>0]
  ll2 <- ll2[sapply(ll2,FUN=length)>0]
  if(addSetSize){
    names(ll) <- paste0(names(ll),"\n(",sapply(ll,FUN=length),")")
    names(ll2) <- paste0(names(ll2)," (",sapply(ll2,FUN=length),")")
  }
  if(!is.null(breakNames)){
    names(ll) <- breakStrings(names(ll),breakNames)
    names(ll2) <- breakStrings(names(ll2),breakNames)
  }
  n <- length(ll)
  j <- length(ll2)
  m <- matrix(0,nrow=n,ncol=j)
  colnames(m) <- names(ll2)
  rownames(m) <- names(ll)
  prob <- m
  enr <- m
  jacc <- m
  for(i in 1:n){
    for(j in 1:length(ll2)){
      m[i,j] <- length(intersect(ll[[i]],ll2[[j]]))
      if(names(ll)[i]==names(ll2)[j]){
        prob[i,j] <- NA
        enr[i,j] <- NA
        jacc[i,j] <- NA
      }else{
        enr[i,j] <- getEnrichment(ll[[i]],ll2[[j]],universe)
        prob[i,j] <- overlap.prob(ll[[i]],ll2[[j]],universe,lower=enr[i,j]<1)
        jacc[i,j] <- m[i,j]/length(unique(c(ll[[i]],ll2[[j]])))
      }
    }
  }
  x <- switch(keyCol,
              prob=prob,
              enrichment=enr,
              log2Enrichment=log2(enr),
              log10Prob=-log10(prob),
              overlap=m,
              jaccard=jacc)
  y <- switch(keyWrite,
              prob=format(prob,digits=1),
              enrichment=format(enr,digits=2,trim=T,drop0trailing=T),
              log2Enrichment=round(log2(enr),1),
              log10Prob=round(-log10(prob),1),
              overlap=m,
              jaccard=round(jacc,3))
  y[is.infinite(y)] <- NA
  x[is.infinite(x) & x>0] <- max(x[which(!is.infinite(x))],na.rm=T)
  x[is.infinite(x) & x<0] <- min(x[which(!is.infinite(x))],na.rm=T)
  lab <- matrix(paste0( rep( gsub("\n"," ",row.names(x),fixed=T),ncol(x) ), "\n",
                        rep( gsub("\n"," ",colnames(x), fixed=T),each=nrow(x)), "\n",
                        "overlap: ", as.numeric(m), "\n",
                        as.character(format(enr,digits=2,trim=T,drop0trailing=T)),"-fold ",
                        sapply(enr,FUN=function(x){if(x>=0) return("enrichment"); return("depletion") }), "\n",
                        "p~",as.character(format(prob,digits=2))
  ),nrow=nrow(x),ncol=ncol(x))

  if(symm){
    for(i in 1:ncol(x)) lab[i,i] <- row.names(x)[i]
  }
  if(length(margin)==1){
    margin <- list(l=margin, r=5, b=margin, t=5, pad=4)
  }
  p <- plot_ly(x=colnames(x), y=row.names(x), z=x, type="heatmap", text=lab, hoverinfo = 'text') %>% colorbar(title = keyCol) %>%
    layout(xaxis = list(showgrid = FALSE), yaxis=list(showgrid = FALSE), margin = margin)
  p <- p %>% add_annotations(x = rep(colnames(x),each=nrow(x)), y = rep(row.names(x),ncol(x)), text = y, xref="x", yref="y", showarrow=FALSE)
  if(returnTables) return(list(overlaps=m, probability=prob, enrichment=enr, jaccard=jacc, plot=p))
  p
}

plEA <- function(allGenes, deGenes, sets, fdr.cutoff = 0.2, minSetLength=5){
  sets <- lapply(sets,y=allGenes,FUN=intersect)
  sets <- sets[which(sapply(sets,length)>minSetLength)]
  deGenes <- intersect(deGenes,allGenes)
  d <- data.frame(row.names = names(sets), nbInSet=sapply(sets,length),
                  overlap=as.numeric(sapply(sets,y=deGenes,FUN=function(x,y){ length(intersect(x,y)) })))
  w <- which(d$overlap>0)
  d <- d[w,]
  sets <- sets[w]
  d$enrichment <- sapply(sets, deGenes=deGenes, allGenes=allGenes, FUN = function(x, deGenes, allGenes){
      getEnrichment(deGenes, x, allGenes)
  })
  d$p.value <- sapply(sets, deGenes=deGenes, allGenes=allGenes, FUN = function(x, deGenes, allGenes) {
      overlap.prob(deGenes, x, allGenes)
  })
  d$FDR <- p.adjust(d$p.value, method = "fdr")
  d <- d[order(d$FDR, d$p.value),]
  if (!any(d$FDR < fdr.cutoff)) {
    message("Nothing found!")
  } else {
      d <- d[which(d$FDR < fdr.cutoff),]
      d$genes <- sapply(sets[row.names(d)], y=deGenes, FUN = function(x,y) {
        paste(intersect(x, y), collapse = ", ")
      })
  }
  return(d)
}



#' dround
#'
#' Trim to a certain number of digits (equivalent to `format(...,digits=digits)`, except that the output is numeric)-
#'
#' @param x A vector of numeric values
#' @param digits The number of digits to keep
#' @param roundGreaterThan1 Whether to trim also numbers greater than 1 (default TRUE)
#'
#' @return A numeric vector of the same length as `x`
#' @export
#'
#' @examples
#' dround( c(0.00002345, 554356, 12.56) )
dround <- function(x, digits=3, roundGreaterThan1=TRUE){
  if(is.matrix(x) || is.data.frame(x)){
    for(i in 1:ncol(x)){
      if(is.numeric(x[,i])){
        tryCatch(x[,i] <- dround(x[,i], digits, roundGreaterThan1), error=function(e) warning(e))
      }
    }
    return(x)
  }
  if(roundGreaterThan1){
    w <- 1:length(x)
  }else{
    w <- which(abs(x)<1)
  }
  if(length(w)==0) return(x)
  e <- ceiling(-log10(abs(x[w])))
  x[w] <- round(10^e*x[w],digits-1)/10^e
  x
}

annoColors <- function(){
  list(genotype=c("WBS"="red", "AtWBS"="orange", "CTRL"="grey", "DUP"="blue"))
}


#' TFA
#'
#' @param se An object of class SummarizedExperiment
#' @param dea A DEA table
#' @param design A formula or model matrix
#' @param regulon A regulon object, or the path to a 3-column network format
#' @param pleiotropy Logical; passed to `viper`
#' @param testCoef Coefficient(s) of `design` to test
#' @param assayName name of the assay to use.
#'
#' @return A `SummarizedExperiment` with TF-level activity
#' @export
TFA <- function(se, dea, design, regulon, pleiotropy=TRUE, testCoef=NULL, assayName=NULL){
  suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(viper)
    library(edgeR)
  })
  mm <- switch(class(design),
              formula=model.matrix(design, data=as.data.frame(colData(se))),
              data.frame=design,
              matrix=design,
              character=model.matrix(as.formula(paste0("~",design)), data=as.data.frame(colData(se))),
              stop("Unknown `design`") )
  if(is.null(testCoef)){
    if(is.character(design)){
      testCoef <- design
    }else{
      testCoef <- colnames(mm)[ncol(mm)]
    }
    message("Testing ", testCoef)
  }
  se <- se[rowSums(assays(se)$counts>20)>=3,]
  if(is.null(assayName)) assayName <- intersect(assayNames(se), c("logcpm", "lognorm", "corrected", "imputed"))[1]
  if(is.character(regulon) & length(regulon)==1){
    regulon <- aracne2regulon(regulon, assays(se)[[assayName]], verbose=FALSE)
  }
  regulon <- regulon[intersect(names(regulon), row.names(se))]
  vi1 <- viper(assays(se)[[assayName]], regulon, pleiotropy=pleiotropy, verbose=FALSE)
  vi2 <- SummarizedExperiment( list(logFC=vi1-rowMeans(vi1), viper=vi1),
                                    colData=colData(se) )
  res1 <- topTable(eBayes(lmFit(vi1, mm)), testCoef, Inf)
  vi2 <- vi2[row.names(res1),]
  rd <- dround(res1[,c(1,4,5)])
  colnames(rd) <- c("activity.logFC", "activity.PValue", "activity.FDR")
  rd2 <- dea[row.names(rd),intersect(c("baseMean","logFC","PValue","FDR"),colnames(dea))]
  rd2 <- dround(rd2)
  colnames(rd2) <- paste("expression",colnames(rd2),sep=".")

  tfs <- lapply(regulon, FUN=function(x){ names(x$tfmode) })
  degs <- getDEGs(subset_mDEA(dea))
  bg <- unique(c(intersect(unlist(tfs),row.names(dea)),degs))
  rowData(vi2) <- tryCatch({
    ea <- plEA(bg,degs,tfs,fdr.cutoff=1)
    ea[,3] <- round(ea[,3],2)
    ea <- ea[row.names(rd),c(2:3,5:6)]
    colnames(ea) <- c("nbTargetsInDEGs", "targetEnrichment", "targetEnrFDR", "DEtargets")
    cbind(rd,ea,rd2)
  }, error=function(e){
    cbind(rd,rd2)
  })
  vi2
}

#' hm.targets
#'
#' @param se An object of class `SummarizedExperiment`
#' @param tf A string specifying the transcription factor
#' @param regulon A regulon object.
#' @param TFA An optional SE object produced by the `TFA` function on the same `se`
#' @param subset A vector of gene names (to restrict what is plotted)
#' @param anno_columns see `?sehm`
#' @param anno_rows see `?sehm`
#' @param assayName see `?sehm`
#' @param ... passed to `?sehm`
#'
#' @return plot
#' @export
hm.targets <- function(se, tf, regulon, TFA=NULL, subset=NULL, anno_columns=c("dataset","batch","genotype"), anno_rows=c(), assayName=NULL, ...){
  ar <- data.frame(row.names=names(regulon[[tf]]$tfmode),
                   tfmode=regulon[[tf]]$tfmode,
                   likelihood=regulon[[tf]]$likelihood)
  g <- intersect(row.names(ar), row.names(se))
  if(!is.null(subset)) g <- intersect(g, subset)
  ar <- ar[g,,drop=FALSE]
  if(tf %in% row.names(se)){
    if(is.null(assayName)) assayName <- intersect(assayNames(se), c("logcpm", "lognorm", "corrected", "imputed"))[1]
    se$TF.expression <- as.numeric(assays(se)[[assayName]][tf,])
    anno_columns <- c(anno_columns, "TF.expression")
  }
  if(!is.null(TFA)){
    se$TF.activity <- as.numeric(assay(TFA)[tf,])
    anno_columns <- c(anno_columns, "TF.activity")
  }
  se <- se[g,]
  rowData(se) <- cbind(rowData(se), ar)
  anno_rows <- c(anno_rows, c("tfmode", "likelihood"))
  sehm(se, anno_rows=anno_rows, anno_columns=anno_columns, assayName=assayName, ...)
}


getWordsFromString <- function (ss){
  ss <- strsplit(gsub(" ", ",", ss, fixed = T), ",", fixed = T)[[1]]
  ss[which(ss != "")]
}

breakStrings <- function (x, minSizeForBreak = 20, lb = "\n", allow3Lines=FALSE){
  sapply(x, minSizeForBreak = minSizeForBreak, lb = lb, FUN=function(x, minSizeForBreak, lb) {
    if ((nchar(x) <= minSizeForBreak))
      return(x)
    g <- gregexpr(" ", x)[[1]]
    if (length(g) == 0)
      return(x)
    if (length(g) == 1 & all(g == -1))
      return(x)
    mid <- nchar(x)/2
    mid <- g[order(abs(g - mid))[1]]
    if(allow3Lines && (mid>minSizeForBreak | (nchar(x)-mid) > minSizeForBreak)){
      m1 <- nchar(x)/3
      m1 <- g[order(abs(g - mid))[1]]
      m2 <- (nchar(x)-m1)/2
      m2 <- g[order(abs(g - (nchar(x)-m2)))[1]]
      substr(x, m1, m1) <- lb
      substr(x, m2, m2) <- lb
    }else{
      substr(x, mid, mid) <- lb
    }
    return(x)
  })
}

cdfplot <- function(ll, by=NULL, k=5, ...){
  if(!is.list(ll)){
    if(is.null(by)) stop("If `ll` is not already a list, `by` should be given.")
    if(length(by)!=length(ll)) stop("Lengths of ll and by differ.")
    if(is.logical(by) || is.character(by)) by <- as.factor(by)
    if(is.factor(by)){
      ll <- split(ll, droplevels(by))
    }else{
      q <- quantile(by, (0:(k+1))/(k+1))
      if(!any(duplicated(round(q))))
        q <- unique(c(floor(q[1]),round(q[c(-1,-length(q))]),ceiling(q[length(q)])))
      ll <- split(ll, cut(by, q))
    }
  }
  p <- format(suppressWarnings(ks.test(ll[[1]], rev(ll)[[1]])$p.value), digits=2)
  message("KS p-value between first and last sets:\n", p)
  d <- dplyr::bind_rows(lapply(ll, FUN=function(x){
    x <- x[!is.na(x) & !is.infinite(x)]
    data.frame( y=seq_along(x)/length(x),
                x=sort(x) )
  }), .id="Genesets")
  d$Genesets <- factor(d$Genesets, levels=unique(d$Genesets))
  p <- ggplot(d, aes(x,y,colour=Genesets)) +
    geom_vline(xintercept=0, linetype="dashed") + geom_line(...)
  p + ylab("Cumulative proportion")
}



compileSE <- function(se, dea=NULL, eas=NULL){
  isREST <- FALSE
  if(is.character(se)){
    isREST <- grepl("RNA\\.rest",se)
    if(file.exists(se)){
      se <- readRDS(se)
    }else if(file.exists(paste0(se,".SE.rds"))){
      if(is.null(dea) && file.exists(f <- paste0(se,".DEA.rds")))
        dea <- readRDS(f)
      if(is.null(eas) && file.exists(f <- paste0(se,".DEA.enrichments.rds")))
        eas <- readRDS(f)
      se <- readRDS(paste0(se,".SE.rds"))
    }
  }
  if(!is.null(dea)){
    if(is.character(dea)) dea <- readRDS(dea)
    a <- grep("^FDR|^logFC", colnames(dea), value=TRUE)
    a <- table(gsub("^FDR\\.|^FDR|^logFC\\.|^logFC","",a))
    a <- names(a)[a==2]
    deas <- lapply(setNames(a, paste0("DEA.",a)), FUN=function(n){
      if(n==""){
        x <- dea[,c("logFC","PValue","FDR")]
      }else if(grepl("^[a-zA-Z0-9_]+$",n)){
        x <- dea[,grep(paste0("\\.",n,"$"), colnames(dea))]
        colnames(x) <- gsub(paste0("\\.",n,"$"), "", colnames(x))
      }else{
        x <- dea[,grep(n, colnames(dea), fixed=TRUE)]
        colnames(x) <- gsub(paste0(".",n), "", colnames(x), fixed=TRUE)
      }
      x <- dround(x, roundGreaterThan1=1)
      if("logCPM" %in% colnames(dea)) x$MeanExpr <- dea$logCPM
      if("baseMean" %in% colnames(dea)) x$MeanExpr <- log1p(dea$baseMean)
      x
    })
    names(deas)[names(deas)=="DEA."] <- ifelse(isREST, "DEA.RestInh", "DEA.regCopyNumber")
    for(f in names(deas)) rowData(se)[[f]] <- deas[[f]][row.names(se),]
  }
  if(!is.null(eas)){
    if(is.character(eas)) eas <- readRDS(eas)
    eas <- lapply(eas, FUN=function(x){
      x <- x[unlist(sapply(x, FUN=function(y) is.data.frame(y) || is(y, "DFrame")))]
      lapply(x, FUN=function(y) dround(y, roundGreaterThan1=TRUE))
    })
    metadata(se)$EA <- eas
  }

  if(!is.null(se$genotype))
    se$genotype <- droplevels(factor(as.character(se$genotype),
                              c("WBS","AtWBS","CTRL","CTL","DUP","7Dup")))
  if(is.null(metadata(se)$default_view))
    metadata(se)$default_view <- list(
      assay=head(intersect(c("log2FC","logFC","corrected","logcpm","vst","lognorm"),
                           assayNames(se)), 1),
      groupvar=head(intersect(c("genotype","treatment","condition"), colnames(colData(se))), 1),
      colvar=head(intersect(c("genotype","treatment","condition"), colnames(colData(se))), 1),
      gridvar=head(intersect(c("Dataset","system"), colnames(colData(se))), 1)
    )
  if(is.null(metadata(se)$anno_colors$genotype))
    metadata(se)$anno_colors$genotype <- annoColors()$genotype
  if(isREST) metadata(se)$anno_colors$treatment <-
      c(DMSO="lightgrey", RestInhibitor="darkgreen")
  se
}


stage <- function(what, against, comp=cor){
  i <- intersect(row.names(what), row.names(against))
  what <- what[i,]
  against <- against[i,]
  diff <- apply( against-rowMeans(what), 1, FUN=function(x){
    x[order(abs(x))[1]]
  })
  t(apply(against-diff, 2, FUN=function(x){
    apply(what, 2, y=x, FUN=comp)
  }))
}
