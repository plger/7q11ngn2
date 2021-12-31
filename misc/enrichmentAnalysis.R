#' runAllEA
#'
#' @param m A merged DEA table, as produced by `cbindDEAs`
#' @param methods The EA methods to use
#' @param addSets An eventual list of additional sets to include.
#' @param dea.thres The FDR threshold for differential expression, default 0.05
#'
#' @return A list.
runAllEA <- function(m, methods=c("camera","goseq"), addSets=NULL, dea.thres=0.05){
  suf <- c("reg"=NA,"WBS"="WBS","DUP"="DUP","DUPvsWBS"="DUPvsWBS")
  res <- suppressWarnings(lapply(suf, m=m, methods=methods, FUN=function(x,m,methods){
    if(is.na(x)){
      m2 <- m
    }else{
      m2 <- try(subset_mDEA(m, x),silent=T)
      if(is(m2,"try-error")) return(NULL)
    }
    m2 <- toSymbol(m2)
    ll <- list()
    if("camera" %in% methods) ll[["camera"]] <- cameraWrapper(m2, dea.thres=dea.thres, addSets=addSets)
    # if("gsea" %in% methods) ll[["gsea"]] <-
    if("goseq" %in% methods) ll[["goseq"]] <- goseq.enrichment(row.names(m2), row.names(m2)[which(m2$FDR<dea.thres)], gotype = "GO:BP", cutoff = 0.01, cutoff.onFDR = F)
    return(ll)
  }))
  res
}

#' goseq.enrichment
#'
#' Gets GO terms enriched in a given set of genes using the `goseq` package, and thereby
#' correcting for the transcript length bias of RNA-seq differential expression data.
#'
#' @param allGenes A character vector containing all the gene symbols in the background (i.e. all tested genes). Set it to NA to use all genes.
#' @param deGenes A character vector containing the dysregulated genes.
#' @param gotype A character vector indicating the GO ontologies to fetch. Can be any combination of 'GO:BP', 'GO:MF', 'GO:CC'.
#' @param cutoff A number between 0 and 1, indicating the FDR threshold below which terms will be retained. Default 0.1.
#' @param cutoff.onFDR Whether to apply the cutoff on FDR (default), otherwise it will be applied on PValue.
#' @param org The organism, e.g. "hg19" (human) or "mm9" (mouse), although other genomes versions should work too.
#' @param minCatSize The minimum size of categories (nb of genes with the annotation) to be considered. Default 10.
#' @param maxCatSize The maximum size of categories (nb of genes with the annotation) to be considered. Default 1000.
#' @param maxResults The maximum number of categories to return. Default 200.
#'
#' @return A data.frame.
#'
#' @export
goseq.enrichment <- function(allGenes, deGenes, gotype=c("GO:BP","GO:MF","GO:CC"), cutoff=0.1, cutoff.onFDR=TRUE, org="hg19", minCatSize=10, maxCatSize=1000, maxResults=200){
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("GO.db"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("AnnotationDbi"))
  gotype <- match.arg(gotype, c("GO:CC","GO:MF","GO:BP"), several.ok=T)
  org <- match.arg(org, c("hg19","mm9"))
  emptyRes <- data.frame( term=vector(mode = "character", length=0),
                          enrichment=vector(mode = "numeric", length=0),
                          PValue=vector(mode = "numeric", length=0),
                          FDR=vector(mode = "numeric", length=0),
                          genes=vector(mode = "character", length=0) )
  if(length(deGenes)<3) return(emptyRes)
  if(all(is.na(allGenes))){
    message("Using all (annotated) genes as a background.")
    if(org=="hg19"){
        data("hg38.refseq")
    }else{
        data("mm10.refseq")
    }
    allGenes <- unique(as.character(gtf$gene_id))
  }
  genes <- as.integer(allGenes %in% as.character(deGenes))
  names(genes) <- allGenes
  pwf=nullp(genes,org,"geneSymbol",plot.fit=F)
  go.all <- suppressMessages(goseq(pwf,org,"geneSymbol",test.cats=gotype))
  go.all <- go.all[which(go.all$numInCat >= minCatSize & go.all$numInCat < maxCatSize),]
  go.all$FDR <- p.adjust(go.all$over_represented_pvalue,method="BH")
  if(cutoff.onFDR){
    go.en <- go.all[go.all$FDR < cutoff,]
  }else{
    go.en <- go.all[go.all$over_represented_pvalue < cutoff,]
  }
  go.en <- go.en[order(go.en$over_represented_pvalue),]
  go.en$Term <- sapply(go.en$category,FUN=function(x){
        if(is.null(x) | is.na(x) | x == ""){
                return("")
        }else{
                return(Term(x))
        }
    })
  go.en$Expected <- length(deGenes)*(go.en$numInCat/length(allGenes))
  go.en$Enrichment <- go.en$numDEInCat/go.en$Expected
  go.en <- go.en[,c("category","Term","numInCat","numDEInCat","Enrichment","over_represented_pvalue","FDR")]
  names(go.en) <- c("GO.ID","term","Annotated","Significant","enrichment","PValue","FDR")

  if(!(nrow(go.en)>=0)) return(emptyRes)

  if(nrow(go.en)>maxResults) go.en <- go.en[1:maxResults,]
  species <- ifelse(substr(org,0,2)=="hg","Hs","Mm")
  db <- paste0('org.',species,'.eg')
  library(package=paste0(db,'.db'), character.only = T)
  eg <- AnnotationDbi::mget(as.character(go.en$GO.ID), get(paste0(db,'GO2ALLEGS')), ifnotfound=NA)
  go.en$genes <- sapply(eg, db=db, sig=as.character(deGenes), FUN=function(x,db,sig){
    x <- x[which(!is.na(x))]
    if(length(x)==0) return("")
    x <- unique(as.character(unlist(AnnotationDbi::mget(as.character(x), get(paste0(db,'SYMBOL'))))))
    x <- intersect(x,sig)
    paste(sort(x),collapse=", ")
  })

  return(go.en)
}



getMsigSets <- function(coll=c("H","C2","C5")){
  library(msigdbr)
  m <- msigdbr(species="Homo sapiens")
  m <- m[which(m$gs_cat %in% coll),]
  m$name2 <- paste0(m$gs_cat,":",m$gs_subcat, ":", m$gs_name)
  split(m$gene_symbol, m$name2)
}
cameraWrapper <- function(dea, gsets=NULL, addSets=NULL, addDEgenes=TRUE, dea.thres=0.05){
  library(limma)
  if(is.null(gsets)){
    gsets <- getMsigSets()
  }
  if(!is.null(addSets)) addSets <- c(gsets, addSets)
  dea <- homogenizeDEAresults(dea)
  dea <- dea[which(!is.na(dea$PValue)),]
  gsets <- lapply(gsets, y=row.names(dea), FUN=intersect)
  gsets <- gsets[which(sapply(gsets,length)>4)]
  if("stat" %in% colnames(dea)){
    gs <- dea$stat
  }else{
    gs <- sign(dea$logFC)*-log10(dea$PValue)
  }
  names(gs) <- row.names(dea)
  gs <- gs[which(!is.na(gs))]
  Cres <- cameraPR(gs, gsets)
  Cres <- Cres[which(Cres$PValue < 0.01),,drop=F]
  if(nrow(Cres)==0) return(NULL)
  dg <- t(sapply(gsets[row.names(Cres)], sig=row.names(dea)[which(dea$FDR<dea.thres)], FUN=function(x,sig){
    c(sum(x %in% sig), paste(intersect(x,sig),collapse=", "))
  }))
  Cres <- cbind(enrichment=rep(NA_real_, nrow(dg)), nDE=as.numeric(dg[,1]), Cres, genes=dg[,2])
  df <- t(sapply(strsplit(row.names(Cres),":",fixed=T), FUN=function(x){
    if(length(x)==2) return(c("custom",x[[2]],paste(x[3:length(x)],collapse=":")))
    if(length(x)==1) return(c("custom",NA,paste(x[3:length(x)],collapse=":")))
    c(x[[1]],x[[2]],paste(x[3:length(x)],collapse=":"))
  }))
  colnames(df) <- c("collection","subcat","term")
  Cres <- as.data.frame(cbind(df,Cres),stringsAsFactors=F)
  row.names(Cres) <- NULL
  Cres
}


#' goDesc
#'
#' Describes one or multiple sets of genes in terms of informative (though not
#' necessarily significantly enriched) GO terms
#'
#' @param allgenes
#' @param genesInSet
#' @param species
#' @param gotype
#' @param maxStartSize
#' @param minSize
#' @param k Number of terms to return
goDesc <- function(allgenes, genesets, species="hg19", gotype=c("GO:CC", "GO:MF", "GO:BP"), maxStartSize=1500, minSize=4, k=3){
  library(GO.db)
  library(goseq)
  gotype <- match.arg(gotype)
  ll <- getgo(allgenes, species, "geneSymbol", gotype)
  catsize <- table(unlist(ll))
  catsize <- catsize[which(catsize<maxStartSize)]
  if(!is.list(genesets)) genesets <- list(genesets)
  res <- lapply(genesets, FUN=function(x){
    m <- .goDescInt(x, ll, catsize, minSize=minSize, k=k)
    if(nrow(m)==0) return(list(summary=NULL, genes=NULL, overlaps=NULL))
    row.names(m) <- Term(row.names(m))
    list( summary=paste0(row.names(m), " (", rowSums(m),"/",length(x),")"),
          genes=apply(m,1,FUN=function(x) colnames(m)[x]),
          overlaps=apply(m,1,FUN=function(x) apply(m,1,FUN=function(y) sum(x & y))) )
  })
  if(length(res)==1 && is.null(names(res))) res <- res[[1]]
  res
}


#' goDesc2
#'
#' Describes multiple sets of genes in terms of informative (though not
#' necessarily significantly enriched) GO terms (one global grid for all sets).
#'
#' @param allgenes
#' @param genesInSet
#' @param species
#' @param gotype
#' @param maxStartSize
#' @param minSize
#' @param k Number of terms to return
goDesc2 <- function(allgenes, genesets, species="hg19", gotype=c("GO:CC", "GO:MF", "GO:BP"), maxStartSize=1500, minSize=4, k=3){
  library(GO.db)
  library(goseq)
  gotype <- match.arg(gotype)
  ll <- getgo(allgenes, species, "geneSymbol", gotype)
  catsize <- table(unlist(ll))
  catsize <- catsize[which(catsize<maxStartSize)]
  ll <- lapply(ll, y=names(catsize), FUN=intersect)
  if(!is.list(genesets)) genesets <- list(genesets)
  mat <- lapply(genesets, ll=ll, catsize=catsize, minSize=minSize, k=k, FUN=.goDescInt)
  cats <- unique(unlist(lapply(mat, row.names)))
  g <- unique(unlist(genesets))
  ll <- lapply(ll[g], y=cats, FUN=intersect)
  catsize <- catsize[intersect(names(catsize), cats)]
  m <- matrix(FALSE, nrow = length(names(catsize)), ncol=length(ll))
  rownames(m) <- names(catsize)
  colnames(m) <- names(ll)
  for (i in 1:length(ll)) {
    m[, i] <- rownames(m) %in% ll[[i]]
  }
  m <- m[order(catsize[row.names(m)]),,drop=FALSE]
  m <- m[!duplicated(m),,drop=FALSE]
  if(nrow(m)==0) return(list(summary=NULL, genes=NULL, overlaps=NULL))
  row.names(m) <- Term(row.names(m))
  m <- as.data.frame(t(m[,!is.na(colnames(m))]))
  m2 <- t(sapply(genesets, FUN=function(x) colSums(m[x,], na.rm=TRUE)))
  list( perGene=m, perSet=m2 )
}

.goDescInt <- function(genes, ll, catsize, minSize=4, k=3){
  ll <- ll[intersect(names(ll),genes)]
  m <- matrix(FALSE, nrow = length(names(catsize)), ncol=length(ll))
  rownames(m) <- names(catsize)
  colnames(m) <- names(ll)
  print(length(ll))
  for (i in 1:length(ll)) {
    m[, i] <- rownames(m) %in% ll[[i]]
  }
  m <- m[order(catsize[row.names(m)]),,drop=FALSE]
  m <- m[!duplicated(m),,drop=FALSE]
  m <- m[rowSums(m)>=minSize,,drop=FALSE]
  if(nrow(m)>k){
    tmp <- cutree(hclust(dist(m, method="manhattan")),k = k)
    tmp <- tmp[!is.na(names(tmp))]
    tmp2 <- sapply(split(names(tmp),tmp), FUN=function(x){ x[order(rowSums(m[x,,drop=FALSE]), decreasing=TRUE)[1]] })
    m <- m[tmp2,]
  }
  m
}

clusterGO <- function(go, k=2:min(20,nrow(go)-1), plot=TRUE, fontsize=4){
  library(ggplot2)
  library(ggrepel)
  go <- as.data.frame(go)
  if(!("FDR" %in% colnames(go))) go$FDR <- go$qvalue
  if(!("enrichment" %in% colnames(go))){
    f <- function(x){ x <- as.numeric(x); x[1]/x[2] }
    go$enrichment <- sapply(strsplit(test2$GeneRatio, split = "/"), FUN=f)/
      sapply(strsplit(test2$BgRatio, split = "/"), FUN=f)
  }
  if(!("Significant" %in% colnames(go))) go$Significant <- go$Count
  if(!("term" %in% colnames(go))) go$term <- go$Description
  if(!("genes" %in% colnames(go))) go$genes <- go$geneID
  bm <- lapply(gsub("/"," ",(go$genes)), getWordsFromString)
  g <- unique(unlist(bm))
  bm <- sapply(bm, FUN=function(x) g %in% x)
  d <- dist(t(bm),method = "binary")
  cc <- lapply(k, FUN=function(k) kmeans(d,k, nstart=3))
  ve <- sapply(cc, FUN=function(x) x$betweenss/x$totss)
  k <- k[farthestPoint(ve)]
  cc <- kmeans(d, k, nstart = 3)$cluster
  go$cluster <- factor(cc, levels=unique(cc))
  #go2 <- go[order(go$cluster, go$FDR),]
  #levels(go$cluster) <- go2$term[!duplicated(go2$cluster)]
  if(!plot) return(go)
  go$Term <- breakStrings(go$term)
  go2 <- go[!duplicated(go$cluster),]
  ggplot(go, aes(log2(enrichment),-log10(FDR), colour=cluster)) +
    geom_point(aes(size=Significant), alpha=0.65) + geom_point(data=go2) +
    geom_label_repel(data=go2, aes(label=Term), fontface="bold", size=fontsize,
                     fill="#FFFFFFC8", min.segment.length=0) +
    scale_color_discrete(guide=FALSE)
}


getWordsFromString <- function(ss){
  for (i in c(" ", "\n", "\r", ";")) ss <- gsub(i, ",", ss, 
                                                fixed = T)
  ss <- strsplit(ss, ",", fixed = T)[[1]]
  ss[which(ss != "")]
}

#' farthestPoint
#'
#' Identifies the point farthest from a line passing through by the first and
#' last points. Used for automatization of the elbow method.
#'
#' @param y Monotonically inscreasing or decreasing values
#' @param x Optional x coordinates corresponding to `y` (defaults to seq)
#'
#' @return The value of `x` farthest from the diagonal.
#'
#' @examples
#' y <- 2^(10:1)
#' plot(y)
#' x <- farthestPoint(y)
#' points(x,y[x],pch=16)
farthestPoint <- function(y, x=NULL){
  if(is.null(x)) x <- seq_len(length(y))
  d <- apply( cbind(x,y), 1,
              a=c(1,y[1]), b=c(length(y),rev(y)[1]),
              FUN=function(y, a, b){
                v1 <- a-b
                v2 <- y-a
                abs(det(cbind(v1,v2)))/sqrt(sum(v1*v1))
              })
  order(d,decreasing=TRUE)[1]
}


CDplot2 <- function(ll, by=NULL, k=5, breaks=NULL, sameFreq=FALSE, addN=FALSE, dig.lab=NULL, ...){
  library(ggplot2)
  if(!is.list(ll)){
    if(is.null(by)) stop("If `ll` is not already a list, `by` should be given.")
    if(length(by)!=length(ll)) stop("Lengths of ll and by differ.")
    w <- which(!is.na(by) & !is.na(ll))
    by <- by[w]
    ll <- ll[w]
    if(is.null(dig.lab)) dig.lab <- max(c(2,3-ceiling(log10(abs(mean(by))))))
    if(is.null(breaks)) breaks <- k
    if(sameFreq) breaks <- quantile(by, prob=seq(from=0, to=1, length.out=k+1), na.rm=TRUE)
    ll <- split(ll, cut(by, breaks, dig.lab=dig.lab))
  }else{
    ll <- lapply(ll, FUN=function(x) x[!is.na(x)])
  }
  p <- format(suppressWarnings(ks.test(ll[[1]], rev(ll)[[1]])$p.value), digits=2)
  message("KS p-value between first and last sets:\n", p)
  d <- dplyr::bind_rows(lapply(ll, FUN=function(x){
    data.frame( y=(seq_along(x)-1)/(length(x)-1),
                x=sort(x) )
  }), .id="Genesets")
  d$Genesets <- factor(d$Genesets, levels=unique(d$Genesets))
  if(addN) levels(d$Genesets) <- paste0(levels(d$Genesets), " (n=",as.numeric(table(d$Genesets)),")")
  p <- ggplot(d, aes(x,y,colour=Genesets)) +
    geom_vline(xintercept=0, linetype="dashed") + geom_line(...)
  p + ylab("Cumulative proportion")
}
