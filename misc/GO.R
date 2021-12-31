#' goseq.enrichment
#'
#' Gets GO terms enriched in a given set of genes using the `goseq` package, and thereby
#' correcting for the transcript length bias of RNA-seq differential expression data.
#'
#' @param allGenes A character vector containing all the gene symbols in the background (i.e. all tested genes). Set it to NA to use all genes.
#' @param deGenes A character vector containing the dysregulated genes.
#' @param gotype A character vector indicating the GO ontologies to fetch. Can be any combination of 'GO:BP', 'GO:MF', 'GO:CC'.
#' @param cutoff A number between 0 and 1, indicating the FDR threshold below which terms will be retained. Default 0.1.
#' @param org The organism, e.g. "hg19" (human) or "mm9" (mouse), although other genomes versions should work too.
#' @param minCatSize The minimum size of categories (nb of genes with the annotation) to be considered. Default 10.
#' @param maxCatSize The maximum size of categories (nb of genes with the annotation) to be considered. Default 1000.
#' @param maxResults The maximum number of categories to return. Default 200.
#'
#' @return A data.frame of class 'GOresults'
#'
#' @export
goseq.enrichment <- function(allGenes, deGenes, gotype=c("GO:BP","GO:MF","GO:CC"), cutoff=0.1, org="hg19", minCatSize=10, maxCatSize=1000, maxResults=200){
  suppressPackageStartupMessages(library("goseq"))
  suppressPackageStartupMessages(library("GO.db"))
  suppressPackageStartupMessages(library("org.Hs.eg.db"))
  suppressPackageStartupMessages(library("AnnotationDbi"))
  gotype <- match.arg(gotype, c("GO:CC","GO:MF","GO:BP"), several.ok=T)
  org <- match.arg(org, c("hg19","mm9"))
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
  go.en <- go.all[go.all$FDR < cutoff,]
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
  names(go.en) <- c("GO.ID","Term","Annotated","Significant","Enrichment","Fisher.pvalue","FDR")
  class(go.en) <- c("data.frame","GOresults")
  if(nrow(go.en)>0){
    if(nrow(go.en)>maxResults) go.en <- go.en[1:maxResults,]
    species <- ifelse(substr(org,0,2)=="hg","Hs","Mm")
    db <- paste0('org.',species,'.eg')
    library(package=paste0(db,'.db'), character.only = T)
    eg <- mget(as.character(go.en$GO.ID), get(paste0(db,'GO2ALLEGS')), ifnotfound=NA)
    go.en$genes <- sapply(eg, db=db, FUN=function(x,db){
      x <- x[which(!is.na(x))]
      if(length(x)==0) return("")
      paste(sort(unique(as.character(unlist(mget(as.character(x), get(paste0(db,'SYMBOL'))))))),collapse=", ")
    })
    return(go.en)
  }
  message("Nothing found!")
}


#' camPlot
#'
#' Plot camera results
#'
#' @param cam Camera results
#' @param dea Optional DEA table
#' @param n Max number of top categories to plot
#' @param cat.fdr.thres Significance threshold for categories
#' @param breakAt Max label length for breaks
#' @param genes.fdr.thres Significance threshold for genes
#' @param ... passed to `plot_grid`
#'
#' @return a `ggplot`
camPlot <- function(cam, dea=NULL, n=15, cat.fdr.thres=0.1, breakAt=30, allow3Lines=TRUE, genes.fdr.thres=0.25, labsize=3.1, ...){
  cam <- cam[cam$collection=="C5" & cam$subcat!="HPO" & cam$FDR<cat.fdr.thres,]
  cam$y <- sign(2*(cam$Direction=="Down")-1)*log10(cam$FDR)
  cam$term <- gsub("GOMF|GOBP|GOCC","GO",cam$term)
  cam$term2 <- breakStrings(gsub("^GO ","",gsub("_"," ",cam$term)), breakAt, allow3Lines=allow3Lines)
  cam <- head(cam,n)
  cam <- cam[order(cam$FDR),]
  cam$term2 <- factor(cam$term2, cam$term2)
  if(is.null(dea)){
    return(ggplot(cam, aes(reorder(term2,y), y, fill=Direction)) + geom_col() +
             labs(y="Direction*-log10(FDR)", x="") + coord_flip() +
             geom_label(aes(y=0, label=nDE), size=4, fontface="bold"))
  }
  cats <- setNames(cam$term, cam$term2)
  d <- dplyr::bind_rows(lapply(cats, FUN=function(x){
    dea[intersect(gsets[[paste0("C5:GO:",x)]],row.names(dea)),]
  }), .id="set")
  y <- -log10(cam$FDR)
  d$set <- factor(d$set, cam$term2)
  d$Direction <- cam$Direction[as.integer(d$set)]

  # p1 <- ggplot(cam, aes(term2, FDR)) + geom_col(width=0.7) + labs(x="") +
  #   scale_y_continuous(trans = ggforce::trans_reverser('log10')) +
  #   coord_flip() + theme(axis.text.y=element_text(hjust=0.5, vjust=0.5)) +
  #   theme( axis.line.y=element_blank(), axis.ticks.y=element_blank())
  p1 <- ggplot(cam, aes(term2, FDR)) + geom_col(width=0.9, fill="grey") +
    geom_text(aes(y=1,label=term2), size=labsize, hjust=0, nudge_y=0.05, lineheight=0.85) +
    labs(x="") + scale_y_continuous(trans = ggforce::trans_reverser('log10')) +
    coord_flip() + theme(axis.text.y=element_text(hjust=0.5, vjust=0.5)) +
    theme( axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.text.y=element_blank())
  ds <- d[d$FDR<genes.fdr.thres,]
  p2 <- ggplot(d, aes(y=set, x=logFC)) +
    geom_violin(aes(fill=Direction, colour=Direction), scale="width", trim=TRUE) +
    geom_vline(xintercept=0, linetype="dashed") + xlim(range(ds$logFC)) +
    geom_point(data=ds, size=3, shape="I") +
    theme( axis.text.y=element_blank(), axis.line.y=element_blank(),
           axis.ticks.y=element_blank(), axis.title.y=element_blank(),
           legend.position="none")
  plot_grid(p2, p1, ...)
}

topgo <- function(all, sig, ontology="BP", mapping="org.Hs.eg.db"){
  library(topGO)
  x2 <- setNames(factor(as.integer(all %in% sig)), all)
  god <- new("topGOdata",  annot=annFUN.org, mapping=mapping, nodeSize=10,
             ontology = ontology, allGenes=x2, ID="symbol")
  res <- runTest(god, statistic = "fisher")
  res <- GenTable(god, fisher=res, topNodes = 200)
  row.names(res) <- res[,1]
  res$enrichment <- res$Significant/res$Expected
  res
}
