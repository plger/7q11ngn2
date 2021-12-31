fn <- function(suffix=".DUP", tissue="Neurons", type="isogenic", ylim=c(-2,2)){

  rpf <- readRDS(paste0("clean/",tissue,"/RPF/",tissue,".RPF.",type,".DEA.rds"))
  prot <- readRDS(paste0("clean/",tissue,"/Protein/",tissue,".Protein.",type,".DEA.rds"))
  te <- readRDS(paste0("clean/",tissue,"/TE/",tissue,".TE.",type,".DEA.rds"))

  rna <- readRDS(paste0("clean/",tissue,"/RNA/",tissue,".RNA.",type,".SE.rds"))
  tpm <- log(rowMeans(assays(toSymbol(rna))$tpm)+0.1)
  gl0 <- plag(as.numeric(as.character(rowData(rna)[,5])), rowData(rna)[,3], mean)
  gl <- log1p(gl0$x)
  names(gl) <- row.names(gl0)
  
  fcf <- paste0("logFC", suffix)
  ylab <- ifelse(suffix=="", "log2FC with copynumber", paste("log2FC", gsub("\\.","",suffix), "vs CTRL"))
    
  layout(matrix(1:4,nrow=2))
  
  LSD::heatscatter(prot$AveExpr, prot[[fcf]], ylab=ylab, xlab="Average log-intensity", main="Protein", ylim=ylim)
  ss <- smooth.spline(prot$AveExpr,prot[[fcf]])
  lines(ss, lwd=2)
  
  te <- te[which(row.names(te) %in% names(tpm) & !is.na(te[[fcf]])),]
  LSD::heatscatter(tpm[row.names(te)], te[[fcf]], ylab=ylab, xlab="log(TPM+0.1)", main="Translation efficiency")
  abline(lm(te[[fcf]]~tpm[row.names(te)]), lwd=2)
  
  rpf <- rpf[intersect(row.names(rpf), names(tpm)),]
  LSD::heatscatter(tpm[row.names(rpf)], rpf[[fcf]], ylab=ylab, xlab="log(TPM+0.1)", main="RPF", ylim=ylim)
  abline(lm(rpf[[fcf]]~tpm[row.names(rpf)]), lwd=2)
  
  # LSD::heatscatter(gl[row.names(rpf)], rpf[[fcf]], ylab=ylab, xlab="Effective length", main="RPF", ylim=ylim)
  # ss <- smooth.spline(tpm[row.names(rpf)], rpf[[fcf]])
  # lines(ss, lwd=2)
  
  LSD::heatscatter(gl[row.names(te)], te[[fcf]], ylab=ylab, xlab="Effective length", main="Translation efficiency",ylim=ylim)
  abline(lm(te[[fcf]]~gl[row.names(te)]), lwd=2)

}