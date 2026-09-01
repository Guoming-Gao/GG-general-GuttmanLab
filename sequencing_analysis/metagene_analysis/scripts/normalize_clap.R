#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) stop("matrix metadata context pseudocpm output")
suppressPackageStartupMessages(library(edgeR))
d <- read.delim(args[1], check.names=FALSE)
m <- read.delim(args[2], check.names=FALSE)
context <- args[3]; pseudo <- as.numeric(args[4]); output <- args[5]
m <- m[m$context == context,,drop=FALSE]
id_cols <- intersect(c("gene_id","gene_name","chrom","strand","gene_length","intron_count","bin","feature"), names(d))
cts <- as.matrix(d[,m$sample,drop=FALSE]); storage.mode(cts) <- "numeric"
y <- DGEList(cts); y <- calcNormFactors(y, method="TMM")
cp <- cpm(y, normalized.lib.sizes=TRUE, log=FALSE)
reps <- unique(m$replicate); enr <- matrix(NA_real_,nrow(cp),length(reps))
for (i in seq_along(reps)) {
  cl <- m$sample[m$replicate==reps[i] & m$assay=="CLAP"]
  inp <- m$sample[m$replicate==reps[i] & m$assay=="Input"]
  if(length(cl)!=1 || length(inp)!=1) stop("Each replicate needs one CLAP and one Input")
  enr[,i] <- log2((cp[,cl]+pseudo)/(cp[,inp]+pseudo))
}
design_meta <- m; design_meta$replicate <- factor(design_meta$replicate); design_meta$assay <- relevel(factor(design_meta$assay),"Input")
design <- model.matrix(~replicate+assay,design_meta)
ql <- tryCatch({
  yy <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(yy,design,robust=TRUE)
  coef <- grep("^assayCLAP$",colnames(design))
  topTags(glmQLFTest(fit,coef=coef),n=Inf,sort.by="none")$table
}, error=function(e) {
  warning(paste("edgeR paired test failed:",conditionMessage(e)))
  data.frame(logFC=rep(NA_real_,nrow(d)),PValue=NA_real_,FDR=NA_real_)
})
out <- cbind(d[,id_cols,drop=FALSE], clap_cpm=rowMeans(cp[,m$sample[m$assay=="CLAP"],drop=FALSE]),
             input_cpm=rowMeans(cp[,m$sample[m$assay=="Input"],drop=FALSE]),
             enrichment=rowMeans(enr), enrichment_repA=enr[,1], enrichment_repB=enr[,2],
             edgeR_logFC=ql$logFC, edgeR_PValue=ql$PValue, edgeR_FDR=ql$FDR)
write.table(out,output,sep="\t",quote=FALSE,row.names=FALSE)
sf <- data.frame(sample=m$sample,lib_size=y$samples$lib.size,norm_factor=y$samples$norm.factors,
                 effective_lib_size=y$samples$lib.size*y$samples$norm.factors)
write.table(sf,sub("\\.tsv$",".size_factors.tsv",output),sep="\t",quote=FALSE,row.names=FALSE)
