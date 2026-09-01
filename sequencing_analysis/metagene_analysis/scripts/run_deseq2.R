#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 8) stop("counts metadata treated reference output alpha normalized session")
counts_file <- args[1]; metadata_file <- args[2]; treated <- args[3]; reference <- args[4]
output <- args[5]; alpha <- as.numeric(args[6]); normalized_out <- args[7]; session_out <- args[8]
suppressPackageStartupMessages({library(DESeq2); library(apeglm)})
cts <- read.delim(counts_file, check.names=FALSE, row.names=1)
meta <- read.delim(metadata_file, check.names=FALSE, row.names=1)
cts <- round(as.matrix(cts[, rownames(meta), drop=FALSE])); storage.mode(cts) <- "integer"
meta$replicate <- factor(meta$replicate)
meta$condition <- relevel(factor(meta$condition), reference)
keep <- rowSums(cts) >= 10
dds <- DESeqDataSetFromMatrix(cts[keep,,drop=FALSE], meta, design=~replicate+condition)
dds <- tryCatch(
  DESeq(dds, quiet=TRUE),
  error=function(e) {
    if (!grepl("all gene-wise dispersion estimates", conditionMessage(e), fixed=TRUE)) stop(e)
    message("Small/synthetic matrix: using gene-wise dispersions for smoke-test compatibility")
    x <- estimateSizeFactors(dds)
    x <- estimateDispersionsGeneEst(x, quiet=TRUE)
    dispersions(x) <- mcols(x)$dispGeneEst
    nbinomWaldTest(x, quiet=TRUE)
  }
)
res <- results(dds, contrast=c("condition",treated,reference), alpha=alpha)
coef_pattern <- paste0("condition_", make.names(treated), "_vs_", make.names(reference))
coef_name <- grep(coef_pattern, resultsNames(dds), value=TRUE, fixed=TRUE)
if (length(coef_name) != 1) {
  coef_name <- grep(paste0("condition_",treated,"_vs_",reference), resultsNames(dds), value=TRUE, fixed=TRUE)
}
if (length(coef_name) != 1) stop(paste("Could not resolve shrinkage coefficient; names:",paste(resultsNames(dds),collapse=",")))
shr <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)
out <- data.frame(gene_id=rownames(res), baseMean=res$baseMean,
                  log2FoldChange_mle=res$log2FoldChange,
                  log2FoldChange_apeglm=shr$log2FoldChange,
                  log2FoldChange=res$log2FoldChange,
                  log2FoldChange_shrunk=shr$log2FoldChange,
                  lfcSE=res$lfcSE, stat=res$stat, pvalue=res$pvalue, padj=res$padj,
                  independent_filtered=is.na(res$padj) & !is.na(res$pvalue),
                  cooks_or_untestable=is.na(res$pvalue), prefilter_excluded=FALSE,
                  check.names=FALSE)
if (any(!keep)) {
  excluded_norm <- sweep(cts[!keep,,drop=FALSE],2,sizeFactors(dds),"/")
  excluded <- data.frame(gene_id=rownames(cts)[!keep],baseMean=rowMeans(excluded_norm),
                         log2FoldChange_mle=NA_real_,log2FoldChange_apeglm=NA_real_,
                         log2FoldChange=NA_real_,log2FoldChange_shrunk=NA_real_,
                         lfcSE=NA_real_,stat=NA_real_,pvalue=NA_real_,padj=NA_real_,
                         independent_filtered=FALSE,cooks_or_untestable=TRUE,
                         prefilter_excluded=TRUE,check.names=FALSE)
  out <- rbind(out,excluded)
}
out <- out[match(rownames(cts),out$gene_id),,drop=FALSE]
write.table(out, output, sep="\t", quote=FALSE, row.names=FALSE)
norm_counts <- sweep(cts,2,sizeFactors(dds),"/")
norm <- data.frame(gene_id=rownames(cts), norm_counts, check.names=FALSE)
write.table(norm, normalized_out, sep="\t", quote=FALSE, row.names=FALSE)
sink(session_out); print(sessionInfo()); sink()
