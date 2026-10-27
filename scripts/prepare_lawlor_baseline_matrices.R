suppressPackageStartupMessages({
  library(Matrix)
})
audit <- "."
ann <- read.csv(file.path(audit, "CZI.PBMC.cell.annotations.csv"), check.names=FALSE, stringsAsFactors=FALSE)
colnames(ann)[1] <- "barcode"
keep <- ann$HTO_Classification == "Baseline" & ann$Demuxlet_Classification == "SNG"
ann_k <- ann[keep, , drop=FALSE]
message("Baseline SNG cells: ", nrow(ann_k))
write.csv(ann_k, file.path(audit, "baseline_sng_annotations.csv"), row.names=FALSE)

# RNA
message("Loading RNA RDS...")
rna <- readRDS(file.path(audit, "CZI.PBMC.RNA.matrix.Rds"))
message("RNA class=", paste(class(rna), collapse=","), " dim=", paste(dim(rna), collapse="x"))
# subset columns
common <- intersect(colnames(rna), ann_k$barcode)
message("RNA barcodes overlapping annotations: ", length(common))
if (length(common) < 1000) stop("Too few overlapping barcodes")
rna_sub <- rna[, common, drop=FALSE]
# ensure dgCMatrix
if (!inherits(rna_sub, "dgCMatrix")) rna_sub <- as(rna_sub, "dgCMatrix")
dir.create(file.path(audit, "baseline_mm"), showWarnings=FALSE)
# write Matrix Market
writeMM(rna_sub, file.path(audit, "baseline_mm", "rna.mtx"))
writeLines(rownames(rna_sub), file.path(audit, "baseline_mm", "rna_genes.txt"))
writeLines(colnames(rna_sub), file.path(audit, "baseline_mm", "rna_barcodes.txt"))
message("RNA written nnz=", length(rna_sub@x))
rm(rna, rna_sub); gc()

# ADT
message("Loading ADT RDS...")
adt <- readRDS(file.path(audit, "CZI.PBMC.ADT.matrix.Rds"))
# ADT is data.frame features x cells
adt_sub <- adt[, common, drop=FALSE]
# drop QC rows
drop <- c("bad_struct","no_match","total_reads")
rn <- rownames(adt_sub)
keep_feat <- !(rn %in% drop) & !startsWith(rn, "control_")
adt_mat <- as(as.matrix(adt_sub[keep_feat, , drop=FALSE]), "dgCMatrix")
writeMM(adt_mat, file.path(audit, "baseline_mm", "adt.mtx"))
writeLines(rownames(adt_mat), file.path(audit, "baseline_mm", "adt_features.txt"))
writeLines(colnames(adt_mat), file.path(audit, "baseline_mm", "adt_barcodes.txt"))
message("ADT written features=", nrow(adt_mat), " cells=", ncol(adt_mat))
# reorder annotations to common order
ann_ord <- ann_k[match(common, ann_k$barcode), , drop=FALSE]
write.csv(ann_ord, file.path(audit, "baseline_mm", "annotations_ordered.csv"), row.names=FALSE)
message("DONE")
