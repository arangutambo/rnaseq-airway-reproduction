# Compare the 316 genes in Himes et al. 2014 (Table S3, Cuffdiff) with this repository's DESeq2 result.
# Table S3 is from a CC BY paper: https://doi.org/10.1371/journal.pone.0099625 (supplementary file s014).
suppressPackageStartupMessages({library(DESeq2); library(airway); library(data.table)})
paper <- fread("data/himes2014_tableS3_316genes.csv")
setnames(paper, c("gene", "dex_fpkm", "untr_fpkm", "ln_fc", "stat", "p", "q", "locus"))
paper[, paper_dir := ifelse(as.numeric(dex_fpkm) > as.numeric(untr_fpkm), 1, -1)]  # + = up with dexamethasone
# The column is headed "Ln[Fold Change]" but its values equal log2(untreated FPKM / dex FPKM), not the natural log
# (for C7: log2(3.76 / 38.41) = -3.35). Check that on every row instead of assuming it:
paper[, log2_from_fpkm := log2(as.numeric(untr_fpkm) / as.numeric(dex_fpkm))]
col_is_log2 <- mean(abs(paper$log2_from_fpkm - as.numeric(paper$ln_fc)) < 0.02, na.rm = TRUE)
cat(sprintf("Table S3 'Ln' column matches log2(untreated/dex FPKM) within 0.02 on %.0f%% of rows\n", 100 * col_is_log2))

data(airway); dds <- DESeqDataSet(airway, ~ cell + dex); dds$dex <- relevel(dds$dex, "untrt")
dds <- dds[rowSums(counts(dds) >= 10) >= 4, ]; dds <- DESeq(dds, quiet = TRUE)
res <- as.data.table(as.data.frame(results(dds, alpha = 0.05))); res[, gene := rowData(dds)$symbol]
res <- res[!is.na(padj)]; res[, rank := rank(padj, ties.method = "first")]

m <- merge(paper, res, by = "gene", all.x = TRUE)
out <- list(
  paper_genes = nrow(paper),
  matched_in_tested_set = sum(!is.na(m$padj)),
  not_tested_after_filter = sum(is.na(m$padj)),
  significant_here = sum(m$padj < 0.05, na.rm = TRUE),
  direction_agrees = sum(sign(m$log2FoldChange) == m$paper_dir, na.rm = TRUE),
  mine_total_sig = sum(res$padj < 0.05),
  mine_sig_not_in_paper = sum(res$padj < 0.05 & !res$gene %in% paper$gene),
  top316_by_padj_overlap = sum(res[rank <= 316]$gene %in% paper$gene),
  lfc1_sig_overlap = sum(res[padj < 0.05 & abs(log2FoldChange) > 1]$gene %in% paper$gene),
  lfc1_sig_total = sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1),
  lfc15_sig_total = sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1.5),
  lfc15_sig_overlap = sum(res[padj < 0.05 & abs(log2FoldChange) > 1.5]$gene %in% paper$gene),
  table_column_is_log2_share = round(col_is_log2, 3),
  median_paper_abs_log2fc = round(median(abs(paper$log2_from_fpkm)), 2),
  median_mine_sig_abs_log2fc = round(median(abs(res[padj < 0.05]$log2FoldChange)), 2)
)
print(unlist(out))
fwrite(m[, .(gene, paper_dex_fpkm = dex_fpkm, paper_untr_fpkm = untr_fpkm, paper_dir, my_log2FC = log2FoldChange, my_padj = padj, my_rank = rank)], "results/himes_comparison.csv")
jsonlite::write_json(out, "results/himes_comparison_summary.json", auto_unbox = TRUE, pretty = TRUE)
