# Glucocorticoid response in airway smooth muscle: a DESeq2 re-analysis

Re-analysis of public bulk RNA-seq (GEO GSE52778, via the Bioconductor `airway` package): four primary human airway
smooth muscle cell lines, each treated with dexamethasone or left untreated. Original study: Himes et al. 2014,
[doi:10.1371/journal.pone.0099625](https://doi.org/10.1371/journal.pone.0099625).

**Report:** [`report/report.pdf`](report/report.pdf) (source: `report/report.qmd`, every number is computed inline).

## Result in brief
- 4,081 genes differentially expressed at adjusted p < 0.05 (951 with |log2FC| > 1) using a paired design (`~ cell + dex`).
- All five genes the original paper highlighted (*CRISPLD2, DUSP1, KLF15, PER1, TSC22D3*) are significantly up-regulated.
- Ignoring the pairing finds 2,773 genes, so modelling the cell line matters.
- PC1 of the variance-stabilised data is the drug (R² 0.95); PC2/PC3 are the cell lines.
- The paper reported 316 genes (Cuffdiff). Compared with its Table S3 (`R/compare_himes.R`): 301 of the 304 testable genes are significant here, all in the same direction, so the paper's list is nearly a subset of mine; the extra genes have much smaller effects (median |log2FC| 0.58 vs 1.73). I did not run Cuffdiff, so I do not claim to know why the tool gave fewer genes.

## Reproduce
```r
# R >= 4.4 with Bioconductor: DESeq2, airway, apeglm, clusterProfiler, org.Hs.eg.db, pheatmap; CRAN: ggplot2, dplyr, gridExtra, knitr
quarto render report/report.qmd
```

## Provenance
Written and run with the assistance of Claude Code (Anthropic). I am responsible for checking and understanding the results. No part derives from university coursework or assessed material.
`report/apa.csl` is the APA style from the [CSL styles repository](https://github.com/citation-style-language/styles) (CC BY-SA 3.0).
