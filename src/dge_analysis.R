# Load libraries
library(edgeR)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(tibble)

# --- Define and create results directory ---
home_dir <- "/Users/isidorabeslic/workspace/egfrKD-AstroRNAseq"
results_dir <- file.path(home_dir, "results")
data_dir <- file.path(home_dir, "data")

dge_dir <- file.path(results_dir, "dge")
if (!dir.exists(dge_dir)) {
  dir.create(dge_dir, recursive = TRUE)
}

feature_counts_dir <- file.path(results_dir, "feature_counts")
if (!dir.exists(feature_counts_dir)) {
  dir.create(feature_counts_dir, recursive = TRUE)
}

# --- Set working directory (optional) ---
setwd(home_dir)

# --- Load design matrix ---
design_matrix <- read.csv(file.path(data_dir, "design_matrix.csv"), stringsAsFactors = TRUE)
design_matrix$condition <- relevel(factor(design_matrix$condition), ref = "sgRosa26")

# --- Read and process individual count files ---
read_counts <- function(filename) {
  full_path <- file.path(feature_counts_dir, filename)
  counts <- read.delim(full_path, skip = 1, row.names = 1)
  counts <- counts[, ncol(counts), drop = FALSE]  # keep only the counts column
  colnames(counts) <- gsub("\\.counts\\.txt$", "", basename(filename))
  return(counts)
}

count_files <- paste0(design_matrix$sample, ".counts.txt")
counts_list <- lapply(count_files, read_counts)
counts_matrix <- do.call(cbind, counts_list)

stopifnot(all(colnames(counts_matrix) == design_matrix$sample))

# --- Collapse technical replicates ---
bio_reps <- unique(design_matrix$group)
collapsed_counts <- sapply(bio_reps, function(rep) {
  cols <- which(design_matrix$group == rep)
  rowSums(counts_matrix[, cols, drop = FALSE])
})
collapsed_counts <- as.matrix(collapsed_counts)
gene_means <- rowMeans(collapsed_counts)
collapsed_counts <- collapsed_counts[gene_means >= 0.5, ]

# --- Build new design matrix ---
collapsed_design <- design_matrix %>%
  distinct(group, condition) %>%
  arrange(group)
collapsed_design$condition <- relevel(factor(collapsed_design$condition), ref = "sgRosa26")

# --- Create DGEList object ---
dge <- DGEList(counts = collapsed_counts, group = collapsed_design$condition)

# --- Normalize and estimate dispersion ---
dge <- calcNormFactors(dge)
design <- model.matrix(~collapsed_design$condition)
dge <- estimateDisp(dge, design)

# --- Differential expression ---
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)
de_genes <- topTags(lrt, n = Inf)
de_genes <- as.data.frame(de_genes)

# --- Annotate genes using biomaRt ---
ensembl_ids <- rownames(de_genes)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")

annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

gene_map <- annotations %>%
  distinct(ensembl_gene_id, external_gene_name) %>%
  deframe()

de_genes$gene <- gene_map[ensembl_ids]
de_genes$gene <- ifelse(is.na(de_genes$gene) | de_genes$gene == "", ensembl_ids, de_genes$gene)

# --- Volcano plot ---
de_genes$logFDR <- -log10(de_genes$PValue)
de_genes$Significant <- de_genes$PValue < 0.05 & abs(de_genes$logFC) > 1

genes_to_label <- c("Nptx1", "Fos", "Cx3cl1", "Ccl7", "Csf3")

volcano_plot <- ggplot(de_genes, aes(x = logFC, y = logFDR)) +
  geom_point(aes(color = Significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text_repel(data = subset(de_genes, gene %in% genes_to_label),
                  aes(label = gene),
                  size = 4,
                  box.padding = 0.3,
                  point.padding = 0.2) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot: sgEGFR vs sgRosa26") +
  theme_minimal() +
  xlim(-12, 12) + 
  ylim(0, 7)

# --- Save plot and results ---
ggsave(file.path(dge_dir, "volcano_plot.pdf"), volcano_plot, width = 10, height = 8)
print(volcano_plot)

write.csv(de_genes, file.path(dge_dir, "dge_results.csv"))
