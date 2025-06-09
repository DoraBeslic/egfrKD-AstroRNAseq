# Install packages if needed:
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "msigdbr", "enrichplot"))

# Load libraries
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(dplyr)
library(ggplot2)

# Set directories
home_dir <- "/Users/isidorabeslic/workspace/egfrKD-AstroRNAseq"
results_dir <- file.path(home_dir, "results")
dge_dir <- file.path(results_dir, "dge")
plot_dir <- file.path(results_dir, "gsea")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Load differential expression results
dge_results <- read.csv(file.path(dge_dir, "dge_results.csv"), row.names = 1)

# Prepare ranked gene list (logFC named by Ensembl ID)
gene_list <- dge_results$logFC
names(gene_list) <- rownames(dge_results)
gene_list <- sort(gene_list, decreasing = TRUE)

### HALLMARK GSEA ----------------------------------------------------------

hallmark_df <- msigdbr(species = "Mus musculus", category = "H") %>%
  select(gs_name, ensembl_gene)

gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_df,
  pvalueCutoff = 1
)

# Target hallmark gene set
hallmark_id <- "HALLMARK_INFLAMMATORY_RESPONSE"
hallmark_row <- as.data.frame(gsea_hallmark) %>%
  filter(ID == hallmark_id)

# Plot and save
hallmark_plot <- gseaplot2(
  gsea_hallmark,
  geneSetID = hallmark_id,
  title = paste0("Hallmark: Inflammatory Response")
)

ggsave(
  filename = file.path(plot_dir, "hallmark_inflammatory_response.png"),
  plot = hallmark_plot,
  width = 8, height = 5
)

### REACTOME GSEA ----------------------------------------------------------

reactome_df <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>%
  select(gs_name, ensembl_gene)

gsea_reactome <- GSEA(
  geneList = gene_list,
  TERM2GENE = reactome_df,
  pvalueCutoff = 1
)

# Target reactome gene set
target_reactome_id <- "REACTOME_RECEPTOR_TYPE_TYROSINE_PROTEIN_PHOSPHATASES"
reactome_row <- as.data.frame(gsea_reactome) %>%
  filter(ID == target_reactome_id)

reactome_plot <- gseaplot2(
  gsea_reactome,
  geneSetID = target_reactome_id,
  title = paste0(gsub("_", " ", target_reactome_id))
  )
  
ggsave(
  filename = file.path(plot_dir, "reactome_tyrosine_phosphatases.png"),
  plot = reactome_plot,
  width = 8, height = 5
  )

### Print p-values ----------------------------------------------------------
cat("\n--- Hallmark Result ---\n")
print(hallmark_row[, c("ID", "Description", "pvalue", "NES")])

cat("\n--- Reactome Result ---\n")
print(reactome_row[, c("ID", "Description", "pvalue", "NES")])
