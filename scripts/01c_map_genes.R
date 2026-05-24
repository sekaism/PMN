# ============================================================
# 1c: Fix gene symbols using NCBI Entrez API
# ============================================================

library(GEOquery)
library(limma)
library(httr)
library(jsonlite)

set.seed(42)
base_dir <- "/root/sources/RTX_MN_NetworkMedicine"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")

# Load previously saved raw DE data
de <- read.csv(file.path(results_dir, "DE_MGN_vs_Normal.csv"))

# Extract numeric Entrez IDs from probe IDs
probes <- as.character(de$Probe)
entrez_ids <- gsub("_at$", "", probes)
entrez_ids <- gsub("_s_at$", "", entrez_ids) 
entrez_ids <- gsub("_x_at$", "", entrez_ids)
numeric_ids <- grep("^[0-9]+$", entrez_ids, value = TRUE)
cat(sprintf("Numeric Entrez IDs: %d\n", length(numeric_ids)))

# Batch query using NCBI efetch
# Use the fData data instead - it already has ENTREZ_GENE_ID
cat("Loading GSE108113 feature data...\n")
suppressMessages({
  gset <- getGEO("GSE108113", GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
})
fdat <- fData(gset[[1]])

# fdat has ID and ENTREZ_GENE_ID columns
# Map probe IDs to Entrez IDs
probe_to_entrez <- setNames(as.character(fdat$ENTREZ_GENE_ID), as.character(fdat$ID))

# Update DE results with Entrez IDs
de$Entrez <- probe_to_entrez[as.character(de$Probe)]
de$Entrez[is.na(de$Entrez) | de$Entrez == ""] <- "0"

# Download gene_info from NCBI for symbol lookup
cat("Downloading NCBI gene info for Homo sapiens...\n")
gene_info_url <- "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
gene_info_file <- file.path(data_dir, "Homo_sapiens.gene_info.gz")

if (!file.exists(gene_info_file)) {
  download.file(gene_info_url, gene_info_file, method = "auto", quiet = TRUE)
}
cat("Loading gene info...\n")
gene_info <- read.delim(gzfile(gene_info_file), stringsAsFactors = FALSE)
cat(sprintf("Loaded %d human genes.\n", nrow(gene_info)))

# Create lookup: Entrez ID -> Symbol
entrez_to_symbol <- setNames(gene_info$Symbol, as.character(gene_info$GeneID))

# Map DE results
unique_entrez <- unique(de$Entrez[de$Entrez != "0"])
cat(sprintf("Unique Entrez IDs in DE results: %d\n", length(unique_entrez)))
mapped <- sum(unique_entrez %in% names(entrez_to_symbol))
cat(sprintf("Mapped to symbols: %d\n", mapped))

de$Gene <- entrez_to_symbol[de$Entrez]
# For unmapped probes, use the probe ID as fallback
de$Gene[is.na(de$Gene) | de$Gene == ""] <- de$Probe[is.na(de$Gene) | de$Gene == ""]

# Remove rows with NA genes
de <- de[!is.na(de$Gene) & de$Gene != "", ]

# Keep most significant per gene
de <- de[order(de$adj.P.Val), ]
de <- de[!duplicated(de$Gene), ]

n_sig <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
n_up <- sum(de$adj.P.Val < 0.05 & de$logFC > 0, na.rm = TRUE)
n_down <- sum(de$adj.P.Val < 0.05 & de$logFC < 0, na.rm = TRUE)

cat(sprintf("\nFinal results:\n"))
cat(sprintf("Total unique genes: %d\n", nrow(de)))
cat(sprintf("Significant DEGs (adj.P < 0.05): %d\n", n_sig))
cat(sprintf("  Up in MN: %d\n", n_up))
cat(sprintf("  Down in MN: %d\n", n_down))

# Save
write.csv(de, file.path(results_dir, "DE_MGN_vs_Normal.csv"), row.names = FALSE)
sig_genes <- de$Gene[de$adj.P.Val < 0.05]
writeLines(sig_genes, file.path(results_dir, "sig_deg_genes.txt"))
cat(sprintf("Saved %d significant DEGs to file.\n", length(sig_genes)))

cat("\nDone! Gene symbols properly mapped.\n")
