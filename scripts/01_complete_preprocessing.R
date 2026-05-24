# ============================================================
# Complete preprocessing: DE analysis + gene symbol mapping
# ============================================================

library(GEOquery)
library(limma)

set.seed(42)
base_dir <- "/root/sources/RTX_MN_NetworkMedicine"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Download data with platform annotation
# ============================================================
cat("Downloading GSE108113...\n")
suppressMessages({
  gset <- getGEO("GSE108113", GSEMatrix = TRUE, getGPL = TRUE, destdir = data_dir)
})

expr_data <- exprs(gset[[1]])
pheno <- pData(gset[[1]])
fdat <- fData(gset[[1]])

cat(sprintf("Expression data: %d probes x %d samples\n", nrow(expr_data), ncol(expr_data)))

# Build probe-to-gene mapping
probe_to_entrez <- setNames(as.character(fdat$ENTREZ_GENE_ID), as.character(fdat$ID))
entrez_to_symbol_file <- file.path(data_dir, "entrez_to_symbol.rds")

if (!file.exists(entrez_to_symbol_file)) {
  cat("Downloading NCBI gene info...\n")
  download.file(
    "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    file.path(data_dir, "Homo_sapiens.gene_info.gz"),
    method = "auto", quiet = TRUE
  )
  gene_info <- read.delim(gzfile(file.path(data_dir, "Homo_sapiens.gene_info.gz")),
                          stringsAsFactors = FALSE)
  entrez_to_symbol <- setNames(gene_info$Symbol, as.character(gene_info$GeneID))
  saveRDS(entrez_to_symbol, entrez_to_symbol_file)
} else {
  entrez_to_symbol <- readRDS(entrez_to_symbol_file)
}
cat(sprintf("NCBI gene info loaded: %d entries\n", length(entrez_to_symbol)))

# ============================================================
# 2. Annotate samples with disease labels
# ============================================================
titles <- as.character(pheno$title)
disease_label <- rep(NA, length(titles))
for (i in seq_along(titles)) {
  t <- titles[i]
  if (grepl("MGN|Membranous", t, ignore.case = TRUE)) disease_label[i] <- "MGN"
  else if (grepl("Donor|Living|Normal|TN|Tumor|Nephrectom", t, ignore.case = TRUE)) disease_label[i] <- "Normal"
  else if (grepl("MCD|Minimal", t, ignore.case = TRUE)) disease_label[i] <- "MCD"
  else if (grepl("FSGS", t, ignore.case = TRUE)) disease_label[i] <- "FSGS"
  else disease_label[i] <- "Other"
}
cat("Sample distribution:\n")
print(table(disease_label))

# ============================================================
# 3. DE: MGN vs Normal
# ============================================================
cat("\n=== DE: MGN vs Normal ===\n")
keep <- disease_label %in% c("MGN", "Normal")
expr_sub <- expr_data[, keep]
group <- factor(disease_label[keep], levels = c("Normal", "MGN"))
cat(sprintf("MGN: %d, Normal: %d\n", sum(group == "MGN"), sum(group == "Normal")))

design <- model.matrix(~ group)
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)
de <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

# Add gene info
de$Probe <- rownames(de)
de$Entrez <- probe_to_entrez[de$Probe]
de$Entrez[is.na(de$Entrez)] <- ""

# Map to gene symbols
de$Gene <- entrez_to_symbol[de$Entrez]
de$Gene[is.na(de$Gene) | de$Gene == ""] <- de$Probe[is.na(de$Gene) | de$Gene == ""]

# Remove non-unique entries
de <- de[!duplicated(de$Gene), ]

n_sig <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
n_up <- sum(de$adj.P.Val < 0.05 & de$logFC > 0, na.rm = TRUE)
n_down <- sum(de$adj.P.Val < 0.05 & de$logFC < 0, na.rm = TRUE)
cat(sprintf("Significant DEGs: %d (Up: %d, Down: %d)\n", n_sig, n_up, n_down))

write.csv(de, file.path(results_dir, "DE_MGN_vs_Normal.csv"), row.names = FALSE)

# Save clean gene list
sig_genes <- de$Gene[de$adj.P.Val < 0.05]
writeLines(sig_genes, file.path(results_dir, "sig_deg_genes.txt"))
cat(sprintf("Saved %d significant DEGs to sig_deg_genes.txt\n", length(sig_genes)))

# Save full data for next steps
saveRDS(list(
  de = de,
  disease_label = disease_label,
  entrez_to_symbol = entrez_to_symbol
), file.path(results_dir, "preprocessed_data.rds"))

cat("\n✓ Preprocessing complete!\n")
