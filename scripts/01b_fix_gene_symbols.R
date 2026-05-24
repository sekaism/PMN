# ============================================================
# 1b: Re-run DE analysis with proper gene symbol annotation
# ============================================================
# The original analysis used probe IDs as "Gene" names.
# This script re-downloads with platform annotation and maps
# probes to official gene symbols.
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
# Download GSE108113 WITH platform annotation
# ============================================================
cat("Downloading GSE108113 with platform annotation...\n")
suppressMessages({
  gset <- getGEO("GSE108113", GSEMatrix = TRUE, getGPL = TRUE, destdir = data_dir)
})

expr_set <- gset[[1]]
expr_data <- exprs(expr_set)
pheno <- pData(expr_set)

# Get feature data (probe-to-gene mapping)
fdat <- fData(expr_set)
cat("Feature data columns:", paste(colnames(fdat), collapse = ", "), "\n")

# Find the Gene Symbol column
gene_col <- NULL
for (col in colnames(fdat)) {
  if (grepl("gene.*sym|symbol|gene_assignment", col, ignore.case = TRUE)) {
    gene_col <- col
    break
  }
}

if (is.null(gene_col)) {
  cat("Looking in fData for gene info...\n")
  # Try common column names
  for (try_col in c("Gene Symbol", "Symbol", "Gene", "GENE_SYMBOL", 
                    "gene_symbol", "GeneSymbol", "genesymbol")) {
    if (try_col %in% colnames(fdat)) {
      gene_col <- try_col
      break
    }
  }
}

cat("Using gene column:", gene_col, "\n")

# Extract gene symbols from the appropriate column
gene_info <- as.character(fdat[, gene_col])

# Some platforms have format "GeneSymbol // Description" - extract first part
gene_symbols <- gsub(" // .*$", "", gene_info)
gene_symbols <- gsub("---", NA, gene_symbols)
gene_symbols <- trimws(gene_symbols)

# Count mapped probes
mapped_count <- sum(!is.na(gene_symbols) & gene_symbols != "")
cat(sprintf("Probes with gene symbols: %d / %d\n", mapped_count, length(gene_symbols)))

# ============================================================
# Re-run DE analysis with gene symbols
# ============================================================
titles <- as.character(pheno$title)

disease_label <- rep(NA, length(titles))
for (i in seq_along(titles)) {
  t <- titles[i]
  if (grepl("MGN|Membranous", t, ignore.case = TRUE)) disease_label[i] <- "MGN"
  else if (grepl("Donor|Living|Normal|TN|Tumor|Nephrectom", t, ignore.case = TRUE)) disease_label[i] <- "Normal"
  else if (grepl("MCD|Minimal", t, ignore.case = TRUE)) disease_label[i] <- "MCD"
  else if (grepl("FSGS", t, ignore.case = TRUE)) disease_label[i] <- "FSGS"
  else if (grepl("IgA", t, ignore.case = TRUE)) disease_label[i] <- "IgAN"
  else if (grepl("RPGN|Rapidly|vasculitis", t, ignore.case = TRUE)) disease_label[i] <- "RPGN"
  else disease_label[i] <- "Other"
}

cat("Disease distribution:\n")
print(table(disease_label))

# ============================================================
# DE: MGN vs Normal (Living Donors)
# ============================================================
cat("\n--- DE: MGN vs Normal ---\n")
keep <- disease_label %in% c("MGN", "Normal")
expr_sub <- expr_data[, keep]
group <- factor(disease_label[keep], levels = c("Normal", "MGN"))
cat(sprintf("MGN: %d, Normal: %d\n", sum(group == "MGN"), sum(group == "Normal")))

design <- model.matrix(~ group)
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)
de <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
de$Probe <- rownames(de)
de$Gene <- gene_symbols[match(de$Probe, rownames(expr_data))]
de$Gene[is.na(de$Gene) | de$Gene == ""] <- de$Probe[is.na(de$Gene) | de$Gene == ""]

# Remove probes without gene symbols
de <- de[!is.na(de$Gene) & de$Gene != "", ]

# For probes mapping to same gene, keep the most significant one
de <- de[order(de$adj.P.Val), ]
de <- de[!duplicated(de$Gene), ]

n_sig <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
cat(sprintf("Significant DEGs (adj.P < 0.05): %d\n", n_sig))
cat(sprintf("  Up: %d, Down: %d\n", 
    sum(de$adj.P.Val < 0.05 & de$logFC > 0, na.rm = TRUE),
    sum(de$adj.P.Val < 0.05 & de$logFC < 0, na.rm = TRUE)))

write.csv(de, file.path(results_dir, "DE_MGN_vs_Normal.csv"), row.names = FALSE)

# Also create a clean gene list
sig_genes <- de$Gene[de$adj.P.Val < 0.05]
directions <- ifelse(de$logFC[de$adj.P.Val < 0.05] > 0, "UP", "DOWN")
gene_list_df <- data.frame(Gene = sig_genes, Direction = directions)
write.csv(gene_list_df, file.path(results_dir, "sig_deg_genes.csv"), row.names = FALSE)
cat(sprintf("Saved %d significant DEG genes.\n", length(sig_genes)))

# ============================================================
# DE: MGN vs All Others (broader signature)
# ============================================================
cat("\n--- DE: MGN vs All Others ---\n")
group_all <- factor(ifelse(disease_label == "MGN", "MGN", "Other"), 
                    levels = c("Other", "MGN"))
design2 <- model.matrix(~ group_all)
fit2 <- lmFit(expr_data, design2)
fit2 <- eBayes(fit2)
de_all <- topTable(fit2, coef = 2, number = Inf, sort.by = "P")
de_all$Probe <- rownames(de_all)
de_all$Gene <- gene_symbols[match(de_all$Probe, rownames(expr_data))]
de_all$Gene[is.na(de_all$Gene) | de_all$Gene == ""] <- de_all$Probe[is.na(de_all$Gene) | de_all$Gene == ""]
de_all <- de_all[!duplicated(de_all$Gene), ]

write.csv(de_all, file.path(results_dir, "DE_MGN_vs_All.csv"), row.names = FALSE)

cat("\nStep 1b complete! Gene symbols properly mapped.\n")
cat(sprintf("Final significant DEGs with gene symbols: %d\n", 
    sum(de$adj.P.Val < 0.05, na.rm = TRUE)))
