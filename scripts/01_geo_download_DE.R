# ============================================================
# Step 1: Download GEO dataset and run differential expression
# ============================================================
# Dataset: GSE108113 - Glomerular transcriptomes from nephrotic
#           syndrome patients (44 MGN + 6 normal + 107 other)
# Also downloads GSE216841 if available
# ============================================================

library(GEOquery)
library(limma)

set.seed(42)
data_dir <- "../data"
results_dir <- "../results"
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Dataset 1: GSE108113 (Affymetrix ST 2.1) - PRIMARY
# 157 glomerular samples: MGN(44), MCD, FSGS, Normal(6), etc.
# ============================================================
cat("=== Downloading GSE108113 (primary dataset) ===\n")
suppressMessages({
  gset <- getGEO("GSE108113", GSEMatrix = TRUE, getGPL = FALSE, destdir = data_dir)
})

expr_set <- gset[[1]]
expr_data <- exprs(expr_set)
pheno <- pData(expr_set)
titles <- as.character(pheno$title)

cat("Total samples:", ncol(expr_data), "\n")

# Extract disease labels and compartment from titles
disease_label <- rep(NA, length(titles))
for (i in seq_along(titles)) {
  t <- titles[i]
  if (grepl("MGN|Membranous", t, ignore.case = TRUE)) disease_label[i] <- "MGN"
  else if (grepl("MCD|Minimal Change", t, ignore.case = TRUE)) disease_label[i] <- "MCD"
  else if (grepl("FSGS", t, ignore.case = TRUE)) disease_label[i] <- "FSGS"
  else if (grepl("Donor|Living donor|Normal", t, ignore.case = TRUE)) disease_label[i] <- "Normal"
  else if (grepl("IgA", t, ignore.case = TRUE)) disease_label[i] <- "IgAN"
  else if (grepl("SLE|Lupus", t, ignore.case = TRUE)) disease_label[i] <- "SLE"
  else if (grepl("RPGN|Rapidly|vasculitis", t, ignore.case = TRUE)) disease_label[i] <- "RPGN"
  else if (grepl("DN|Diabetic", t, ignore.case = TRUE)) disease_label[i] <- "DN"
  else if (grepl("HT|Hypertensive", t, ignore.case = TRUE)) disease_label[i] <- "HT"
  else if (grepl("TN|Tumor|Nephrectom", t, ignore.case = TRUE)) disease_label[i] <- "Normal"
  else disease_label[i] <- "Other"
}

cat("Disease distribution:\n")
print(table(disease_label))

# ============================================================
# DE Analysis 1: MGN vs Normal (Living Donors + Tumor Nephrectomies)
# ============================================================
cat("\n--- DE: MGN vs Normal ---\n")
keep_mn <- disease_label %in% c("MGN", "Normal")
expr_mn <- expr_data[, keep_mn]
group_mn <- factor(disease_label[keep_mn], levels = c("Normal", "MGN"))
cat("MGN:", sum(group_mn == "MGN"), ", Normal:", sum(group_mn == "Normal"), "\n")

design <- model.matrix(~ group_mn)
fit <- lmFit(expr_mn, design)
fit <- eBayes(fit)
de_mn <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
de_mn$Gene <- rownames(de_mn)

n_sig <- sum(de_mn$adj.P.Val < 0.05, na.rm = TRUE)
cat("Significant DEGs (adj.P < 0.05):", n_sig, "\n")
cat("Up in MGN:", sum(de_mn$adj.P.Val < 0.05 & de_mn$logFC > 0, na.rm = TRUE), "\n")
cat("Down in MGN:", sum(de_mn$adj.P.Val < 0.05 & de_mn$logFC < 0, na.rm = TRUE), "\n")

write.csv(de_mn, file.path(results_dir, "DE_MGN_vs_Normal.csv"), row.names = FALSE)
saveRDS(list(expr = expr_mn, group = group_mn, de = de_mn),
        file = file.path(results_dir, "GSE108113_MGNvsNormal.rds"))

# ============================================================
# DE Analysis 2: MGN vs All Others (for broader disease signature)
# ============================================================
cat("\n--- DE: MGN vs All Others ---\n")
group_all <- ifelse(disease_label == "MGN", "MGN", "Other")
group_all <- factor(group_all, levels = c("Other", "MGN"))
design2 <- model.matrix(~ group_all)
fit2 <- lmFit(expr_data, design2)
fit2 <- eBayes(fit2)
de_all <- topTable(fit2, coef = 2, number = Inf, sort.by = "P")
de_all$Gene <- rownames(de_all)

n_sig2 <- sum(de_all$adj.P.Val < 0.05, na.rm = TRUE)
cat("Significant DEGs (adj.P < 0.05):", n_sig2, "\n")
cat("Up in MGN:", sum(de_all$adj.P.Val < 0.05 & de_all$logFC > 0, na.rm = TRUE), "\n")
cat("Down in MGN:", sum(de_all$adj.P.Val < 0.05 & de_all$logFC < 0, na.rm = TRUE), "\n")

write.csv(de_all, file.path(results_dir, "DE_MGN_vs_All.csv"), row.names = FALSE)
saveRDS(list(expr = expr_data, group = group_all, de = de_all, disease_label = disease_label),
        file = file.path(results_dir, "GSE108113_all.rds"))

# ============================================================
# Enrichment analysis of top DEGs using enrichR
# ============================================================
cat("\n--- Running GO enrichment on MGN vs Normal DEGs ---\n")
sig_genes <- de_mn$Gene[which(de_mn$adj.P.Val < 0.05)]
if (length(sig_genes) >= 10) {
  library(enrichR)
  dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
  enriched <- enrichr(sig_genes, dbs)
  
  for (db in names(enriched)) {
    if (nrow(enriched[[db]]) > 0) {
      outfile <- file.path(results_dir, paste0("enrichment_", gsub(" ", "_", db), ".csv"))
      write.csv(enriched[[db]][1:min(50, nrow(enriched[[db]])), ], outfile, row.names = FALSE)
      cat("  Saved:", basename(outfile), "\n")
    }
  }
}

cat("\nStep 1 complete!\n")
cat("Output files:\n")
cat("  - DE_MGN_vs_Normal.csv (primary DEG list)\n")
cat("  - DE_MGN_vs_All.csv (MGN vs others DEG list)\n")
cat("  - enrichment_GO*.csv, enrichment_KEGG*.csv\n")
