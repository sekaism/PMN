#!/usr/bin/env Rscript
# ============================================================
# External validation preprocessing and differential expression
# Dataset: GSE99340-GPL19184 glomerular MGN vs control samples
# ============================================================

suppressPackageStartupMessages({
  library(limma)
})

set.seed(20260519)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE))
} else {
  script_dir <- getwd()
}
base_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
if (!dir.exists(file.path(base_dir, "data"))) {
  base_dir <- "/workspaces/platform/MN"
}
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results", "external")
tables_dir <- file.path(base_dir, "tables")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

matrix_path <- file.path(data_dir, "GSE99340-GPL19184_series_matrix.txt.gz")
if (!file.exists(matrix_path)) {
  stop("Missing input file: ", matrix_path)
}

cat("=== External validation: GSE99340-GPL19184 ===\n")
cat("Input:", matrix_path, "\n")

# ------------------------------------------------------------
# Parse GEO series matrix metadata and expression table locally.
# ------------------------------------------------------------
con <- gzfile(matrix_path, open = "rt")
lines <- readLines(con, warn = FALSE)
close(con)

get_field <- function(prefix) {
  idx <- grep(paste0("^", prefix, "\\t"), lines)
  if (length(idx) == 0) return(NULL)
  strsplit(lines[idx[1]], "\\t", fixed = FALSE)[[1]][-1]
}
get_fields <- function(prefix) {
  idx <- grep(paste0("^", prefix, "\\t"), lines)
  if (length(idx) == 0) return(list())
  lapply(idx, function(i) strsplit(lines[i], "\\t", fixed = FALSE)[[1]][-1])
}
clean <- function(x) gsub('^"|"$', '', x)

accession <- clean(get_field("!Sample_geo_accession"))
title <- clean(get_field("!Sample_title"))
source <- clean(get_field("!Sample_source_name_ch1"))
platform <- clean(get_field("!Sample_platform_id"))
chars <- lapply(get_fields("!Sample_characteristics_ch1"), clean)

if (length(accession) == 0 || length(title) == 0) {
  stop("Could not parse sample metadata from ", matrix_path)
}

tissue <- if (length(chars) >= 1) chars[[1]] else rep(NA_character_, length(accession))
batch <- if (length(chars) >= 2) sub("^batch:\\s*", "", chars[[2]], ignore.case = TRUE) else rep(NA_character_, length(accession))
tissue <- sub("^tissue:\\s*", "", tissue, ignore.case = TRUE)

sample_info <- data.frame(
  geo_accession = accession,
  title = title,
  source_name = source,
  tissue = tissue,
  batch = batch,
  platform = platform,
  stringsAsFactors = FALSE
)

sample_info$disease_group <- "Other"
sample_info$disease_group[grepl("Membranous|MGN", sample_info$title, ignore.case = TRUE)] <- "MGN"
sample_info$disease_group[grepl("Tumor Nephrectomy|Cadaveric Donor|Living|Normal|control", sample_info$title, ignore.case = TRUE)] <- "Control"
sample_info$disease_group[grepl("Minimal|MCD", sample_info$title, ignore.case = TRUE)] <- "MCD"
sample_info$disease_group[grepl("FSGS", sample_info$title, ignore.case = TRUE)] <- "FSGS"
sample_info$disease_group[grepl("Diabetic|DN", sample_info$title, ignore.case = TRUE)] <- "DN"
sample_info$disease_group[grepl("IgA|IgAN", sample_info$title, ignore.case = TRUE)] <- "IgAN"
sample_info$disease_group[grepl("SLE|Lupus", sample_info$title, ignore.case = TRUE)] <- "SLE"

sample_info$compartment <- "Other"
sample_info$compartment[grepl("glomeruli", sample_info$tissue, ignore.case = TRUE) | grepl("-Glom-", sample_info$title, ignore.case = TRUE)] <- "Glomeruli"
sample_info$compartment[grepl("tubulo|interstitium", sample_info$tissue, ignore.case = TRUE) | grepl("-Tub-", sample_info$title, ignore.case = TRUE)] <- "Tubulointerstitium"

write.csv(sample_info, file.path(results_dir, "GSE99340_GPL19184_sample_metadata.csv"), row.names = FALSE)
cat("Sample disease distribution:\n")
print(table(sample_info$compartment, sample_info$disease_group))

begin <- grep("^!series_matrix_table_begin", lines)
end <- grep("^!series_matrix_table_end", lines)
if (length(begin) != 1 || length(end) != 1 || end <= begin) {
  stop("Could not locate expression matrix table in ", matrix_path)
}

table_text <- paste(lines[(begin + 1):(end - 1)], collapse = "\n")
expr_df <- read.delim(textConnection(table_text), check.names = FALSE, stringsAsFactors = FALSE)
colnames(expr_df) <- clean(colnames(expr_df))
expr_df[[1]] <- clean(expr_df[[1]])
rownames(expr_df) <- expr_df[[1]]
expr_df[[1]] <- NULL
expr <- as.matrix(expr_df)
storage.mode(expr) <- "numeric"

# Ensure expression sample order follows metadata order.
missing_samples <- setdiff(sample_info$geo_accession, colnames(expr))
if (length(missing_samples) > 0) {
  stop("Expression matrix missing samples: ", paste(missing_samples, collapse = ", "))
}
expr <- expr[, sample_info$geo_accession, drop = FALSE]

keep <- sample_info$compartment == "Glomeruli" & sample_info$disease_group %in% c("Control", "MGN")
if (sum(keep & sample_info$disease_group == "MGN") < 3 || sum(keep & sample_info$disease_group == "Control") < 3) {
  stop("Insufficient glomerular MGN/control samples after filtering")
}

expr_sub <- expr[, keep, drop = FALSE]
meta_sub <- sample_info[keep, , drop = FALSE]
group <- factor(meta_sub$disease_group, levels = c("Control", "MGN"))
cat(sprintf("Glomerular comparison: MGN=%d, Control=%d\n", sum(group == "MGN"), sum(group == "Control")))
write.csv(meta_sub, file.path(results_dir, "GSE99340_GPL19184_glom_MGN_vs_control_metadata.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Differential expression with limma.
# Values in the GEO matrix are already log2-scale normalized values.
# ------------------------------------------------------------
design <- model.matrix(~ group)
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)
de <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
de$Probe <- rownames(de)

# Probe IDs on GPL19184 are Entrez-like identifiers with suffixes such as _at.
de$Entrez <- sub("_.*$", "", de$Probe)
de$Entrez[!grepl("^[0-9]+$", de$Entrez)] <- ""

entrez_to_symbol_file <- file.path(data_dir, "entrez_to_symbol.rds")
if (file.exists(entrez_to_symbol_file)) {
  entrez_to_symbol <- readRDS(entrez_to_symbol_file)
} else {
  gene_info_path <- file.path(data_dir, "Homo_sapiens.gene_info.gz")
  if (!file.exists(gene_info_path)) {
    stop("Missing gene annotation file: ", gene_info_path)
  }
  gene_info <- read.delim(gzfile(gene_info_path), stringsAsFactors = FALSE)
  entrez_to_symbol <- setNames(gene_info$Symbol, as.character(gene_info$GeneID))
  saveRDS(entrez_to_symbol, entrez_to_symbol_file)
}

de$Gene <- unname(entrez_to_symbol[de$Entrez])
de$Gene[is.na(de$Gene) | de$Gene == ""] <- de$Probe[is.na(de$Gene) | de$Gene == ""]
de$Gene <- trimws(as.character(de$Gene))
de <- de[!is.na(de$Gene) & de$Gene != "" & de$Gene != "NA", ]

# Keep the most significant probe per gene and use this gene-level table for
# all downstream external validation outputs.
de <- de[order(de$P.Value), ]
de <- de[!duplicated(de$Gene), ]
de <- de[order(de$P.Value), ]

n_sig_05 <- sum(de$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_10 <- sum(de$adj.P.Val < 0.10, na.rm = TRUE)
n_nom <- sum(de$P.Value < 0.05, na.rm = TRUE)
cat(sprintf("External DE genes: FDR<0.05=%d, FDR<0.10=%d, nominal P<0.05=%d\n", n_sig_05, n_sig_10, n_nom))
cat(sprintf("External DE direction at FDR<0.05: up=%d, down=%d\n", sum(de$adj.P.Val < 0.05 & de$logFC > 0, na.rm = TRUE), sum(de$adj.P.Val < 0.05 & de$logFC < 0, na.rm = TRUE)))

write.csv(de, file.path(results_dir, "GSE99340_GPL19184_glom_MGN_vs_control_DE.csv"), row.names = FALSE)
write.csv(de[de$adj.P.Val < 0.05, ], file.path(results_dir, "GSE99340_GPL19184_glom_MGN_vs_control_sig_genes.csv"), row.names = FALSE)

saveRDS(
  list(expr = expr_sub, metadata = meta_sub, group = group, de = de),
  file = file.path(results_dir, "GSE99340_GPL19184_glom_MGN_vs_control.rds")
)

cat("Saved external validation preprocessing outputs to ", results_dir, "\n", sep = "")
