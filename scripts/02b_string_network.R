# ============================================================
# Step 2b: Build full PPI network using STRINGdb
# ============================================================
# Downloads the human STRING network (full) via the STRINGdb
# R package. Extracts subnetwork for MN DEGs + drug targets.
# Output: graphml file for Python NetMedPy analysis
# ============================================================

library(STRINGdb)

set.seed(42)
base_dir <- "/root/sources/RTX_MN_NetworkMedicine"
data_dir <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load STRING database (human, version 11.5)
# ============================================================
cat("Loading STRINGdb...\n")
string_db <- STRINGdb$new(version = "11.5", species = 9606, 
                          score_threshold = 400, 
                          input_directory = data_dir)
cat("STRINGdb loaded.\n")

# ============================================================
# Load MN DEGs
# ============================================================
cat("Loading MN DEGs...\n")
de_mn <- read.csv(file.path(results_dir, "DE_MGN_vs_Normal.csv"))
sig_de <- de_mn[de_mn$adj.P.Val < 0.05, ]
deg_genes <- unique(as.character(sig_de$Gene))
deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != ""]
cat(sprintf("  %d significant DEGs loaded.\n", length(deg_genes)))

# GWAS risk genes
gwas_genes <- c("PLA2R1", "NFKB1", "IRF4", "HLA-DQA1", "HLA-DRB1")

# Map gene symbols to STRING IDs
all_seed_genes <- unique(c(deg_genes, gwas_genes))
cat(sprintf("  %d total seed genes.\n", length(all_seed_genes)))

# Map to STRING identifiers
mapped <- string_db$map(data.frame(gene = all_seed_genes), 
                        "gene", removeUnmappedRows = TRUE)
cat(sprintf("  Mapped %d / %d genes to STRING IDs.\n", 
            nrow(mapped), length(all_seed_genes)))

# Save mapped genes
write.csv(mapped, file.path(results_dir, "mapped_genes.csv"), row.names = FALSE)

# ============================================================
# Get interaction partners and build network
# ============================================================
cat("Retrieving interaction network...\n")
# Get the subnetwork for our proteins
# Use get_subnetwork which returns all interactions between mapped proteins
string_ids <- mapped$STRING_id[1:min(500, nrow(mapped))]  # Top 500

# Get interactions between these proteins
interactions <- string_db$get_interactions(string_ids)
cat(sprintf("  Retrieved %d interactions.\n", nrow(interactions)))

# Get additional neighbors (one step) to expand the network
# This makes the network richer for network medicine proximity analysis
neighbors <- string_db$get_neighbors(string_ids)
cat(sprintf("  Retrieved %d neighboring proteins.\n", length(neighbors)))

# Combine to get expanded network
all_string_ids <- unique(c(string_ids, neighbors))
cat(sprintf("  Total proteins in expanded network: %d\n", length(all_string_ids)))

# Get all interactions among the expanded set
if (length(all_string_ids) > 500) {
  # Process in batches
  all_interactions <- data.frame()
  for (i in seq(1, length(all_string_ids), 200)) {
    batch <- all_string_ids[i:min(i+199, length(all_string_ids))]
    rows <- string_db$get_interactions(batch)
    all_interactions <- rbind(all_interactions, rows)
    cat(sprintf("  Batch %d: %d interactions\n", i, nrow(rows)))
  }
} else {
  all_interactions <- string_db$get_interactions(all_string_ids)
}

cat(sprintf("Total interactions: %d\n", nrow(all_interactions)))

# ============================================================
# Build network and save as edge list for Python
# ============================================================
cat("Building network...\n")

# Map STRING IDs back to gene symbols
# Create a mapping dictionary
hits <- string_db$mp(mapped$gene)
gene_map <- setNames(mapped$gene, mapped$STRING_id)

# Build edge list with gene symbols
edge_list <- data.frame(
  from = character(),
  to = character(),
  weight = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(all_interactions)) {
  p1 <- as.character(all_interactions$protein1[i])
  p2 <- as.character(all_interactions$protein2[i])
  combined_score <- all_interactions$combined_score[i] / 1000  # Normalize to 0-1
  
  # Map to gene symbols if available
  g1 <- ifelse(p1 %in% names(gene_map), gene_map[p1], p1)
  g2 <- ifelse(p2 %in% names(gene_map), gene_map[p2], p2)
  
  edge_list <- rbind(edge_list, data.frame(
    from = g1, to = g2, weight = combined_score,
    stringsAsFactors = FALSE
  ))
}

# Remove self-loops
edge_list <- edge_list[edge_list$from != edge_list$to, ]

# Filter by combined_score >= 0.4 (standard STRING threshold)
edge_list <- edge_list[edge_list$weight >= 0.4, ]

# Keep unique edges (undirected)
edge_list <- edge_list[!duplicated(t(apply(edge_list[,1:2], 1, sort))), ]

cat(sprintf("Edge list: %d edges, %d unique nodes\n", 
            nrow(edge_list),
            length(unique(c(edge_list$from, edge_list$to)))))

# Save edge list
write.csv(edge_list, file.path(data_dir, "string_edge_list.csv"), row.names = FALSE)

# Also save as a simple format for NetMedPy
# NetMedPy expects a graph with node names
# We'll save as edgelist format that networkx can read
write.table(edge_list[, c("from", "to", "weight")],
            file.path(data_dir, "string_network_edgelist.txt"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Save gene mapping
write.csv(data.frame(STRING_id = names(gene_map), Gene = gene_map),
          file.path(results_dir, "string_gene_mapping.csv"), row.names = FALSE)

cat("\nStep 2b complete!\n")
cat(sprintf("Network saved with %d nodes and %d edges.\n",
            length(unique(c(edge_list$from, edge_list$to))),
            nrow(edge_list)))
cat("Ready for Step 3: NetMedPy network medicine analysis.\n")
