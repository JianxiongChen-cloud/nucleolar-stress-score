setwd("/home/xxm_xxm/CJX_workspace/geneset")
# Literature-derived upregulated genes
literature_up_genes <- c(
  "AATF", "AKT1", "ATF3", "ATM", "ATR", "BAX", "BBC3", "CCDC137",
  "CDKN1A", "CDKN1B", "CDKN2A", "CHEK1", "CNOT2", "CORO2B", "CTNNB1",
  "E2F1", "EEF2", "EEF2K", "EMG1", "FCN3", "GADD45A", "HIF1A",
  "HNRNPK", "MAP3K5", "MAP3K8", "MAPK14", "MDM2", "MDM4", "NEAT1",
  "NCL", "NOP53", "NPM1", "KGD4", "MYBBP1A", "MYC", "MYCN", "NEDD8",
  "PAK1IP1", "PIM1", "PML", "PRKAA1", "PRDM1", "RGS6", "RPL11",
  "RPL22", "RPL23", "RPL26", "RPL3", "RPL4", "RPL5", "RPL6", "RPS14",
  "RPS2", "RPS25", "RPS27", "RPS27A", "RPS27L", "RPS7", "RPS9",
  "SIRT7", "SNORA13", "SOCS3", "SRSF1", "TOPBP1", "TP53", "WNT4"
)

# Literature-derived downregulated genes
literature_down_genes <- c(
  "AGTPBP1", "AKT1", "AKT1S1", "BCL2", "BOP1", "BRIX1", "CCDC137",
  "CCND1", "CDK4", "CDK6", "DHODH", "DHX16", "DHX33", "E2F1", "EGR1",
  "EIF2AK4", "EXOSC8", "FANCA", "FBXO7", "FTO", "GNL3", "GRK5", "GRWD1",
  "HEATR1", "IGF1", "IMPDH2", "KMT5A", "LAS1L", "LIN28A", "MDM2", "MDM4",
  "MRPL44", "MRPS16", "MTREX", "MYC", "NBN", "NCL", "NDC80", "NEDD8",
  "NFE2L2", "NIFK", "NOD2", "NOL12", "NOL7", "NOLC1", "NOP2", "NOP53",
  "NPM1", "NUF2", "NUMA1", "PELP1", "PIM1", "POLR1A", "POLR1G", "PPAN",
  "PSMD9", "PYGO2", "RAP1GDS1", "RBM28", "RIOK1", "RIOK2", "RPL11",
  "RPL13", "RPL18", "RPL23", "RPL27A", "RPL32", "RPL37", "RPL9", "RPS14",
  "RPS15A", "RPS19", "RPS2", "RPS26", "RPS27A", "RPS27L", "RPS6",
  "RPS6KA1", "RPS7", "RPS8", "RPS9", "RRP12", "RRP15", "RRP8", "RRS1",
  "SBDS", "SCD", "SIRT1", "SPEN", "SURF2", "TAF1B", "TCOF1", "TERF2",
  "TRIM24", "UBTF", "USP47", "UTP11", "UTP14A", "UTP18", "WDR3", "WDR5",
  "WDR75", "YBX1"
)

# Remove duplicates
literature_up_genes <- unique(literature_up_genes)
literature_down_genes <- unique(literature_down_genes)

# Identify conflicting genes
literature_conflict_genes <- intersect(literature_up_genes, literature_down_genes)

# Remove conflicting genes
literature_up_highconf <- setdiff(literature_up_genes, literature_conflict_genes)
literature_down_highconf <- setdiff(literature_down_genes, literature_conflict_genes)

# Construct high-confidence gene table (excluding conflicts)
up_df <- data.frame(
  gene_symbol = literature_up_highconf,
  group = "UP",
  source = "literature",
  stringsAsFactors = FALSE
)

down_df <- data.frame(
  gene_symbol = literature_down_highconf,
  group = "DOWN",
  source = "literature",
  stringsAsFactors = FALSE
)

literature_gene_df <- rbind(up_df, down_df)

# Save final gene table
write.csv(
  literature_gene_df,
  file = "literature_gene_list_highconf.csv",
  row.names = FALSE
)

# Save conflicting genes separately for record
conflict_df <- data.frame(
  gene_symbol = sort(literature_conflict_genes),
  stringsAsFactors = FALSE
)

write.csv(
  conflict_df,
  file = "literature_conflict_genes.csv",
  row.names = FALSE
)

# Print summary statistics
cat("Conflict genes removed:", length(literature_conflict_genes), "\n")
cat("High-confidence UP genes:", length(literature_up_highconf), "\n")
cat("High-confidence DOWN genes:", length(literature_down_highconf), "\n")
cat("Total rows in final table:", nrow(literature_gene_df), "\n")

#########################################################################
########################### Core genes ##################################

# Literature-derived core upregulated genes
literature_core_up_genes <- c(
  "ATM", "ATR", "BAX", "CDKN1A", "CDKN2A", "HIF1A", "MAPK14",
  "MDM2", "MYC", "RPL11", "RPL22", "RPL23", "RPL26", "RPL3",
  "RPL5", "RPS2", "RPS7", "SRSF1", "TP53"
)

# Literature-derived core downregulated genes
literature_core_down_genes <- c(
  "AKT1", "BCL2", "BOP1", "E2F1", "FTO", "GNL3", "GRWD1", "HEATR1",
  "IMPDH2", "MDM2", "MDM4", "MYC", "NEDD8", "NOLC1", "NOP53", "NPM1",
  "PIM1", "POLR1A", "PPAN", "RIOK1", "RPL11", "RPL23", "RPL27A",
  "RPL9", "RPS15A", "RPS19", "RPS6", "RPS7", "RRP15", "SBDS",
  "SIRT1", "TAF1B", "TRIM24", "UBTF"
)

# Remove duplicates
literature_core_up_genes <- unique(literature_core_up_genes)
literature_core_down_genes <- unique(literature_core_down_genes)

# Identify conflicting genes between UP and DOWN groups
literature_core_conflict_genes <- intersect(
  literature_core_up_genes,
  literature_core_down_genes
)

# Remove conflicting genes to obtain high-confidence sets
literature_core_up_highconf <- setdiff(
  literature_core_up_genes,
  literature_core_conflict_genes
)

literature_core_down_highconf <- setdiff(
  literature_core_down_genes,
  literature_core_conflict_genes
)

# Construct final high-confidence table (excluding conflicts)
up_df <- data.frame(
  gene_symbol = sort(literature_core_up_highconf),
  group = "UP",
  dataset = "literature_core",
  stringsAsFactors = FALSE
)

down_df <- data.frame(
  gene_symbol = sort(literature_core_down_highconf),
  group = "DOWN",
  dataset = "literature_core",
  stringsAsFactors = FALSE
)

literature_core_gene_df <- rbind(up_df, down_df)

# Construct conflict gene table
literature_core_conflict_df <- data.frame(
  gene_symbol = sort(literature_core_conflict_genes),
  dataset = "literature_core",
  stringsAsFactors = FALSE
)

# Save CSV files
write.csv(
  literature_core_gene_df,
  file = "literature_core_gene_list_highconf.csv",
  row.names = FALSE
)

write.csv(
  literature_core_conflict_df,
  file = "literature_core_conflict_genes.csv",
  row.names = FALSE
)

# Print summary statistics
cat("Conflict genes removed:", length(literature_core_conflict_genes), "\n")
cat("High-confidence UP genes:", length(literature_core_up_highconf), "\n")
cat("High-confidence DOWN genes:", length(literature_core_down_highconf), "\n")
cat("Total rows in final table:", nrow(literature_core_gene_df), "\n")

cat("\nConflict genes:\n")
print(sort(literature_core_conflict_genes))
