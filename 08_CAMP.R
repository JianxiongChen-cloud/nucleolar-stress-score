###############################################################################
## CMap/LINCS 2020 data processing workflow - GSVA-based NuS analysis
## Author: CJX
## Final workflow:
##   1) read LINCS expression + metadata
##   2) map drug names
##   3) convert gene IDs to symbols
##   4) calculate NuS score by GSVA
##   5) summarize at drug level
##   6) read normality score table
##   7) intersect drugs
##   8) correlation analysis and visualization
###############################################################################

###############################################################################
## 0. Environment
###############################################################################
setwd("/home/xxm_xxm/CJX_workspace/clue_data/")
rm(list = ls())

library(readr)
library(readxl)
library(cmapR)
library(rhdf5)
library(dplyr)
library(stringr)
library(GSVA)
library(pbapply)
library(ggplot2)
library(ggpubr)

###############################################################################
## 1. Read expression matrix
###############################################################################
cat("Reading expression matrix...\n")
gctx_file <- "level5_beta_trt_cp_n720216x12328.gctx"
all_data <- parse_gctx(gctx_file)
expr_raw <- all_data@mat

cat(sprintf("Raw expression matrix dimension: %d genes x %d samples\n",
            nrow(expr_raw), ncol(expr_raw)))
cat(sprintf("Expression value range: [%.2f, %.2f]\n",
            min(expr_raw, na.rm = TRUE), max(expr_raw, na.rm = TRUE)))

###############################################################################
## 2. Read metadata
###############################################################################
cat("Reading sample metadata...\n")
met <- read_delim("instinfo_beta.txt", delim = "\t", show_col_types = FALSE)
save(met, file = "instinfo.Rdata")
cat(sprintf("Metadata dimension: %d rows x %d columns\n", nrow(met), ncol(met)))

cat("Reading gene information...\n")
gene_info <- read_delim("geneinfo_beta.txt", delim = "\t", show_col_types = FALSE)

###############################################################################
## 3. Match sample IDs to drug names
###############################################################################
cat("Matching drug information...\n")

sample_ids <- colnames(expr_raw)
id_extract <- sapply(strsplit(sample_ids, ":"), "[", 2)

sample_info <- data.frame(
  col_id = sample_ids,
  id = id_extract,
  stringsAsFactors = FALSE
)

drug_meta <- data.frame(
  mfcid = met$pert_mfc_id,
  id = met$pert_id,
  name = met$cmap_name,
  stringsAsFactors = FALSE
)

meta_long <- bind_rows(
  drug_meta %>% select(mfcid, name) %>% rename(id_match = mfcid),
  drug_meta %>% select(id, name) %>% rename(id_match = id)
) %>%
  distinct(id_match, .keep_all = TRUE)

sample_info <- sample_info %>%
  left_join(meta_long, by = c("id" = "id_match"))

drug_filtered <- sample_info %>%
  filter(!is.na(name) & !str_detect(name, "^BRD-"))

cat(sprintf("Sample count before filtering: %d\n", nrow(sample_info)))
cat(sprintf("Sample count after filtering: %d\n", nrow(drug_filtered)))
cat(sprintf("Number of unique drugs: %d\n", length(unique(drug_filtered$name))))

col_keep <- colnames(expr_raw) %in% drug_filtered$col_id
expr_filtered <- expr_raw[, col_keep, drop = FALSE]

cat(sprintf("Filtered expression matrix dimension: %d genes x %d samples\n",
            nrow(expr_filtered), ncol(expr_filtered)))

save(expr_filtered, drug_filtered, file = "Step1_filtered_data.Rdata")

###############################################################################
## 4. Convert gene IDs to gene symbols
###############################################################################
cat("Converting gene IDs to gene symbols...\n")

gene_symbols <- gene_info$gene_symbol[match(rownames(expr_filtered), gene_info$gene_id)]

cat(sprintf("Original gene count: %d\n", nrow(expr_filtered)))
cat(sprintf("Successfully matched gene count: %d\n", sum(!is.na(gene_symbols))))
cat(sprintf("Unmatched gene count: %d\n", sum(is.na(gene_symbols))))

rownames(expr_filtered) <- gene_symbols
expr_filtered <- expr_filtered[!is.na(rownames(expr_filtered)) & rownames(expr_filtered) != "", , drop = FALSE]

# Merge duplicated symbols by mean
if (any(duplicated(rownames(expr_filtered)))) {
  cat("Duplicated gene symbols detected; merging by mean...\n")
  expr_filtered <- rowsum(expr_filtered, group = rownames(expr_filtered), reorder = FALSE) /
    as.vector(table(rownames(expr_filtered)))
}

cat(sprintf("Final expression matrix dimension: %d genes x %d samples\n",
            nrow(expr_filtered), ncol(expr_filtered)))

save(expr_filtered, file = "Step2_gene_symbol_matrix.Rdata")

###############################################################################
## 5. Load NuS gene sets
###############################################################################
cat("Loading NuS gene sets...\n")

nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"
loaded_obj_names <- load(nus_rdata)
cat("Objects in Rdata: ", paste(loaded_obj_names, collapse = ", "), "\n")

if (!"geneSets_final" %in% loaded_obj_names) {
  stop("Load failed: object 'geneSets_final' was not found in the Rdata file")
}

if (is.null(geneSets_final$NuStress_UP) || is.null(geneSets_final$NuStress_DOWN)) {
  stop("Load failed: NuStress_UP or NuStress_DOWN was not found in 'geneSets_final'")
}

geneSets_all <- list(
  NuStress_UP   = unique(as.character(geneSets_final$NuStress_UP)),
  NuStress_DOWN = unique(as.character(geneSets_final$NuStress_DOWN))
)

geneSets_all <- lapply(geneSets_all, function(x) {
  x <- x[!is.na(x)]
  x <- x[x != ""]
  unique(x)
})

cat(sprintf("Original NuStress_UP gene count: %d\n", length(geneSets_all$NuStress_UP)))
cat(sprintf("Original NuStress_DOWN gene count: %d\n", length(geneSets_all$NuStress_DOWN)))

# Intersect with the expression matrix
geneSets_use <- lapply(geneSets_all, function(gs) intersect(gs, rownames(expr_filtered)))

cat(sprintf("NuStress_UP matched gene count: %d\n", length(geneSets_use$NuStress_UP)))
cat(sprintf("NuStress_DOWN matched gene count: %d\n", length(geneSets_use$NuStress_DOWN)))

if (length(geneSets_use$NuStress_UP) < 5 || length(geneSets_use$NuStress_DOWN) < 5) {
  stop("Too few genes overlap between gene sets and the expression matrix for stable GSVA")
}

###############################################################################
## 6. GSVA-based NuS scoring
###############################################################################
cat("Starting GSVA-based NuS scoring...\n")

expr_matrix <- as.matrix(expr_filtered)
n_samples <- ncol(expr_matrix)
cat(sprintf("Total sample count: %d\n", n_samples))

chunk_size <- 1000
chunks <- split(seq_len(n_samples), ceiling(seq_len(n_samples) / chunk_size))

process_chunk <- function(col_indices) {
  expr_chunk <- expr_matrix[, col_indices, drop = FALSE]
  
  gsva_param <- GSVA::gsvaParam(
    exprData = expr_chunk,
    geneSets = geneSets_use,
    minSize = 5,
    maxSize = 500,
    kcdf = "Gaussian"
  )
  
  gsva_res <- gsva(gsva_param)
  
  data.frame(
    Sample = colnames(expr_chunk),
    NuS_up_score = as.numeric(gsva_res["NuStress_UP", ]),
    NuS_down_score = as.numeric(gsva_res["NuStress_DOWN", ]),
    NuS_score = as.numeric(gsva_res["NuStress_UP", ] - gsva_res["NuStress_DOWN", ]),
    stringsAsFactors = FALSE
  )
}

cat("Calculating GSVA scores by chunks...\n")
nus_score_list <- pbapply::pblapply(chunks, process_chunk)
nus_score_df <- do.call(rbind, nus_score_list)

cat(sprintf("Successfully obtained NuS scores for %d samples\n", nrow(nus_score_df)))
save(nus_score_df, file = "Step3_nus_score_by_GSVA.Rdata")

###############################################################################
## 7. Merge sample-level NuS scores with drug info
###############################################################################
drug_meta_with_score <- merge(
  drug_filtered,
  nus_score_df,
  by.x = "col_id",
  by.y = "Sample",
  all.x = TRUE
)

cat(sprintf("Merged data dimension: %d x %d\n",
            nrow(drug_meta_with_score), ncol(drug_meta_with_score)))

write.csv(drug_meta_with_score, "drug_meta_with_GSVA_NuS_score.csv", row.names = FALSE)
save(expr_filtered, drug_meta_with_score, file = "Step4_final_results.Rdata")

###############################################################################
## 8. Summarize at drug level
###############################################################################
cat("Summarizing drug-level NuS scores...\n")

drug_summary <- drug_meta_with_score %>%
  group_by(name) %>%
  summarise(
    mean_score   = mean(NuS_score, na.rm = TRUE),
    median_score = median(NuS_score, na.rm = TRUE),
    max_score    = max(NuS_score, na.rm = TRUE),
    min_score    = min(NuS_score, na.rm = TRUE),
    sd_score     = sd(NuS_score, na.rm = TRUE),
    n_samples    = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))

write.csv(drug_summary, "drug_summary_GSVA_NuS_score.csv", row.names = FALSE)

cat("\nTop 10 drugs by mean NuS_score:\n")
print(head(drug_summary, 10))

###############################################################################
## 9. Basic plots for NuS screening
###############################################################################
cat("Generating NuS drug screening plots...\n")

p1 <- ggplot(drug_summary, aes(x = mean_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.8) +
  theme_minimal() +
  labs(
    x = "Mean NuS score per drug",
    y = "Count",
    title = "Distribution of drug-induced NuS scores"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

top_drugs_mean <- drug_summary %>%
  slice_max(order_by = mean_score, n = 20)

p2 <- ggplot(top_drugs_mean, aes(x = reorder(name, mean_score), y = mean_score)) +
  geom_col(fill = "#7DC18C", alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Drug",
    y = "Mean NuS score",
    title = "Top 20 drugs by mean NuS score"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

top_drugs_max <- drug_summary %>%
  slice_max(order_by = max_score, n = 20)

p3 <- ggplot(top_drugs_max, aes(x = reorder(name, max_score), y = max_score)) +
  geom_col(fill = "#47957F", alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Drug",
    y = "Max NuS score",
    title = "Top 20 drugs by max NuS score"
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

ggsave("figure1_drug_score_distribution.pdf", p1, width = 8, height = 4, bg = "white")
ggsave("figure2_top20_mean_score.pdf", p2, width = 4, height = 8, bg = "white")
ggsave("figure3_top20_max_score.pdf", p3, width = 4, height = 8, bg = "white")

###############################################################################
## 10. Read normality score Excel
###############################################################################
cat("Reading normality score table...\n")
library(readxl)
# Change this to the actual path of your Excel file
normality_xlsx <- "/home/xxm_xxm/CJX_workspace/clue_data/517150_file02.xlsx"

norm_raw <- read_excel(normality_xlsx)
cat("Excel column names:\n")
print(colnames(norm_raw))

norm_df <- norm_raw %>%
  dplyr::select(
    Product_Name = `Product Name`,
    Cat_No = `Cat.#`,
    normality_10uM = `10um normality score, normalized`,
    normality_1uM  = `1um normality score, normalized`,
    Target = Target,
    Pathway = Pathway
  ) %>%
  filter(!is.na(Product_Name))

# Remove group headers and other non-drug rows
norm_df <- norm_df %>%
  filter(!grepl("^\\s*1 and 10uM hits\\s*$", Product_Name, ignore.case = TRUE))

###############################################################################
## 11. Clean drug names and intersect
###############################################################################
clean_drug_name <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("hydrochloride| hcl$|hcl| sodium$|potassium|mesylate|maleate|tartrate|phosphate|acetate|bromide", "") %>%
    str_replace_all("\\(.*?\\)", "") %>%
    str_replace_all("hydrobromide", "") %>%
    str_replace_all("dihydrochloride", "") %>%
    str_replace_all("monohydrate", "") %>%
    str_replace_all("[^a-z0-9]", "") %>%
    str_trim()
}

norm_df <- norm_df %>%
  mutate(drug_clean = clean_drug_name(Product_Name))

nus_df <- drug_summary %>%
  dplyr::select(
    drug_name = name,
    NuS_score = mean_score,
    median_score,
    max_score,
    sd_score,
    n_samples
  ) %>%
  filter(!is.na(drug_name)) %>%
  mutate(drug_clean = clean_drug_name(drug_name))

merged_df <- inner_join(nus_df, norm_df, by = "drug_clean")

cat(sprintf("Number of intersected drugs: %d\n", nrow(merged_df)))

write.csv(merged_df, "NuS_normality_intersection.csv", row.names = FALSE)

###############################################################################
## 12. Correlation analysis
###############################################################################
cat("Running correlation analysis between NuS and normality scores...\n")

# All intersected drugs
cor_10 <- cor.test(
  merged_df$NuS_score,
  merged_df$normality_10uM,
  method = "spearman",
  use = "complete.obs"
)

cor_1 <- cor.test(
  merged_df$NuS_score,
  merged_df$normality_1uM,
  method = "spearman",
  use = "complete.obs"
)

cat("\n===== Spearman correlation: NuS vs 10uM normality =====\n")
print(cor_10)

cat("\n===== Spearman correlation: NuS vs 1uM normality =====\n")
print(cor_1)

# More robust subset with n_samples >= 3
merged_df_n3 <- merged_df %>%
  filter(n_samples >= 3)

cat(sprintf("\nNumber of intersected drugs with n_samples >= 3: %d\n", nrow(merged_df_n3)))

if (nrow(merged_df_n3) >= 5) {
  cor_10_n3 <- cor.test(
    merged_df_n3$NuS_score,
    merged_df_n3$normality_10uM,
    method = "spearman",
    use = "complete.obs"
  )
  
  cor_1_n3 <- cor.test(
    merged_df_n3$NuS_score,
    merged_df_n3$normality_1uM,
    method = "spearman",
    use = "complete.obs"
  )
  
  cat("\n===== n_samples >= 3: NuS vs 10uM =====\n")
  print(cor_10_n3)
  
  cat("\n===== n_samples >= 3: NuS vs 1uM =====\n")
  print(cor_1_n3)
}

###############################################################################
## 13. Visualization of correlation
###############################################################################
cat("Generating correlation plots...\n")

p_cor_10 <- ggplot(merged_df, aes(x = NuS_score, y = normality_10uM)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(
    x = "Mean NuS score (GSVA)",
    y = "10uM normality score",
    title = "Correlation between NuS score and 10uM normality score"
  )

p_cor_1 <- ggplot(merged_df, aes(x = NuS_score, y = normality_1uM)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  ggpubr::stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_classic() +
  labs(
    x = "Mean NuS score (GSVA)",
    y = "1uM normality score",
    title = "Correlation between NuS score and 1uM normality score"
  )

ggsave("figure4_NuS_vs_10uM_normality.pdf", p_cor_10, width = 6, height = 5, bg = "white")
ggsave("figure5_NuS_vs_1uM_normality.pdf", p_cor_1, width = 6, height = 5, bg = "white")

final_tbl_pub <- merged_df %>%
  transmute(
    Compound = drug_name,
    NuS = NuS_score,
    Normality_score_10uM = normality_10uM,
    Normality_score_1uM = normality_1uM,
    NuS_median = median_score,
    NuS_max = max_score,
    NuS_sd = sd_score,
    Profile_number = n_samples,
    Target = Target,
    Pathway = Pathway
  ) %>%
  arrange(desc(NuS))

write.csv(final_tbl_pub, "Supplementary_Table_NuS_normality_shared_compounds.csv", row.names = FALSE)

print(head(final_tbl_pub, 20))


###############################################################################
## 14. Optional: highlight key drugs
###############################################################################
key_pattern <- "cx-5461|actinomycin|bmh-21|oxaliplatin|cisplatin|carboplatin"

key_drugs_df <- merged_df %>%
  filter(grepl(key_pattern, Product_Name, ignore.case = TRUE) |
           grepl(key_pattern, drug_name, ignore.case = TRUE))

write.csv(key_drugs_df, "key_drugs_in_intersection.csv", row.names = FALSE)

###############################################################################
## 15. Save summary
###############################################################################
cat(paste0("\n", strrep("=", 70), "\n"))
cat("Final analysis completed. Result summary:\n")
cat(paste0(strrep("-", 70), "\n"))
cat(sprintf("1. Raw expression matrix: %d genes x %d samples\n", nrow(expr_raw), ncol(expr_raw)))
cat(sprintf("2. Filtered expression matrix: %d genes x %d samples\n", nrow(expr_filtered), ncol(expr_filtered)))
cat(sprintf("3. NuStress_UP matched gene count: %d\n", length(geneSets_use$NuStress_UP)))
cat(sprintf("4. NuStress_DOWN matched gene count: %d\n", length(geneSets_use$NuStress_DOWN)))
cat(sprintf("5. Samples with NuS scores: %d\n", nrow(nus_score_df)))
cat(sprintf("6. Number of unique drugs: %d\n", nrow(drug_summary)))
cat(sprintf("7. Drugs intersecting with normality score: %d\n", nrow(merged_df)))
cat(sprintf("8. NuS score range: [%.3f, %.3f]\n",
            min(drug_summary$mean_score, na.rm = TRUE),
            max(drug_summary$mean_score, na.rm = TRUE)))
cat("9. Output figures: figure1-5.pdf\n")
cat("10. Output tables: drug_summary_GSVA_NuS_score.csv, NuS_normality_intersection.csv\n")
cat(paste0(strrep("=", 70), "\n"))

###############################################################################
## 16. Keep useful objects
###############################################################################
rm(list = setdiff(ls(), c(
  "expr_filtered",
  "drug_meta_with_score",
  "drug_summary",
  "merged_df",
  "merged_df_n3",
  "p1", "p2", "p3", "p_cor_10", "p_cor_1"
)))

cat("GSVA-based NuS analysis workflow finished.\n")
