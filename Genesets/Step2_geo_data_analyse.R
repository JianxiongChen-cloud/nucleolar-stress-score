rm(list = ls())
options(stringsAsFactors = FALSE)
setwd("/home/xxm_xxm/CJX_workspace/geneset")
############################################################
##Step2 ### 0. 路径设置
############################################################
input_dir <- "/home/xxm_xxm/CJX_workspace/Dataset"
output_base_dir <- getwd()

cat("Input directory:", input_dir, "\n")
cat("Output base directory:", output_base_dir, "\n")

summary_dir <- file.path(output_base_dir, "summary_results")
if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)

final_dir <- file.path(output_base_dir, "final_nucleolar_gene_sets")
if (!dir.exists(final_dir)) dir.create(final_dir, recursive = TRUE)

############################################################
## 1. 数据库文件信息
## subdir 相对于 input_dir
############################################################
file_info <- data.frame(
  dataset = c(
    "GSE6400",
    "GSE12666",
    "GSE33417_ips",
    "GSE33417_hff",
    "GSE62593",
    "GSE62963",
    "GSE118565",
    "GSE198178",
    "GSE255898",
    "GSE261563_603",
    "GSE261563_543",
    "GSE267499",
    "GSE267501",
    "GSE282212",
    "GSE282214",
    "GSE298220_BMH21",
    "GSE298220_CX5461"
  ),
  subdir = c(
    "GSE6400",
    "GSE12666",
    "GSE33417",
    "GSE33417",
    "GSE62593",
    "GSE62963",
    "GSE118565",
    "GSE198178",
    "GSE255898",
    "GSE261563",
    "GSE261563",
    "GSE267499",
    "GSE267501",
    "GSE282212",
    "GSE282214",
    "GSE298220",
    "GSE298220"
  ),
  file = c(
    "GSE6400_diffgenesymbol_long.csv",
    "GSE12666_diffgenesymbol_long.csv",
    "GSE33417_ips_diffgenesymbol_long.csv",
    "GSE33417_hff_diffgenesymbol_long.csv",
    "GSE62593_diffgenesymbol_long.csv",
    "GSE62963_common_strict_DEGs_diffgenesymbol_long.csv",
    "GSE118565_diffgenesymbol_long.csv",
    "GSE198178_common_4h_2h_1h_diffgenesymbol_long.csv",
    "GSE255898_diffgenesymbol_long.csv",
    "GSE261563_603_diffgenesymbol_long.csv",
    "GSE261563_543_diffgenesymbol_long.csv",
    "GSE267499_diffgenesymbol_long.csv",
    "GSE267501_diffgenesymbol_long.csv",
    "GSE282212_diffgenesymbol_long.csv",
    "GSE282214_diffgenesymbol_long.csv",
    "GSE298220_BMH21_diffgenesymbol_long.csv",
    "GSE298220_CX5461_diffgenesymbol_long.csv"
  ),
  stringsAsFactors = FALSE
)

############################################################
## 2. 读取并清洗单个数据集
############################################################
read_one_dataset <- function(file_path, dataset_name) {
  x <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!("group" %in% colnames(x))) {
    stop(paste("Missing 'group' column in:", file_path))
  }
  if (!("gene_symbol" %in% colnames(x))) {
    stop(paste("Missing 'gene_symbol' column in:", file_path))
  }
  
  x <- x[, c("gene_symbol", "group")]
  x$gene_symbol <- toupper(trimws(x$gene_symbol))
  x$group <- toupper(trimws(x$group))
  
  x <- x[!is.na(x$gene_symbol) & x$gene_symbol != "", ]
  x <- x[!is.na(x$group) & x$group %in% c("UP", "DOWN"), ]
  x <- unique(x)
  
  gene_group_n <- aggregate(group ~ gene_symbol, data = x, FUN = function(z) length(unique(z)))
  inner_conflict_genes <- gene_group_n$gene_symbol[gene_group_n$group > 1]
  
  inner_conflict_df <- data.frame(
    dataset = rep(dataset_name, length(inner_conflict_genes)),
    gene_symbol = inner_conflict_genes,
    stringsAsFactors = FALSE
  )
  
  if (length(inner_conflict_genes) > 0) {
    x <- x[!(x$gene_symbol %in% inner_conflict_genes), ]
  }
  
  x$dataset <- dataset_name
  x <- x[, c("dataset", "group", "gene_symbol")]
  
  return(list(
    data = x,
    inner_conflict = inner_conflict_df
  ))
}

############################################################
## 3. 读取所有数据库文件
############################################################
all_data_list <- list()
all_inner_conflicts <- list()

for (i in seq_len(nrow(file_info))) {
  dataset_name <- file_info$dataset[i]
  file_path <- file.path(input_dir, file_info$subdir[i], file_info$file[i])
  
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  cat("Reading:", dataset_name, "->", file_path, "\n")
  res <- read_one_dataset(file_path, dataset_name)
  
  all_data_list[[dataset_name]] <- res$data
  all_inner_conflicts[[dataset_name]] <- res$inner_conflict
}

all_long <- do.call(rbind, all_data_list)
all_long <- unique(all_long)

inner_conflict_table <- do.call(rbind, all_inner_conflicts)
if (is.null(inner_conflict_table) || nrow(inner_conflict_table) == 0) {
  inner_conflict_table <- data.frame(
    dataset = character(0),
    gene_symbol = character(0),
    stringsAsFactors = FALSE
  )
}

write.csv(
  all_long,
  file.path(summary_dir, "all_datasets_diffgenesymbol_long_raw.csv"),
  row.names = FALSE
)

write.csv(
  inner_conflict_table,
  file.path(summary_dir, "inner_conflict_genes_within_dataset.csv"),
  row.names = FALSE
)

cat("\nStep 3 finished.\n")
cat("Raw long table rows:", nrow(all_long), "\n")
cat("Inner conflict genes:", nrow(inner_conflict_table), "\n\n")

############################################################
## 4. 原始频率统计
############################################################
gene_freq_raw <- aggregate(
  dataset ~ gene_symbol + group,
  data = all_long,
  FUN = function(x) length(unique(x))
)
colnames(gene_freq_raw)[3] <- "freq"

gene_freq_raw <- gene_freq_raw[order(
  gene_freq_raw$group,
  -gene_freq_raw$freq,
  gene_freq_raw$gene_symbol
), ]

write.csv(
  gene_freq_raw,
  file.path(summary_dir, "gene_up_down_frequency_raw.csv"),
  row.names = FALSE
)

############################################################
## 5. 跨数据库冲突基因
############################################################
gene_trend_summary <- aggregate(
  group ~ gene_symbol,
  data = all_long,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)
colnames(gene_trend_summary)[2] <- "trend_type"

conflict_genes <- gene_trend_summary$gene_symbol[gene_trend_summary$trend_type == "DOWN;UP"]

conflict_genes_df <- data.frame(
  gene_symbol = sort(unique(conflict_genes)),
  stringsAsFactors = FALSE
)

write.csv(
  conflict_genes_df,
  file.path(summary_dir, "conflict_genes_across_datasets.csv"),
  row.names = FALSE
)

############################################################
## 6. 去掉跨数据库冲突基因，保留一致趋势基因
############################################################
all_long_consistent <- all_long[!(all_long$gene_symbol %in% conflict_genes), ]
all_long_consistent <- unique(all_long_consistent)

write.csv(
  all_long_consistent,
  file.path(summary_dir, "all_datasets_diffgenesymbol_long_consistent.csv"),
  row.names = FALSE
)

############################################################
## 7. 一致趋势基因频率统计
############################################################
gene_freq_consistent <- aggregate(
  dataset ~ gene_symbol + group,
  data = all_long_consistent,
  FUN = function(x) length(unique(x))
)
colnames(gene_freq_consistent)[3] <- "freq"

gene_freq_consistent <- gene_freq_consistent[order(
  gene_freq_consistent$group,
  -gene_freq_consistent$freq,
  gene_freq_consistent$gene_symbol
), ]

write.csv(
  gene_freq_consistent,
  file.path(summary_dir, "gene_up_down_frequency_consistent.csv"),
  row.names = FALSE
)

############################################################
##补充. 分别给上调基因和下调基因绘制频率分布直方图
############################################################
up_freq <- subset(gene_freq_consistent, group == "UP")
down_freq <- subset(gene_freq_consistent, group == "DOWN")

# 统一横坐标范围，方便比较
max_freq <- max(gene_freq_consistent$freq, na.rm = TRUE)
breaks_seq <- seq(0.5, max_freq + 0.5, by = 1)

# 上调基因频率直方图
pdf(
  file.path(summary_dir, "histogram_up_gene_frequency.pdf"),
  width = 8, height = 5
)
hist(
  up_freq$freq,
  breaks = breaks_seq,
  main = "Frequency distribution of upregulated genes",
  xlab = "Frequency across datasets",
  ylab = "Number of genes",
  xaxt = "n"
)
axis(1, at = 1:max_freq)
dev.off()

# 下调基因频率直方图
pdf(
  file.path(summary_dir, "histogram_down_gene_frequency.pdf"),
  width = 8, height = 5
)
hist(
  down_freq$freq,
  breaks = breaks_seq,
  main = "Frequency distribution of downregulated genes",
  xlab = "Frequency across datasets",
  ylab = "Number of genes",
  xaxt = "n"
)
axis(1, at = 1:max_freq)
dev.off()


############################################################
## 补充. 更适合论文展示的离散频率柱状图
############################################################
up_tab <- table(factor(up_freq$freq, levels = 1:max_freq))
down_tab <- table(factor(down_freq$freq, levels = 1:max_freq))

pdf(
  file.path(summary_dir, "barplot_up_gene_frequency.pdf"),
  width = 8, height = 5
)
barplot(
  up_tab,
  main = "Upregulated genes",
  xlab = "Frequency across datasets",
  ylab = "Number of genes"
)
dev.off()

pdf(
  file.path(summary_dir, "barplot_down_gene_frequency.pdf"),
  width = 8, height = 5
)
barplot(
  down_tab,
  main = "Downregulated genes",
  xlab = "Frequency across datasets",
  ylab = "Number of genes"
)
dev.off()



############################################################
## 8. 一致上调/下调总表与统计
############################################################
final_up_genes <- sort(unique(gene_freq_consistent$gene_symbol[gene_freq_consistent$group == "UP"]))
final_down_genes <- sort(unique(gene_freq_consistent$gene_symbol[gene_freq_consistent$group == "DOWN"]))

write.csv(
  data.frame(gene_symbol = final_up_genes, stringsAsFactors = FALSE),
  file.path(summary_dir, "final_consistent_up_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = final_down_genes, stringsAsFactors = FALSE),
  file.path(summary_dir, "final_consistent_down_genes.csv"),
  row.names = FALSE
)

dataset_group_count <- aggregate(
  gene_symbol ~ dataset + group,
  data = all_long_consistent,
  FUN = function(x) length(unique(x))
)
colnames(dataset_group_count)[3] <- "gene_count"

dataset_group_count <- dataset_group_count[order(dataset_group_count$dataset, dataset_group_count$group), ]

write.csv(
  dataset_group_count,
  file.path(summary_dir, "dataset_up_down_gene_count_consistent.csv"),
  row.names = FALSE
)

cat("Step 8 finished.\n")
cat("Consistent UP genes:", length(final_up_genes), "\n")
cat("Consistent DOWN genes:", length(final_down_genes), "\n")
cat("Cross-dataset conflict genes removed:", length(conflict_genes), "\n\n")

############################################################
## 9. 数据库来源高信度基因：至少 4 个数据库
############################################################
freq_cutoff_highconf <- 4

gene_freq_ge4 <- gene_freq_consistent[gene_freq_consistent$freq >= freq_cutoff_highconf, ]
gene_freq_ge4 <- gene_freq_ge4[order(
  gene_freq_ge4$group,
  -gene_freq_ge4$freq,
  gene_freq_ge4$gene_symbol
), ]

database_up_highconf <- sort(unique(gene_freq_ge4$gene_symbol[gene_freq_ge4$group == "UP"]))
database_down_highconf <- sort(unique(gene_freq_ge4$gene_symbol[gene_freq_ge4$group == "DOWN"]))

write.csv(
  gene_freq_ge4,
  file.path(summary_dir, "database_highconf_freq_ge4_detail.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = database_up_highconf, stringsAsFactors = FALSE),
  file.path(summary_dir, "database_up_highconf_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = database_down_highconf, stringsAsFactors = FALSE),
  file.path(summary_dir, "database_down_highconf_genes.csv"),
  row.names = FALSE
)

cat("Step 9 finished.\n")
cat("Database highconf UP (>=4 datasets):", length(database_up_highconf), "\n")
cat("Database highconf DOWN (>=4 datasets):", length(database_down_highconf), "\n\n")

############################################################
## 10. 数据集-处理方式注释表
############################################################
dataset_treatment_info <- data.frame(
  dataset = c(
    "GSE6400",
    "GSE12666",
    "GSE33417_ips",
    "GSE33417_hff",
    "GSE62593",
    "GSE62963",
    "GSE118565",
    "GSE198178",
    "GSE255898",
    "GSE261563_603",
    "GSE261563_543",
    "GSE267499",
    "GSE267501",
    "GSE282212",
    "GSE282214",
    "GSE298220_BMH21",
    "GSE298220_CX5461"
  ),
  treatment = c(
    "ActD",
    "BMH21",
    "ActD",
    "ActD",
    "ActD",
    "ActD",
    "CX5461",
    "ActD",
    "CX5461",
    "CX5461",
    "CX5461",
    "BMH21",
    "BMH21",
    "BMH21",
    "BMH21",
    "BMH21",
    "CX5461"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  dataset_treatment_info,
  file.path(summary_dir, "dataset_treatment_info.csv"),
  row.names = FALSE
)

all_long_consistent2 <- merge(
  all_long_consistent,
  dataset_treatment_info,
  by = "dataset",
  all.x = TRUE
)

if (any(is.na(all_long_consistent2$treatment))) {
  warning("Some datasets do not have treatment annotation.")
  print(unique(all_long_consistent2$dataset[is.na(all_long_consistent2$treatment)]))
}

write.csv(
  all_long_consistent2,
  file.path(summary_dir, "all_datasets_diffgenesymbol_long_consistent_with_treatment.csv"),
  row.names = FALSE
)

############################################################
## 11. 数据库来源核心基因
## 条件：至少 5 个数据库 + 至少 2 种不同抑制剂
############################################################
freq_cutoff_core <- 5
treatment_cutoff_core <- 2

db_core_candidate <- gene_freq_consistent[gene_freq_consistent$freq >= freq_cutoff_core, ]

db_core_support <- aggregate(
  cbind(dataset, treatment) ~ gene_symbol + group,
  data = all_long_consistent2,
  FUN = function(x) length(unique(x))
)
colnames(db_core_support) <- c("gene_symbol", "group", "dataset_n", "treatment_n")

db_core_detail <- merge(
  db_core_candidate,
  db_core_support,
  by = c("gene_symbol", "group"),
  all.x = TRUE
)

db_core_detail <- db_core_detail[
  db_core_detail$freq >= freq_cutoff_core &
    db_core_detail$treatment_n >= treatment_cutoff_core,
]

db_core_detail <- db_core_detail[order(
  db_core_detail$group,
  -db_core_detail$freq,
  -db_core_detail$treatment_n,
  db_core_detail$gene_symbol
), ]

database_core_up_genes <- sort(unique(db_core_detail$gene_symbol[db_core_detail$group == "UP"]))
database_core_down_genes <- sort(unique(db_core_detail$gene_symbol[db_core_detail$group == "DOWN"]))

write.csv(
  db_core_detail,
  file.path(summary_dir, "database_core_gene_candidates_detail.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = database_core_up_genes, stringsAsFactors = FALSE),
  file.path(summary_dir, "database_core_up_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = database_core_down_genes, stringsAsFactors = FALSE),
  file.path(summary_dir, "database_core_down_genes.csv"),
  row.names = FALSE
)

cat("Step 11 finished.\n")
cat("Database core UP (>=5 datasets & >=2 treatments):", length(database_core_up_genes), "\n")
cat("Database core DOWN (>=5 datasets & >=2 treatments):", length(database_core_down_genes), "\n\n")

############################################################
## 12. 数据库核心基因的支持表
############################################################
get_gene_support_table <- function(gene_vec, group_name, dat) {
  if (length(gene_vec) == 0) return(data.frame())
  
  sub <- dat[dat$gene_symbol %in% gene_vec & dat$group == group_name, ]
  sub <- unique(sub[, c("gene_symbol", "group", "dataset", "treatment")])
  sub <- sub[order(sub$gene_symbol, sub$dataset), ]
  rownames(sub) <- NULL
  return(sub)
}

database_core_up_support <- get_gene_support_table(database_core_up_genes, "UP", all_long_consistent2)
database_core_down_support <- get_gene_support_table(database_core_down_genes, "DOWN", all_long_consistent2)

write.csv(
  database_core_up_support,
  file.path(summary_dir, "database_core_up_support_table.csv"),
  row.names = FALSE
)

write.csv(
  database_core_down_support,
  file.path(summary_dir, "database_core_down_support_table.csv"),
  row.names = FALSE
)

############################################################
## 13. 读取文献来源高信度基因和核心基因
## 这两个文件在当前工作目录 geneset 下
############################################################
literature_highconf_file <- file.path(output_base_dir, "literature_gene_list_highconf.csv")
literature_core_file <- file.path(output_base_dir, "literature_core_gene_list_highconf.csv")

if (!file.exists(literature_highconf_file)) {
  stop("File not found: ", literature_highconf_file)
}
if (!file.exists(literature_core_file)) {
  stop("File not found: ", literature_core_file)
}

literature_gene_df <- read.csv(literature_highconf_file, stringsAsFactors = FALSE, check.names = FALSE)
literature_core_gene_df <- read.csv(literature_core_file, stringsAsFactors = FALSE, check.names = FALSE)

literature_gene_df$gene_symbol <- toupper(trimws(literature_gene_df$gene_symbol))
literature_gene_df$group <- toupper(trimws(literature_gene_df$group))

literature_core_gene_df$gene_symbol <- toupper(trimws(literature_core_gene_df$gene_symbol))
literature_core_gene_df$group <- toupper(trimws(literature_core_gene_df$group))

literature_up_highconf <- sort(unique(literature_gene_df$gene_symbol[literature_gene_df$group == "UP"]))
literature_down_highconf <- sort(unique(literature_gene_df$gene_symbol[literature_gene_df$group == "DOWN"]))

literature_core_up_highconf <- sort(unique(literature_core_gene_df$gene_symbol[literature_core_gene_df$group == "UP"]))
literature_core_down_highconf <- sort(unique(literature_core_gene_df$gene_symbol[literature_core_gene_df$group == "DOWN"]))

cat("Literature highconf UP:", length(literature_up_highconf), "\n")
cat("Literature highconf DOWN:", length(literature_down_highconf), "\n")
cat("Literature core UP:", length(literature_core_up_highconf), "\n")
cat("Literature core DOWN:", length(literature_core_down_highconf), "\n\n")

############################################################
## 14. 合并函数：文献 + 数据库，再去除冲突
############################################################
merge_gene_sets_clean <- function(up1, down1, up2, down2,
                                  source1 = "literature",
                                  source2 = "database") {
  up1 <- unique(toupper(trimws(up1)))
  down1 <- unique(toupper(trimws(down1)))
  up2 <- unique(toupper(trimws(up2)))
  down2 <- unique(toupper(trimws(down2)))
  
  long1 <- rbind(
    data.frame(source = source1, group = "UP", gene_symbol = up1, stringsAsFactors = FALSE),
    data.frame(source = source1, group = "DOWN", gene_symbol = down1, stringsAsFactors = FALSE)
  )
  
  long2 <- rbind(
    data.frame(source = source2, group = "UP", gene_symbol = up2, stringsAsFactors = FALSE),
    data.frame(source = source2, group = "DOWN", gene_symbol = down2, stringsAsFactors = FALSE)
  )
  
  merged_raw <- unique(rbind(long1, long2))
  
  trend_summary <- aggregate(
    group ~ gene_symbol,
    data = merged_raw,
    FUN = function(x) paste(sort(unique(x)), collapse = ";")
  )
  colnames(trend_summary)[2] <- "trend_type"
  
  conflict_genes <- sort(unique(trend_summary$gene_symbol[trend_summary$trend_type == "DOWN;UP"]))
  
  merged_final <- merged_raw[!(merged_raw$gene_symbol %in% conflict_genes), ]
  merged_final <- unique(merged_final)
  
  final_up <- sort(unique(merged_final$gene_symbol[merged_final$group == "UP"]))
  final_down <- sort(unique(merged_final$gene_symbol[merged_final$group == "DOWN"]))
  
  list(
    merged_raw = merged_raw,
    merged_final = merged_final,
    conflict_genes = conflict_genes,
    final_up = final_up,
    final_down = final_down
  )
}

############################################################
## 15. 生成核仁应激基因总集 NuStress
############################################################
nustress_res <- merge_gene_sets_clean(
  up1 = literature_up_highconf,
  down1 = literature_down_highconf,
  up2 = database_up_highconf,
  down2 = database_down_highconf,
  source1 = "literature_highconf",
  source2 = "database_highconf"
)

write.csv(
  nustress_res$merged_raw,
  file.path(final_dir, "nucleolar_stress_gene_set_merged_raw.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_res$conflict_genes, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_gene_set_conflict_genes.csv"),
  row.names = FALSE
)

write.csv(
  nustress_res$merged_final,
  file.path(final_dir, "nucleolar_stress_gene_set_long_table.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_res$final_up, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_up_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_res$final_down, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_down_genes.csv"),
  row.names = FALSE
)

cat("Step 15 finished.\n")
cat("NuStress conflict genes removed:", length(nustress_res$conflict_genes), "\n")
cat("NuStress final UP:", length(nustress_res$final_up), "\n")
cat("NuStress final DOWN:", length(nustress_res$final_down), "\n\n")


############################################################
## 16. 生成核仁应激核心基因集 NuStressCore
############################################################
nustress_core_res <- merge_gene_sets_clean(
  up1 = literature_core_up_highconf,
  down1 = literature_core_down_highconf,
  up2 = database_core_up_genes,
  down2 = database_core_down_genes,
  source1 = "literature_core",
  source2 = "database_core"
)

write.csv(
  nustress_core_res$merged_raw,
  file.path(final_dir, "nucleolar_stress_core_gene_set_merged_raw.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_core_res$conflict_genes, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_core_gene_set_conflict_genes.csv"),
  row.names = FALSE
)

write.csv(
  nustress_core_res$merged_final,
  file.path(final_dir, "nucleolar_stress_core_gene_set_long_table.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_core_res$final_up, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_core_up_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = nustress_core_res$final_down, stringsAsFactors = FALSE),
  file.path(final_dir, "nucleolar_stress_core_down_genes.csv"),
  row.names = FALSE
)

cat("Step 16 finished.\n")
cat("NuStressCore conflict genes removed:", length(nustress_core_res$conflict_genes), "\n")
cat("NuStressCore final UP:", length(nustress_core_res$final_up), "\n")
cat("NuStressCore final DOWN:", length(nustress_core_res$final_down), "\n\n")


############################################################
## 补充. 绘制 NuStress 和 NuStressCore 的 UpSet 图
## 目的：展示“最终保留基因集”中各基因的来源构成
############################################################
library(UpSetR)

# 将基因列表转换为 UpSetR 所需的二进制矩阵
make_upset_input <- function(gene_lists) {
  gene_lists <- lapply(gene_lists, function(x) sort(unique(na.omit(x))))
  all_genes <- sort(unique(unlist(gene_lists)))
  binary_matrix <- sapply(gene_lists, function(x) as.integer(all_genes %in% x))
  binary_matrix <- as.data.frame(binary_matrix)
  rownames(binary_matrix) <- all_genes
  return(binary_matrix)
}

# 生成来源长表，便于检查每个最终基因来自哪些来源
build_source_long_table <- function(final_up, final_down,
                                    lit_up, lit_down,
                                    db_up, db_down) {
  final_up <- sort(unique(final_up))
  final_down <- sort(unique(final_down))
  
  up_df <- data.frame(
    gene_symbol = final_up,
    final_group = "UP",
    literature_up = as.integer(final_up %in% lit_up),
    literature_down = 0,
    database_up = as.integer(final_up %in% db_up),
    database_down = 0,
    stringsAsFactors = FALSE
  )
  
  down_df <- data.frame(
    gene_symbol = final_down,
    final_group = "DOWN",
    literature_up = 0,
    literature_down = as.integer(final_down %in% lit_down),
    database_up = 0,
    database_down = as.integer(final_down %in% db_down),
    stringsAsFactors = FALSE
  )
  
  out <- rbind(up_df, down_df)
  
  source_label <- apply(
    out[, c("literature_up", "literature_down", "database_up", "database_down")],
    1,
    function(z) {
      labs <- c("Literature_Up", "Literature_Down", "Database_Up", "Database_Down")[z == 1]
      paste(labs, collapse = ";")
    }
  )
  
  out$source_pattern <- source_label
  out <- out[order(out$final_group, out$gene_symbol), ]
  rownames(out) <- NULL
  return(out)
}

############################################################
## 补充.1 NuStress 总集：仅基于最终保留基因集绘制来源 UpSet
############################################################
nustress_gene_lists_finalonly <- list(
  Literature_Up = intersect(literature_up_highconf, nustress_res$final_up),
  Literature_Down = intersect(literature_down_highconf, nustress_res$final_down),
  Database_Up = intersect(database_up_highconf, nustress_res$final_up),
  Database_Down = intersect(database_down_highconf, nustress_res$final_down)
)

nustress_binary <- make_upset_input(nustress_gene_lists_finalonly)

pdf(
  file.path(final_dir, "NuStress_UpSet_final_source.pdf"),
  width = 6, height = 5
)
upset(
  nustress_binary,
  nsets = 4,
  nintersects = 20,
  order.by = "freq",
  keep.order = TRUE,
  mainbar.y.label = "Gene count",
  sets.x.label = "Set size",
  main.bar.color = "steelblue",
  sets.bar.color = "darkred",
  matrix.color = "black"
)
dev.off()

# 输出来源明细表
nustress_source_table <- build_source_long_table(
  final_up = nustress_res$final_up,
  final_down = nustress_res$final_down,
  lit_up = literature_up_highconf,
  lit_down = literature_down_highconf,
  db_up = database_up_highconf,
  db_down = database_down_highconf
)

write.csv(
  nustress_source_table,
  file.path(final_dir, "NuStress_final_gene_source_table.csv"),
  row.names = FALSE
)

############################################################
## 补充.2 NuStressCore 核心基因集：仅基于最终保留基因集绘制来源 UpSet
############################################################
nustress_core_gene_lists_finalonly <- list(
  LiteratureCore_Up = intersect(literature_core_up_highconf, nustress_core_res$final_up),
  LiteratureCore_Down = intersect(literature_core_down_highconf, nustress_core_res$final_down),
  DatabaseCore_Up = intersect(database_core_up_genes, nustress_core_res$final_up),
  DatabaseCore_Down = intersect(database_core_down_genes, nustress_core_res$final_down)
)

nustress_core_binary <- make_upset_input(nustress_core_gene_lists_finalonly)

pdf(
  file.path(final_dir, "NuStressCore_UpSet_final_source.pdf"),
  width = 6, height = 5
)
upset(
  nustress_core_binary,
  nsets = 4,
  nintersects = 20,
  order.by = "freq",
  keep.order = TRUE,
  mainbar.y.label = "Gene count",
  sets.x.label = "Set size",
  main.bar.color = "steelblue",
  sets.bar.color = "darkred",
  matrix.color = "black"
)
dev.off()

# 输出来源明细表
nustress_core_source_table <- build_source_long_table(
  final_up = nustress_core_res$final_up,
  final_down = nustress_core_res$final_down,
  lit_up = literature_core_up_highconf,
  lit_down = literature_core_down_highconf,
  db_up = database_core_up_genes,
  db_down = database_core_down_genes
)

write.csv(
  nustress_core_source_table,
  file.path(final_dir, "NuStressCore_final_gene_source_table.csv"),
  row.names = FALSE
)


############################################################
## 17. 保存最终 gene sets
############################################################
geneSets_final <- list(
  NuStress_UP = nustress_res$final_up,
  NuStress_DOWN = nustress_res$final_down,
  NuStressCore_UP = nustress_core_res$final_up,
  NuStressCore_DOWN = nustress_core_res$final_down
)

save(
  geneSets_final,
  file = file.path(final_dir, "NuStress_geneSets_final.Rdata")
)

write.csv(
  data.frame(gene_symbol = geneSets_final$NuStress_UP, stringsAsFactors = FALSE),
  file.path(final_dir, "NuStress_UP.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = geneSets_final$NuStress_DOWN, stringsAsFactors = FALSE),
  file.path(final_dir, "NuStress_DOWN.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = geneSets_final$NuStressCore_UP, stringsAsFactors = FALSE),
  file.path(final_dir, "NuStressCore_UP.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(gene_symbol = geneSets_final$NuStressCore_DOWN, stringsAsFactors = FALSE),
  file.path(final_dir, "NuStressCore_DOWN.csv"),
  row.names = FALSE
)

############################################################
## 18. 汇总统计表
############################################################
final_summary <- data.frame(
  category = c(
    "database_highconf_up_freq_ge4",
    "database_highconf_down_freq_ge4",
    "database_core_up_freq_ge5_treat_ge2",
    "database_core_down_freq_ge5_treat_ge2",
    "literature_highconf_up",
    "literature_highconf_down",
    "literature_core_up",
    "literature_core_down",
    "NuStress_conflict_removed",
    "NuStress_up_final",
    "NuStress_down_final",
    "NuStressCore_conflict_removed",
    "NuStressCore_up_final",
    "NuStressCore_down_final"
  ),
  count = c(
    length(database_up_highconf),
    length(database_down_highconf),
    length(database_core_up_genes),
    length(database_core_down_genes),
    length(literature_up_highconf),
    length(literature_down_highconf),
    length(literature_core_up_highconf),
    length(literature_core_down_highconf),
    length(nustress_res$conflict_genes),
    length(nustress_res$final_up),
    length(nustress_res$final_down),
    length(nustress_core_res$conflict_genes),
    length(nustress_core_res$final_up),
    length(nustress_core_res$final_down)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  final_summary,
  file.path(final_dir, "final_gene_set_summary_statistics.csv"),
  row.names = FALSE
)

############################################################
## 19. 控制台输出
############################################################
cat("====================================================\n")
cat("summary_results saved to:\n", summary_dir, "\n\n")
cat("final_nucleolar_gene_sets saved to:\n", final_dir, "\n\n")

cat("[Database]\n")
cat("Raw all_long rows:", nrow(all_long), "\n")
cat("Consistent all_long rows:", nrow(all_long_consistent), "\n")
cat("Cross-dataset conflict genes:", length(conflict_genes), "\n")
cat("Database highconf UP:", length(database_up_highconf), "\n")
cat("Database highconf DOWN:", length(database_down_highconf), "\n")
cat("Database core UP:", length(database_core_up_genes), "\n")
cat("Database core DOWN:", length(database_core_down_genes), "\n\n")

cat("[Final gene sets]\n")
cat("NuStress UP:", length(nustress_res$final_up), "\n")
cat("NuStress DOWN:", length(nustress_res$final_down), "\n")
cat("NuStressCore UP:", length(nustress_core_res$final_up), "\n")
cat("NuStressCore DOWN:", length(nustress_core_res$final_down), "\n")
cat("====================================================\n")