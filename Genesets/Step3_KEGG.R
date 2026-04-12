rm(list = ls())
options(stringsAsFactors = FALSE)

#===============================
# 0. 初始化
#===============================
setwd("/home/xxm_xxm/CJX_workspace/geneset/")

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)

# step2 输出目录（按你的 step2 文件命名）
input_dir <- file.path(getwd(), "final_nucleolar_gene_sets")
out_dir <- file.path(getwd(), "step3_KEGG_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#===============================
# 1. 读取基因集
#===============================
read_gene_symbols <- function(file_path) {
  if (!file.exists(file_path)) stop(paste("File not found:", file_path))
  
  x <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("gene_symbol" %in% colnames(x))) {
    stop(paste("Missing 'gene_symbol' column in:", file_path))
  }
  
  g <- unique(trimws(toupper(x$gene_symbol)))
  g <- g[!is.na(g) & g != ""]
  return(sort(g))
}

NuStress_UP <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_up_genes.csv"))
NuStress_DOWN <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_down_genes.csv"))
NuStressCore_UP <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_core_up_genes.csv"))
NuStressCore_DOWN <- read_gene_symbols(file.path(input_dir, "nucleolar_stress_core_down_genes.csv"))

NuStress_ALL <- sort(unique(c(NuStress_UP, NuStress_DOWN)))
NuStressCore_ALL <- sort(unique(c(NuStressCore_UP, NuStressCore_DOWN)))

cat("NuStress UP:", length(NuStress_UP), "\n")
cat("NuStress DOWN:", length(NuStress_DOWN), "\n")
cat("NuStress ALL:", length(NuStress_ALL), "\n\n")

cat("NuStressCore UP:", length(NuStressCore_UP), "\n")
cat("NuStressCore DOWN:", length(NuStressCore_DOWN), "\n")
cat("NuStressCore ALL:", length(NuStressCore_ALL), "\n\n")

#===============================
# 2. SYMBOL 转 ENTREZID
#===============================
convert_to_entrez_df <- function(gene_symbols) {
  gene_symbols <- unique(trimws(toupper(gene_symbols)))
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  
  mapping <- mapping[!is.na(mapping$ENTREZID) & mapping$ENTREZID != "", ]
  mapping <- mapping[!duplicated(mapping$SYMBOL), ]
  rownames(mapping) <- NULL
  
  return(mapping)
}

#===============================
# 3. 统一绘图函数
#    配色严格沿用你的参考代码
#===============================
green_blue_gradient <- colorRampPalette(
  c("#C4E7C1", "#93D4BC", "#51B3D1", "#08589E")
)(100)

plot_kegg_bar <- function(result_df, set_name, out_subdir, top_n = 10) {
  plot_data <- result_df %>%
    head(top_n) %>%
    mutate(
      log10_padj = -log10(p.adjust),
      Description = factor(Description, levels = rev(Description))
    )
  
  p1 <- ggplot(plot_data, aes(x = log10_padj, y = Description, fill = log10_padj)) +
    geom_bar(stat = "identity", width = 0.8, alpha = 0.9) +
    scale_fill_gradientn(
      colors = green_blue_gradient,
      name = "-log10(Adj.p)"
    ) +
    labs(
      title = "KEGG Pathway Enrichment Analysis",
      subtitle = paste0("Top ", min(top_n, nrow(result_df)), " Significant Pathways (", set_name, ")"),
      x = "-log10(Adjusted p-value)",
      y = "Pathway"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", margin = margin(b = 15)),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.text.x = element_text(margin = margin(t = 5)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)))
  
  ggsave(
    file.path(out_subdir, paste0("KEGG_", set_name, "_barplot.png")),
    p1, width = 8, height = 10, dpi = 600
  )
  ggsave(
    file.path(out_subdir, paste0("KEGG_", set_name, "_barplot.pdf")),
    p1, width = 8, height = 10
  )
  
  print(p1)
}

plot_kegg_dot <- function(result_df, set_name, out_subdir, top_n = 10) {
  plot_data <- result_df %>%
    head(top_n) %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        xx <- strsplit(x, "/")[[1]]
        as.numeric(xx[1]) / as.numeric(xx[2])
      }),
      Description = factor(Description, levels = rev(Description)),
      log10_padj = -log10(p.adjust)
    )
  
  p2 <- ggplot(plot_data, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, fill = log10_padj), shape = 21, colour = "black", alpha = 0.9) +
    scale_fill_gradientn(
      colors = green_blue_gradient,
      name = "-log10(Adj.p)"
    ) +
    labs(
      title = "KEGG Pathway Enrichment Analysis",
      subtitle = paste0("Top ", min(top_n, nrow(result_df)), " Significant Pathways (", set_name, ")"),
      x = "GeneRatio",
      y = "Pathway",
      size = "Count"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40", margin = margin(b = 15)),
      panel.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.text.x = element_text(margin = margin(t = 5)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(
    file.path(out_subdir, paste0("KEGG_", set_name, "_dotplot.png")),
    p2, width = 8, height = 10, dpi = 600
  )
  ggsave(
    file.path(out_subdir, paste0("KEGG_", set_name, "_dotplot.pdf")),
    p2, width = 8, height = 10
  )
  
  print(p2)
}

#===============================
# 4. 单个基因集 KEGG 分析函数
#===============================
run_kegg <- function(gene_symbols, set_name, out_dir) {
  cat("====================================\n")
  cat("Running KEGG for:", set_name, "\n")
  cat("Input genes:", length(gene_symbols), "\n")
  
  out_subdir <- file.path(out_dir, set_name)
  if (!dir.exists(out_subdir)) dir.create(out_subdir, recursive = TRUE)
  
  write.csv(
    data.frame(gene_symbol = gene_symbols, stringsAsFactors = FALSE),
    file.path(out_subdir, paste0(set_name, "_merged_genes.csv")),
    row.names = FALSE
  )
  
  mapping_df <- convert_to_entrez_df(gene_symbols)
  entrez_ids <- unique(mapping_df$ENTREZID)
  
  write.csv(
    mapping_df,
    file.path(out_subdir, paste0(set_name, "_SYMBOL_to_ENTREZID.csv")),
    row.names = FALSE
  )
  
  cat("成功转换的 ENTREZID 数量:", length(entrez_ids), "/", length(gene_symbols), "\n")
  
  if (length(entrez_ids) < 5) {
    cat("可用于富集分析的 ENTREZID 数量太少，跳过:", set_name, "\n\n")
    return(NULL)
  }
  
  cat("正在进行 KEGG 富集分析...\n")
  kegg_result <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
    cat("没有找到显著富集的 KEGG 通路:", set_name, "\n\n")
    return(NULL)
  }
  
  result_df <- as.data.frame(kegg_result) %>%
    arrange(p.adjust)
  
  cat("找到", nrow(result_df), "个显著富集的 KEGG 通路\n")
  
  write.csv(
    result_df,
    file.path(out_subdir, paste0("KEGG_Enrichment_", set_name, ".csv")),
    row.names = FALSE
  )
  
  plot_kegg_bar(result_df, set_name, out_subdir, top_n = 10)
  plot_kegg_dot(result_df, set_name, out_subdir, top_n = 10)
  
  cat("\n=== 最显著的前10个 KEGG 通路 ===\n")
  print(result_df[1:min(10, nrow(result_df)), c("Description", "p.adjust", "Count", "GeneRatio")])
  cat("\n")
  
  return(result_df)
}

#===============================
# 5. 分别分析 NuStress 和 NuStressCore
#===============================
res_nustress <- run_kegg(
  gene_symbols = NuStress_ALL,
  set_name = "NuStress_ALL",
  out_dir = out_dir
)

res_nustress_core <- run_kegg(
  gene_symbols = NuStressCore_ALL,
  set_name = "NuStressCore_ALL",
  out_dir = out_dir
)

#===============================
# 6. 汇总统计
#===============================
summary_df <- data.frame(
  gene_set = c("NuStress_ALL", "NuStressCore_ALL"),
  input_gene_n = c(length(NuStress_ALL), length(NuStressCore_ALL)),
  sig_pathway_n = c(
    ifelse(is.null(res_nustress), 0, nrow(res_nustress)),
    ifelse(is.null(res_nustress_core), 0, nrow(res_nustress_core))
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  file.path(out_dir, "KEGG_summary_statistics.csv"),
  row.names = FALSE
)

cat("Step 3 finished.\n")
print(summary_df)