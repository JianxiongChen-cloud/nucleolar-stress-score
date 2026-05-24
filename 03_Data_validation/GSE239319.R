############################################################
## GSE239319
## Fully optimized pipeline:
## mouse-to-human gene conversion +
## expression matrix construction +
## ssGSEA + statistical analysis + visualization
## Full gene set only
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

setwd("/home/xxm_xxm/CJX_workspace/Data validation/GSE239319/")

library(readxl)
library(dplyr)
library(babelgene)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)

geo_id <- "GSE239319"
input_file <- "GSE239319_3RNAseq_IAO_results.xlsx"
geneset_file <- "/home/xxm_xxm/CJX_workspace/geneset/genesets.Rdata"

############################################################
## Step 1. Read raw input table
############################################################
dat_raw <- read_excel(input_file)
dat_raw <- as.data.frame(dat_raw)

cat("===== Raw table dimension =====\n")
print(dim(dat_raw))
cat("===== Raw column names =====\n")
print(colnames(dat_raw))

############################################################
## Step 2. Extract expression matrix
## Columns 10:13 were used in the original pipeline
############################################################
expr_raw <- dat_raw[, 10:13, drop = FALSE]
colnames(expr_raw) <- as.character(unlist(expr_raw[2, ]))
expr_raw <- expr_raw[-c(1, 2), , drop = FALSE]

gene_mouse <- dat_raw$...2
gene_mouse <- gene_mouse[-c(1, 2)]

expr_raw <- as.matrix(expr_raw)
mode(expr_raw) <- "numeric"
rownames(expr_raw) <- gene_mouse

cat("===== Mouse expression matrix dimension =====\n")
print(dim(expr_raw))
cat("===== Mouse expression range =====\n")
print(range(expr_raw, na.rm = TRUE))
cat("===== Sample names =====\n")
print(colnames(expr_raw))

############################################################
## Step 3. Convert mouse genes to human orthologs
############################################################
converted_genes <- orthologs(
  genes = rownames(expr_raw),
  species = "mouse",
  human = FALSE
)

cat("===== Ortholog mapping rows =====\n")
print(nrow(converted_genes))

## Retain genes successfully matched to the expression matrix
match_idx <- match(converted_genes$symbol, rownames(expr_raw))
keep <- !is.na(match_idx)

converted_genes <- converted_genes[keep, , drop = FALSE]
expr_human <- expr_raw[match_idx[keep], , drop = FALSE]
rownames(expr_human) <- converted_genes$human_symbol

cat("===== Human-mapped expression dimension =====\n")
print(dim(expr_human))

############################################################
## Step 4. Remove empty row names / NA / duplicated genes
############################################################
valid_symbol <- !is.na(rownames(expr_human)) & rownames(expr_human) != ""
expr_human <- expr_human[valid_symbol, , drop = FALSE]

expr_human_df <- as.data.frame(expr_human, check.names = FALSE)
expr_human_df$median_expr <- apply(expr_human, 1, median, na.rm = TRUE)
expr_human_df$symbol <- rownames(expr_human)

expr_human_df <- expr_human_df[order(expr_human_df$symbol, -expr_human_df$median_expr), ]
expr_human_df <- expr_human_df[!duplicated(expr_human_df$symbol), ]

expr_gene <- as.matrix(expr_human_df[, colnames(expr_human), drop = FALSE])
mode(expr_gene) <- "numeric"
rownames(expr_gene) <- expr_human_df$symbol

cat("===== Final gene-level matrix dimension =====\n")
print(dim(expr_gene))
cat("===== Final gene-level range =====\n")
print(range(expr_gene, na.rm = TRUE))

save(expr_gene, file = paste0(geo_id, "_step1_workspace.Rdata"))
write.csv(expr_gene, file = "expressionmetrix_GSE.csv")

############################################################
## Step 5. Read expression matrix and perform log2 transformation
############################################################
expr <- expr_gene
mode(expr) <- "numeric"

expr <- log2(expr + 1)
expr <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]
rownames(expr) <- toupper(trimws(rownames(expr)))
expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

cat("===== Expression matrix after log2(x+1) =====\n")
print(dim(expr))
print(range(expr, na.rm = TRUE))

############################################################
## Step 6. Group information
############################################################
group_info <- data.frame(
  Sample = colnames(expr),
  Group = rep(c("inducible", "non-inducible"), each = 2),
  stringsAsFactors = FALSE
)

control_group <- "non-inducible"

group_info$Group <- factor(
  group_info$Group,
  levels = c("inducible", "non-inducible")
)

cat("===== Group info =====\n")
print(group_info)

############################################################
## Step 7. Read the full gene set
############################################################
load(geneset_file)

geneSets_all <- list(
  NuS_Up = unique(toupper(NuS_up_all_unique)),
  NuS_Down = unique(toupper(NuS_down_all_unique))
)

geneSets_all$NuS_Up <- intersect(geneSets_all$NuS_Up, rownames(expr))
geneSets_all$NuS_Down <- intersect(geneSets_all$NuS_Down, rownames(expr))

cat("===== Gene set overlap =====\n")
cat("NuS_Up:", length(geneSets_all$NuS_Up), "\n")
cat("NuS_Down:", length(geneSets_all$NuS_Down), "\n")

if (length(geneSets_all$NuS_Up) < 5 || length(geneSets_all$NuS_Down) < 5) {
  stop("Too few overlapping genes for ssGSEA.")
}

############################################################
## Step 8. ssGSEA
############################################################
ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr,
  geneSets = geneSets_all,
  minSize = 5,
  maxSize = 500,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_res <- gsva(ssgsea_param)

if (!all(c("NuS_Up", "NuS_Down") %in% rownames(ssgsea_res))) {
  stop("ssGSEA result does not contain NuS_Up or NuS_Down.")
}

nus_score <- ssgsea_res["NuS_Up", ] - ssgsea_res["NuS_Down", ]

scores <- data.frame(
  Sample = colnames(expr),
  NuS_Up = as.numeric(ssgsea_res["NuS_Up", ]),
  NuS_Down = as.numeric(ssgsea_res["NuS_Down", ]),
  NuS_Score = as.numeric(nus_score),
  stringsAsFactors = FALSE
)

scores <- merge(scores, group_info, by = "Sample")
scores$Group <- factor(scores$Group, levels = levels(group_info$Group))

plot_data <- data.frame(
  Score = scores$NuS_Score,
  Group = scores$Group
)

cat("===== NuS scores =====\n")
print(scores)

write.csv(scores, paste0(geo_id, "_NuS_all_ssGSEA_scores.csv"), row.names = FALSE)
write.csv(as.data.frame(ssgsea_res), paste0(geo_id, "_NuS_all_ssGSEA_matrix.csv"), row.names = TRUE)

############################################################
## Step 9. Statistical analysis
############################################################
n_groups <- length(unique(plot_data$Group))
cat("===== Number of groups =====\n")
print(n_groups)

if (n_groups == 2) {
  overall_test <- t.test(Score ~ Group, data = plot_data)
  overall_p <- overall_test$p.value
  
  stat_summary <- data.frame(
    method = "t.test",
    p.value = overall_p,
    stringsAsFactors = FALSE
  )
} else {
  overall_test <- aov(Score ~ Group, data = plot_data)
  overall_p <- summary(overall_test)[[1]][["Pr(>F)"]][1]
  
  tukey_res <- TukeyHSD(overall_test)
  tukey_df <- as.data.frame(tukey_res$Group)
  tukey_df$comparison <- rownames(tukey_df)
  rownames(tukey_df) <- NULL
  write.csv(tukey_df, paste0(geo_id, "_NuS_all_TukeyHSD.csv"), row.names = FALSE)
  
  stat_summary <- data.frame(
    method = "ANOVA",
    p.value = overall_p,
    stringsAsFactors = FALSE
  )
}

write.csv(stat_summary, paste0(geo_id, "_NuS_all_overall_stats.csv"), row.names = FALSE)

## vs control
control_data <- plot_data$Score[plot_data$Group == control_group]
vs_control_res <- data.frame()

for (grp in unique(plot_data$Group)) {
  if (grp != control_group) {
    grp_data <- plot_data$Score[plot_data$Group == grp]
    pval <- t.test(grp_data, control_data)$p.value
    
    vs_control_res <- rbind(
      vs_control_res,
      data.frame(
        comparison = paste(grp, "vs", control_group),
        p.value = pval,
        stringsAsFactors = FALSE
      )
    )
  }
}

write.csv(vs_control_res, paste0(geo_id, "_NuS_all_vsControl_stats.csv"), row.names = FALSE)

############################################################
## Step 10. Violin plot
## Use the predefined fixed visualization style
############################################################
fill_colors <- c(
  "inducible" = "#93D4BC",
  "non-inducible" = "#08589E"
)

p1 <- ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(
    alpha = 0.85,
    trim = TRUE,
    scale = "width",
    width = 0.7,
    color = "black",
    linewidth = 1.0,
    adjust = 1.2
  ) +
  geom_boxplot(
    width = 0.12,
    alpha = 0.9,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = paste0("NuS Score (ssGSEA) - ", geo_id, " - all"),
    x = NULL,
    y = "NuS Score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(
      hjust = 0.5, size = 16, face = "bold", color = "#13646B",
      margin = margin(b = 15)
    ),
    panel.background = element_rect(fill = "white", colour = "black", linewidth = 1),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    label.x.npc = "center",
    size = 5
  )

print(p1)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.pdf"),
  plot = p1,
  width = 2.5,
  height = 4,
  bg = "white"
)

ggsave(
  filename = paste0(geo_id, "_NuS_all_groups_violin.png"),
  plot = p1,
  width = 2.5,
  height = 4,
  dpi = 600,
  bg = "white"
)

############################################################
## Step 11. Save workspace
############################################################
save(
  expr,
  expr_gene,
  group_info,
  geneSets_all,
  ssgsea_res,
  scores,
  plot_data,
  stat_summary,
  vs_control_res,
  file = paste0(geo_id, "_NuS_all_workspace.Rdata")
)

cat("Analysis completed.\n")
cat("Output files saved in:\n")
cat(getwd(), "\n")