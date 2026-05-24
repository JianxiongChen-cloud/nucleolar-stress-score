rm(list = ls())
options(stringsAsFactors = FALSE)

############################################################
## 0. Path settings
############################################################
out_dir <- "/home/xxm_xxm/CJX_workspace/Data validation/GSE6400"
data_dir <- "/home/xxm_xxm/CJX_workspace/Dataset/GSE6400"
geneset_dir <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
setwd(out_dir)

############################################################
## 1. Load R packages
############################################################
library(readr)
library(AnnotationDbi)
library(hgu133plus2.db)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(multcomp)

############################################################
## 2. Read GSE6400 expression matrix
############################################################
expr_probe <- read_table(
  file.path(data_dir, "GSE6400_series_matrix.txt.gz"),
  comment = "!",
  show_col_types = FALSE
)

expr_probe <- as.data.frame(expr_probe, check.names = FALSE)
rownames(expr_probe) <- expr_probe[, 1]
expr_probe <- expr_probe[, -1, drop = FALSE]

rownames(expr_probe) <- gsub('"', "", rownames(expr_probe))
colnames(expr_probe) <- gsub('"', "", colnames(expr_probe))

expr_probe <- as.matrix(expr_probe)
mode(expr_probe) <- "numeric"

cat("===== Raw matrix dimension =====\n")
print(dim(expr_probe))
cat("===== First few sample names =====\n")
print(colnames(expr_probe)[1:min(10, ncol(expr_probe))])

############################################################
## 3. Probe annotation
############################################################
ids <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr_probe),
  columns = c("SYMBOL"),
  keytype = "PROBEID"
)

colnames(ids) <- c("probe_id", "symbol")
ids <- ids[!is.na(ids$symbol) & ids$symbol != "", ]
ids <- ids[ids$probe_id %in% rownames(expr_probe), ]
ids <- unique(ids)

cat("===== Annotated probes =====\n")
print(nrow(ids))

############################################################
## 4. Keep matched probes
############################################################
expr_probe_annot <- expr_probe[ids$probe_id, , drop = FALSE]
ids <- ids[match(rownames(expr_probe_annot), ids$probe_id), ]

cat("===== Matched annotated matrix dimension =====\n")
print(dim(expr_probe_annot))

############################################################
## 5. Collapse duplicated probes
## Strategy: keep the probe with highest median expression
############################################################
probe_median <- apply(expr_probe_annot, 1, median, na.rm = TRUE)
ids$median_expr <- probe_median

ids <- ids[order(ids$symbol, -ids$median_expr), ]
ids <- ids[!duplicated(ids$symbol), ]

expr_gene <- expr_probe_annot[ids$probe_id, , drop = FALSE]
rownames(expr_gene) <- ids$symbol

cat("===== Gene-level matrix dimension =====\n")
print(dim(expr_gene))

############################################################
## 6. Select samples
############################################################
expr_gene <- expr_gene[, 1:12, drop = FALSE]

cat("===== Selected samples =====\n")
print(colnames(expr_gene))

############################################################
## 7. Group information
############################################################
group_info <- data.frame(
  Sample = colnames(expr_gene),
  Group = factor(
    c(
      rep("ActD", 3),
      rep("Control", 3),
      rep("1.25uM sapphyrin", 3),
      rep("2.5uM sapphyrin", 3)
    ),
    levels = c("Control", "ActD", "1.25uM sapphyrin", "2.5uM sapphyrin")
  ),
  stringsAsFactors = FALSE
)

print(group_info)

############################################################
## 8. Read gene sets
############################################################
read_gene_vector <- function(file) {
  x <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  x <- x[, 1, drop = TRUE]
  x <- unique(toupper(trimws(x)))
  x <- x[!is.na(x) & x != ""]
  return(x)
}

full_up   <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_up_genes.csv"))
full_down <- read_gene_vector(file.path(geneset_dir, "nucleolar_stress_down_genes.csv"))

core_up_file   <- file.path(geneset_dir, "nucleolar_stress_core_up_genes.csv")
core_down_file <- file.path(geneset_dir, "nucleolar_stress_core_down_genes.csv")

has_core <- file.exists(core_up_file) && file.exists(core_down_file)

if (has_core) {
  core_up   <- read_gene_vector(core_up_file)
  core_down <- read_gene_vector(core_down_file)
} else {
  warning("Core gene set files not found. Only full gene set will be analyzed.")
}

############################################################
## 9. Intersect gene sets with the expression matrix
############################################################
full_up_in   <- intersect(full_up, rownames(expr_gene))
full_down_in <- intersect(full_down, rownames(expr_gene))

cat("Full gene set overlap:\n")
cat("NuS_Up:", length(full_up_in), "\n")
cat("NuS_Down:", length(full_down_in), "\n")

if (length(full_up_in) < 5 || length(full_down_in) < 5) {
  stop("Too few overlapping genes for the full gene set.")
}

geneSets_full <- list(
  NuS_Up = full_up_in,
  NuS_Down = full_down_in
)

if (has_core) {
  core_up_in   <- intersect(core_up, rownames(expr_gene))
  core_down_in <- intersect(core_down, rownames(expr_gene))
  
  cat("Core gene set overlap:\n")
  cat("NuSCore_Up:", length(core_up_in), "\n")
  cat("NuSCore_Down:", length(core_down_in), "\n")
  
  if (length(core_up_in) >= 5 && length(core_down_in) >= 5) {
    geneSets_core <- list(
      NuSCore_Up = core_up_in,
      NuSCore_Down = core_down_in
    )
  } else {
    warning("Too few overlapping genes for the core gene set. Core analysis skipped.")
    has_core <- FALSE
  }
}

############################################################
## 10. NuS score calculation function
############################################################
calc_nus_score <- function(expr_gene, geneSets, method = c("gsva", "ssgsea"), set_name = "Full") {
  method <- match.arg(method)
  
  if (method == "gsva") {
    param <- gsvaParam(
      exprData = expr_gene,
      geneSets = geneSets,
      minSize = 5,
      maxSize = 500
    )
  } else {
    param <- ssgseaParam(
      exprData = expr_gene,
      geneSets = geneSets,
      minSize = 5,
      maxSize = 500,
      alpha = 0.25,
      normalize = TRUE
    )
  }
  
  res <- gsva(param)
  
  if (all(c("NuS_Up", "NuS_Down") %in% rownames(res))) {
    score <- res["NuS_Up", ] - res["NuS_Down", ]
  } else if (all(c("NuSCore_Up", "NuSCore_Down") %in% rownames(res))) {
    score <- res["NuSCore_Up", ] - res["NuSCore_Down", ]
  } else {
    stop("Up/Down score rows not found in GSVA output.")
  }
  
  out <- data.frame(
    Sample = colnames(expr_gene),
    Score = as.numeric(score),
    Method = method,
    GeneSet = set_name,
    stringsAsFactors = FALSE
  )
  return(out)
}

############################################################
## 11. Calculate GSVA and ssGSEA scores
############################################################
score_list <- list()

score_list[["full_gsva"]]   <- calc_nus_score(expr_gene, geneSets_full, "gsva", "Full")
score_list[["full_ssgsea"]] <- calc_nus_score(expr_gene, geneSets_full, "ssgsea", "Full")

if (has_core) {
  score_list[["core_gsva"]]   <- calc_nus_score(expr_gene, geneSets_core, "gsva", "Core")
  score_list[["core_ssgsea"]] <- calc_nus_score(expr_gene, geneSets_core, "ssgsea", "Core")
}

scores_all <- bind_rows(score_list)
scores_all <- left_join(scores_all, group_info, by = "Sample")
scores_all$Group <- factor(
  scores_all$Group,
  levels = c("Control", "ActD", "1.25uM sapphyrin", "2.5uM sapphyrin")
)
scores_all$Method <- factor(scores_all$Method, levels = c("gsva", "ssgsea"))
scores_all$GeneSet <- factor(scores_all$GeneSet, levels = c("Full", "Core"))

write.csv(scores_all, "GSE6400_NuS_scores_GSVA_ssGSEA.csv", row.names = FALSE)

############################################################
## 12. Statistical testing: one-way ANOVA + Dunnett's test
############################################################
test_res <- scores_all %>%
  group_by(Method, GeneSet) %>%
  group_modify(~{
    dat <- .x
    dat$Group <- droplevels(dat$Group)
    
    if (!("Control" %in% levels(dat$Group))) {
      stop("Control group not found in Group levels.")
    }
    
    dat$Group <- relevel(dat$Group, ref = "Control")
    
    ## one-way ANOVA
    fit <- aov(Score ~ Group, data = dat)
    anova_p <- tryCatch(
      summary(fit)[[1]][["Pr(>F)"]][1],
      error = function(e) NA_real_
    )
    
    ## group summary
    mean_sd_df <- dat %>%
      group_by(Group) %>%
      summarise(
        Mean = mean(Score, na.rm = TRUE),
        SD = sd(Score, na.rm = TRUE),
        N = sum(!is.na(Score)),
        .groups = "drop"
      )
    
    control_row <- mean_sd_df %>% filter(Group == "Control")
    if (nrow(control_row) != 1) {
      stop("Control group summary is not unique.")
    }
    
    ## Dunnett
    dunnett_res <- tryCatch(
      summary(glht(fit, linfct = mcp(Group = "Dunnett"))),
      error = function(e) NULL
    )
    
    if (is.null(dunnett_res)) {
      out <- data.frame(
        comparison = NA_character_,
        treatment = NA_character_,
        p_anova = anova_p,
        p_dunnett = NA_real_,
        mean_control = control_row$Mean,
        sd_control = control_row$SD,
        n_control = control_row$N,
        mean_treat = NA_real_,
        sd_treat = NA_real_,
        n_treat = NA_real_,
        diff_mean = NA_real_,
        stringsAsFactors = FALSE
      )
      return(out)
    }
    
    comp_names <- names(dunnett_res$test$coefficients)
    comp_diff  <- as.numeric(dunnett_res$test$coefficients)
    comp_p     <- as.numeric(dunnett_res$test$pvalues)
    
    treat_names <- sub(" - Control$", "", comp_names)
    
    treat_df <- mean_sd_df %>%
      filter(Group != "Control") %>%
      as.data.frame()
    
    match_idx <- match(treat_names, as.character(treat_df$Group))
    
    out <- data.frame(
      comparison = comp_names,
      treatment = treat_names,
      p_anova = rep(anova_p, length(comp_names)),
      p_dunnett = comp_p,
      mean_control = rep(control_row$Mean, length(comp_names)),
      sd_control = rep(control_row$SD, length(comp_names)),
      n_control = rep(control_row$N, length(comp_names)),
      mean_treat = treat_df$Mean[match_idx],
      sd_treat = treat_df$SD[match_idx],
      n_treat = treat_df$N[match_idx],
      diff_mean = comp_diff,
      stringsAsFactors = FALSE
    )
    
    out
  }) %>%
  ungroup()

write.csv(test_res, "GSE6400_NuS_score_statistics_ANOVA_Dunnett.csv", row.names = FALSE)
print(test_res)

############################################################
## 13. Plot violin plots only
############################################################
fill_colors <- c(
  "Control" = "#C4E7C1",
  "ActD" = "#93D4BC",
  "1.25uM sapphyrin" = "#51B3D1",
  "2.5uM sapphyrin" = "#08589E"
)

p <- ggplot(scores_all, aes(x = Group, y = Score, fill = Group)) +
  geom_violin(trim = FALSE, color = "black", linewidth = 0.8, alpha = 0.85) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.7
  ) +
  facet_grid(GeneSet ~ Method, scales = "free_y") +
  scale_fill_manual(values = fill_colors) +
  labs(
    title = "NuS scores in GSE6400",
    x = NULL,
    y = "NuS score"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("GSE6400_NuS_violinplot.pdf", p, width = 8, height = 5.5, bg = "white")
ggsave("GSE6400_NuS_violinplot.png", p, width = 8, height = 5.5, dpi = 300, bg = "white")

############################################################
## 14. Save workspace
############################################################
save(
  expr_gene,
  group_info,
  scores_all,
  file = "GSE6400_NuS_validation_workspace.Rdata"
)