# Step4_liver 手动注释 + RiboSis/NuS 可视化与导出（最终版）
# ======================================================================
# 1. 创建环境
# ======================================================================
getwd()
setwd("/home/xxm_xxm/CJX_workspace/SC-ST/")
rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 1e6)

# 加载必要包
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)
library(dplyr)

# ======================================================================
# 2. 参数设置
# ======================================================================
# 基因集文件绝对路径
nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

# 细胞类型颜色
cell_type_colors <- c(
  "Fibroblast" = "#1f77b4",
  "Tumor" = "#d62728",
  "Hepatocytes" = "#00A087FF",
  "Endothelial" = "#F39B7FFF"
)

# 各样本的人工注释映射（来自原 Step4）
cluster_mapping_list <- list(
  "ST-liver1" = c(
    "Hepatocytes",
    "Tumor",
    "Tumor",
    "Hepatocytes",
    "Endothelial",
    "Fibroblast",
    "Tumor",
    "Endothelial",
    "Hepatocytes",
    "Hepatocytes"
  ),
  
  "ST-liver2" = c(
    "Hepatocytes",
    "Tumor",
    "Hepatocytes",
    "Fibroblast",
    "Tumor",
    "Hepatocytes",
    "Tumor",
    "Tumor",
    "Tumor",
    "Hepatocytes",
    "Endothelial",
    "Fibroblast",
    "Tumor",
    "Endothelial",
    "Fibroblast"
  ),
  
  "ST-liver3" = c(
    "Hepatocytes",
    "Hepatocytes",
    "Hepatocytes",
    "Hepatocytes",
    "Hepatocytes",
    "Tumor",
    "Hepatocytes",
    "Fibroblast",
    "Hepatocytes",
    "Hepatocytes"
  ),
  
  "ST-liver4" = c(
    "Hepatocytes",
    "Tumor",
    "Hepatocytes",
    "Fibroblast",
    "Fibroblast",
    "Tumor",
    "Hepatocytes",
    "Tumor",
    "Hepatocytes",
    "Endothelial"
  )
)

# 要批处理的样本
sample_names <- c("ST-liver1", "ST-liver2", "ST-liver3", "ST-liver4")

# ======================================================================
# 3. 公共函数
# ======================================================================

# ------------------------------
# 3.1 空间聚类绘图函数
# ------------------------------
MySpatialClusterPlot <- function(obj, group_by, point_size = 2, alpha = 0.8,
                                 custom_colors = NULL, show_legend = TRUE) {
  
  if (length(obj@images) < 1) stop("对象中没有 images 信息，请检查空间对象。")
  if (!group_by %in% colnames(obj@meta.data)) {
    stop("group_by 列不存在于 meta.data 中。可用列: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  }
  
  img <- obj@images[[1]]@image
  coords <- GetTissueCoordinates(obj)
  
  df <- data.frame(
    imagecol = coords$imagecol,
    imagerow = coords$imagerow,
    cluster = obj@meta.data[[group_by]]
  )
  
  cell_types <- unique(df$cluster)
  
  if (is.null(custom_colors)) {
    n_types <- length(cell_types)
    if (n_types <= 9) {
      custom_colors <- RColorBrewer::brewer.pal(n_types, "Set1")
    } else {
      custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_types)
    }
    names(custom_colors) <- cell_types
  } else {
    missing_types <- setdiff(cell_types, names(custom_colors))
    if (length(missing_types) > 0) {
      warning("以下细胞类型没有定义颜色，将使用默认颜色: ",
              paste(missing_types, collapse = ", "))
      n_missing <- length(missing_types)
      default_colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_missing)
      names(default_colors) <- missing_types
      custom_colors <- c(custom_colors, default_colors)
    }
  }
  
  img_height <- nrow(img)
  img_width <- ncol(img)
  
  p <- ggplot() +
    annotation_raster(
      img,
      xmin = 0, xmax = img_width,
      ymin = 0, ymax = img_height
    ) +
    geom_point(
      data = df,
      aes(x = imagecol, y = imagerow, color = cluster),
      size = point_size, alpha = alpha
    ) +
    scale_color_manual(values = custom_colors, name = group_by) +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none"
    ) +
    ggtitle(paste("Spatial distribution of", group_by))
  
  return(p)
}
# ------------------------------
# 3.2 连续变量空间绘图函数
# ------------------------------
MySpatialPlot2 <- function(obj, feature, point_size = 2, alpha = 0.8,
                           discrete = FALSE, color_option = "plasma",
                           limits = NULL, na.value = "grey80") {
  
  if (length(obj@images) < 1) stop("对象中没有 images 信息。")
  if (!feature %in% colnames(obj@meta.data)) {
    stop("特征列不存在于 meta.data 中。可用列: ",
         paste(colnames(obj@meta.data), collapse = ", "))
  }
  
  img <- obj@images[[1]]@image
  coords <- GetTissueCoordinates(obj)
  
  df <- data.frame(
    imagecol = coords$imagecol,
    imagerow = coords$imagerow,
    feature_value = obj@meta.data[[feature]]
  )
  
  p <- ggplot() +
    annotation_raster(
      img,
      xmin = 0, xmax = ncol(img),
      ymin = 0, ymax = nrow(img)
    ) +
    geom_point(
      data = df,
      aes(x = imagecol, y = imagerow, color = feature_value),
      size = point_size,
      alpha = alpha
    ) +
    scale_y_reverse() +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(feature)
  
  if (discrete) {
    p <- p + scale_color_brewer(
      palette = "Set2",
      name = feature,
      na.value = na.value
    )
  } else {
    if (is.null(limits)) {
      p <- p + scale_color_viridis_c(
        option = color_option,
        name = feature,
        na.value = na.value
      )
    } else {
      p <- p + scale_color_viridis_c(
        option = color_option,
        name = feature,
        limits = limits,
        oob = scales::squish,
        na.value = na.value
      )
    }
  }
  
  return(p)
}
# ------------------------------
# 3.3 NuS 评分函数
# ------------------------------
CalculateNuSScore <- function(obj, up_genes, down_genes, score_name = "NuS_Score") {
  
  if (length(up_genes) == 0) stop("up_genes 不能为空")
  if (length(down_genes) == 0) stop("down_genes 不能为空")
  
  available_genes <- rownames(obj)
  up_genes <- up_genes[up_genes %in% available_genes]
  down_genes <- down_genes[down_genes %in% available_genes]
  
  cat("UP基因:", length(up_genes), "个可用\n")
  cat("DOWN基因:", length(down_genes), "个可用\n")
  
  if (length(up_genes) == 0 || length(down_genes) == 0) {
    stop("没有足够的基因存在于数据中")
  }
  
  obj <- AddModuleScore(obj, features = list(up_genes), name = "UP_tmp")
  up_score <- obj@meta.data[["UP_tmp1"]]
  
  obj <- AddModuleScore(obj, features = list(down_genes), name = "DOWN_tmp")
  down_score <- obj@meta.data[["DOWN_tmp1"]]
  
  obj@meta.data[[score_name]] <- up_score - down_score
  
  obj@meta.data[["UP_tmp1"]] <- NULL
  obj@meta.data[["DOWN_tmp1"]] <- NULL
  
  cat("评分计算完成，添加到列:", score_name, "\n")
  return(obj)
}

# ======================================================================
# 4. 主循环：逐个样本处理
# ======================================================================
for (name in sample_names) {
  
  cat("\n=====================================================\n")
  cat("开始处理样本:", name, "\n")
  cat("=====================================================\n")
  
  rda_file <- paste0(name, ".rda")
  if (!file.exists(rda_file)) {
    warning("找不到数据文件: ", rda_file, "，跳过。")
    next
  }
  
  load(rda_file)
  
  if (!exists("obj_liver")) {
    warning("load 后未发现 obj_liver，对样本 ", name, " 跳过。")
    next
  }
  
  cat("对象类型:", class(obj_liver), "\n")
  cat("默认 Assay:", DefaultAssay(obj_liver), "\n")
  cat("图像数量:", length(obj_liver@images), "\n")
  cat("元数据列名:\n")
  print(colnames(obj_liver@meta.data))
  
  # --------------------------------------------------
  # 4.1 保存 HE 图
  # --------------------------------------------------
  p_he <- SpatialFeaturePlot(
    obj_liver,
    features = "nCount_Spatial",
    slot = "data",
    alpha = 0
  ) +
    theme(legend.position = "right") +
    NoLegend()
  
  ggsave(
    filename = paste0(name, "_HE.pdf"),
    plot = p_he,
    width = 4,
    height = 4
  )
  
  # --------------------------------------------------
  # 4.2 可视化 seurat_clusters
  # --------------------------------------------------
  cat("seurat_clusters 分布:\n")
  print(table(obj_liver$seurat_clusters))
  
  p_cluster <- MySpatialClusterPlot(
    obj_liver,
    "seurat_clusters",
    point_size = 2,
    alpha = 0.8
  )
  
  ggsave(
    filename = paste0(name, "_seurat_clusters.pdf"),
    plot = p_cluster,
    width = 8,
    height = 4
  )
  
  # --------------------------------------------------
  # 4.3 手动注释
  # --------------------------------------------------
  Idents(obj_liver) <- "seurat_clusters"
  
  if (!name %in% names(cluster_mapping_list)) {
    warning("未找到样本 ", name, " 的手动注释映射，跳过 manual.annot。")
  } else {
    new.cluster.ids <- cluster_mapping_list[[name]]
    
    current_levels <- levels(obj_liver)
    cat("当前 levels:\n")
    print(current_levels)
    
    if (length(new.cluster.ids) != length(current_levels)) {
      warning("样本 ", name, " 的 new.cluster.ids 长度与 levels(obj_liver) 不一致，跳过 RenameIdents。")
    } else {
      names(new.cluster.ids) <- current_levels
      obj_liver <- RenameIdents(obj_liver, new.cluster.ids)
      obj_liver$manual.annot <- Idents(obj_liver)
      
      p_manual <- MySpatialClusterPlot(
        obj_liver,
        "manual.annot",
        point_size = 2,
        alpha = 0.9,
        custom_colors = cell_type_colors
      )
      
      ggsave(
        filename = paste0(name, "_celltype.pdf"),
        plot = p_manual,
        width = 8,
        height = 4
      )
    }
  }
  
  # --------------------------------------------------
  # 4.4 计算 RiboSis
  # --------------------------------------------------
  cat("\n====================\n")
  cat("开始计算 RiboSis score\n")
  cat("====================\n")
  
  ribo_limits <- c(-0.15, 0.68)
  
  if (file.exists(ribosis_rdata)) {
    env_ribo <- new.env()
    load(ribosis_rdata, envir = env_ribo)
    
    cat("RiboSis 文件中的对象:\n")
    print(ls(env_ribo))
    
    if (exists("ribosis", envir = env_ribo)) {
      ribosis_raw <- get("ribosis", envir = env_ribo)
    } else {
      obj_names <- ls(env_ribo)
      ribo_candidate <- obj_names[grepl("ribo", obj_names, ignore.case = TRUE)][1]
      
      if (is.na(ribo_candidate)) {
        stop("RiboSis 文件中未找到名为 ribosis 的对象，也未自动匹配到相关对象。")
      } else {
        ribosis_raw <- get(ribo_candidate, envir = env_ribo)
        cat("自动使用对象:", ribo_candidate, "\n")
      }
    }
    
    cat("ribosis_raw class:\n")
    print(class(ribosis_raw))
    cat("ribosis_raw length:\n")
    print(length(ribosis_raw))
    
    if (is.data.frame(ribosis_raw) || is.matrix(ribosis_raw)) {
      ribosis_genes <- unique(as.character(unlist(ribosis_raw)))
    } else if (is.list(ribosis_raw)) {
      ribosis_genes <- unique(as.character(unlist(ribosis_raw)))
    } else {
      ribosis_genes <- unique(as.character(ribosis_raw))
    }
    
    ribosis_genes <- ribosis_genes[!is.na(ribosis_genes)]
    ribosis_genes <- ribosis_genes[ribosis_genes != ""]
    
    cat("整理后 ribosis_genes 数量:", length(ribosis_genes), "\n")
    cat("前10个基因:\n")
    print(head(ribosis_genes, 10))
    
    obj_genes <- rownames(obj_liver)
    ribosis_genes_present <- ribosis_genes[ribosis_genes %in% obj_genes]
    
    if (length(ribosis_genes_present) < 10) {
      cat("直接匹配到的基因太少，尝试大小写不敏感匹配...\n")
      obj_genes_upper <- toupper(obj_genes)
      ribosis_genes_upper <- toupper(ribosis_genes)
      match_idx <- match(ribosis_genes_upper, obj_genes_upper)
      ribosis_genes_present <- unique(obj_genes[na.omit(match_idx)])
    }
    
    cat("obj_liver 中存在的 RiboSis 基因数:", length(ribosis_genes_present), "\n")
    
    if (length(ribosis_genes_present) < 5) {
      warning("RiboSis 可用基因太少，跳过 ", name, " 的 RiboSis 评分。", immediate. = TRUE)
    } else {
      obj_liver <- AddModuleScore(
        object = obj_liver,
        features = list(ribosis_genes_present),
        name = "Ribosis"
      )
      
      if ("Ribosis1" %in% colnames(obj_liver@meta.data)) {
        p_ribo <- MySpatialPlot2(
          obj = obj_liver,
          feature = "Ribosis1",
          point_size = 2,
          alpha = 0.9,
          discrete = FALSE,
          color_option = "plasma",
          limits = ribo_limits
        )
        
        ggsave(
          filename = paste0(name, "_Ribo.pdf"),
          plot = p_ribo,
          width = 8,
          height = 4
        )
      } else {
        warning("Ribosis1 未出现在 meta.data 中，RiboSis 作图失败。")
      }
    }
    
  } else {
    warning("文件不存在: ", ribosis_rdata)
  }
  
  # --------------------------------------------------
  # 4.5 计算 NuS
  # --------------------------------------------------
  cat("\n====================\n")
  cat("开始计算 NuS score\n")
  cat("====================\n")
  
  nus_limits <- c(-0.15, 0.35)
  
  if (file.exists(nus_rdata)) {
    env_nus <- new.env()
    load(nus_rdata, envir = env_nus)
    
    cat("NuStress 文件中的对象:\n")
    print(ls(env_nus))
    
    if (!exists("geneSets_final", envir = env_nus)) {
      stop("NuStress 文件中未找到 geneSets_final 对象。")
    }
    
    geneSets_final <- get("geneSets_final", envir = env_nus)
    
    cat("geneSets_final 的名字:\n")
    print(names(geneSets_final))
    
    up_genes <- geneSets_final$NuStress_UP
    down_genes <- geneSets_final$NuStress_DOWN
    
    up_genes <- unique(as.character(unlist(up_genes)))
    down_genes <- unique(as.character(unlist(down_genes)))
    
    up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
    down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]
    
    obj_genes <- rownames(obj_liver)
    up_present <- up_genes[up_genes %in% obj_genes]
    down_present <- down_genes[down_genes %in% obj_genes]
    
    if (length(up_present) < 10) {
      up_idx <- match(toupper(up_genes), toupper(obj_genes))
      up_present <- unique(obj_genes[na.omit(up_idx)])
    }
    
    if (length(down_present) < 10) {
      down_idx <- match(toupper(down_genes), toupper(obj_genes))
      down_present <- unique(obj_genes[na.omit(down_idx)])
    }
    
    cat("obj_liver 中可用的 UP 基因数:", length(up_present), "\n")
    cat("obj_liver 中可用的 DOWN 基因数:", length(down_present), "\n")
    
    if (length(up_present) < 5 || length(down_present) < 5) {
      warning("NuS 上/下调基因集在当前对象中匹配到的基因太少，跳过 ", name, " 的 NuS 评分。", immediate. = TRUE)
    } else {
      obj_liver <- CalculateNuSScore(
        obj = obj_liver,
        up_genes = up_present,
        down_genes = down_present,
        score_name = "NuS_Score"
      )
      
      if ("NuS_Score" %in% colnames(obj_liver@meta.data)) {
        p_nus <- MySpatialPlot2(
          obj = obj_liver,
          feature = "NuS_Score",
          point_size = 2,
          alpha = 0.9,
          discrete = FALSE,
          color_option = "plasma",
          limits = nus_limits
        )
        
        ggsave(
          filename = paste0(name, "_NuS.pdf"),
          plot = p_nus,
          width = 8,
          height = 4
        )
      } else {
        warning("NuS_Score 未出现在 meta.data 中，NuS 作图失败。")
      }
    }
    
  } else {
    warning("文件不存在: ", nus_rdata)
  }
  
  # --------------------------------------------------
  # 4.6 保存结果
  # --------------------------------------------------
  cat("\n=== 分析完成:", name, " ===\n")
  cat("元数据列:\n")
  print(colnames(obj_liver@meta.data))
  
  if (all(c("Ribosis1", "NuS_Score", "manual.annot") %in% colnames(obj_liver@meta.data))) {
    data4 <- data.frame(
      Ribosis = obj_liver@meta.data$Ribosis1,
      NuS = obj_liver@meta.data$NuS_Score,
      name = obj_liver@meta.data$manual.annot
    )
    
    save(data4, file = paste0(name, "_data.rda"))
    cat("结果已保存:", paste0(name, "_data.rda"), "\n")
  } else {
    warning("缺少 Ribosis1 / NuS_Score / manual.annot 中的部分列，未保存 ", name, "_data.rda")
  }
  
  save(obj_liver, file = rda_file)
  rm(obj_liver)
}