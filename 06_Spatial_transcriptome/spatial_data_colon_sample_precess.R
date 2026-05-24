# Step3_manual annotation + RiboSis/NuS visualization and export (final version)
# ======================================================================
# 1. Create environment
# ======================================================================
getwd()
setwd("/home/xxm_xxm/CJX_workspace/SC-ST/")
rm(list = ls())
options(stringsAsFactors = FALSE)
options(timeout = 1e6)

# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(scales)
library(dplyr)

# ======================================================================
# 2. Parameter settings
# ======================================================================
# Absolute paths to gene set files
nus_rdata <- "/home/xxm_xxm/CJX_workspace/geneset/final_nucleolar_gene_sets/NuStress_geneSets_final.Rdata"
ribosis_rdata <- "/home/xxm_xxm/CJX_workspace/proteinomics/RiboSis activity.Rdata"

# Cell type colors (consistent with the final Step1 style)
cell_type_colors <- c(
  "Fibroblast" = "#1f77b4",
  "lamina propria" = "#ff7f0e",
  "Normal Epithelium" = "#2ca02c",
  "Tumor" = "#d62728",
  "Subserosa" = "#9467bd",
  "Smooth Muscle" = "#8491B4FF",
  "Endothelial" = "#F39B7FFF",
  "Hepatocytes" = "#00A087FF",
  "Erythrocyte" = "#7A52A3"
)

# Manual annotation mapping for each sample (from the original Step3)
cluster_mapping_list <- list(
  "ST-colon1" = c(
    "Tumor",
    "Tumor",
    "Tumor",
    "Fibroblast",
    "Normal Epithelium",
    "lamina propria",
    "Tumor",
    "lamina propria",
    "Tumor"
  ),
  
  "ST-colon2" = c(
    "Smooth Muscle",
    "Normal Epithelium",
    "Subserosa",
    "Smooth Muscle",
    "lamina propria",
    "Smooth Muscle",
    "lamina propria",
    "Tumor",
    "Smooth Muscle",
    "lamina propria",
    "Normal Epithelium",
    "Normal Epithelium",
    "Tumor",
    "Smooth Muscle",
    "Smooth Muscle",
    "Smooth Muscle"
  ),
  
  "ST-colon3" = c(
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Tumor",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Tumor",
    "Fibroblast"
  ),
  
  "ST-colon4" = c(
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Tumor",
    "Fibroblast",
    "Fibroblast",
    "Fibroblast",
    "Tumor",
    "Fibroblast"
  )
)

# Samples to process in batch
sample_names <- c("ST-colon1", "ST-colon2", "ST-colon3", "ST-colon4")

# ======================================================================
# 3. Shared functions
# ======================================================================

# ------------------------------
# 3.1 Spatial cluster plotting function
# ------------------------------
MySpatialClusterPlot <- function(obj, group_by, point_size = 2, alpha = 0.8,
                                 custom_colors = NULL, show_legend = TRUE) {
  
  if (length(obj@images) < 1) stop("No image information found in the object. Please check the spatial object.")
  if (!group_by %in% colnames(obj@meta.data)) {
    stop("The group_by column does not exist in meta.data. Available columns: ",
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
      warning("The following cell types do not have defined colors and will use default colors: ",
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
# 3.2 Spatial plotting function for continuous variables
# ------------------------------
MySpatialPlot2 <- function(obj, feature, point_size = 2, alpha = 0.8,
                           discrete = FALSE, color_option = "plasma",
                           limits = NULL, na.value = "grey80") {
  
  if (length(obj@images) < 1) stop("No image information found in the object.")
  if (!feature %in% colnames(obj@meta.data)) {
    stop("The feature column does not exist in meta.data. Available columns: ",
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
# 3.3 NuS scoring function
# ------------------------------
CalculateNuSScore <- function(obj, up_genes, down_genes, score_name = "NuS_Score") {
  
  if (length(up_genes) == 0) stop("up_genes cannot be empty")
  if (length(down_genes) == 0) stop("down_genes cannot be empty")
  
  available_genes <- rownames(obj)
  up_genes <- up_genes[up_genes %in% available_genes]
  down_genes <- down_genes[down_genes %in% available_genes]
  
  cat("UP genes:", length(up_genes), " available\n")
  cat("DOWN genes:", length(down_genes), " available\n")
  
  if (length(up_genes) == 0 || length(down_genes) == 0) {
    stop("Not enough genes are present in the data")
  }
  
  obj <- AddModuleScore(obj, features = list(up_genes), name = "UP_tmp")
  up_score <- obj@meta.data[["UP_tmp1"]]
  
  obj <- AddModuleScore(obj, features = list(down_genes), name = "DOWN_tmp")
  down_score <- obj@meta.data[["DOWN_tmp1"]]
  
  obj@meta.data[[score_name]] <- up_score - down_score
  
  obj@meta.data[["UP_tmp1"]] <- NULL
  obj@meta.data[["DOWN_tmp1"]] <- NULL
  
  cat("Score calculation completed and added to column:", score_name, "\n")
  return(obj)
}

# ======================================================================
# 4. Main loop: process samples one by one
# ======================================================================
for (name in sample_names) {
  
  cat("\n=====================================================\n")
  cat("Start processing sample:", name, "\n")
  cat("=====================================================\n")
  
  rda_file <- paste0(name, ".rda")
  if (!file.exists(rda_file)) {
    warning("Data file not found: ", rda_file, "; skipped.")
    next
  }
  
  load(rda_file)
  
  if (!exists("obj_colon")) {
    warning("obj_colon was not found after load; sample ", name, " skipped.")
    next
  }
  
  cat("Object type:", class(obj_colon), "\n")
  cat("Default assay:", DefaultAssay(obj_colon), "\n")
  cat("Number of images:", length(obj_colon@images), "\n")
  cat("Metadata column names:\n")
  print(colnames(obj_colon@meta.data))
  
  # --------------------------------------------------
  # 4.1 Save HE image
  # --------------------------------------------------
  p_he <- SpatialFeaturePlot(
    obj_colon,
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
  # 4.2 Visualize seurat_clusters
  # --------------------------------------------------
  cat("seurat_clusters distribution:\n")
  print(table(obj_colon$seurat_clusters))
  
  p_cluster <- MySpatialClusterPlot(
    obj_colon,
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
  # 4.3 Run additional marker analysis for ST-colon1 (keeps the original Step3 logic)
  # --------------------------------------------------
  Idents(obj_colon) <- "seurat_clusters"
  
  if (name == "ST-colon1") {
    cat("Run ", name, " FindAllMarkers...\n", sep = "")
    
    markers <- FindAllMarkers(
      obj_colon,
      assay = "SCT",
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    write.csv(markers, paste0(name, "_FindAllMarkers.csv"), row.names = FALSE)
    
    top20_markers <- markers %>%
      group_by(cluster) %>%
      top_n(20, avg_log2FC) %>%
      arrange(cluster, desc(avg_log2FC))
    
    write.csv(top20_markers, paste0(name, "_top20_markers.csv"), row.names = FALSE)
  }
  
  # --------------------------------------------------
  # 4.4 Manual annotation
  # --------------------------------------------------
  if (!name %in% names(cluster_mapping_list)) {
    warning("Manual annotation mapping was not found for sample ", name, "; skipped manual.annot.")
  } else {
    new.cluster.ids <- cluster_mapping_list[[name]]
    
    current_levels <- levels(obj_colon)
    cat("Current levels:\n")
    print(current_levels)
    
    if (length(new.cluster.ids) != length(current_levels)) {
      warning("Sample ", name, " has a new.cluster.ids length that does not match levels(obj_colon); skipped RenameIdents.")
    } else {
      names(new.cluster.ids) <- current_levels
      obj_colon <- RenameIdents(obj_colon, new.cluster.ids)
      obj_colon$manual.annot <- Idents(obj_colon)
      
      p_manual <- MySpatialClusterPlot(
        obj_colon,
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
  # 4.5 Calculate RiboSis
  # --------------------------------------------------
  cat("\n====================\n")
  cat("Start calculating RiboSis score\n")
  cat("====================\n")
  
  ribo_limits <- c(-0.15, 0.68)
  
  if (file.exists(ribosis_rdata)) {
    env_ribo <- new.env()
    load(ribosis_rdata, envir = env_ribo)
    
    cat("Objects in the RiboSis file:\n")
    print(ls(env_ribo))
    
    if (exists("ribosis", envir = env_ribo)) {
      ribosis_raw <- get("ribosis", envir = env_ribo)
    } else {
      obj_names <- ls(env_ribo)
      ribo_candidate <- obj_names[grepl("ribo", obj_names, ignore.case = TRUE)][1]
      
      if (is.na(ribo_candidate)) {
        stop("No object named ribosis was found in the RiboSis file, and no related object was matched automatically.")
      } else {
        ribosis_raw <- get(ribo_candidate, envir = env_ribo)
        cat("Automatically using object:", ribo_candidate, "\n")
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
    
    cat("Number of ribosis_genes after cleaning:", length(ribosis_genes), "\n")
    cat("First 10 genes:\n")
    print(head(ribosis_genes, 10))
    
    obj_genes <- rownames(obj_colon)
    ribosis_genes_present <- ribosis_genes[ribosis_genes %in% obj_genes]
    
    if (length(ribosis_genes_present) < 10) {
      cat("Too few genes matched directly; trying case-insensitive matching...\n")
      obj_genes_upper <- toupper(obj_genes)
      ribosis_genes_upper <- toupper(ribosis_genes)
      match_idx <- match(ribosis_genes_upper, obj_genes_upper)
      ribosis_genes_present <- unique(obj_genes[na.omit(match_idx)])
    }
    
    cat("Number of RiboSis genes present in obj_colon:", length(ribosis_genes_present), "\n")
    
    if (length(ribosis_genes_present) < 5) {
      warning("Too few available RiboSis genes; skipped RiboSis scoring for ", name, ".", immediate. = TRUE)
    } else {
      obj_colon <- AddModuleScore(
        object = obj_colon,
        features = list(ribosis_genes_present),
        name = "Ribosis"
      )
      
      if ("Ribosis1" %in% colnames(obj_colon@meta.data)) {
        p_ribo <- MySpatialPlot2(
          obj = obj_colon,
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
        warning("Ribosis1 was not found in meta.data; RiboSis plotting failed.")
      }
    }
    
  } else {
    warning("File does not exist: ", ribosis_rdata)
  }
  
  # --------------------------------------------------
  # 4.6 Calculate NuS
  # --------------------------------------------------
  cat("\n====================\n")
  cat("Start calculating NuS score\n")
  cat("====================\n")
  
  nus_limits <- c(-0.15, 0.35)
  
  if (file.exists(nus_rdata)) {
    env_nus <- new.env()
    load(nus_rdata, envir = env_nus)
    
    cat("Objects in the NuStress file:\n")
    print(ls(env_nus))
    
    if (!exists("geneSets_final", envir = env_nus)) {
      stop("Object geneSets_final was not found in the NuStress file.")
    }
    
    geneSets_final <- get("geneSets_final", envir = env_nus)
    
    cat("Names in geneSets_final:\n")
    print(names(geneSets_final))
    
    up_genes <- geneSets_final$NuStress_UP
    down_genes <- geneSets_final$NuStress_DOWN
    
    up_genes <- unique(as.character(unlist(up_genes)))
    down_genes <- unique(as.character(unlist(down_genes)))
    
    up_genes <- up_genes[!is.na(up_genes) & up_genes != ""]
    down_genes <- down_genes[!is.na(down_genes) & down_genes != ""]
    
    obj_genes <- rownames(obj_colon)
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
    
    cat("Number of available UP genes in obj_colon:", length(up_present), "\n")
    cat("Number of available DOWN genes in obj_colon:", length(down_present), "\n")
    
    if (length(up_present) < 5 || length(down_present) < 5) {
      warning("Too few NuS up/down genes matched in the current object; skipped NuS scoring for ", name, ".", immediate. = TRUE)
    } else {
      obj_colon <- CalculateNuSScore(
        obj = obj_colon,
        up_genes = up_present,
        down_genes = down_present,
        score_name = "NuS_Score"
      )
      
      if ("NuS_Score" %in% colnames(obj_colon@meta.data)) {
        p_nus <- MySpatialPlot2(
          obj = obj_colon,
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
        warning("NuS_Score was not found in meta.data; NuS plotting failed.")
      }
    }
    
  } else {
    warning("File does not exist: ", nus_rdata)
  }
  
  # --------------------------------------------------
  # 4.7 Save results
  # --------------------------------------------------
  cat("\n=== Analysis completed:", name, " ===\n")
  cat("Metadata columns:\n")
  print(colnames(obj_colon@meta.data))
  
  if (all(c("Ribosis1", "NuS_Score", "manual.annot") %in% colnames(obj_colon@meta.data))) {
    data4 <- data.frame(
      Ribosis = obj_colon@meta.data$Ribosis1,
      NuS = obj_colon@meta.data$NuS_Score,
      name = obj_colon@meta.data$manual.annot
    )
    
    save(data4, file = paste0(name, "_data.rda"))
    cat("Results saved:", paste0(name, "_data.rda"), "\n")
  } else {
    warning("Some columns among Ribosis1 / NuS_Score / manual.annot are missing; did not save ", name, "_data.rda")
  }
  
  # Save the new object with annotation and scores
  save(obj_colon, file = rda_file)
  
  # Clean objects to avoid contamination between loop iterations
  rm(obj_colon)
}
