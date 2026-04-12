###############################################################################
## TCGA + GTEx pan-cancer project
## Final integrated preprocessing pipeline
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 1e6)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(UCSCXenaTools)
})

###############################################################################
## 0. Paths and setup
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"

data_dir        <- file.path(base_dir, "raw_data", "xena_toil")
out_dir         <- file.path(base_dir, "processed_data")
clinical_base   <- file.path(base_dir, "clinical_xena")
clinical_rawdir <- file.path(clinical_base, "raw_download")
clinical_rdsdir <- file.path(clinical_base, "rds_by_cancer")
clinical_tabdir <- file.path(clinical_base, "merged_tables")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_base, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_rawdir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_rdsdir, recursive = TRUE, showWarnings = FALSE)
dir.create(clinical_tabdir, recursive = TRUE, showWarnings = FALSE)

expr_file       <- file.path(data_dir, "TcgaTargetGtex_rsem_gene_tpm.gz")
pheno_file      <- file.path(data_dir, "TcgaTargetGTEX_phenotype.txt.gz")
tcga_pheno_file <- file.path(data_dir, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
mut_file        <- file.path(data_dir, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz")
surv_file       <- file.path(clinical_rawdir, "Survival_SupplementalTable_S1_20171025_xena_sp")

stopifnot(file.exists(expr_file))
stopifnot(file.exists(pheno_file))
stopifnot(file.exists(tcga_pheno_file))
stopifnot(file.exists(mut_file))

###############################################################################
## 1. Helper functions
###############################################################################
pick_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
  df[[hit[1]]]
}

safe_numeric_matrix <- function(df) {
  mat <- as.matrix(df)
  mode(mat) <- "numeric"
  mat
}

qc_tab <- function(x, name = deparse(substitute(x))) {
  cat("\n====", name, "====\n")
  print(table(x, useNA = "ifany"))
}

tumor_pattern  <- regex("tumor|tumour|primary|metastatic|recurrent", ignore_case = TRUE)
normal_pattern <- regex("normal", ignore_case = TRUE)

abbr_from_cohort <- function(x) {
  sub(".*\\(([^()]*)\\).*", "\\1", x)
}

get_patient_id <- function(x) {
  ifelse(str_detect(x, "^TCGA-"), substr(x, 1, 12), NA_character_)
}

get_sample_id_short <- function(x) {
  ifelse(str_detect(x, "^TCGA-"), substr(x, 1, 16), x)
}

find_sample_col <- function(df) {
  cn <- colnames(df)
  hit <- cn[tolower(cn) %in% c(
    "sample", "sampleid", "sample_id",
    "submitter_id", "cases.submitter_id",
    "bcr_patient_barcode", "bcr_sample_barcode"
  )]
  if (length(hit) > 0) return(hit[1])
  
  hit2 <- cn[str_detect(tolower(cn), "sample|submitter|barcode")]
  if (length(hit2) > 0) return(hit2[1])
  
  stop("No sample-like column found.")
}

normalize_sample_col <- function(df) {
  sample_col <- find_sample_col(df)
  df %>%
    rename(sample_id = all_of(sample_col)) %>%
    mutate(
      sample_id = as.character(sample_id),
      patient_id = get_patient_id(sample_id),
      sample_id_short = get_sample_id_short(sample_id)
    )
}

find_treatment_columns <- function(df) {
  cn <- colnames(df)
  patt <- paste(
    c(
      "treat", "therapy", "therap", "drug", "pharm",
      "radiat", "neoadjuvant", "adjuvant", "prior",
      "postop", "preop", "targeted", "immuno",
      "chemotherapy", "chemother", "rx"
    ),
    collapse = "|"
  )
  cn[str_detect(tolower(cn), patt)]
}

classify_treated_value <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  case_when(
    is.na(x0) ~ NA,
    x0 %in% c("", "na", "n/a", "not reported", "not available", "unknown", "[not available]") ~ NA,
    str_detect(x0, "yes|true|treated|received|prior|completed|ongoing|administered") ~ TRUE,
    str_detect(x0, "no|false|untreated|none|not received|no prior") ~ FALSE,
    TRUE ~ NA
  )
}

collapse_treatment_flags <- function(df_flag) {
  apply(df_flag, 1, function(z) {
    z <- as.logical(z)
    if (all(is.na(z))) return(NA)
    if (any(z %in% TRUE, na.rm = TRUE)) return(TRUE)
    if (all(z %in% FALSE | is.na(z))) return(FALSE)
    NA
  })
}

safe_select_dataset <- function(df) {
  if (nrow(df) == 0) return(df[0, , drop = FALSE])
  
  out <- df %>%
    filter(
      str_detect(tolower(Label), "phenotype|clinical"),
      str_detect(tolower(Type), "clinicalmatrix|phenotype")
    )
  
  if (nrow(out) == 0) return(out)
  
  if ("XenaHostNames" %in% colnames(out)) {
    out <- out %>%
      mutate(host_rank = case_when(
        str_detect(tolower(XenaHostNames), "gdc") ~ 1,
        str_detect(tolower(XenaHostNames), "tcga") ~ 2,
        TRUE ~ 9
      ))
  } else {
    out$host_rank <- 9
  }
  
  if ("SampleCount" %in% colnames(out)) {
    out$SampleCount <- suppressWarnings(as.numeric(out$SampleCount))
  } else {
    out$SampleCount <- NA_real_
  }
  
  out %>% arrange(host_rank, desc(SampleCount))
}

download_one_xena_dataset <- function(dataset_id, destdir) {
  query <- XenaGenerate(subset = XenaDatasets == dataset_id)
  XenaQuery(query) %>%
    XenaDownload(destdir = destdir, download_probeMap = FALSE) %>%
    XenaPrepare()
}

choose_first_non_na_safe <- function(x) {
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x_chr <- x_chr[!is.na(x_chr)]
  if (length(x_chr) == 0) return(NA_character_)
  x_chr[1]
}

###############################################################################
## 2. Expression preprocessing
###############################################################################
cat("=== Step 1: Read expression matrix ===\n")
expr_dt <- fread(expr_file, data.table = FALSE)
colnames(expr_dt)[1] <- "Gene"

expr_dt <- dplyr::filter(expr_dt, !is.na(Gene), Gene != "")

expr_num <- expr_dt[, -1, drop = FALSE]
expr_num[] <- lapply(expr_num, function(x) suppressWarnings(as.numeric(as.character(x))))

expr_dt$mean_expr <- rowMeans(expr_num, na.rm = TRUE)

expr_dt <- dplyr::arrange(expr_dt, Gene, dplyr::desc(mean_expr))
expr_dt <- dplyr::distinct(expr_dt, Gene, .keep_all = TRUE)
expr_dt <- dplyr::select(expr_dt, -mean_expr)

expr_mat <- safe_numeric_matrix(expr_dt[, -1, drop = FALSE])
rownames(expr_mat) <- expr_dt$Gene

expr_log2 <- log2(expr_mat + 1)
saveRDS(expr_log2, file.path(out_dir, "expr_log2_TcgaTargetGtex_rsem_tpm.rds"))
###############################################################################
## 3. Combined phenotype preprocessing
###############################################################################
cat("\n=== Step 2: Read combined phenotype ===\n")
pheno <- fread(pheno_file, data.table = FALSE)
colnames(pheno)[1] <- "sample"

sample_annot <- data.frame(
  sample = pheno$sample,
  study  = pick_first_existing(pheno, c("_study", "study", "cohort")),
  sample_type_raw = pick_first_existing(pheno, c("_sample_type", "sample_type")),
  primary_disease = pick_first_existing(pheno, c("_primary_disease", "primary_disease")),
  primary_site    = pick_first_existing(pheno, c("_primary_site", "primary_site")),
  tissue          = pick_first_existing(pheno, c("SMTSD", "SMTS", "_primary_site")),
  stringsAsFactors = FALSE
) %>%
  mutate(
    source = case_when(
      str_detect(sample, "^TCGA-") ~ "TCGA",
      str_detect(sample, "^GTEX-") ~ "GTEx",
      TRUE ~ ifelse(!is.na(study), study, "Other")
    ),
    patient_id = case_when(
      source == "TCGA" ~ substr(sample, 1, 12),
      source == "GTEx" ~ str_extract(sample, "^[^-]+-[^-]+"),
      TRUE ~ sample
    ),
    sample_id_short = case_when(
      source == "TCGA" ~ substr(sample, 1, 16),
      TRUE ~ sample
    ),
    group = case_when(
      source == "GTEx" ~ "Normal",
      str_detect(coalesce(sample_type_raw, ""), normal_pattern) ~ "Normal",
      str_detect(coalesce(sample_type_raw, ""), tumor_pattern) ~ "Tumor",
      TRUE ~ NA_character_
    )
  )

write.csv(sample_annot, file.path(out_dir, "sample_annotation_basic.csv"), row.names = FALSE)

###############################################################################
## 4. TCGA phenotype merge
###############################################################################
cat("\n=== Step 3: Read TCGA phenotype ===\n")
tcga_pheno <- fread(tcga_pheno_file, data.table = FALSE)
colnames(tcga_pheno)[1] <- "sample"

tcga_annot <- data.frame(
  sample = tcga_pheno$sample,
  tcga_sample_type     = pick_first_existing(tcga_pheno, c("_sample_type", "sample_type")),
  tcga_primary_disease = pick_first_existing(tcga_pheno, c("_primary_disease", "primary_disease")),
  tcga_primary_site    = pick_first_existing(tcga_pheno, c("_primary_site", "primary_site")),
  stringsAsFactors = FALSE
)

sample_annot2 <- sample_annot %>%
  left_join(tcga_annot, by = "sample") %>%
  mutate(
    cancer_type = case_when(
      source == "TCGA" ~ coalesce(tcga_primary_disease, primary_disease, primary_site),
      source == "GTEx" ~ tissue,
      TRUE ~ primary_disease
    ),
    sample_type_final = case_when(
      source == "TCGA" ~ coalesce(tcga_sample_type, sample_type_raw),
      TRUE ~ sample_type_raw
    ),
    group = case_when(
      source == "GTEx" ~ "Normal",
      str_detect(coalesce(sample_type_final, ""), normal_pattern) ~ "Normal",
      str_detect(coalesce(sample_type_final, ""), tumor_pattern) ~ "Tumor",
      TRUE ~ group
    )
  )

write.csv(sample_annot2, file.path(out_dir, "sample_annotation_merged.csv"), row.names = FALSE)

###############################################################################
## 5. TP53 mutation extraction
###############################################################################
cat("\n=== Step 4: Read mutation matrix and extract TP53 ===\n")
mut_dt <- fread(mut_file, data.table = FALSE)
colnames(mut_dt)[1] <- "Gene"

mut_dt <- mut_dt %>%
  filter(!is.na(Gene), Gene != "")

if (!"TP53" %in% mut_dt$Gene) stop("TP53 not found in Gene column.")

tp53_rows <- mut_dt[mut_dt$Gene == "TP53", , drop = FALSE]
tp53_mat  <- safe_numeric_matrix(tp53_rows[, -1, drop = FALSE])
tp53_vec  <- apply(tp53_mat, 2, max, na.rm = TRUE)

tp53_status <- data.frame(
  sample = colnames(tp53_rows)[-1],
  TP53 = as.numeric(tp53_vec),
  stringsAsFactors = FALSE
) %>%
  mutate(
    TP53_mut = case_when(
      TP53 == 1 ~ "Mut",
      TP53 == 0 ~ "WT",
      TRUE ~ NA_character_
    )
  )

sample_annot3 <- sample_annot2 %>%
  left_join(tp53_status, by = "sample") %>%
  mutate(TP53_call_method = "mc3_nonsilentGene_xena_collapsed_by_max_across_duplicated_TP53_rows")

write.csv(sample_annot3, file.path(out_dir, "sample_annotation_with_TP53.csv"), row.names = FALSE)
saveRDS(sample_annot3, file.path(out_dir, "sample_annotation_with_TP53.rds"))

###############################################################################
## 6. Match expression and annotation
###############################################################################
cat("\n=== Step 5: Match expression and annotation ===\n")
common_samples <- intersect(colnames(expr_log2), sample_annot3$sample)

expr_log2_sub <- expr_log2[, common_samples, drop = FALSE]
annot_final <- sample_annot3 %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

stopifnot(identical(colnames(expr_log2_sub), annot_final$sample))

saveRDS(expr_log2_sub, file.path(out_dir, "expr_log2_matched.rds"))
saveRDS(annot_final,   file.path(out_dir, "sample_annotation_final.rds"))
write.csv(annot_final, file.path(out_dir, "sample_annotation_final.csv"), row.names = FALSE)

###############################################################################
## 7. Final cleanup for TCGA + GTEx
###############################################################################
cat("\n=== Step 6: Final cleanup for downstream analysis ===\n")
annot_clean <- annot_final %>%
  mutate(
    source = case_when(
      source %in% c("GTEx", "GTEX") ~ "GTEx",
      TRUE ~ source
    )
  )

annot_tcga_gtex <- annot_clean %>%
  filter(source %in% c("TCGA", "GTEx"))

expr_tcga_gtex <- expr_log2_sub[, annot_tcga_gtex$sample, drop = FALSE]

annot_main <- annot_tcga_gtex %>%
  filter((source == "TCGA" & group == "Tumor") | (source == "GTEx" & group == "Normal"))

expr_main <- expr_log2_sub[, annot_main$sample, drop = FALSE]

annot_tcga_tumor2 <- annot_tcga_gtex %>%
  filter(source == "TCGA", group == "Tumor")

expr_tcga_tumor2 <- expr_log2_sub[, annot_tcga_tumor2$sample, drop = FALSE]

saveRDS(annot_tcga_gtex, file.path(out_dir, "annot_tcga_gtex_only.rds"))
saveRDS(expr_tcga_gtex,  file.path(out_dir, "expr_tcga_gtex_only.rds"))
saveRDS(annot_main,      file.path(out_dir, "annot_main_tcga_tumor_gtex_normal.rds"))
saveRDS(expr_main,       file.path(out_dir, "expr_main_tcga_tumor_gtex_normal.rds"))
saveRDS(annot_tcga_tumor2, file.path(out_dir, "annot_tcga_tumor_only_final.rds"))
saveRDS(expr_tcga_tumor2,  file.path(out_dir, "expr_tcga_tumor_only_final.rds"))

write.csv(annot_tcga_gtex, file.path(out_dir, "annot_tcga_gtex_only.csv"), row.names = FALSE)
write.csv(annot_main, file.path(out_dir, "annot_main_tcga_tumor_gtex_normal.csv"), row.names = FALSE)
write.csv(annot_tcga_tumor2, file.path(out_dir, "annot_tcga_tumor_only_final.csv"), row.names = FALSE)
###############################################################################
## 8. Download TCGA clinical phenotype via UCSCXenaTools
###############################################################################
cat("\n=== Step 7: Load UCSC Xena dataset catalog ===\n")

all_datasets <- XenaData %>%
  dplyr::filter(stringr::str_detect(XenaCohorts, "TCGA")) %>%
  as.data.frame()

write.csv(
  all_datasets,
  file.path(clinical_tabdir, "UCSCXena_TCGA_dataset_catalog.csv"),
  row.names = FALSE
)

cancer_cohorts <- c(
  "GDC TCGA Mesothelioma (MESO)",
  "GDC TCGA Stomach Cancer (STAD)",
  "GDC TCGA Ocular melanomas (UVM)",
  "GDC TCGA Pancreatic Cancer (PAAD)",
  "GDC TCGA Colon Cancer (COAD)",
  "GDC TCGA Thymoma (THYM)",
  "GDC TCGA Ovarian Cancer (OV)",
  "GDC TCGA Kidney Clear Cell Carcinoma (KIRC)",
  "GDC TCGA Lower Grade Glioma (LGG)",
  "GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)",
  "GDC TCGA Kidney Chromophobe (KICH)",
  "GDC TCGA Thyroid Cancer (THCA)",
  "GDC TCGA Pheochromocytoma & Paraganglioma (PCPG)",
  "GDC TCGA Rectal Cancer (READ)",
  "GDC TCGA Cervical Cancer (CESC)",
  "GDC TCGA Sarcoma (SARC)",
  "GDC TCGA Acute Myeloid Leukemia (LAML)",
  "GDC TCGA Bladder Cancer (BLCA)",
  "GDC TCGA Lung Squamous Cell Carcinoma (LUSC)",
  "GDC TCGA Esophageal Cancer (ESCA)",
  "GDC TCGA Lung Adenocarcinoma (LUAD)",
  "GDC TCGA Testicular Cancer (TGCT)",
  "GDC TCGA Glioblastoma (GBM)",
  "GDC TCGA Breast Cancer (BRCA)",
  "GDC TCGA Uterine Carcinosarcoma (UCS)",
  "GDC TCGA Endometrioid Cancer (UCEC)",
  "GDC TCGA Bile Duct Cancer (CHOL)",
  "GDC TCGA Adrenocortical Cancer (ACC)",
  "GDC TCGA Liver Cancer (LIHC)",
  "GDC TCGA Head and Neck Cancer (HNSC)",
  "GDC TCGA Prostate Cancer (PRAD)",
  "GDC TCGA Large B-cell Lymphoma (DLBC)",
  "GDC TCGA Melanoma (SKCM)"
)

clinical_list <- list()
download_log  <- list()

for (cohort_name in cancer_cohorts) {
  cancer_abbr <- abbr_from_cohort(cohort_name)
  rds_path <- file.path(clinical_rdsdir, paste0(cancer_abbr, "_clinical.rds"))
  csv_path <- file.path(clinical_rdsdir, paste0(cancer_abbr, "_clinical.csv"))
  
  cat("\n====================================================\n")
  cat("Processing clinical phenotype:", cohort_name, "[", cancer_abbr, "]\n")
  
  ## 优先复用已存在文件
  if (file.exists(rds_path)) {
    cat("Existing clinical RDS found, reusing:", rds_path, "\n")
    clinical_df <- readRDS(rds_path)
    clinical_list[[cancer_abbr]] <- clinical_df
    
    download_log[[cancer_abbr]] <- tibble::tibble(
      cancer_abbr = cancer_abbr,
      cohort_name = cohort_name,
      clinical_dataset = NA_character_,
      clinical_status = "reused_existing_rds"
    )
    next
  }
  
  one <- all_datasets %>%
    dplyr::filter(XenaCohorts == cohort_name) %>%
    as.data.frame()
  
  if (nrow(one) == 0) {
    download_log[[cancer_abbr]] <- tibble::tibble(
      cancer_abbr = cancer_abbr,
      cohort_name = cohort_name,
      clinical_dataset = NA_character_,
      clinical_status = "not_found"
    )
    next
  }
  
  clinical_candidates <- safe_select_dataset(one)
  clinical_status <- "not_found"
  clinical_dataset_id <- NA_character_
  clinical_df <- NULL
  
  if (nrow(clinical_candidates) > 0) {
    clinical_dataset_id <- clinical_candidates$XenaDatasets[1]
    cat("Clinical dataset:", clinical_dataset_id, "\n")
    
    clinical_df <- tryCatch(
      download_one_xena_dataset(clinical_dataset_id, clinical_rawdir),
      error = function(e) {
        message("Clinical download failed: ", e$message)
        NULL
      }
    )
    
    if (!is.null(clinical_df)) {
      clinical_df <- tryCatch({
        normalize_sample_col(as.data.frame(clinical_df)) %>%
          dplyr::mutate(
            cancer_abbr = cancer_abbr,
            cohort_name = cohort_name
          )
      }, error = function(e) {
        message("Clinical normalization failed: ", e$message)
        NULL
      })
      
      if (!is.null(clinical_df)) {
        saveRDS(clinical_df, rds_path)
        write.csv(clinical_df, csv_path, row.names = FALSE)
        clinical_list[[cancer_abbr]] <- clinical_df
        clinical_status <- "downloaded"
      }
    }
  }
  
  download_log[[cancer_abbr]] <- tibble::tibble(
    cancer_abbr = cancer_abbr,
    cohort_name = cohort_name,
    clinical_dataset = clinical_dataset_id,
    clinical_status = clinical_status
  )
}

download_log_df <- dplyr::bind_rows(download_log)

write.csv(
  download_log_df,
  file.path(clinical_tabdir, "download_log_clinical.csv"),
  row.names = FALSE
)

###############################################################################
## 9. Download / read PanCanAtlas survival
###############################################################################
cat("\n=== Step 8: Download PanCanAtlas survival directly ===\n")

surv_url <- "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp"

if (!file.exists(surv_file)) {
  cat("Downloading PanCanAtlas survival file...\n")
  tryCatch(
    {
      download.file(
        url = surv_url,
        destfile = surv_file,
        mode = "wb",
        method = "libcurl"
      )
    },
    error = function(e) {
      stop(
        "Failed to download survival file.\n",
        "Please manually download:\n",
        surv_url, "\n",
        "and save it to:\n",
        surv_file, "\n\n",
        "Original error: ", e$message
      )
    }
  )
} else {
  cat("Survival file already exists:", surv_file, "\n")
}

survival_all <- fread(surv_file, data.table = FALSE)

## 不再用 rename，避免包冲突
if ("sample" %in% colnames(survival_all)) {
  colnames(survival_all)[colnames(survival_all) == "sample"] <- "sample_id"
} else {
  colnames(survival_all)[1] <- "sample_id"
}

survival_all <- survival_all %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    sample_id_short = get_sample_id_short(sample_id),
    patient_id = if ("_PATIENT" %in% colnames(.)) as.character(`_PATIENT`) else get_patient_id(sample_id),
    cancer_abbr = if ("cancer type abbreviation" %in% colnames(.)) as.character(`cancer type abbreviation`) else NA_character_,
    cohort_name = ifelse(!is.na(cancer_abbr), paste0("TCGA-", cancer_abbr), NA_character_)
  )

saveRDS(
  survival_all,
  file.path(clinical_tabdir, "TCGA_survival_all_samples_pancan.rds")
)

write.csv(
  survival_all,
  file.path(clinical_tabdir, "TCGA_survival_all_samples_pancan.csv"),
  row.names = FALSE
)

###############################################################################
## 10. Merge all clinical phenotype tables
###############################################################################
cat("\n=== Step 9: Merge all clinical phenotype tables ===\n")

if (length(clinical_list) == 0) {
  existing_rds <- list.files(
    path = clinical_rdsdir,
    pattern = "_clinical\\.rds$",
    full.names = TRUE
  )
  
  if (length(existing_rds) > 0) {
    cat("No newly downloaded clinical data, loading existing RDS files from disk...\n")
    clinical_list <- setNames(lapply(existing_rds, readRDS), basename(existing_rds))
  }
}

if (length(clinical_list) == 0) {
  stop("No clinical phenotype data available. Check downloads or existing *_clinical.rds files.")
}

clinical_all <- dplyr::bind_rows(clinical_list)

saveRDS(
  clinical_all,
  file.path(clinical_tabdir, "TCGA_clinical_all_samples.rds")
)

write.csv(
  clinical_all,
  file.path(clinical_tabdir, "TCGA_clinical_all_samples.csv"),
  row.names = FALSE
)
###############################################################################
## 11. Infer treatment-related variables
###############################################################################
cat("\n=== Step 10: Infer treatment-related variables ===\n")

treat_cols <- find_treatment_columns(clinical_all)

write.csv(
  data.frame(treatment_candidate_columns = treat_cols),
  file.path(clinical_tabdir, "treatment_candidate_columns.csv"),
  row.names = FALSE
)

if (length(treat_cols) > 0) {
  
  ## 不再用 count()，避免包冲突
  treat_value_summary_list <- lapply(treat_cols, function(cc) {
    raw_value <- as.character(clinical_all[[cc]])
    raw_value <- ifelse(is.na(raw_value), NA_character_, trimws(raw_value))
    
    tb <- table(raw_value, useNA = "ifany")
    out <- data.frame(
      column = cc,
      raw_value = names(tb),
      n = as.integer(tb),
      stringsAsFactors = FALSE
    )
    
    out$raw_value[is.na(out$raw_value) | out$raw_value == "<NA>"] <- NA_character_
    out <- out[order(-out$n), , drop = FALSE]
    rownames(out) <- NULL
    out
  })
  
  treat_value_summary <- do.call(rbind, treat_value_summary_list)
  
  write.csv(
    treat_value_summary,
    file.path(clinical_tabdir, "treatment_candidate_value_summary.csv"),
    row.names = FALSE
  )
  
  treat_flags_list <- lapply(treat_cols, function(cc) {
    classify_treated_value(clinical_all[[cc]])
  })
  
  treat_flags_df <- as.data.frame(treat_flags_list, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(treat_flags_df) <- paste0("flag__", treat_cols)
  
  clinical_all_treat <- dplyr::bind_cols(clinical_all, treat_flags_df) %>%
    dplyr::mutate(
      treated_flag_sample = collapse_treatment_flags(treat_flags_df),
      treated_flag_detail = dplyr::case_when(
        treated_flag_sample %in% TRUE ~ "Treated",
        treated_flag_sample %in% FALSE ~ "Untreated_or_no_evidence_of_treatment",
        TRUE ~ "Unknown"
      )
    )
  
} else {
  
  clinical_all_treat <- clinical_all %>%
    dplyr::mutate(
      treated_flag_sample = NA,
      treated_flag_detail = "Unknown"
    )
}

saveRDS(
  clinical_all_treat,
  file.path(clinical_tabdir, "TCGA_clinical_all_samples_with_treatment.rds")
)

write.csv(
  clinical_all_treat,
  file.path(clinical_tabdir, "TCGA_clinical_all_samples_with_treatment.csv"),
  row.names = FALSE
)

###############################################################################
## 12. Merge clinical phenotype with PanCanAtlas survival
###############################################################################
cat("\n=== Step 11: Build sample-level merged clinical table ===\n")

join_keys <- c("sample_id", "patient_id", "sample_id_short", "cancer_abbr")
missing_in_clin <- setdiff(join_keys, colnames(clinical_all_treat))
missing_in_surv <- setdiff(join_keys, colnames(survival_all))

if (length(missing_in_clin) > 0) {
  stop("Missing join keys in clinical_all_treat: ", paste(missing_in_clin, collapse = ", "))
}
if (length(missing_in_surv) > 0) {
  stop("Missing join keys in survival_all: ", paste(missing_in_surv, collapse = ", "))
}

clinical_surv_sample <- dplyr::full_join(
  clinical_all_treat,
  survival_all,
  by = join_keys,
  suffix = c(".clinical", ".survival")
)

saveRDS(
  clinical_surv_sample,
  file.path(clinical_tabdir, "TCGA_clinical_survival_merged_sample_level.rds")
)

write.csv(
  clinical_surv_sample,
  file.path(clinical_tabdir, "TCGA_clinical_survival_merged_sample_level.csv"),
  row.names = FALSE
)

###############################################################################
## 13. Build patient-level merged clinical table
###############################################################################
cat("\n=== Step 12: Build patient-level merged clinical table ===\n")

collapse_patient_df <- function(df) {
  data.frame(
    n_samples = nrow(df),
    sample_ids = paste(unique(df$sample_id), collapse = ";"),
    sample_id_short_all = paste(unique(df$sample_id_short), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

non_id_cols <- setdiff(
  colnames(clinical_surv_sample),
  c("sample_id", "sample_id_short", "patient_id", "cancer_abbr", "cohort_name", "cohort")
)

patient_core <- clinical_surv_sample %>%
  dplyr::group_by(cancer_abbr, patient_id) %>%
  dplyr::group_modify(~ collapse_patient_df(.x)) %>%
  dplyr::ungroup()

patient_rest <- clinical_surv_sample %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(non_id_cols), as.character)) %>%
  dplyr::group_by(cancer_abbr, patient_id) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(non_id_cols), choose_first_non_na_safe),
    .groups = "drop"
  )

clinical_surv_patient <- dplyr::left_join(
  patient_core,
  patient_rest,
  by = c("cancer_abbr", "patient_id")
)

saveRDS(
  clinical_surv_patient,
  file.path(clinical_tabdir, "TCGA_clinical_survival_merged_patient_level.rds")
)

write.csv(
  clinical_surv_patient,
  file.path(clinical_tabdir, "TCGA_clinical_survival_merged_patient_level.csv"),
  row.names = FALSE
)

###############################################################################
## 14. Build simplified patient-level analysis-ready clinical summary
###############################################################################
cat("\n=== Step 13: Build simplified patient-level clinical summary ===\n")

analysis_clinical_patient <- clinical_surv_patient %>%
  dplyr::mutate(
    age = suppressWarnings(as.numeric(age_at_index.demographic)),
    sex = as.character(gender.demographic),
    
    stage_pathologic = as.character(ajcc_pathologic_stage.diagnoses),
    stage_clinical   = as.character(ajcc_clinical_stage.diagnoses),
    stage = dplyr::coalesce(stage_pathologic, stage_clinical),
    
    vital_status = as.character(vital_status.demographic),
    days_to_death = suppressWarnings(as.numeric(days_to_death.demographic)),
    days_to_last_follow_up = suppressWarnings(as.numeric(days_to_last_follow_up.diagnoses)),
    
    os = suppressWarnings(as.numeric(OS)),
    os_time = suppressWarnings(as.numeric(OS.time)),
    dss = suppressWarnings(as.numeric(DSS)),
    dss_time = suppressWarnings(as.numeric(DSS.time)),
    dfi = suppressWarnings(as.numeric(DFI)),
    dfi_time = suppressWarnings(as.numeric(DFI.time)),
    pfi = suppressWarnings(as.numeric(PFI)),
    pfi_time = suppressWarnings(as.numeric(PFI.time)),
    
    treated_flag_sample = as.character(treated_flag_sample),
    treated_flag_detail = as.character(treated_flag_detail),
    
    prior_treatment_raw = if ("prior_treatment.diagnoses" %in% colnames(clinical_surv_patient)) {
      as.character(prior_treatment.diagnoses)
    } else {
      NA_character_
    }
  ) %>%
  dplyr::select(
    dplyr::any_of(c(
      "cancer_abbr",
      "patient_id",
      "n_samples",
      "sample_ids",
      "sample_id_short_all",
      "primary_site",
      "disease_type",
      "primary_diagnosis.diagnoses",
      "tissue_or_organ_of_origin.diagnoses",
      "age",
      "sex",
      "race.demographic",
      "ethnicity.demographic",
      "stage",
      "stage_pathologic",
      "stage_clinical",
      "vital_status",
      "os",
      "os_time",
      "dss",
      "dss_time",
      "dfi",
      "dfi_time",
      "pfi",
      "pfi_time",
      "days_to_death",
      "days_to_last_follow_up",
      "treated_flag_sample",
      "treated_flag_detail",
      "prior_treatment_raw"
    )),
    dplyr::everything()
  )

saveRDS(
  analysis_clinical_patient,
  file.path(clinical_tabdir, "TCGA_patient_clinical_analysis_ready.rds")
)

write.csv(
  analysis_clinical_patient,
  file.path(clinical_tabdir, "TCGA_patient_clinical_analysis_ready.csv"),
  row.names = FALSE
)

###############################################################################
## 15. Merge clinical info back to expression annotation
###############################################################################
cat("\n=== Step 14: Merge clinical info back to expression annotation ===\n")

patient_merge_df <- analysis_clinical_patient %>%
  dplyr::select(
    dplyr::any_of(c(
      "patient_id",
      "cancer_abbr",
      "primary_site",
      "disease_type",
      "primary_diagnosis.diagnoses",
      "age",
      "sex",
      "stage",
      "vital_status",
      "os",
      "os_time",
      "dss",
      "dss_time",
      "dfi",
      "dfi_time",
      "pfi",
      "pfi_time",
      "treated_flag_sample",
      "treated_flag_detail",
      "prior_treatment_raw"
    ))
  ) %>%
  dplyr::distinct(patient_id, .keep_all = TRUE)

annot_tcga_tumor2_clin_patient <- dplyr::left_join(
  annot_tcga_tumor2,
  patient_merge_df,
  by = "patient_id"
)

saveRDS(
  annot_tcga_tumor2_clin_patient,
  file.path(out_dir, "annot_tcga_tumor_only_with_patient_clinical.rds")
)

write.csv(
  annot_tcga_tumor2_clin_patient,
  file.path(out_dir, "annot_tcga_tumor_only_with_patient_clinical.csv"),
  row.names = FALSE
)

annot_main_with_clinical <- dplyr::left_join(
  annot_main,
  patient_merge_df,
  by = "patient_id"
)

saveRDS(
  annot_main_with_clinical,
  file.path(out_dir, "annot_main_tcga_tumor_gtex_normal_with_clinical.rds")
)

write.csv(
  annot_main_with_clinical,
  file.path(out_dir, "annot_main_tcga_tumor_gtex_normal_with_clinical.csv"),
  row.names = FALSE
)

###############################################################################
## 16. QC summaries
###############################################################################
cat("\n=== Final QC summaries ===\n")

qc_tab(annot_clean$source, "annot_clean$source")
qc_tab(annot_main$group, "annot_main$group")
qc_tab(annot_tcga_tumor2$TP53_mut, "annot_tcga_tumor2$TP53_mut")

if ("treated_flag_detail" %in% colnames(annot_tcga_tumor2_clin_patient)) {
  qc_tab(
    annot_tcga_tumor2_clin_patient$treated_flag_detail,
    "annot_tcga_tumor2_clin_patient$treated_flag_detail"
  )
}

cat("\nExpression dimensions:\n")
cat("expr_main:", nrow(expr_main), "genes x", ncol(expr_main), "samples\n")
cat("expr_tcga_tumor2:", nrow(expr_tcga_tumor2), "genes x", ncol(expr_tcga_tumor2), "samples\n")

cat("\nClinical dimensions:\n")
cat("clinical_all:", dim(clinical_all)[1], "x", dim(clinical_all)[2], "\n")
cat("survival_all:", dim(survival_all)[1], "x", dim(survival_all)[2], "\n")
cat("clinical_surv_sample:", dim(clinical_surv_sample)[1], "x", dim(clinical_surv_sample)[2], "\n")
cat("clinical_surv_patient:", dim(clinical_surv_patient)[1], "x", dim(clinical_surv_patient)[2], "\n")
cat("analysis_clinical_patient:", dim(analysis_clinical_patient)[1], "x", dim(analysis_clinical_patient)[2], "\n")

cat("\nOutput directories:\n")
cat("Processed expression/annotation:", out_dir, "\n")
cat("Clinical download/merge:", clinical_base, "\n")
cat("=== Integrated preprocessing completed successfully ===\n")

print(dim(analysis_clinical_patient))
print(table(analysis_clinical_patient$treated_flag_detail, useNA = "ifany"))
print(sum(!is.na(analysis_clinical_patient$os)))
print(sum(!is.na(analysis_clinical_patient$os_time)))
print(dim(annot_tcga_tumor2_clin_patient))
print(table(annot_tcga_tumor2_clin_patient$treated_flag_detail, useNA = "ifany"))