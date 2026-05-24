###############################################################################
## Step 1. Download TCGA + GTEx raw data from UCSC Xena
## Final robust version
###############################################################################

rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 1e6)

suppressPackageStartupMessages({
  library(data.table)
})

###############################################################################
## 0. Paths
###############################################################################
base_dir <- "/home/xxm_xxm/CJX_workspace/TCGA_GTEx_project"

raw_dir <- file.path(base_dir, "raw_data", "xena_toil")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 1. Helper functions
###############################################################################
timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

msg <- function(...) {
  cat(sprintf("[%s] ", timestamp()), ..., "\n", sep = "")
  flush.console()
}

download_one_file <- function(url,
                              destfile,
                              overwrite = FALSE,
                              min_size_bytes = 1000,
                              retries = 2) {
  if (file.exists(destfile) && !overwrite) {
    finfo <- file.info(destfile)
    if (!is.na(finfo$size) && finfo$size >= min_size_bytes) {
      msg("Skip existing file: ", basename(destfile),
          " (", format(finfo$size, big.mark = ","), " bytes)")
      return(invisible(TRUE))
    } else {
      msg("Existing file seems incomplete, removing: ", basename(destfile))
      unlink(destfile)
    }
  }
  
  success <- FALSE
  
  for (k in seq_len(retries + 1)) {
    msg("Downloading [attempt ", k, "]: ", basename(destfile))
    
    try({
      download.file(
        url = url,
        destfile = destfile,
        mode = "wb",
        method = "libcurl",
        quiet = FALSE
      )
      success <- TRUE
    }, silent = TRUE)
    
    if (success && file.exists(destfile)) {
      finfo <- file.info(destfile)
      if (!is.na(finfo$size) && finfo$size >= min_size_bytes) {
        msg("Finished: ", basename(destfile),
            " (", format(finfo$size, big.mark = ","), " bytes)")
        return(invisible(TRUE))
      }
    }
    
    msg("Download attempt failed for: ", basename(destfile))
    if (file.exists(destfile)) unlink(destfile)
  }
  
  stop("Download failed after retries: ", url)
}

quick_read_test <- function(filepath, nrows = 5) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  if (grepl("\\.gz$", filepath, ignore.case = TRUE)) {
    cmd_txt <- paste("gzip -dc", shQuote(filepath), "| head -n", nrows)
    out <- fread(
      cmd = cmd_txt,
      sep = "\t",
      header = TRUE,
      data.table = FALSE
    )
  } else {
    out <- fread(
      filepath,
      sep = "\t",
      header = TRUE,
      nrows = nrows,
      data.table = FALSE
    )
  }
  
  return(out)
}

###############################################################################
## 2. Download table
###############################################################################
## Notes:
## 1) TCGA-GTEx-TARGET expression TPM matrix
## 2) phenotype / clinical / survival / mutation / CNV files
## 3) MC3_PUBLIC.maf.gz is intentionally not included here for now
##    because the current workflow mainly needs gene-level nonsilent mutation
##    annotation rather than the very large raw MAF file.

download_tbl <- data.frame(
  name = c(
    "TcgaTargetGtex_rsem_gene_tpm.gz",
    "TcgaTargetGTEX_phenotype.txt.gz",
    "TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
    "Survival_SupplementalTable_S1_20171025_xena_sp",
    "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",
    "gtex_RSEM_gene_tpm.gz",
    "GTEX_phenotype.gz",
    "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
  ),
  url = c(
    "https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz",
    "https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
    "https://tcga.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp",
    "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz",
    "https://toil.xenahubs.net/download/gtex_RSEM_gene_tpm.gz",
    "https://toil.xenahubs.net/download/GTEX_phenotype.gz",
    "https://tcga.xenahubs.net/download/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
  ),
  stringsAsFactors = FALSE
)

download_tbl$destfile <- file.path(raw_dir, download_tbl$name)

###############################################################################
## 3. Remove incomplete files if present
###############################################################################
msg("Checking incomplete files before download ...")

for (i in seq_len(nrow(download_tbl))) {
  f <- download_tbl$destfile[i]
  if (file.exists(f)) {
    fs <- file.info(f)$size
    min_size_now <- if (grepl("Survival_SupplementalTable_S1_20171025_xena_sp$", basename(f))) 100 else 1000
    if (is.na(fs) || fs < min_size_now) {
      msg("Removing incomplete file: ", basename(f))
      unlink(f)
    }
  }
}

###############################################################################
## 4. Download files
###############################################################################
msg("Start downloading TCGA/GTEx raw files ...")

for (i in seq_len(nrow(download_tbl))) {
  min_size_now <- if (grepl("Survival_SupplementalTable_S1_20171025_xena_sp$", download_tbl$name[i])) {
    100
  } else {
    1000
  }
  
  download_one_file(
    url = download_tbl$url[i],
    destfile = download_tbl$destfile[i],
    overwrite = FALSE,
    min_size_bytes = min_size_now,
    retries = 2
  )
}

msg("All downloads completed.")

###############################################################################
## 5. File summary
###############################################################################
file_summary <- data.frame(
  file = basename(download_tbl$destfile),
  exists = file.exists(download_tbl$destfile),
  size_bytes = file.info(download_tbl$destfile)$size,
  stringsAsFactors = FALSE
)

print(file_summary)

fwrite(
  file_summary,
  file = file.path(raw_dir, "download_file_summary.tsv"),
  sep = "\t"
)

###############################################################################
## 6. Quick read test
###############################################################################
msg("Running quick read tests ...")

expr_test <- quick_read_test(file.path(raw_dir, "TcgaTargetGtex_rsem_gene_tpm.gz"), nrows = 5)
pheno_test <- quick_read_test(file.path(raw_dir, "TcgaTargetGTEX_phenotype.txt.gz"), nrows = 5)
tcga_pheno_test <- quick_read_test(file.path(raw_dir, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"), nrows = 5)
surv_test <- quick_read_test(file.path(raw_dir, "Survival_SupplementalTable_S1_20171025_xena_sp"), nrows = 5)
mut_test <- quick_read_test(file.path(raw_dir, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"), nrows = 5)
gtex_pheno_test <- quick_read_test(file.path(raw_dir, "GTEX_phenotype.gz"), nrows = 5)
cnv_test <- quick_read_test(file.path(raw_dir, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"), nrows = 5)

msg("Expression test dim: ", nrow(expr_test), " x ", ncol(expr_test))
msg("Combined phenotype test dim: ", nrow(pheno_test), " x ", ncol(pheno_test))
msg("TCGA phenotype test dim: ", nrow(tcga_pheno_test), " x ", ncol(tcga_pheno_test))
msg("Survival test dim: ", nrow(surv_test), " x ", ncol(surv_test))
msg("Mutation test dim: ", nrow(mut_test), " x ", ncol(mut_test))
msg("GTEx phenotype test dim: ", nrow(gtex_pheno_test), " x ", ncol(gtex_pheno_test))
msg("CNV test dim: ", nrow(cnv_test), " x ", ncol(cnv_test))

###############################################################################
## 7. Save session info and key file paths
###############################################################################
writeLines(capture.output(sessionInfo()),
           con = file.path(raw_dir, "Step1_sessionInfo.txt"))

key_files <- data.frame(
  object = c(
    "expr_file",
    "pheno_file",
    "tcga_pheno_file",
    "surv_file",
    "mut_gene_file",
    "gtex_expr_file",
    "gtex_pheno_file",
    "cnv_file"
  ),
  path = c(
    file.path(raw_dir, "TcgaTargetGtex_rsem_gene_tpm.gz"),
    file.path(raw_dir, "TcgaTargetGTEX_phenotype.txt.gz"),
    file.path(raw_dir, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
    file.path(raw_dir, "Survival_SupplementalTable_S1_20171025_xena_sp"),
    file.path(raw_dir, "mc3.v0.2.8.PUBLIC.nonsilentGene.xena.gz"),
    file.path(raw_dir, "gtex_RSEM_gene_tpm.gz"),
    file.path(raw_dir, "GTEX_phenotype.gz"),
    file.path(raw_dir, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz")
  ),
  stringsAsFactors = FALSE
)

fwrite(
  key_files,
  file = file.path(raw_dir, "Step1_key_files.tsv"),
  sep = "\t"
)

msg("Step 1 finished successfully.")