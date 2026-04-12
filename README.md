# Nucleolar Stress Score (NuS)

This repository contains the code used in the study:

**"A nucleolar stress gene signature for quantitative scoring across multi-omics contexts"**

---

## Overview

NuS (Nucleolar Stress Score) is a gene expression–based metric designed to quantify nucleolar stress across multiple data modalities, including:

- Bulk transcriptomics
- Single-cell RNA-seq
- Spatial transcriptomics
- Proteomics

---

## Repository Structure

- `Genesets/`  
  Gene sets used for NuS calculation (upregulated and downregulated signatures)

- `TCGA/`  
  TCGA pan-cancer analysis scripts

- `Spatial_transcriptome/`  
  Spatial transcriptomics analysis

- `Dataset_analyse/`  
  Dataset-level processing and analysis

- `Data_validation/`  
  Validation datasets and analyses

- `scRNA-seq.R`  
  Single-cell analysis pipeline

- `Proteomics.R`  
  Proteomics-based NuS analysis

- `CAMP.R`  
  Drug screening and perturbation analysis (CMap)

---

## Requirements

R (>= 4.0 recommended)

Key packages:
- GSVA
- Seurat
- limma
- tidyverse

---

## Usage

NuS is calculated as:
NuS = Score_up − Score_down
where enrichment scores are computed using GSVA or ssGSEA.

Users may need to modify file paths and adjust parameters depending on their datasets.

---

## Notes

- This repository contains analysis scripts rather than a packaged pipeline.
- Large datasets (e.g., TCGA) are not included and should be downloaded separately.
- Scripts may require adaptation to local environments.

---

## Code Availability

All code used in this study is publicly available at:
https://github.com/JianxiongChen-cloud/nucleolar-stress-score
