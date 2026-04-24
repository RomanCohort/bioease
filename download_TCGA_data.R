# ==============================================================================
# TCGA 数据下载脚本
# 功能：下载TCGA RNA-seq counts 数据和临床生存数据
# 输出：counts_matrix.rds, metadata.rds, clinical_data.rds
# ==============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || args[1] != "--cancer") {
  stop("Usage: Rscript download_TCGA_data.R --cancer <CANCER_TYPE>")
}
cancer <- args[2]

# 创建结果目录
dir.create("results", recursive = TRUE, showWarnings = FALSE)

cat("Starting download for project:", cancer, "\n")

# ==============================================================================
# 1. 下载RNA-seq数据
# ==============================================================================

tryCatch({
  # 查询TCGA RNA-seq数据
  cat("Querying data...\n")
  query <- GDCquery(project = cancer,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  cat("Query completed. Found", nrow(getResults(query)), "files.\n")
  
  # 下载数据
  cat("Downloading data...\n")
  GDCdownload(query)
  cat("Download completed.\n")
  
  # 准备数据
  cat("Preparing data...\n")
  brca_data <- GDCprepare(query)
  cat("Data preparation completed.\n")
  
  # 提取counts矩阵
  counts <- assay(brca_data, "unstranded")
  gene_info <- rowData(brca_data)
  
  # 样本信息
  col_data <- colData(brca_data)
  
  # 保存counts矩阵
  saveRDS(list(counts = counts, gene_name = gene_info$gene_name), "results/counts_matrix.rds")
  cat("Counts matrix saved.\n")
  
  # 保存元数据
  metadata <- data.frame(
    sampleID = colnames(counts),
    sample_file = paste0(colnames(counts), ".tsv"),
    sample_path = NA,  # 不适用
    condition = ifelse(grepl("-01A-", colnames(counts)), "Tumor", "Normal")
  )
  saveRDS(metadata, "results/metadata.rds")
  cat("Metadata saved.\n")
  
  cat("RNA-seq data downloaded and saved.\n")
}, error = function(e) {
  cat("Error in RNA-seq download:", e$message, "\n")
  stop(e$message)
})

# ==============================================================================
# 2. 下载临床数据
# ==============================================================================

tryCatch({
  cat("Querying clinical data...\n")
  # 查询临床数据
  clinical_query <- GDCquery(project = cancer,
                             data.category = "Clinical",
                             data.type = "Clinical Supplement",
                             data.format = "BCR Biotab")
  cat("Clinical query completed.\n")
  
  # 下载临床数据
  cat("Downloading clinical data...\n")
  GDCdownload(clinical_query)
  cat("Clinical download completed.\n")
  
  # 准备临床数据
  cat("Preparing clinical data...\n")
  clinical_data <- GDCprepare(clinical_query)
  cat("Clinical data preparation completed.\n")
  
  # 提取生存信息
  cancer_lower <- tolower(gsub("TCGA-", "", cancer))
  patient_file <- list.files(path = "GDCdata", pattern = paste0("clinical_patient_", cancer_lower, ".txt"), recursive = TRUE, full.names = TRUE)
  if (length(patient_file) > 0) {
    patient_data <- read.delim(patient_file[1], stringsAsFactors = FALSE)
    survival_data <- patient_data %>%
      select(bcr_patient_barcode,
             days_to_death = days_to_death,
             days_to_last_followup = days_to_last_followup,
             vital_status = vital_status) %>%
      mutate(OS.time = ifelse(!is.na(days_to_death), days_to_death,
                             ifelse(!is.na(days_to_last_followup), days_to_last_followup, NA)),
             OS = ifelse(vital_status == "Dead", 1, 0)) %>%
      filter(!is.na(OS.time) & OS.time > 0) %>%
      select(sample = bcr_patient_barcode, OS.time, OS)
    
    # 保存临床数据
    saveRDS(survival_data, "results/clinical_survival.rds")
    cat("Clinical survival data saved.\n")
  } else {
    stop("Failed to find clinical data file")
  }
  
  cat("Clinical data downloaded and saved.\n")
  cat("Data download completed. You can now run analysis_script.R and upregulated_TAA_survival_analysis.R\n")
}, error = function(e) {
  cat("Error in clinical download:", e$message, "\n")
  stop(e$message)
})