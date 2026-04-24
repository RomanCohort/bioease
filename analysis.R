# ==============================================================================
# TCGA BRCA 差异表达基因分析脚本
# 功能：从转录组数据中筛选差异表达基因，绘制火山图并标注Top 10显著基因
# 输入：TCGA_BRCA_Cancer 和 TCGA_BRCA_Normal 文件夹中的 .tsv 文件
# 输出：volcano_plot.png (火山图), upregulated_TAA_list.txt (上调基因列表)
# 依赖包：DESeq2, ggplot2, ggrepel
# ==============================================================================

# ==============================================================================
# 包管理部分：安装和加载所需的R包
# ==============================================================================

# 检查并安装 BiocManager (用于安装 Bioconductor 包)
# BiocManager 是 Bioconductor 项目的包管理器，类似于 CRAN 的 install.packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 检查并安装 DESeq2 (差异表达分析的核心包)
# DESeq2 用于 RNA-seq 数据的标准化和差异表达分析
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

# 检查并安装 ggplot2 (绘图包)
# ggplot2 用于创建高质量的统计图形，包括火山图
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

# 检查并安装 ggrepel (用于标注基因名，避免重叠)
# ggrepel 在ggplot2图表中智能放置文本标签，避免重叠
if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel")

# 用于 ID->基因名 映射的包（按需安装）
# 这些包用于将基因ID转换为基因符号
if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

# 加载包到环境
# 将所有必需的包加载到R会话中
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(httr)
library(jsonlite)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)

# ==============================================================================
# 1. 数据读取与预处理
# ==============================================================================

# ------------------------- 配置与命令行参数 -------------------------
# 脚本支持灵活的命令行参数配置
# 支持通过命令行传入 base_dir（例如：Rscript script_analysis.R "D:/颜/生信文件夹"）
# 如果不传参则使用脚本中默认路径。
args <- commandArgs(trailingOnly = TRUE)

# 初始化候选文件变量
candidate_file <- NULL

# 解析命令行参数
# 支持两种参数格式：
# 1. --base_dir=PATH 和 --candidate_file=PATH (标准格式)
# 2. 位置参数 (向后兼容)
if (length(args) >= 1) {
  # 先解析形如 --key=val 的参数
  for (a in args) {
    if (grepl('^--candidate_file=', a)) {
      candidate_file <- sub('^--candidate_file=', '', a)
    } else if (grepl('^--base_dir=', a)) {
      base_dir <- sub('^--base_dir=', '', a)
    }
  }
  # 兼容旧式位置参数：若第一个参数不是以 -- 开头，则视为 base_dir；若第二参数存在且不是 -- 开头，则视为 candidate_file
  if (length(args) >= 1 && !grepl('^--', args[1]) && !nzchar(base_dir)) {
    base_dir <- args[1]
  }
  if (length(args) >= 2 && !grepl('^--', args[2]) && is.null(candidate_file)) {
    candidate_file <- args[2]
  }
}

# 若未通过参数提供 base_dir，则使用脚本内默认路径
if (!exists("base_dir") || !nzchar(base_dir)) base_dir <- "D:/颜/生信文件夹"

# 如果提供了 candidate_file，则规范化路径
# candidate_file 包含基因ID到基因名的映射表
if (!is.null(candidate_file) && nzchar(candidate_file)) {
  candidate_file <- tryCatch({ normalizePath(candidate_file, winslash = "/", mustWork = FALSE) }, error = function(e) candidate_file)
  cat("Using candidate_file:", candidate_file, "\n")
}

# 结果输出目录（将所有产物写入此目录，便于管理）
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
if (length(args) >= 1) {
  # 先解析形如 --key=val 的参数
  for (a in args) {
    if (grepl('^--candidate_file=', a)) {
      candidate_file <- sub('^--candidate_file=', '', a)
    } else if (grepl('^--base_dir=', a)) {
      base_dir <- sub('^--base_dir=', '', a)
    }
  }
  # 兼容旧式位置参数：若第一个参数不是以 -- 开头，则视为 base_dir；若第二参数存在且不是 -- 开头，则视为 candidate_file
  if (length(args) >= 1 && !grepl('^--', args[1]) && !nzchar(base_dir)) {
    base_dir <- args[1]
  }
  if (length(args) >= 2 && !grepl('^--', args[2]) && is.null(candidate_file)) {
    candidate_file <- args[2]
  }
}

# 若未通过参数提供 base_dir，则使用脚本内默认路径
if (!exists("base_dir") || !nzchar(base_dir)) base_dir <- "D:/颜/生信文件夹"

# 如果提供了 candidate_file，则规范化路径
if (!is.null(candidate_file) && nzchar(candidate_file)) {
  candidate_file <- tryCatch({ normalizePath(candidate_file, winslash = "/", mustWork = FALSE) }, error = function(e) candidate_file)
  cat("Using candidate_file:", candidate_file, "\n")
}

# 结果输出目录（将所有产物写入此目录，便于管理）
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 1.1 文件路径与元数据构建
# 本节负责自动发现和配置数据文件路径，以及构建样本元数据

# 定义癌症和正常组织数据文件夹路径（在 base_dir 下）
# 支持常见命名变体并递归搜索子目录中的 .tsv 文件，增强容错性
find_existing_dir <- function(base, candidates) {
  # 策略1: 先尝试 base 下的直接子目录
  # 例如：base/TCGA_BRCA_Cancer
  for (d in candidates) {
    path <- file.path(base, d)
    if (dir.exists(path)) {
      return(path)
    }
  }

  # 策略2: 如果未找到，递归搜索 base 下的所有子目录，匹配目录名（basename）等于候选名
  # 这处理数据被放在某个子文件夹（例如 Candidate_Test_Data）的情况
  if (dir.exists(base)) {
    all_dirs <- list.dirs(base, recursive = TRUE, full.names = TRUE)
    # 优先严格匹配 basename
    for (d in candidates) {
      matches <- all_dirs[basename(all_dirs) == d]
      if (length(matches) > 0) return(matches[1])
    }
    # 然后尝试不区分大小写的匹配
    for (d in candidates) {
      matches <- all_dirs[tolower(basename(all_dirs)) == tolower(d)]
      if (length(matches) > 0) return(matches[1])
    }
    # 最后尝试包含关系（目录名或路径中包含候选字符串）
    for (d in candidates) {
      idx <- grep(tolower(d), tolower(all_dirs), fixed = TRUE)
      if (length(idx) > 0) return(all_dirs[idx[1]])
    }
  }

  # 若仍然未找到，返回第一个候选路径（原有行为），后续会抛出友好错误
  return(file.path(base, candidates[1]))
}

# 候选目录名称（可根据需要扩展）
# 支持多种可能的目录命名方式，提高脚本的鲁棒性
cancer_dir_candidates <- c("TCGA_BRCA_Cancer", "TCGA_BRCA_Tumor")
normal_dir_candidates <- c("TCGA_BRCA_Normal", "TCGA_BRCA_Norm", "TCGA_BRCA_Control")

# 使用智能搜索函数找到实际的数据目录
cancer_dir <- find_existing_dir(base_dir, cancer_dir_candidates)
normal_dir <- find_existing_dir(base_dir, normal_dir_candidates)

cat("Using cancer_dir:", cancer_dir, "\n")
cat("Using normal_dir:", normal_dir, "\n")

# 使用 list.files() 获取所有相关计数文件的完整路径（递归搜索子目录）
# 支持常见扩展：.tsv, .tsv.gz, .txt；且忽略大小写以提高兼容性
# 这样可以处理压缩文件和不同命名方式
file_pattern <- "\\.(tsv|tsv\\.gz|txt)$"
cancer_files <- if (dir.exists(cancer_dir)) list.files(cancer_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE) else character(0)
normal_files <- if (dir.exists(normal_dir)) list.files(normal_dir, pattern = file_pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE) else character(0)

cat(sprintf("Found %d cancer files and %d normal files matching pattern '%s'\n", length(cancer_files), length(normal_files), file_pattern))

# 输出发现的文件样本，便于用户验证
if (length(cancer_files) > 0) {
  cat("Sample cancer files:\n", paste(head(cancer_files, 10), collapse = "\n"), "\n")
}
if (length(normal_files) > 0) {
  cat("Sample normal files:\n", paste(head(normal_files, 10), collapse = "\n"), "\n")
}

# 写出诊断日志到 results 目录，包含找到的文件和每组各前两个文件的前 20 行预览
# 这有助于调试和验证数据文件的内容和格式
diag_file <- file.path(results_dir, "file_discovery_log.txt")
diag_lines <- character()
diag_lines <- c(diag_lines, sprintf("Found %d cancer files and %d normal files (pattern='%s')", length(cancer_files), length(normal_files), file_pattern))
if (length(cancer_files) > 0) {
  diag_lines <- c(diag_lines, "\nSample cancer files:")
  diag_lines <- c(diag_lines, head(cancer_files, 20))
  # 预览前两个文件的前 20 行，帮助用户了解数据格式
  for (f in head(cancer_files, 2)) {
    diag_lines <- c(diag_lines, sprintf("\n--- Preview: %s ---", f))
    preview <- tryCatch({
      if (grepl("\\.gz$", f, ignore.case = TRUE)) readLines(gzfile(f), n = 20) else readLines(f, n = 20)
    }, error = function(e) { paste0("<error reading file: ", e$message, ">") })
    diag_lines <- c(diag_lines, preview)
  }
}
if (length(normal_files) > 0) {
  diag_lines <- c(diag_lines, "\nSample normal files:")
  diag_lines <- c(diag_lines, head(normal_files, 20))
  for (f in head(normal_files, 2)) {
    diag_lines <- c(diag_lines, sprintf("\n--- Preview: %s ---", f))
    preview <- tryCatch({
      if (grepl("\\.gz$", f, ignore.case = TRUE)) readLines(gzfile(f), n = 20) else readLines(f, n = 20)
    }, error = function(e) { paste0("<error reading file: ", e$message, ">") })
    diag_lines <- c(diag_lines, preview)
  }
}
tryCatch({
  writeLines(diag_lines, con = diag_file)
  cat("Wrote discovery diagnostics to:", diag_file, "\n")
}, error = function(e) {
  cat("Failed to write discovery diagnostics:", e$message, "\n")
})

# 检查是否找到任何数据文件，如果没有则提供详细的错误信息和诊断建议
if (length(cancer_files) == 0 && length(normal_files) == 0) {
  stop(paste0("未在 cancer_dir='", cancer_dir, "' 或 normal_dir='", normal_dir, "' 中找到任何匹配的文件。",
              "可能原因：目录路径不正确，文件扩展名不是 .tsv（或为 .tsv.gz/.txt），",
              "或文件位于更深的子目录中/有权限问题。\n",
              "建议诊断命令（在 PowerShell 中运行）：\n",
              "  Get-ChildItem -Path '", cancer_dir, "' -Recurse -Include *.tsv,*.tsv.gz,*.txt -File | Select-Object -First 20 FullName\n",
              "  Get-ChildItem -Path '", normal_dir, "' -Recurse -Include *.tsv,*.tsv.gz,*.txt -File | Select-Object -First 20 FullName\n",
              "或在 R 中检查：\n",
              "  list.files('", cancer_dir, "', pattern='\\\\.(tsv|tsv\\\\.gz|txt)$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)\n",
              "  list.files('", normal_dir, "', pattern='\\\\.(tsv|tsv\\\\.gz|txt)$', recursive=TRUE, full.names=TRUE, ignore.case=TRUE)\n"
  ))
}

# 合并所有文件路径 (癌症文件在前，正常文件在后)
# 这样确保了数据矩阵的列顺序与分组一致
all_files <- c(cancer_files, normal_files)

# 从文件路径生成样本ID (移除路径和扩展名)
# 样本ID将用作数据矩阵的列名
sample_files <- basename(all_files)  # 获取文件名
sample_ids <- sub("\\.tsv$", "", sample_files)  # 移除 .tsv 扩展名
sample_ids <- sub("\\.txt$", "", sample_ids)  # 也移除 .txt 扩展名（以防万一）

# 生成条件标签 (Tumor/Normal) —— 基于实际文件数量动态生成，避免固定长度导致的 row.names 错误
# 动态生成确保与实际文件数量匹配
n_cancer <- length(cancer_files)
n_normal <- length(normal_files)
conditions <- c(rep("Tumor", n_cancer), rep("Normal", n_normal))

# 检查 sample_ids 与 conditions 长度是否匹配
# 这是一个重要的完整性检查，确保数据一致性
if (length(conditions) != length(sample_ids)) {
  stop(sprintf("样本数不匹配：从文件检测到 %d 个样本，但生成的 conditions 长度为 %d，请检查文件夹路径和文件数量。",
               length(sample_ids), length(conditions)))
}

# 处理样本名重复的情况（若存在重复则使其唯一并给出提示）
# 重复的样本名可能导致DESeq2分析错误
if (any(duplicated(sample_ids))) {
  warning("检测到重复的样本名（文件名去掉扩展名后）。将对重复名使用 make.unique() 以确保唯一性。")
  sample_ids <- make.unique(sample_ids)
}

# 创建元数据数据框 (colData)
# colData 是DESeq2必需的样本信息数据框，包含实验设计信息
colData <- data.frame(
  sampleID = sample_ids,                    # 样本唯一标识符
  sample_file = sample_files,               # 原始文件名
  sample_path = normalizePath(all_files, winslash = "/", mustWork = FALSE),  # 文件完整路径
  condition = factor(conditions, levels = c("Normal", "Tumor")),  # 实验条件，Normal设为对照组
  stringsAsFactors = FALSE
)

# 设置行名为样本ID，确保与 countData 列名匹配
# 这是DESeq2的要求：colData的行名必须与countData的列名完全匹配
rownames(colData) <- colData$sampleID

# 在分析开始前保存元数据（TSV 与 RDS），方便复现与检查
# 保存为两种格式：TSV用于人工查看，RDS用于R程序读取
metadata_tsv <- file.path(results_dir, "metadata.tsv")
metadata_rds <- file.path(results_dir, "metadata.rds")
write.table(colData, file = metadata_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(colData, metadata_rds)
cat("Saved metadata to:\n", metadata_tsv, "\n", metadata_rds, "\n")

# 1.2 数据读取与预处理
# 本节负责读取所有RNA-seq计数文件，合并成统一的计数矩阵

# 定义读取单个RNA-seq计数文件的函数
# 这个函数需要处理各种文件格式和潜在的数据质量问题
read_count_file <- function(file_path) {
  # 函数功能：读取单个TCGA RNA-seq计数文件，返回基因表达计数向量
  # 支持 gz 压缩文件并智能识别计数列
  # 处理步骤：
  # 1) 读取文件（支持 .gz）；2) 过滤明显的注释行；3) 解析为表格；
  # 4) 自动检测 gene id 列与 count 列（按列名匹配或数值检测）

  # 步骤1: 读取原始文本（支持 gz 压缩文件）
  raw_lines <- NULL
  if (grepl("\\.gz$", file_path, ignore.case = TRUE)) {
    con <- gzfile(file_path, open = "rt")
    raw_lines <- readLines(con)
    close(con)
  } else {
    raw_lines <- readLines(file_path)
  }

  # 步骤2: 过滤掉明显的注释/元信息行
  # TCGA文件可能包含以 "N_" 开头的统计行，这些不是基因表达数据
  # 使用正则表达式过滤，但避免删除以 N 为基因名的正常基因
  filtered <- raw_lines[!grepl("^(N_|N[[:space:]]|N$)", raw_lines)]

  # 如果过滤后为空，则回退使用原始内容（以防过度过滤误删数据）
  if (length(filtered) == 0) filtered <- raw_lines

  # 步骤3: 解析过滤后的文本为数据框
  # 先尝试标准解析方法
  df <- NULL
  parse_error <- NULL
  try({
    df <- read.delim(text = filtered, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  }, silent = TRUE)

  if (is.null(df)) {
    # 记录解析错误并尝试更宽松的解析方法
    parse_error <- "read.delim failed"
    warning(sprintf("标准解析失败（%s）；尝试更宽松的解析方式。", file_path))

    # 备选方案1: 使用 data.table::fread (如果可用)
    if (requireNamespace("data.table", quietly = TRUE)) {
      try({
        df <- data.table::fread(text = paste(filtered, collapse = "\n"), sep = "\t", header = TRUE, fill = TRUE, data.table = FALSE, showProgress = FALSE)
      }, silent = TRUE)
    }

    # 备选方案2: 使用 base::read.table 进行更宽松的解析
    if (is.null(df)) {
      try({
        df <- read.table(text = paste(filtered, collapse = "\n"), header = TRUE, sep = "\t", fill = TRUE, quote = "\"", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
      }, silent = TRUE)
    }
  }

  # 如果仍未成功解析为多列表格，生成诊断并返回 NULL（该文件将被跳过）
  if (is.null(df) || ncol(df) < 2) {
    # 诊断信息：记录出现问题的文件路径、前 200 行样例、每行的字段数分布
    diag_path <- file.path(results_dir, "problematic_files.txt")
    max_preview <- min(length(filtered), 200)
    lines_preview <- filtered[1:max_preview]
    # 计算每行字段数
    field_counts <- sapply(filtered, function(l) length(strsplit(l, "\t", fixed = TRUE)[[1]]))
    txt <- c(sprintf("--- Problem parsing file: %s ---", file_path),
             sprintf("parse methods attempted: %s", paste(c("read.delim","data.table::fread","read.table"), collapse = ", ")),
             sprintf("max lines previewed: %d", max_preview),
             "\n--- Lines (first preview) ---",
             lines_preview,
             "\n--- Field counts summary ---",
             sprintf("max fields in a line: %d", max(field_counts, na.rm = TRUE)),
             paste0("field counts (first 100 lines): ", paste(head(field_counts, 100), collapse = ", ")),
             "\n")
    # 将诊断信息追加到诊断文件
    tryCatch({
      cat(paste0(txt, collapse = "\n"), file = diag_path, append = TRUE)
      warning(sprintf("无法解析文件 '%s'（已记录到 %s）。尝试使用手工解析策略。", file_path, diag_path))
    }, error = function(e) {
      warning(sprintf("无法写入诊断文件 '%s'：%s", diag_path, e$message))
    })

    # 手工解析回退：逐行按制表符分割，构建矩阵并尝试找到基因列与计数列
    pieces <- strsplit(filtered, "\t", fixed = TRUE)
    nfields <- sapply(pieces, length)
    maxf <- max(nfields, na.rm = TRUE)
    # 构建字符矩阵并填充 NA
    mat <- matrix(NA_character_, nrow = length(pieces), ncol = maxf)
    for (i in seq_along(pieces)) mat[i, seq_len(length(pieces[[i]]))] <- pieces[[i]]
    df_manual <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
    # 尝试把第一行为 header 的可能性：如果第一行中包含任何非空且非数字元素，则作为 header
    first_row <- df_manual[1, ]
  # 判断第一行是否像 header：包含非空且非数字的元素
  is_header_like <- any(!grepl('^[[:space:]]*$', first_row) & is.na(suppressWarnings(as.numeric(first_row))))
    if (is_header_like) {
      colnames(df_manual) <- as.character(unlist(df_manual[1, ]))
      df_manual <- df_manual[-1, , drop = FALSE]
      # 如果列名有空则填充默认名
      empty_names <- which(colnames(df_manual) == "" | is.na(colnames(df_manual)))
      if (length(empty_names) > 0) colnames(df_manual)[empty_names] <- paste0("V", empty_names)
    }

    # 现在尝试识别计数列：选择那些大多数值可转为数值的列
    numeric_prop <- sapply(df_manual, function(col) {
      vals <- col[!is.na(col) & nzchar(col)]
      if (length(vals) == 0) return(0)
      mean(!is.na(suppressWarnings(as.numeric(vals))))
    })
    # 选择 numeric_prop > 0.6 的列作为候选计数列
    candidate_idx <- which(numeric_prop > 0.6)
    if (length(candidate_idx) == 0) {
      # 退化策略：选第一列之外的列中数值比例最高的列
      candidate_idx <- which.max(numeric_prop)
    }
    # 选择第一个 candidate_idx
    count_col_idx <- candidate_idx[1]
    # 基因列优先使用第一列；如果第一列也是高度数值，则尝试使用第二列作为基因列
    gene_col_idx <- 1
    if (!is_header_like) {
      # 如果第一列看起来像数值（多数为数字），则尝试把非数字最多的列作为 gene
      gene_numeric_prop <- numeric_prop[gene_col_idx]
      if (gene_numeric_prop > 0.5) {
        # 选择数值比例最低的前3列中的第一个作为 gene 列
        ord <- order(numeric_prop)
        gene_col_idx <- ord[1]
      }
    }

    # 构建最终 data.frame，尝试命名列
    gene_col_name <- if (!is.null(colnames(df_manual))) colnames(df_manual)[gene_col_idx] else paste0("V", gene_col_idx)
    count_col_name <- if (!is.null(colnames(df_manual))) colnames(df_manual)[count_col_idx] else paste0("V", count_col_idx)
    # 提取并清理
    gene_vals <- as.character(df_manual[[gene_col_idx]])
    count_vals <- suppressWarnings(as.integer(as.numeric(df_manual[[count_col_idx]])))
    # 如果 count_vals 都为 NA，则放弃手工解析
    if (all(is.na(count_vals))) {
      warning(sprintf("手工解析失败：在文件 '%s' 未能找到有效的计数列，已跳过。", file_path))
      return(NULL)
    }
    counts <- data.frame(count = count_vals, stringsAsFactors = FALSE)
    # 去除行名为 NA 或 空字符串的基因
    valid_idx <- which(!is.na(gene_vals) & nzchar(gene_vals))
    if (length(valid_idx) == 0) {
      warning(sprintf("手工解析失败：在文件 '%s' 未能找到有效的基因ID列，已跳过。", file_path))
      return(NULL)
    }
    counts <- counts[valid_idx, , drop = FALSE]
    rownames(counts) <- gene_vals[valid_idx]
    # 设置样本列名
    samp_name <- basename(file_path)
    samp_name <- sub("\\.gz$", "", samp_name, ignore.case = TRUE)
    samp_name <- sub("\\.rna_seq.*$", "", samp_name, ignore.case = TRUE)
    samp_name <- sub("\\.tsv$", "", samp_name, ignore.case = TRUE)
    colnames(counts) <- samp_name
    return(counts)
  }
  

  if (is.null(df) || ncol(df) < 2) {
    stop(sprintf("文件 '%s' 列数少于2，无法识别 gene id 与 count 列。", file_path))
  }

  # 自动识别 gene id 列：优先选择名为 gene/ensembl/gene_id/gene_name 的列，若无则取第一列
  gene_col_candidates <- c("gene", "gene_id", "gene_name", "ensembl_gene_id", "Name")
  gene_col <- NULL
  for (nm in gene_col_candidates) {
    if (nm %in% colnames(df)) { gene_col <- nm; break }
  }
  if (is.null(gene_col)) gene_col <- colnames(df)[1]

  # 自动识别计数列：优先按常见名称匹配，否则选中第一个数值列（排除基因列）
  count_col_candidates <- c("count", "counts", "read_count", "expression", "raw_count")
  count_col <- NULL
  for (nm in count_col_candidates) {
    if (nm %in% colnames(df)) { count_col <- nm; break }
  }
  if (is.null(count_col)) {
    # 在剩余列中寻找数值型列
    other_cols <- setdiff(colnames(df), gene_col)
    numeric_cols <- sapply(df[other_cols], function(x) all(suppressWarnings(!is.na(as.numeric(x[!is.na(x)])))) )
    if (any(numeric_cols)) {
      count_col <- other_cols[which(numeric_cols)[1]]
    } else {
      # 退化为第二列
      count_col <- colnames(df)[2]
    }
  }

  # 构建 counts 数据框：基因ID 为行名，计数为单列
  counts <- data.frame(count = as.integer(as.numeric(df[[count_col]])), stringsAsFactors = FALSE)
  rownames(counts) <- as.character(df[[gene_col]])

  # 如果转换产生 NA（例如列中包含非数值），尽量保留原始并尝试强制为整数
  if (all(is.na(counts$count))) {
    stop(sprintf("在文件 '%s' 中未能解析出有效的计数列（尝试的列：%s）。", file_path, count_col))
  }

  # 设置列名为样本ID：若文件名包含 ".rna_seq" 这样的后缀，去掉以得到更短的样本名
  samp_name <- basename(file_path)
  samp_name <- sub("\\.gz$", "", samp_name, ignore.case = TRUE)
  samp_name <- sub("\\.rna_seq.*$", "", samp_name, ignore.case = TRUE)
  samp_name <- sub("\\.tsv$", "", samp_name, ignore.case = TRUE)
  colnames(counts) <- samp_name

  return(counts)
}

# 使用for循环逐一读取所有文件
# 初始化空列表存储结果
count_list <- list()
for (file_path in all_files) {
  # 调用读取函数，获取计数矩阵；对单个文件读取添加 tryCatch，以便记录具体出错文件并继续处理其他文件
  counts <- tryCatch({
    read_count_file(file_path)
  }, error = function(e) {
    warning(sprintf("Error while reading file '%s': %s", file_path, e$message))
    # 将错误信息追加写入诊断文件，便于排查
    err_log <- file.path(results_dir, "read_errors.log")
    tryCatch({
      cat(sprintf("%s\tError reading '%s' : %s\n", Sys.time(), file_path, e$message), file = err_log, append = TRUE)
    }, error = function(e2) {})
    return(NULL)
  })
  # 将结果添加到列表末尾（NULL 会在后续被过滤掉）
  count_list[[length(count_list) + 1]] <- counts
}

# 移除可能为 NULL 的元素（如果某些文件读取失败）并在没有读取到任何样本时给出友好错误
count_list <- Filter(Negate(is.null), count_list)
if (length(count_list) == 0) {
  stop(sprintf("未读取到任何计数文件，请检查目录和文件：\n  cancer_dir='%s'\n  normal_dir='%s'\n检测到的文件总数: %d (cancer=%d, normal=%d)。",
               cancer_dir, normal_dir, length(all_files), n_cancer, n_normal))
}

# 1.3 数据合并与格式化
# 使用 do.call(cbind, ...) 将所有单样本矩阵按列合并为统一计数矩阵
# 结果：行名为基因ID，列名为样本ID，值为整数计数值
countData <- do.call(cbind, count_list)

# 将数据框转换为矩阵，并确保存储模式为整数 (DESeq2 要求)
countData <- as.matrix(countData)

# ---------------------------- NA 检查与修复 ----------------------------
# DESeq2 要求计数矩阵中不能含有 NA。下面的逻辑：
# 1) 如果检测到 NA，记录诊断（前 100 个位置写到 results 目录）；
# 2) 将 NA 替换为 0（常见做法，但请确认对你数据合理）；
# 3) 若你希望改为删除含 NA 的基因/样本，可手动调整此处逻辑。
na_count <- sum(is.na(countData))
if (na_count > 0) {
  warning(sprintf("检测到 %d 个 NA 值在 countData 中。将写出诊断并把 NA 替换为 0。", na_count))
  na_pos <- which(is.na(countData), arr.ind = TRUE)
  # 构建诊断表（行名、列名、位置索引），只保存前 100 条以避免文件过大
  na_diag <- data.frame(
    gene_id = rownames(countData)[na_pos[, 1]],
    sample_id = colnames(countData)[na_pos[, 2]],
    row = na_pos[, 1],
    col = na_pos[, 2],
    stringsAsFactors = FALSE
  )
  na_diag_file <- file.path(results_dir, "count_NA_positions.tsv")
  tryCatch({
    write.table(head(na_diag, 100), file = na_diag_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("Wrote NA diagnostics to: ", na_diag_file)
  }, error = function(e) {
    warning(sprintf("无法写出 NA 诊断文件：%s", e$message))
  })

  # 默认策略：将 NA 当作 0（通常表示未检测到的计数）
  countData[is.na(countData)] <- 0L
  message(sprintf("已将 %d 个 NA 值替换为 0。", na_count))
}

# ----------------------------------------------------------------------

# 确保为整数（DESeq2 期望原始整数计数）
countData <- round(countData)
storage.mode(countData) <- "integer"

# ==============================================================================
# 2. 数据预处理
# ==============================================================================

# 2.1 colData 与 countData 比对
# 确保 countData 的列名 (样本ID) 与 colData 的行名一致且顺序相同
# 如果不匹配，尝试智能对齐；对齐失败则给出清晰诊断并停止
simple_name <- function(x) {
  x2 <- sub("\\.gz$", "", x, ignore.case = TRUE)
  x2 <- sub("\\.rna_seq.*$", "", x2, ignore.case = TRUE)
  x2 <- sub("\\.tsv$", "", x2, ignore.case = TRUE)
  x2 <- basename(x2)
  return(x2)
}

if (!all(colnames(countData) == rownames(colData))) {
  # 报告差异
  missing_in_colData <- setdiff(colnames(countData), rownames(colData))
  missing_in_countData <- setdiff(rownames(colData), colnames(countData))
  warning(sprintf("样本名不匹配：%d 列在 countData 中不存在于 colData，%d 行在 colData 中不存在于 countData。",
                  length(missing_in_colData), length(missing_in_countData)))

  # 尝试按简化名字对齐（去掉扩展与 .rna_seq 后缀）
  col_simple <- simple_name(colnames(countData))
  row_simple <- simple_name(rownames(colData))

  if (all(col_simple %in% row_simple)) {
    # 将 colData 重排以匹配 countData 的列顺序
    ord <- match(col_simple, row_simple)
    colData <- colData[ord, , drop = FALSE]
    # 使 rownames(colData) 与 colnames(countData) 完全一致（保留原始 sampleID 在列中）
    rownames(colData) <- colnames(countData)
    message("已根据简化样本名自动重排 colData 与 countData。")
  } else if (all(row_simple %in% col_simple)) {
    # 将 countData 列重排以匹配 colData
    ord <- match(row_simple, col_simple)
    countData <- countData[, ord, drop = FALSE]
    colnames(countData) <- rownames(colData)
    message("已根据简化样本名自动重排 countData 与 colData。")
  } else {
    # 尝试使用 colData$sample_file（去掉扩展）匹配 countData 列名
    if ("sample_file" %in% colnames(colData)) {
      sample_file_simple <- simple_name(colData$sample_file)
      if (all(colnames(countData) %in% sample_file_simple)) {
        ord <- match(colnames(countData), sample_file_simple)
        colData <- colData[ord, , drop = FALSE]
        rownames(colData) <- colnames(countData)
        message("已根据 colData$sample_file 自动重排 colData 与 countData。")
      } else if (all(sample_file_simple %in% colnames(countData))) {
        ord <- match(sample_file_simple, colnames(countData))
        countData <- countData[, ord, drop = FALSE]
        colnames(countData) <- rownames(colData)
        message("已根据 colData$sample_file 自动重排 countData 与 colData。")
      } else {
        stop(sprintf("样本名仍不匹配。countData 中有 %d 个样本未找到于 colData，colData 中有 %d 个样本未找到于 countData。\n缺失示例 (countData->colData): %s\n缺失示例 (colData->countData): %s\n请检查样本命名或传入正确的 base_dir。",
                     length(missing_in_colData), length(missing_in_countData),
                     paste(head(missing_in_colData, 10), collapse = ", "),
                     paste(head(missing_in_countData, 10), collapse = ", ")))
      }
    } else {
      stop(sprintf("样本名不匹配，且 colData 中无 sample_file 字段以供匹配。请检查: %d vs %d。",
                   length(missing_in_colData), length(missing_in_countData)))
    }
  }
}

# 2.2 低计数基因过滤
# 过滤逻辑：至少在 1 个样本中计数 > 1 的基因
# rowSums(countData > 1) 计算每个基因在多少样本中计数 > 1
# >= 1 确保基因在至少一个样本中有表达
keep <- rowSums(countData > 1) >= 1
countData <- countData[keep, ]

# ==============================================================================
# 注：下面我们在写出计数矩阵与后续结果时，优先输出基因 symbol（gene_name），并保留原始 gene_id 列
# ==============================================================================

# ID -> SYMBOL 映射函数，能处理 Ensembl 带版本号情形（如 ENSG000001234.5）
map_ids_to_symbol <- function(ids) {
  # 使用 biomaRt 优先进行 ID -> SYMBOL 的映射（支持 Ensembl 带版本号、Entrez）
  # 逻辑：去版本号 -> 保留看起来像符号的 id -> biomaRt 批量查询 ENSG/Entrez -> MyGene 回退 -> 写出未映射列表

  ids_orig <- as.character(ids)
  ids2 <- sub("\\.\\d+$", "", ids_orig)

  sym <- rep(NA_character_, length(ids2))
  names(sym) <- ids_orig

  # 若看起来像 gene symbol（字母开头且非 ENSG），先保留
  likely_symbol <- grepl("^[A-Za-z]", ids_orig) & !grepl("^ENS", ids2)
  sym[likely_symbol] <- ids_orig[likely_symbol]

  # 尝试用 biomaRt 批量查询（优先 Ensembl id -> hgnc_symbol）
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    try({
      ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

      # 1) ENSG 类型的 id
      ens_idx <- which(is.na(sym) & grepl("^ENSG", ids2))
      if (length(ens_idx) > 0) {
        keys <- unique(ids2[ens_idx])
        bm <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = keys, mart = ensembl)
        if (nrow(bm) > 0) {
          bm_map <- setNames(as.character(bm$hgnc_symbol), as.character(bm$ensembl_gene_id))
          for (i in ens_idx) {
            k <- ids2[i]
            if (!is.na(k) && k %in% names(bm_map) && nzchar(bm_map[[k]])) sym[i] <- bm_map[[k]]
          }
        }
      }

      # 2) Entrez ID（纯数字）的情况
      entrez_idx <- which(is.na(sym) & grepl("^\\d+$", ids2))
      if (length(entrez_idx) > 0) {
        keys_e <- unique(ids2[entrez_idx])
        # 有些 biomaRt 版本使用 entrezgene_id 或 entrezgene; 尝试常见的 filters
        bm2 <- NULL
        try({ bm2 <- biomaRt::getBM(attributes = c("entrezgene_id", "hgnc_symbol"), filters = "entrezgene_id", values = keys_e, mart = ensembl) }, silent = TRUE)
        if (is.null(bm2) || nrow(bm2) == 0) {
          try({ bm2 <- biomaRt::getBM(attributes = c("entrezgene", "hgnc_symbol"), filters = "entrezgene", values = keys_e, mart = ensembl) }, silent = TRUE)
        }
        if (!is.null(bm2) && nrow(bm2) > 0) {
          # 取第一个 symbol
          keycol <- intersect(c("entrezgene_id", "entrezgene"), colnames(bm2))[1]
          bm_map2 <- setNames(as.character(bm2$hgnc_symbol), as.character(bm2[[keycol]]))
          for (i in entrez_idx) {
            k <- ids2[i]
            if (!is.na(k) && k %in% names(bm_map2) && nzchar(bm_map2[[k]])) sym[i] <- bm_map2[[k]]
          }
        }
      }

      # 3) 若还有剩余未映射，可尝试用 external_gene_name 作为 filters（当输入实际为 symbol 的不同形式）
      remain_idx <- which(is.na(sym))
      if (length(remain_idx) > 0) {
        keys_remain <- unique(ids2[remain_idx])
        bm3 <- tryCatch({ biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "external_gene_name", values = keys_remain, mart = ensembl) }, error = function(e) NULL)
        if (!is.null(bm3) && nrow(bm3) > 0) {
          bm_map3 <- setNames(as.character(bm3$hgnc_symbol), as.character(bm3$ensembl_gene_id))
          # 对于匹配到的 external_gene_name，getBM 返回 ensembl id 与 hgnc_symbol; 我们用 hgnc_symbol 来填充
          for (i in remain_idx) {
            k <- ids2[i]
            # find by hgnc symbol match
            matched <- bm3$hgnc_symbol[bm3$hgnc_symbol == k]
            if (length(matched) > 0 && nzchar(matched[1])) sym[i] <- matched[1]
          }
        }
      }
    }, silent = TRUE)
  }

  # MyGene.info 回退（保持可选）：对仍未映射的尝试在线查询
  remain_idx <- which(is.na(sym))
  if (length(remain_idx) > 0 && requireNamespace("httr", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE)) {
    try({
      keys_need <- unique(ids2[remain_idx])
      chunk_size <- 400
      chunks <- split(keys_need, ceiling(seq_along(keys_need)/chunk_size))
      mg_map <- list()
      for (ch in chunks) {
        ids_chunk <- paste(ch, collapse = ",")
        url <- sprintf("https://mygene.info/v3/genes/%s?fields=symbol&species=human", ids_chunk)
        resp <- httr::GET(url, httr::user_agent("R (mygene fallback)"), httr::timeout(30))
        if (httr::status_code(resp) != 200) next
        txt <- httr::content(resp, as = "text", encoding = "UTF-8")
        res_json <- jsonlite::fromJSON(txt)
        if (is.data.frame(res_json)) {
          for (i in seq_len(nrow(res_json))) {
            row <- res_json[i, , drop = TRUE]
            qid <- as.character(row[["_id"]])
            symb <- if ("symbol" %in% names(row)) as.character(row[["symbol"]]) else NA_character_
            if (!is.null(qid) && !is.na(qid) && nzchar(qid) && !is.null(symb) && !is.na(symb) && nzchar(symb)) mg_map[[qid]] <- symb
          }
        } else if (is.list(res_json)) {
          for (elem in res_json) {
            if (is.null(elem)) next
            qid <- NULL
            if (!is.null(elem[["_id"]])) qid <- as.character(elem[["_id"]])
            if (is.null(qid) && !is.null(elem[["query"]])) qid <- as.character(elem[["query"]])
            symb <- if (!is.null(elem[["symbol"]])) as.character(elem[["symbol"]]) else NA_character_
            if (!is.null(qid) && !is.na(qid) && nzchar(qid) && !is.null(symb) && !is.na(symb) && nzchar(symb)) mg_map[[qid]] <- symb
          }
        }
      }
      if (length(mg_map) > 0) {
        for (i in remain_idx) {
          key <- ids2[i]
          if (!is.na(key) && key %in% names(mg_map)) sym[i] <- mg_map[[key]]
        }
      }
    }, silent = TRUE)
  }

  # 最终仍未映射的，保留原始 id（避免丢失），并写出 unmapped 列表便于诊断
  remain_idx <- which(is.na(sym))
  if (length(remain_idx) > 0) {
    sym[remain_idx] <- ids_orig[remain_idx]
    unmapped_file <- file.path(results_dir, "counts_unmapped_ids.txt")
    tryCatch({
      writeLines(unique(ids_orig[remain_idx]), con = unmapped_file)
      message("Wrote unmapped IDs to: ", unmapped_file)
    }, error = function(e) {
      warning(sprintf("无法写出 unmapped id 文件：%s", e$message))
    })
  }

  names(sym) <- ids_orig
  return(as.character(sym))
}

# ------------------------------------------------------------------------------
# 使用原始 candidate 数据将 id 映射为 name（优先使用 candidate 表中的注释）
# 参数:
#  ids: 向量，基因 id（如 ENSG00000 或 Entrez）
#  base_dir: 项目根目录，用于在其下搜索 candidate 文件
#  results_dir: 输出诊断文件的目录
#  patterns: 在文件名中匹配 candidate 文件的关键词
map_ids_to_names_using_candidate <- function(ids, base_dir, results_dir, patterns = c("candidate","candidate_list","candidates","candidate_genes"), candidate_file = NULL) {
  ids_in <- as.character(ids)
  ids_strip <- sub("\\.\\d+$", "", ids_in)

  # 如果指定了显式的 candidate_file，则优先使用它
  cand_files <- character(0)
  if (!is.null(candidate_file) && nzchar(candidate_file)) {
    if (file.exists(candidate_file)) {
      cand_files <- candidate_file
    } else {
      warning(sprintf("指定的 candidate_file '%s' 不存在；将回退为在 base_dir 下搜索。", candidate_file))
    }
  }
  if (length(cand_files) == 0) {
    # 搜索 base_dir 下可能的 candidate 文件（递归）
    all_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
    cand_files <- all_files[grepl(paste(patterns, collapse = "|"), basename(all_files), ignore.case = TRUE)]
  }

  mapping <- data.frame(id = character(0), name = character(0), source = character(0), stringsAsFactors = FALSE)
  for (f in cand_files) {
    df <- NULL
    # 尝试常见的表格读取方法
    try({
      df <- tryCatch(read.delim(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (is.null(df)) df <- tryCatch(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    }, silent = TRUE)
    if (is.null(df) || ncol(df) < 2) next

    cols_low <- tolower(colnames(df))
    id_cols <- c("ensembl_gene_id","ensembl","gene_id","id","entrezgene","entrez","ensembl_id")
    name_cols <- c("gene_name","symbol","name","gene","hgnc_symbol")
    id_col <- colnames(df)[which(cols_low %in% id_cols)[1]]
    name_col <- colnames(df)[which(cols_low %in% name_cols)[1]]

    # 如果没有检测到常见列，则尝试基于列的数据类型推断（数值比例高的列更可能是 id/entrez）
    if (is.na(id_col) || is.na(name_col)) {
      if (ncol(df) >= 2) {
        numeric_prop <- sapply(df, function(col) {
          vals <- as.character(col)
          vals <- vals[!is.na(vals) & nzchar(vals)]
          if (length(vals) == 0) return(0)
          mean(!is.na(suppressWarnings(as.numeric(vals))))
        })
        col_id_idx <- which.max(numeric_prop)
        col_name_idx <- setdiff(seq_len(ncol(df)), col_id_idx)[1]
        id_col <- colnames(df)[col_id_idx]
        name_col <- colnames(df)[col_name_idx]
      }
    }

    if (!is.na(id_col) && !is.na(name_col)) {
      ids_vals <- as.character(df[[id_col]])
      names_vals <- as.character(df[[name_col]])
      ids_vals <- sub("\\.\\d+$", "", ids_vals)
      tmp <- data.frame(id = ids_vals, name = names_vals, source = basename(f), stringsAsFactors = FALSE)
      tmp <- tmp[!is.na(tmp$id) & nzchar(tmp$id) & !is.na(tmp$name) & nzchar(tmp$name), , drop = FALSE]
      if (nrow(tmp) > 0) mapping <- unique(rbind(mapping, tmp))
    }
  }

  if (nrow(mapping) == 0) {
    msg <- sprintf("No candidate mapping files found under '%s' (patterns: %s) or no usable columns detected.", base_dir, paste(patterns, collapse = ","))
    warning(msg)
    tryCatch(writeLines(msg, con = file.path(results_dir, "candidate_mapping_not_found.txt")), error = function(e) {})
    return(rep(NA_character_, length(ids_in)))
  }

  # 构建 id->name 字典（若重复 id 则保留第一个）
  id2name <- tapply(mapping$name, mapping$id, function(x) x[1])
  out <- vapply(ids_strip, function(k) {
    if (is.na(k) || !nzchar(k)) return(NA_character_)
    if (k %in% names(id2name)) return(as.character(id2name[[k]]))
    return(NA_character_)
  }, character(1))

  # 写出 mapping 诊断文件
  tryCatch(write.table(mapping, file = file.path(results_dir, "candidate_id2name_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE), error = function(e) {})

  return(out)
}

# ------------------------------------------------------------------------------
# 从计数文件中提取 gene_id -> gene_name 映射（如果计数文件包含 gene_name 列）
# 参数:
#  files: 含有计数文件路径的向量（优先从第一个文件开始尝试），
#  results_dir: 写出诊断文件路径
# 返回: 命名字符向量，names 为原始 id（去掉版本号），值为基因名；若未找到返回 NULL
extract_name_map_from_countfiles <- function(files, results_dir) {
  if (length(files) == 0) return(NULL)
  for (f in files) {
    # 尝试读取前 200 行以检测 header
    txt <- NULL
    try({
      if (grepl("\\.gz$", f, ignore.case = TRUE)) {
        con <- gzfile(f, open = "rt")
        txt <- readLines(con, n = 200)
        close(con)
      } else {
        txt <- readLines(f, n = 200)
      }
    }, silent = TRUE)
    if (is.null(txt) || length(txt) == 0) next

    # 尝试用 read.delim on the preview to get header and columns
    df <- NULL
    try({ df <- read.delim(text = txt, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) }, silent = TRUE)
    if (is.null(df) || ncol(df) < 2) {
      # 如果失败，尝试把整个文件读入小范围行做 header 识别
      next
    }
    cols_low <- tolower(colnames(df))
    if ("gene_name" %in% cols_low || "symbol" %in% cols_low || "name" %in% cols_low) {
      # 找到 name 列与 id 列
      name_col <- colnames(df)[which(cols_low %in% c("gene_name","symbol","name"))[1]]
      # id 列优先是 gene_id 或 ensembl_gene_id，否则第一列
      id_col <- colnames(df)[which(cols_low %in% c("gene_id","ensembl_gene_id","ensembl","id"))[1]]
      if (is.na(id_col) || !nzchar(id_col)) id_col <- colnames(df)[1]

      # 读取整个文件的这两列以确保完整性（支持大文件）
      df_full <- NULL
      try({ df_full <- read.delim(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE) }, silent = TRUE)
      if (is.null(df_full) || !(id_col %in% colnames(df_full)) || !(name_col %in% colnames(df_full))) next

      ids_vals <- as.character(df_full[[id_col]])
      names_vals <- as.character(df_full[[name_col]])
      # 清理 ensembl 版本号
      ids_vals2 <- sub("\\.\\d+$", "", ids_vals)
      keep_idx <- which(!is.na(ids_vals2) & nzchar(ids_vals2) & !is.na(names_vals) & nzchar(names_vals))
      if (length(keep_idx) == 0) next
      map_df <- data.frame(id = ids_vals2[keep_idx], name = names_vals[keep_idx], stringsAsFactors = FALSE)
      # 若有重复 id，保留第一个出现的 name
      id2name <- tapply(map_df$name, map_df$id, function(x) x[1])
      # 写出诊断文件
      tryCatch(write.table(map_df, file = file.path(results_dir, "countfile_id2name_mapping_preview.tsv"), sep = "\t", quote = FALSE, row.names = FALSE), error = function(e) {})
      return(id2name)
    }
  }
  return(NULL)
}

# 保存基因 x 样本的计数矩阵（供后续分析或查看）
# 我们将生成一个带有 gene_id 与 gene_name 的表格，便于外部查看
counts_tsv <- file.path(results_dir, "counts_matrix.tsv")
counts_rds <- file.path(results_dir, "counts_matrix.rds")
gene_ids <- rownames(countData)

# 如果不存在 candidate 文件，则从计数文件中自动创建
candidate_file_path <- file.path(base_dir, "candidate.txt")
if (!file.exists(candidate_file_path) && length(all_files) > 0) {
  cat("未找到 candidate 文件，正在从计数文件中自动创建...\n")
  # 选择第一个计数文件来提取映射
  sample_file <- all_files[1]
  tryCatch({
    # 读取文件，跳过以 N_ 开头的行
    df <- read.delim(sample_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
    # 过滤掉以 N_ 开头的行（这些是统计行，不是基因）
    df <- df[!grepl("^N_", df[[1]]), ]
    # 确保有 gene_id 和 gene_name 列
    cols_low <- tolower(colnames(df))
    id_col <- which(cols_low %in% c("gene_id", "ensembl_gene_id", "ensembl", "id"))[1]
    name_col <- which(cols_low %in% c("gene_name", "symbol", "name"))[1]
    if (!is.na(id_col) && !is.na(name_col)) {
      mapping_df <- data.frame(
        gene_id = df[[id_col]],
        gene_name = df[[name_col]],
        stringsAsFactors = FALSE
      )
      # 移除 NA 和空值
      mapping_df <- mapping_df[!is.na(mapping_df$gene_id) & nzchar(mapping_df$gene_id) & 
                               !is.na(mapping_df$gene_name) & nzchar(mapping_df$gene_name), ]
      # 移除重复
      mapping_df <- mapping_df[!duplicated(mapping_df$gene_id), ]
      # 写出文件
      write.table(mapping_df, file = candidate_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
      cat("已创建 candidate 文件：", candidate_file_path, "\n")
    } else {
      warning("无法从计数文件中提取 gene_id 和 gene_name 列")
    }
  }, error = function(e) {
    warning("自动创建 candidate 文件失败：", e$message)
  })
}

# 优先尝试使用原始 candidate 数据将 id -> name 映射（若存在 candidate 文件）
candidate_names_counts <- tryCatch({
  # 优先从计数文件中提取 gene_name（若存在）
  id2name_from_counts <- extract_name_map_from_countfiles(all_files, results_dir)
  if (!is.null(id2name_from_counts)) {
    ids_strip <- sub("\\.\\d+$", "", gene_ids)
    out <- vapply(ids_strip, function(k) if (!is.na(k) && k %in% names(id2name_from_counts)) id2name_from_counts[[k]] else NA_character_, character(1))
    return(out)
  }
  # 若计数文件中未找到，则尝试 candidate 文件（显式或搜索）
  map_ids_to_names_using_candidate(gene_ids, base_dir, results_dir, candidate_file = candidate_file)
}, error = function(e) { warning("candidate 映射失败（counts）：", e$message); return(rep(NA_character_, length(gene_ids))) })

# 退回到 biomaRt 等通用映射以补全未映射的 id
biomart_names_counts <- tryCatch({ map_ids_to_symbol(gene_ids) }, error = function(e) { warning("biomaRt 映射失败：", e$message); return(rep(NA_character_, length(gene_ids))) })

# 最终 gene_name：若 candidate 提供了注释则优先使用，否则使用 biomaRt 的注释（若有）
gene_symbols <- ifelse(!is.na(candidate_names_counts) & nzchar(candidate_names_counts), candidate_names_counts, biomart_names_counts)

# 构建带注释的 data.frame（保留原始计数矩阵用于分析）
# 显式设置 row.names = NULL，避免当 countData 的行名中存在 NA 时触发 "行名里有缺失值" 错误
if (any(is.na(gene_ids))) {
  warning("检测到 gene_ids 中存在 NA，将在输出表中使用默认行号作为行名以避免错误。")
}
counts_df <- data.frame(gene_id = gene_ids, gene_name = gene_symbols, as.data.frame(countData), check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
counts_df$candidate_mapped_name <- candidate_names_counts

# 写出带注释的 counts 矩阵及带 candidate 名称的备份
write.table(counts_df, file = counts_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(list(counts = countData, gene_id = gene_ids, gene_name = gene_symbols, candidate_mapped_name = candidate_names_counts), counts_rds)
counts_tsv2 <- file.path(results_dir, "counts_matrix_with_candidate_names.tsv")
tryCatch({ write.table(counts_df, file = counts_tsv2, sep = "\t", quote = FALSE, row.names = FALSE) }, error = function(e) {})
cat("Saved annotated count matrix to:\n", counts_tsv, "\n", counts_rds, "\n")
cat("Also saved counts with candidate names to:", counts_tsv2, "\n")

# ==============================================================================
# 3. 差异表达分析 (DEA) - 使用 DESeq2
# ==============================================================================

# 3.1 DESeq2 对象构建
# DESeq2 是RNA-seq数据差异表达分析的标准工具
# 使用 DESeqDataSetFromMatrix() 创建 DESeqDataSet 对象
# 参数说明：
# - countData: 计数矩阵 (基因 x 样本)，行是基因，列是样本
# - colData: 元数据数据框 (样本信息)，行名必须与countData列名匹配
# - design = ~ condition: 设计公式，指定实验设计中的变量
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)

# 3.2 质量控制与预过滤
# 过滤低表达基因：至少在最小样本数中有表达的基因
# 这里使用 rowSums(counts(dds) >= 10) >= 3 表示基因至少在3个样本中有>=10的计数
# 这个阈值可以根据实验设计和测序深度调整
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# 3.3 设置对照组并运行分析
# 使用 relevel() 将 Normal 组设为对照组 (参考水平)
# 这样 log2FoldChange 将是 Tumor 相对于 Normal 的变化
dds$condition <- relevel(dds$condition, ref = "Normal")

# 运行 DESeq() 函数进行完整的差异表达分析
# 这个函数执行以下步骤：
# 1. 估计尺寸因子 (size factors) 用于标准化
# 2. 估计分散度 (dispersion)
# 3. 进行负二项分布检验 (Wald test)
dds <- DESeq(dds)

# 3.4 结果提取
# 使用 results() 函数提取分析结果
# 先打印可用的系数名，便于确认模型中哪个系数对应我们感兴趣的对比
cat("DESeq2 resultsNames:", paste(resultsNames(dds), collapse = ", "), "\n")

# 明确指定对比：Tumor vs Normal
# 参数说明：
# - contrast = c("condition", "Tumor", "Normal"): 指定对比的变量名、分子水平、分母水平
# - alpha = 0.05: 用于多重检验校正的显著性水平
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)

# 保存完整的结果表（TSV 与 RDS），便于后续分析或复现
## 延迟写出 DESeq2 结果：先对基因 ID 做注释（gene_name），然后再写出带注释的结果表
## 这样输出中将包含 gene_id 与 gene_name 两列，便于外部查看与下游分析

# ==============================================================================
# 4. 结果筛选与可视化
# ==============================================================================

# 4.1 结果整理与筛选
# 将 DESeq2 结果转换为数据框，便于操作
res_df <- as.data.frame(res)

# 添加基因 id 列 (从行名提取) 并注释为 gene_name
res_df$gene_id <- rownames(res_df)
# 先尝试使用 candidate 数据映射 id -> name，再回退到通用映射
res_candidate_names <- tryCatch({
  # 优先尝试从计数文件中提取映射
  id2name_from_counts <- extract_name_map_from_countfiles(all_files, results_dir)
  if (!is.null(id2name_from_counts)) {
    ids_strip <- sub("\\.\\d+$", "", res_df$gene_id)
    out <- vapply(ids_strip, function(k) if (!is.na(k) && k %in% names(id2name_from_counts)) id2name_from_counts[[k]] else NA_character_, character(1))
    return(out)
  }
  # 否则尝试 candidate 文件
  map_ids_to_names_using_candidate(res_df$gene_id, base_dir, results_dir, candidate_file = candidate_file)
}, error = function(e) { warning("candidate 映射失败（res_df）：", e$message); return(rep(NA_character_, nrow(res_df))) })
res_biomart_names <- tryCatch({ map_ids_to_symbol(res_df$gene_id) }, error = function(e) { warning("biomaRt 映射失败（res_df）：", e$message); return(rep(NA_character_, nrow(res_df))) })
# 优先使用 candidate 名称
res_df$gene_name <- ifelse(!is.na(res_candidate_names) & nzchar(res_candidate_names), res_candidate_names, res_biomart_names)
res_df$candidate_mapped_name <- res_candidate_names

# 为兼容旧代码，保留 gene 列指向 gene_name
res_df$gene <- res_df$gene_name

# 创建 status 列，标记基因状态
# 初始化为 "Not Significant"
res_df$status <- "Not Significant"

# 上调：log2FoldChange > 1 且 padj < 0.05
res_df$status[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Upregulated"

# 下调：log2FoldChange < -1 且 padj < 0.05
res_df$status[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Downregulated"

# 写出带注释的完整 DESeq2 结果表（包含 gene_id 与 gene_name）
res_full_tsv <- file.path(results_dir, "DESeq2_results_Tumor_vs_Normal.annotated.tsv")
res_full_rds <- file.path(results_dir, "DESeq2_results_Tumor_vs_Normal.annotated.rds")
write.table(res_df, file = res_full_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(res_df, res_full_rds)
cat("Saved annotated DESeq2 results to:\n", res_full_tsv, "\n", res_full_rds, "\n")

# 4.2 火山图绘制（基础）
# 火山图是差异表达分析的经典可视化方法
# X 轴: log2FoldChange (表达变化倍数) - 正值表示上调，负值表示下调
# Y 轴: -log10(padj) (显著性，取负对数) - 值越大越显著
# 颜色: 根据表达状态进行分类
p <- ggplot(na.omit(res_df), aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point() +
  # 设置颜色映射：不显著基因用灰色，上调用红色，下调用蓝色
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  # 使用简洁主题
  theme_minimal()

# 4.3 火山图绘制（进阶）- Top 10 基因标注
# 为增强图表的信息量，标注最显著的基因
# 策略：按 padj 升序选择前10个最显著的基因

# 首先过滤掉 padj 为 NA 的基因
res_nonNA <- res_df[!is.na(res_df$padj), ]
if (nrow(res_nonNA) == 0) {
  warning("没有可用于绘图的显著性(padj)值；火山图将为空。")
  top10 <- res_df[0, ]  # 创建空数据框
} else {
  # 按 padj 升序排序，取前10个
  top10 <- head(res_nonNA[order(res_nonNA$padj), ], 10)
}

# 为 top10 基因添加标注信息
if (nrow(top10) > 0) {
  top10$rank <- seq_len(nrow(top10))  # 添加排名
  top10$padj_fmt <- formatC(top10$padj, format = "e", digits = 2)  # 格式化p值

  # 创建标签文本：包含排名、基因名和p值
  top10$label_text <- paste0(top10$rank, ". ", top10$gene_name, "\n(padj=", top10$padj_fmt, ")")

  # 写出详细的 top10 基因表格，便于进一步分析
  top10_out <- top10[, c("rank", "gene_id", "gene_name", "log2FoldChange", "pvalue", "padj", "status")]
  write.table(top10_out, file = file.path(results_dir, "top10_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # 在火山图上突出显示 Top10 基因
  # 使用更大的点和边框来突出
  p <- p +
    geom_point(data = top10, aes(x = log2FoldChange, y = -log10(padj), color = status), size = 3.5, stroke = 0.6)

  # 添加排名数字标签
  p <- p +
    geom_text(data = top10, aes(x = log2FoldChange, y = -log10(padj), label = rank),
              color = "black", size = 3, vjust = -1)

  # 使用 geom_label_repel 添加基因名标签
  # 参数说明：
  # - color/white: 文字颜色
  # - box.padding: 标签周围的填充空间
  # - segment.size: 连接线的粗细
  # - alpha: 透明度
  # - fontface: 字体样式
  p <- p + geom_label_repel(data = top10,
                            aes(x = log2FoldChange, y = -log10(padj), label = gene_name, fill = status),
                            color = "white", size = 3.0, max.overlaps = 30,
                            box.padding = 0.35, segment.size = 0.45, segment.color = "grey40",
                            label.size = 0.25, label.r = 0.15, alpha = 0.95, fontface = "bold", show.legend = FALSE)

  # 为标签设置填充颜色（深色，便于白字显示）
  p <- p + scale_fill_manual(values = c("Not Significant" = "#000000", "Upregulated" = "#9B111E", "Downregulated" = "#1E40AF"))
}

# 保存火山图为 PNG 文件到 results 目录
# 参数说明：
# - width = 8, height = 6: 图表尺寸 (英寸)，适合论文发表
# - 默认分辨率为 300 DPI，适合高质量打印
volcano_file <- file.path(results_dir, "volcano_plot.png")
ggsave(volcano_file, p, width = 8, height = 6)

# ==============================================================================
# 5. 交付成果与总结
# ==============================================================================

# 5.1 结果文件总结
# 脚本生成的主要输出文件：
# 1. metadata.tsv/rds - 样本元数据
# 2. counts_matrix.tsv/rds - 基因表达计数矩阵
# 3. DESeq2_results_Tumor_vs_Normal.annotated.tsv/rds - 差异表达分析结果
# 4. top10_genes.tsv - Top 10 最显著差异表达基因
# 5. volcano_plot.png - 火山图
# 6. upregulated_TAA_list.txt - 上调基因列表

# 5.2 生成上调基因列表
# 创建一个简单的文本文件，只包含上调基因名（每行一个基因）
# 这通常用于下游分析，如通路富集分析或文献检索
upregulated <- res_df$gene[res_df$status == "Upregulated"]
upregulated_file <- file.path(results_dir, "upregulated_TAA_list.txt")
writeLines(upregulated, upregulated_file)

# 5.3 运行完成提示
# 打印结果目录路径，便于用户快速找到输出文件
cat("Results written to:", results_dir, "\n")
cat("Analysis completed successfully!\n")
