# BioEASE 🧬

> 一站式生物信息学数据分析平台 — 从 TCGA 数据下载到生存分析，零门槛开箱即用。

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.0%2B-276DC3?logo=r)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-TCGA%20%7C%20RNA--Seq-9C27B0)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux-blue?logo=windows)

## ✨ 功能特性

### 🔬 数据获取
- **TCGA 数据自动下载** — 一键从 GDC 接口获取 TCGA 肿瘤样本的基因表达矩阵、临床元数据
- **多癌种覆盖** — 支持所有 TCGA 收录的 33+ 种肿瘤类型

### 📊 差异表达分析
- **自动化差异分析** — 对正常组织 vs 肿瘤组织进行差异表达基因（DEG）筛选
- **多重校正** — 内置 FDR 校正，支持 p-value / Log2FC 阈值自定义
- **可视化结果** — Volcano Plot / Heatmap 等图表直接导出

### 🧪 生存分析
- **Kaplan-Meier 曲线** — 按基因表达高低分组，绘制生存曲线
- **Cox 比例风险模型** — 多因素生存回归，输出 HR / 95%CI / p-value
- **TAA 生存分析** — 肿瘤相关抗原（TAA）专项分析，筛选预后biomarker

### 🖥️ 可视化界面
- **Streamlit 前端** — 浏览器直连，数据上传、分析、导出全程可视化
- **无需 R 语言基础** — 所有 R 分析引擎化，一键调用

## 🎯 适用场景

| 场景 | 说明 |
|------|------|
| 🎓 科研入门 | 生信零基础也能完成完整的 TCGA 分析流程 |
| 📝 论文写作 | 直接输出 publication-ready 图表 |
| 🔬 生物课题 | 快速筛选癌症预后相关基因 |
| 🏥 临床转化 | 探索 TAA 作为潜在biomarker |

## 📁 项目结构

```
bioease/
├── download_TCGA_data.R     # TCGA 数据下载脚本
├── analysis.R               # 差异表达 + 生存分析
├── upregulated_TAA_survival_analysis.R  # TAA 专项分析
├── front.py                 # Streamlit Web 界面
├── launcher.py              # 一键启动器
├── pyproject.toml           # Python 依赖配置
└── results/                 # 分析结果输出目录
```

## 🚀 快速开始

### 前置依赖

- Python 3.10+
- R 4.0+
- R 包：`TCGAbiolinks`, `survival`, `ggplot2`, `patchwork`

### 安装 R 依赖

```bash
# 在 R 控制台运行
install.packages(c("TCGAbiolinks", "survival", "ggplot2", "patchwork"))
```

### 安装 Python 依赖

```bash
pip install streamlit requests pandas numpy
```

### 运行

```bash
# 1. 下载 TCGA 数据
Rscript download_TCGA_data.R

# 2. 执行分析
Rscript analysis.R

# 3. 启动可视化界面
python front.py
# 浏览器打开 http://localhost:8501
```

## 📈 输出示例

- `volcano_plot.pdf` — 差异基因火山图
- `heatmap.pdf` — 差异基因热图
- `survival_plot.pdf` — Kaplan-Meier 生存曲线
- `cox_results.csv` — Cox 回归结果表
- `results.rds` — 完整分析结果（支持后续自定义可视化）

## 🛠️ 技术栈

| 层级 | 技术 |
|------|------|
| 数据下载 | GDC API (`TCGAbiolinks`) |
| 统计分析 | R (`survival`, `ggplot2`) |
| 可视化界面 | Python (`Streamlit`) |
| 数据格式 | RDS / CSV / TSV |

## 📖 相关文档

- [TCGA Official](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) — TCGA 数据库官方文档
- [TCGAbiolinks](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html) — Bioconductor TCGA 分析包
- [Survival Analysis in R](https://www.empirica.biz/en/survival-analysis-with-r/) — R 语言生存分析教程

## 📝 License

MIT License &copy; 2024 [RomanCohort](https://github.com/RomanCohort)

---

*如果你觉得 BioEASE 有帮助，请点一个 ⭐ Star！*
