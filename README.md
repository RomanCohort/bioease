# BioEASE

生物信息学数据分析工具，用于 TCGA 数据下载、差异表达分析与生存分析可视化。

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![TCGA](https://img.shields.io/badge/TCGA-Data-yellow)

## 功能

- **TCGA 数据下载**：自动化获取癌症基因组数据
- **差异表达分析**：基因表达量比较
- **生存分析**：Kaplan-Meier 曲线与 Cox 回归
- **可视化**：R ggplot2 图表生成
- **前端界面**：Python GUI 数据分析操作

## 项目结构

```
├── front.py                            # Python 前端界面
├── launcher.py                         # 启动器
├── download_TCGA_data.R                # TCGA 数据下载脚本
├── analysis.R                          # 主分析流程
├── upregulated_TAA_survival_analysis.R # 上调基因生存分析
├── df.rds / results.rds                # 预处理 R 数据对象
├── MANIFEST.txt                        # TCGA 数据清单
└── README.md
```

## 环境要求

- Python 3.10+
- R >= 4.0 及以下包：
  - `TCGAbiolines`
  - `survival`, `survminer`
  - `ggplot2`, `dplyr`

## 运行

```bash
# 前端界面
python front.py

# 或直接运行 R 分析
Rscript analysis.R
```

## 数据

项目包含乳腺癌 (BRCA) TCGA 数据分析示例。

## 许可证

MIT License

## Related Projects

| Project | Description |
|---------|-------------|
| [paper-search-tool](https://github.com/RomanCohort/paper-search-tool) | AI 论文搜索与整理工具 |
| [ai-desktop-pet](https://github.com/RomanCohort/ai-desktop-pet) | AI 桌面宠物 |
| [web-crawler-v2](https://github.com/RomanCohort/web-crawler-v2) | 网站爬取器 |
| [berlin-tank-commander](https://github.com/RomanCohort/berlin-tank-commander) | 柏林车长文字冒险 |
| [bioease](https://github.com/RomanCohort/bioease) | 生物信息学分析 |
| [IGEM-sama](https://github.com/RomanCohort/IGEM-sama) | IGEM AI 虚拟主播 |
| [ppt-agent](https://github.com/RomanCohort/ppt-agent) | PPT 草稿生成器 |
