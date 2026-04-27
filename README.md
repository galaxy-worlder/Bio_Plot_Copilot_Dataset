
# 🧬 Bio-Plot Copilot Dataset
> “流星世界演绎者的生信绘图专属微调数据集” —— 致力于让大语言模型成为顶尖的科研插画师。

![Python](https://img.shields.io/badge/Python-3.11-blue)
![R](https://img.shields.io/badge/R-ggplot2-276DC3)
![LLM](https://img.shields.io/badge/FineTuning-LLaMA--Factory-orange)
![License](https://img.shields.io/badge/License-MIT-green.svg)

📌 简介 (Introduction)
本项目是一个专为**生物信息学与结构生物学数据可视化**打造的高质量指令微调数据集仓库。

通过持续收集顶级期刊 (CNS) 级别的绘图代码（如复杂多层热图、双样本 Logo 图、带显著性检验的差异柱状图等），并将其转化为标准的 Alpaca JSON 格式，本数据集旨在微调开源大语言模型 (如 Qwen2.5)，打造一个精通 R 语言生信绘图的“专属私有化代码副驾驶 (Copilot)”。

📂 仓库结构 (Repository Structure)
采用“本地模块化语料库 + 自动化脚本采集与打包”的双层管理架构：

```text
Bio_Plot_Copilot_Dataset/
├── 00_Raw_Code_Pool/            # 🤖 [自动] GitHub API 爬取的待清洗 R 源码池
├── 01_Heatmaps/                 # 📊 热图类代码标本
├── 02_Scatter_Plots/            # 📉 散点图/火山图代码标本
├── 03_Structural_Biology/       # 🧬 结构生物学特征可视化代码标本
├── github_r_scraper.py          # 🚀 高分 R 代码自动化收集爬虫
├── build_dataset.py             # 🌟 核心 JSON 打包脚本
├── bio_copilot_dataset.json     # 自动生成的 Alpaca 格式微调文件
└── .gitignore                   # 严格防范大文件/隐私数据泄露
```

🚀 工作流 (Workflow)

本仓库采用“自动化采集 + 三步曲”的高效数据累积工作流：

### 阶段零：自动化采集代码池 (Auto-Scraping)
告别手动低效搜索。通过根目录的爬虫脚本，调用 GitHub API 批量抓取高分期刊和顶级实验室的 `.R` 源码。
1. 在本地配置 `.env` 文件填入 `GITHUB_TOKEN`。
2. 运行 `python github_r_scraper.py`。
3. 脚本会根据设定的关键词（如 `"ComplexHeatmap" AND "ggseqlogo"`），自动将优质脚本下载至 `00_Raw_Code_Pool/` 供后续人工筛选。

### 第一步：清洗并存入标本 (Daily Accumulation)
从源码池中筛选出优秀代码，在对应的分类目录（如 `01_Heatmaps`）下新建文件夹，并放入“核心三件套”：
* `meta.md`: 用自然语言描述绘图需求（作为 Instruction）。
* `code.R`: 高质量、可运行且经过格式化的 R 语言绘图脚本（作为 Output）。
* `mock_data.csv`: 极小型的脱敏测试数据表头（作为 Input 上下文，防止模型产生幻觉）。

### 第二步：一键炼丹准备 (Build Dataset)
在根目录下运行自动化构建脚本：
```bash
python build_dataset.py
```
脚本会自动遍历所有标本目录，将其压榨并组装成 LLaMA-Factory 可直接读取的 `bio_copilot_dataset.json`。

### 第三步：微调训练 (Fine-Tuning)
将生成的 JSON 文件挂载至 LLaMA-Factory 环境中，启动 LoRA 微调，将新的生信制图知识注入大模型。

💡 数据集样例 (Data Sample)
自动生成的 Alpaca 数据格式示例：
```json
[
  {
    "instruction": "作为生物信息学代码助手，请编写一个 R 语言函数，使用 ggplot2 绘制带显著性标注的靶点富集差异柱状图。",
    "input": "输入数据格式如下 (CSV 预览):\nposition,amino_acid,freq_diff,p_value,sig_label\n6,Y,0.15,0.01,**\n...",
    "output": "```r\nlibrary(ggplot2)\n# ...高质量绘图代码...\n```"
  }
]
```

🗺️ 未来路线图 (Roadmap)
- [x] 确立项目架构与模块化“三件套”语料库方案
- [x] 编写自动化 GitHub 高质量 R 代码采集脚本
- [ ] 扩充 `ggseqlogo` 序列特征可视化的核心数据集
- [ ] 在数据集中引入常见的“报错-修复”对（Negative Samples）
- [ ] 引入 CoT (思维链)，让模型在输出代码前先输出绘图逻辑

🛡️ 安全与隐私注意事项
**严禁上传真实测序数据或未经脱敏的临床数据！** 本仓库的 `.gitignore` 已配置严格规则，默认拦截所有 `.csv`, `.rds`, `.RData` 等大型文件，仅允许各标本文件夹下名为 `mock_data.csv` 的微型样例数据入库。

🤝 参与贡献 (Contributing)
如果你同样被糟糕的生信作图折磨过，或者手里有引以为傲的 R 语言绘图脚本，欢迎提交 Pull Request，一起充实这个 Copilot 的“知识库”！
```

***

这份更新后的版本不仅保留了你原创的“核心三件套”这种非常出色的工程化表达，还把自动化爬虫无缝衔接了进去，整体的专业度直接拉满了！你可以直接覆盖提交了。
