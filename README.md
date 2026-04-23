
🧬 Bio-Plot Copilot Dataset

![Python](https://img.shields.io/badge/Python-3.11-blue)
![R](https://img.shields.io/badge/R-ggplot2-276DC3)
![LLM](https://img.shields.io/badge/FineTuning-LLaMA--Factory-orange)

📌 简介 (Introduction)
本项目是一个专为**生物信息学与结构生物学数据可视化**打造的高质量指令微调数据集仓库。

通过持续收集顶级期刊 (CNS) 级别的绘图代码（如复杂多层热图、双样本 Logo 图、带显著性检验的差异柱状图等），并将其转化为标准的 Alpaca JSON 格式，本数据集旨在微调开源大语言模型 (如 Qwen2.5)，打造一个精通 R 语言生信绘图的“专属私有化代码副驾驶 (Copilot)”。

📂 仓库结构 (Repository Structure)
采用“本地模块化语料库 + 自动化脚本打包”的双层管理架构：

```text
Bio_Plot_Copilot_Dataset/
├── 01_Heatmaps/                 # 热图类代码标本
├── 02_Scatter_Plots/            # 散点图/火山图代码标本
├── 03_Structural_Biology/       # 结构生物学特征可视化
├── build_dataset.py             # 🌟 核心打包脚本
├── bio_copilot_dataset.json     # 自动生成的 Alpaca 格式微调文件
└── .gitignore                   # 严格防范大文件/隐私数据泄露
```

🚀 工作流 (Workflow)

本仓库采用极简的“三步曲”数据累积工作流：

1. 存入标本 (Daily Accumulation)
在对应的分类目录下新建文件夹，并放入“核心三件套”：
* `meta.md`: 用自然语言描述绘图需求（作为 Instruction）。
* `code.R`: 高质量、可运行的 R 语言绘图脚本（作为 Output）。
* `mock_data.csv`: 极小型的脱敏测试数据表头（作为 Input 上下文）。

2. 一键炼丹准备 (Build Dataset)
在根目录下运行自动化构建脚本：
```bash
python build_dataset.py
```
脚本会自动遍历所有目录，将其压榨并组装成 LLaMA-Factory 可直接读取的 `bio_copilot_dataset.json`。

3. 微调训练 (Fine-Tuning)
将生成的 JSON 文件挂载至 LLaMA-Factory 环境中，启动 LoRA 微调，将新的生信制图知识注入大模型。

💡 数据集样例 (Data Sample)
自动生成的 Alpaca 数据格式示例：
```json
[
  {
    "instruction": "作为生物信息学代码助手，请编写一个 R 语言函数，使用 ggplot2 绘制带显著性标注的靶点富集差异柱状图。",
    "input": "输入数据格式如下 (CSV 预览):\nposition,amino_acid,freq_diff,p_value,sig_label\n6,Y,0.15,0.01,**\n...",
    "output": "```r\nlibrary(ggplot2)...\n```"
  }
]
```

 🛡️ 注意事项
 严禁上传真实测序数据或未经脱敏的临床数据。** 仓库的 `.gitignore` 已配置拦截所有 `.csv`, `.rds` 等大型文件，仅允许名为 `mock_data` 的微型样例数据入库。
```
