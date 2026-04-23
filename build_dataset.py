import os
import json
import pandas as pd

# 配置区
REPO_DIR = "."  # 当前仓库目录
OUTPUT_FILE = "bio_copilot_dataset.json"

def build_dataset():
    dataset = []
    
    # 遍历当前目录下所有的子文件夹
    for root, dirs, files in os.walk(REPO_DIR):
        # 忽略隐藏文件夹（如 .git）
        if "/." in root or root == REPO_DIR:
            continue
            
        # 检查是否包含核心三件套
        if "meta.md" in files and "code.R" in files:
            print(f"正在处理: {root}")
            
            # 1. 提取 Instruction (来自 meta.md)
            with open(os.path.join(root, "meta.md"), "r", encoding="utf-8") as f:
                instruction = f.read().strip()
                
            # 2. 提取 Input
            input_text = ""
            if "mock_data.csv" in files:
                try:
                    df = pd.read_csv(os.path.join(root, "mock_data.csv"), nrows=2)
                    input_text = "输入数据格式如下 (CSV 预览):\n" + df.to_csv(index=False)
                except Exception as e:
                    input_text = "暂无有效输入数据。"
            
            # 3. 提取 Output (来自 code.R)
            with open(os.path.join(root, "code.R"), "r", encoding="utf-8") as f:
                output_code = "```r\n" + f.read().strip() + "\n```"
                
            # 组装成 Alpaca 格式
            dataset.append({
                "instruction": instruction,
                "input": input_text,
                "output": output_code
            })
            
    # 保存为 JSON 文件
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(dataset, f, ensure_ascii=False, indent=2)
    print(f"\n✅ 成功生成数据集！共收录 {len(dataset)} 条微调数据。")
    print(f"文件已保存至: {OUTPUT_FILE}")

if __name__ == "__main__":
    build_dataset()
