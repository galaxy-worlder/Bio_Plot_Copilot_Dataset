import requests
import os
import time

# ================= 配置区 =================
# 替换为你刚才申请的 GitHub Token
GITHUB_TOKEN = "your_token_here"

# 定义你的搜索关键词组合 (这是这个脚本的灵魂)
# 语法：关键词 + 扩展名 + 甚至是指定的科研机构/用户
SEARCH_QUERIES = [
    '"Nature Communications" AND "ggplot2" extension:R',
    '"ComplexHeatmap" AND "ggseqlogo" extension:R',
    '"Cell" AND "Seurat" AND "FeaturePlot" extension:R'
]

# 下载保存的基础目录，刚好契合你项目里的暂存区概念
SAVE_DIR = "00_Raw_Code_Pool"
# ==========================================

HEADERS = {
    "Authorization": f"token {GITHUB_TOKEN}",
    "Accept": "application/vnd.github.v3+json"
}

def create_dir_if_not_exists(path):
    if not os.path.exists(path):
        os.makedirs(path)

def download_raw_file(raw_url, save_path):
    """根据 raw URL 下载文件并保存"""
    try:
        response = requests.get(raw_url, headers=HEADERS, timeout=10)
        if response.status_code == 200:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write(response.text)
            print(f"  [+] 成功下载: {os.path.basename(save_path)}")
        else:
            print(f"  [-] 下载失败, 状态码: {response.status_code}")
    except Exception as e:
        print(f"  [!] 下载异常: {e}")

def convert_to_raw_url(html_url):
    """
    将 GitHub 网页链接转换为 Raw 下载链接
    例如: https://github.com/user/repo/blob/master/script.R 
    转换为: https://raw.githubusercontent.com/user/repo/master/script.R
    """
    return html_url.replace("github.com", "raw.githubusercontent.com").replace("/blob/", "/")

def search_and_download(query):
    print(f"\n🚀 开始执行检索任务: {query}")
    
    # 构建当前 query 的专属保存文件夹，按合法名称清理
    safe_folder_name = query.replace('"', '').replace(':', '_').replace(' ', '_')[:30]
    target_dir = os.path.join(SAVE_DIR, safe_folder_name)
    create_dir_if_not_exists(target_dir)

    search_url = "https://api.github.com/search/code"
    
    # 注意：为了避免一次性拉取太多被封，这里演示拉取前 20 条结果
    params = {
        "q": query,
        "per_page": 20, 
        "page": 1
    }

    try:
        response = requests.get(search_url, headers=HEADERS, params=params)
        
        if response.status_code != 200:
            print(f"❌ 检索失败: {response.json().get('message', '未知错误')}")
            return

        data = response.json()
        items = data.get('items', [])
        print(f"🔍 找到 {len(items)} 个匹配的文件，准备下载...\n")

        for index, item in enumerate(items):
            repo_name = item['repository']['full_name'].replace("/", "_")
            file_name = item['name']
            html_url = item['html_url']
            
            # 构建一个具有辨识度的文件名：仓库名_原文件名
            safe_file_name = f"{repo_name}_{file_name}"
            save_path = os.path.join(target_dir, safe_file_name)

            # 获取 raw 链接并下载
            raw_url = convert_to_raw_url(html_url)
            print(f"[{index+1}/{len(items)}] 正在处理来源: {item['repository']['full_name']}")
            download_raw_file(raw_url, save_path)
            
            # 🚦 礼貌性延时，保护你的 GitHub 账号不被 Rate Limit 触发
            time.sleep(1.5) 

    except Exception as e:
        print(f"任务执行异常: {e}")

if __name__ == "__main__":
    create_dir_if_not_exists(SAVE_DIR)
    
    # 遍历我们定义的查询组合
    for query in SEARCH_QUERIES:
        search_and_download(query)
        # 两个大查询之间休息一下
        time.sleep(5)
        
    print("\n🎉 所有自动抓取任务已完成！快去查看你的语料库吧。")