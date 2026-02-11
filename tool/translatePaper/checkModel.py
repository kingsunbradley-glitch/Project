import google.generativeai as genai
import os

# ================= 1. 代理设置 (必填) =================
# 必须和 translatePaper.py 里保持一致，否则连不上
os.environ['HTTP_PROXY'] = 'http://127.0.0.1:7890'
os.environ['HTTPS_PROXY'] = 'http://127.0.0.1:7890'
# =====================================================

# 2. 配置 API Key
# 记得把这里换成你的真实 Key
genai.configure(api_key="AIzaSyD-b9YY9KP9vDtFC6f5eyjq0XcPPbucjwo") 

print("正在查询可用模型列表...\n")

try:
    # 列出所有模型
    for m in genai.list_models():
        # 我们只关心能“生成内容”的模型 (过滤掉只能做 embedding 的模型)
        if 'generateContent' in m.supported_generation_methods:
            print(f"名称: {m.name}")
            print(f"描述: {m.description}")
            print("-" * 30)
            
except Exception as e:
    print(f"查询失败: {e}")