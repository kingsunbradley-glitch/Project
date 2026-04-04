import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 1. 读取并清理数据
data = []
with open('HsData.txt', 'r', encoding='utf-8') as f:
    for line in f:
        # 忽略注释和空行
        line = line.split('#')[0].strip()
        if not line or 'energy_keV' in line:
            continue
        
        parts = line.split()
        
        # 逐个提取符合 [energy, sigma, source, label] 格式的数据块
        while len(parts) >= 4:
            try:
                energy = float(parts[0])
                sigma = float(parts[1])
                source = parts[2]
                label = parts[3]
                data.append({
                    'energy': energy, 
                    'sigma': sigma, 
                    'source': source, 
                    'label': label
                })
            except ValueError:
                pass
            parts = parts[4:]

df = pd.DataFrame(data)

# 2. 创建画布
plt.figure(figsize=(12, 7))

# X轴为数据点序号（从1开始）
x = np.arange(1, len(df) + 1)
y = df['energy']
yerr = df['sigma']

# 3. 绘制带有误差棒的散点图（按 EVR 和 Decay 分类设定颜色）
colors = {'EVR': '#1f77b4', 'Decay': '#2ca02c'} # 蓝色和绿色
for src in df['source'].unique():
    mask = df['source'] == src
    plt.errorbar(x[mask], y[mask], yerr=yerr[mask], fmt='o', 
                 color=colors.get(src, 'black'), ecolor='gray', 
                 capsize=5, elinewidth=1.5, markeredgewidth=1.5, 
                 label=f'Data ({src})')

# 4. 绘制拟合结果的基准线
# 单峰 (M1) 红色虚线
plt.axhline(9164, color='red', linestyle='--', linewidth=2, label='Single Peak (M1): 9164 keV')
# 双峰 (M2) 红色实线
plt.axhline(9108, color='red', linestyle='-', linewidth=2, label='Double Peak (M2): 9108 keV')
plt.axhline(9199, color='red', linestyle='-', linewidth=2, label='Double Peak (M2): 9199 keV')

# 5. 图表格式化设置
plt.xlabel('Data Point Index', fontsize=12)
plt.ylabel('Energy (keV)', fontsize=12)
plt.title('Experimental Data Distribution vs Fitted Peaks', fontsize=14)
plt.xticks(x) # 让X轴只显示整数序号
plt.legend(loc='best', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.7)

# 6. 调整布局并显示/保存
plt.tight_layout()
plt.show()
# 如果你想直接保存为图片，可以取消下面这行的注释：
# plt.savefig('peak_comparison.png', dpi=300)