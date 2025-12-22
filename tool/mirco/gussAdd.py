import numpy as np
import matplotlib.pyplot as plt

def gaussian(x, amp, mu, sigma):
    """
    定义高斯函数
    x: 自变量数组
    amp: 振幅 (峰的高度)
    mu: 均值 (峰的中心位置)
    sigma: 标准差 (峰的宽度)
    """
    return amp * np.exp(-((x - mu)**2) / (2 * sigma**2))

# ===========================
# 1. 在这里设置您的参数
# ===========================

# 第一个波峰 (Peak 1)
amp1 = 750     # 高度
mu1 =  6112.72     # 中心位置
sigma1 = 40/2.35   # 宽度 (数值越小越尖锐)

# 第二个波峰 (Peak 2)
amp2 = 250    # 高度
mu2 = 6069.43      # 中心位置
sigma2 = 40/2.35   # 宽度

# 设置 X 轴的范围 (从 -6 到 6)
x = np.linspace(6030, 6150, 1000)

# ===========================
# 2. 计算数据
# ===========================
y1 = gaussian(x, amp1, mu1, sigma1)
y2 = gaussian(x, amp2, mu2, sigma2)
y_total = y1 + y2  # 线性叠加

# ===========================
# 3. 绘制图像
# ===========================
plt.figure(figsize=(10, 6))

# 绘制第一个峰（绿色虚线）
plt.plot(x, y1, 'g--', label=f'Peak 1 ($\mu={mu1}$)', alpha=0.6)
plt.fill_between(x, y1, color='green', alpha=0.1) # 填充颜色

# 绘制第二个峰（蓝色虚线）
plt.plot(x, y2, 'b--', label=f'Peak 2 ($\mu={mu2}$)', alpha=0.6)
plt.fill_between(x, y2, color='blue', alpha=0.1) # 填充颜色

# 绘制叠加后的峰（红色实线）
plt.plot(x, y_total, 'r-', linewidth=2.5, label='Sum (Superposition)')

# 图表装饰
plt.title("Two Gaussian Peaks and their Superposition")
plt.xlabel("Position (x)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)

# 显示图表
plt.show()