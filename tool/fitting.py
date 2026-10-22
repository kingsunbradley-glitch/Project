import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import sys

def parse_input(prompt_text):
    """处理输入，支持空格或逗号分隔"""
    try:
        raw_str = input(prompt_text)
        # 兼容中文逗号
        raw_str = raw_str.replace('，', ',')
        if ',' in raw_str:
            data = [float(x) for x in raw_str.split(',')]
        else:
            data = [float(x) for x in raw_str.split()]
        return np.array(data)
    except ValueError:
        print("❌ 输入格式错误，请输入纯数字。")
        sys.exit(1)

def main():
    print("=== 最小二乘法工具 ===")
    
    # 1. 获取输入
    # 注意：如果 Y 使用默认的3个数据，X 也必须输入3个数据，否则会报错
    x = np.sort(parse_input("请输入三个峰的道值 (顺序任意，空格或逗号隔开): "))
    
    # --- 修改开始：默认 Y 值设置 ---
    # 模式 A: 使用默认值 (当前开启)
    print("👉 正在使用默认 Y 数据: [5156.6, 5485.6, 5804.8]")
    y = np.array([5156.6, 5485.6, 5804.8])

    # 模式 B: 手动输入 (若要使用，请取消下面这行的注释，并注释掉上面两行)
    # y = parse_input("请输入 y 组数据 (空格或逗号隔开): ")
    # --- 修改结束 ---

    if len(x) != len(y):
        print(f"❌ 错误：数据长度不一致 (x={len(x)}, y={len(y)})")
        print(f"💡 提示：当前默认 Y 包含 {len(y)} 个数据，请确保 X 也输入 {len(y)} 个数据。")
        return

    # 2. 计算统计量
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_squared = r_value ** 2

    # 3. 控制台输出中文统计报告
    print("-" * 40)
    print(f"📈 拟合方程: y = {slope:.4f}x + {intercept:.4f}")
    print(f"📊 统计结果:")
    print(f"   R-squared (R²): {r_squared:.4f}")
    print(f"   P-value       : {p_value:.4e}")
    print(f"   Std Error     : {std_err:.4f}")
    print("-" * 40)

    # 4. 绘图 (严格使用英文，避免乱码)
    plt.figure(figsize=(10, 6))
    
    # 绘制散点 (Data Points)
    plt.scatter(x, y, label='Original Data', color='blue', alpha=0.6)
    
    # 绘制拟合线 (Fit Line)
    x_fit = np.linspace(min(x), max(x), 100)
    y_fit = slope * x_fit + intercept
    
    # 图例标签使用英文公式
    label_str = f'Fit: $y={slope:.3f}x + {intercept:.3f}$'
    if intercept < 0:
        label_str = f'Fit: $y={slope:.3f}x - {abs(intercept):.3f}$'
        
    plt.plot(x_fit, y_fit, color='red', linewidth=2, label=label_str)

    # 标题和坐标轴全部英文
    plt.title(f'Linear Regression Analysis ($R^2={r_squared:.3f}$)', fontsize=14)
    plt.xlabel('X Variable', fontsize=12)
    plt.ylabel('Y Variable', fontsize=12)
    
    # 显示图例、网格
    plt.legend(fontsize=11)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    #print("🖼️  正在显示图像 (Plot window)...")
    #plt.show()

if __name__ == "__main__":
    main()
