import pandas as pd
import matplotlib.pyplot as plt

def main():
    # 读取数据，假设数据是以制表符(Tab)或空格分隔
    try:
        df = pd.read_csv('data.dat', sep='\s+')
    except FileNotFoundError:
        print("错误：未找到 'data.dat' 文件，请确保该文件与脚本在同一目录下。")
        return

    # 提取需要的列
    N = df['N']
    T_exp = df['T0.5']
    err_plus = df['+']
    err_minus = df['-']
    T_urf = df['URF']

    # 创建图形
    plt.figure(figsize=(10, 6))

    # 绘制实验半衰期数据（带误差棒的点，【不使用线条连接】）
    # 设置 linestyle='none' 即可去掉连接线
    plt.errorbar(N, T_exp, yerr=[err_minus, err_plus], 
                 fmt='o', color='blue', linestyle='none', 
                 capsize=4, label='Experimental $T_{1/2}$')

    # 绘制URF理论预测半衰期数据（空心点，【使用虚线连接】）
    # markerfacecolor='none' 用于实现空心点
    plt.plot(N, T_urf, marker='o', color='red', linestyle='--', linewidth=1.5, 
             markerfacecolor='none', markeredgecolor='red', label='URF Theoretical $T_{1/2}$')

    # 设置Y轴为以10为底的对数坐标
    plt.yscale('log')

    # 设置坐标轴标签和标题
    plt.xlabel('Neutron Number ($N$)', fontsize=14)
    plt.ylabel('Half-life $T_{1/2}$ (s)', fontsize=14)
    plt.title('Evolution of Half-life ($T_{1/2}$) vs Neutron Number ($N$)', fontsize=16)

    # 优化刻度大小
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    # 添加图例
    plt.legend(fontsize=12)

    # 添加网格线以增加可读性
    plt.grid(True, which="both", ls="--", alpha=0.5)

    # 自动调整布局并显示/保存图像
    plt.tight_layout()
    
    # 将图像保存到本地（可选）
    plt.savefig('halflife_vs_N.png', dpi=300)
    print("图像已保存为 'halflife_vs_N.png'")
    
    # 显示图像
    plt.show()

if __name__ == '__main__':
    main()