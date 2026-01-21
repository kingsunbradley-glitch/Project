import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# --- 1. 设置期刊风格的全局参数 ---
def set_style():
    # 优先使用 Arial 字体 (Nature/Science 标准)
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    
    # 基础字号设置
    # 为了让格子里的数字能放下，基础字号不能太大
    mpl.rcParams['font.size'] = 10 
    mpl.rcParams['axes.labelsize'] = 11
    mpl.rcParams['axes.titlesize'] = 12
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    
    # 刻度线朝内
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    
    # 高清输出设置
    mpl.rcParams['savefig.dpi'] = 300
    mpl.rcParams['savefig.bbox'] = 'tight'

def plot_readable_heatmap(input_file):
    set_style()
    
    # 读取数据
    try:
        df = pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        print(f"找不到文件 {input_file}")
        return

    # 创建输出目录
    output_dir = 'heatmaps_readable'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    d1_values = sorted(df['D1'].unique())

    # 使用高对比度的色谱
    # 'viridis' 或 'Blues_r' (深蓝到白) 也是不错的选择
    # 这里继续用 'mako_r' 因为它很高级，且深色背景下白色文字很清晰
    cmap_name = 'mako_r' 

    for d1 in d1_values:
        subset = df[df['D1'] == d1]
        
        # 数据透视
        heatmap_data = subset.pivot(index='Q2', columns='Q1', values='N')
        
        # 排序: Y轴降序 (大数值在上面), X轴升序
        heatmap_data.sort_index(ascending=False, inplace=True) 
        heatmap_data.sort_index(axis=1, ascending=True, inplace=True) 

        # --- 动态调整画布大小 ---
        # 根据格子的数量自动调整图片长宽，防止格子太扁或太挤
        n_rows, n_cols = heatmap_data.shape
        # 基础尺寸 + 每个格子给一点空间
        fig_width = 4 + n_cols * 0.4 
        fig_height = 3 + n_rows * 0.3
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # --- 核心绘图 ---
        sns.heatmap(
            heatmap_data,
            ax=ax,
            annot=True,           # 【开启】显示数字
            fmt=".0f",            # 【格式】只显示整数，不带小数点，节省空间
            annot_kws={"size": 9},# 【字号】格子里的数字字号 (根据需要微调)
            cmap=cmap_name,       
            square=True,          # 正方形格子
            linewidths=0.5,       # 白色分割线
            linecolor='white',
            cbar_kws={
                'label': 'Average N', 
                'shrink': 0.6     # 色条不要太长
            },
            mask=heatmap_data.isnull() # 隐藏空数据
        )
        
        # --- 标签优化 ---
        ax.set_title(f'D1 = {d1:.3f}', pad=15, fontweight='bold')
        ax.set_xlabel('Q1 Parameter', labelpad=10)
        ax.set_ylabel('Q2 Parameter', labelpad=10)

        # 【关键】X轴标签旋转，防止重叠 (你发来的图片里这里重叠了)
        plt.xticks(rotation=45, ha='right') # 45度倾斜，对齐到右边
        plt.yticks(rotation=0)              # Y轴保持水平

        # 保存
        filename = os.path.join(output_dir, f'heatmap_D1_{d1:.3f}.png')
        plt.savefig(filename)
        plt.close(fig)
        print(f"已保存: {filename}")

if __name__ == "__main__":
    plot_readable_heatmap('processed_result.txt')