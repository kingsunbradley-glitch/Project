import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
import pandas as pd

# ==============================================================================
# --- 1. 全局字体与样式设置 (用于发表级绘图) ---
# ==============================================================================
# 这里的设置是为了让图片的字体和 Physical Review 等期刊的风格一致 (Times New Roman 风格)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'  # 公式字体设置 (如上标、希腊字母)

# ==============================================================================
# --- 2. 数据读取与预处理 ---
# ==============================================================================
# 【用户配置】：输入文件名
file_name = 'HsData.txt'

try:
    # pd.read_csv: 读取文本文件
    # sep='\t': 指定“制表符”为分隔符。如果你的数据用空格分隔，改成 sep='\s+'
    df = pd.read_csv(file_name, sep='\t')
    
    # dropna(): 去掉空值 (NaN)。因为 Decay 列可能比 EVR 列短，必须去掉空才能画图
    # .values: 转换成 numpy 数组，方便后续计算
    decay_data = df['Decay'].dropna().values
    evr_data = df['EVR'].dropna().values
    
    # np.concatenate: 将两列数据拼接到一起，用来画第一张 "All" 的总图
    all_data = np.concatenate([decay_data, evr_data])
    
except Exception as e:
    # 错误处理：如果读不到文件，打印错误并生成空数组，防止程序直接崩溃
    print(f"读取文件出错: {e}")
    decay_data, evr_data, all_data = np.array([]), np.array([]), np.array([])

# ==============================================================================
# --- 3. 画布与布局配置 ---
# ==============================================================================
# figsize=(8, 10): 画布大小，宽8英寸，高10英寸。修改这里可以改变图片的整体长宽比
fig = plt.figure(figsize=(8, 10)) 

# GridSpec(3, 1): 将画布切分为 3行 1列
# hspace=0: 【重要】设置子图之间的垂直间距为0，实现“上下紧贴”的效果
gs = GridSpec(3, 1, figure=fig, hspace=0)

# ==============================================================================
# --- 4. 坐标轴范围与 Bins 设置 ---
# ==============================================================================
# 自动计算 X 轴范围，确保所有数据都能显示，并且左右留有余地
if len(all_data) > 0:
    data_min, data_max = np.min(all_data), np.max(all_data)
    # padding: 在最大最小值两边各留 20% 的空白，避免柱子贴着边框
    padding = (data_max - data_min) * 0.2
    energy_min, energy_max = data_min - padding, data_max + padding
else:
    # 默认范围 (如果没有数据)
    energy_min, energy_max = 8.0, 10.0

# np.linspace: 创建直方图的“桶” (bins)。
# 400: 【用户配置】桶的数量。
# - 改大 (如 800): 柱子变细，分辨率更高。
# - 改小 (如 50): 柱子变粗，看起来更平滑。
bins = np.linspace(energy_min, energy_max, 400)

# ==============================================================================
# --- 5. 绘图配置列表 (核心控制区) ---
# ==============================================================================
# 列表中的每个元组代表一张图的配置：(数据变量, 颜色代码, 显示标签)
plots_config = [
    # 第一张图：所有数据，黑色
    (all_data,   "black",   "All"),
    # 第二张图：EVR数据，红色 (#D81525 是 RGB 十六进制颜色)
    (evr_data,  "#126782" , "EVR"),
    # 第三张图：Decay数据，蓝色，标签用了 LaTeX 语法 ($...$) 显示上标，r表示原始字符串，避免Latex的\转义问题
    (decay_data, "#D81525", r"$^{273}$Ds $\alpha$ Decay")
]

axes = [] # 用来存储生成的三个子图对象

# ==============================================================================
# --- 6. 循环绘图逻辑 (Main Loop) ---
# ==============================================================================
# enumerate: 同时获取索引(i) 和 配置内容(data, color, label)
for i, (data, color, label) in enumerate(plots_config):
    
    # gs[i, 0]: 在第 i 行生成一个子图
    ax = fig.add_subplot(gs[i, 0])
    axes.append(ax)
    
    # --- A. 绘制直方图 ---
    if len(data) > 0:
        # ax.hist: 绘制直方图函数
        # alpha=0.9: 透明度 (0-1)，0.9 表示稍有透明
        # edgecolor='white': 【用户配置】柱子边缘颜色。'white' 为白衬线，'none' 为无边框
        # linewidth=0.5: 衬线的粗细
        # 返回值: n(每个桶的计数), bins_out(桶的边界), patches(所有的柱子对象)
        n, bins_out, patches = ax.hist(data, bins=bins, color=color, alpha=0.9, 
                                      label=label, edgecolor='white', linewidth=0.5)
        
        # --- B. 特殊高亮逻辑 (仅针对第3张图 i==2) ---
        if i == 2:
            target_value = 8.915 # 【用户配置】需要高亮的能量值
            
            # zip: 将“柱子对象”和“桶的左右边界”打包在一起遍历
            for patch, left_edge, right_edge in zip(patches, bins_out[:-1], bins_out[1:]):
                # 判断：如果 target_value 落在当前这个柱子的范围内
                if left_edge <= target_value <= right_edge:
                    # set_hatch: 设置填充样式
                    # '///': 斜线。你可以改成 'xxx' (交叉), '...' (点), '|||' (竖线)
                    patch.set_hatch('///')

    # --- C. 添加文字标注 ---
    # transform=ax.transAxes: 【重要】使用相对坐标系。
    # (0,0)是左下角，(1,1)是右上角。这样无论数据是多少，文字位置都固定在框内。
    
    # 1. 添加左上角的大标签 (如 "All", "EVR")
    # 0.05, 0.85: X位置(靠左), Y位置(靠上)
    ax.text(0.05, 0.85, label, transform=ax.transAxes, fontsize=16)
    
    # 2. 添加事件统计 (如 "13 events")
    # 0.05, 0.75: Y位置比上面那个低 0.1，形成两行字的效果
    count_text = f"{len(data)} events"
    ax.text(0.05, 0.75, count_text, transform=ax.transAxes, fontsize=14)
    
    # --- D. 坐标轴美化 ---
    ax.set_xlim(energy_min, energy_max) # 限制 X 轴显示范围
    
    # MaxNLocator: 强制 Y 轴只显示整数刻度 (integer=True)，最多显示 4 个刻度 (nbins=4)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=4))
    # AutoMinorLocator: 自动添加 X 轴的小刻度
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    
    # 显示顶部和右侧的边框 (spines)，构成一个封闭的方框
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    
    # tick_params: 刻度线设置
    # direction='in': 刻度线朝里
    # top=True, right=True: 上边和右边也要有刻度线
    ax.tick_params(direction='in', top=True, right=True, which='both', labelsize=12)
    
    # --- E. 隐藏不需要的 X 轴标签 ---
    # 如果不是最后一张图 (i < 2)，就把 X 轴的数字隐藏掉，避免中间重叠
    if i < 2:
        ax.tick_params(labelbottom=False)
    else:
        # 最后一张图，加上 X 轴名称
        ax.set_xlabel("Energy / MeV", fontsize=14)

# ==============================================================================
# --- 7. 全局标签与保存 ---
# ==============================================================================
# fig.text: 在整张画布上写字，不属于任何一个子图
# 0.10, 0.5: 【用户配置】X, Y 坐标。
# - 如果觉得字离图太近，就把 0.10 改成 0.05 (往左移)
# - 如果觉得字被切掉了，就结合 subplots_adjust 的 left 参数一起改
fig.text(0.10, 0.5, 'Event Counts / 5 keV', va='center', rotation='vertical', fontsize=16)

# subplots_adjust: 调整图表边缘留白
# left=0.15: 左边留 15% 空白 (给 Y 轴标题腾位置)
# hspace=0:再次确保垂直间距为 0
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.08)

# savefig: 保存图片
# dpi=300: 分辨率 (300是出版级标准)
# bbox_inches='tight': 自动裁剪掉多余的白边
plt.savefig('HsData_WithEvents.png', dpi=300, bbox_inches='tight')

# show: 在屏幕上显示出来
plt.show()