import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
from pathlib import Path

# ==============================================================================
# --- 1. 全局字体与样式设置 (用于发表级绘图) ---
# ==============================================================================

# ==============================================================================
# --- 1. 全局字体与样式设置 (用于发表级绘图) ---
# ==============================================================================
# 【用户配置】：统一字体与子图标注（放在前面方便统一调整）
PANEL_LABELS = ['(a)', '(b)', '(c)']
FONT_SIZE_MAIN = 18
FONT_SIZE_COUNT = 16
FONT_SIZE_AXIS = 18
FONT_SIZE_TICK = 14
FONT_WEIGHT = 'bold'

# 将 Mac 的中文字体 (Songti SC 或 Arial Unicode MS) 加入到字体列表中排在前面
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'  # 公式字体设置保持不变，继续渲染上标和希腊字母
plt.rcParams['axes.unicode_minus'] = False # 【关键】防止引入中文字体后，坐标轴的负号 '-' 变成小方块
plt.rcParams['font.weight'] = FONT_WEIGHT
plt.rcParams['axes.labelweight'] = FONT_WEIGHT
plt.rcParams['axes.titleweight'] = FONT_WEIGHT

# ==============================================================================
# --- 2. 数据读取与预处理 ---
# ==============================================================================
# 【用户配置】：输入文件名（相对脚本所在目录，避免受运行目录影响）
file_name = Path(__file__).with_name('HsData.txt')

try:
    decay_values = []
    evr_values = []

    # 逐行解析：兼容“只填 EVR，Decay 留空”的行
    with open(file_name, 'r', encoding='utf-8') as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue  # 跳过表头

            if not line.strip():
                continue

            tokens = line.split()
            if not tokens:
                continue

            # 以空白开头时，视为只有 EVR 数据（如 "\t\t8.92"）
            if line[:1].isspace():
                evr_values.append(float(tokens[-1]))
            else:
                decay_values.append(float(tokens[0]))
                if len(tokens) > 1:
                    evr_values.append(float(tokens[1]))

    decay_data = np.array(decay_values)
    evr_data = np.array(evr_values)
    
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
bins = np.linspace(energy_min, energy_max, 100)

# ==============================================================================
# --- 5. 绘图配置列表 (核心控制区) ---
# ==============================================================================
# 列表中的每个元组代表一张图的配置：(数据变量, 颜色代码, 显示标签)
plots_config = [
    # 第一张图：所有数据，黑色
    (all_data,   "black",   "Total (EVR + Decay)"),
    # 第二张图：EVR数据，红色 (#D81525 是 RGB 十六进制颜色)
    (evr_data,  "black" , "EVR"),
    # 第三张图：Decay数据，蓝色，标签用了 LaTeX 语法 ($...$) 显示上标，r表示原始字符串，避免Latex的\转义问题
    (decay_data, "#D81525", r"From $^{273}$Ds $\alpha$ decay ")
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
    # --- A. 绘制直方图 ---
    if len(data) > 0:
        if i == 0:
            # 针对第一张图 (i==0)：绘制 EVR 和 Decay 的堆叠直方图 (Stacked Histogram)
            # 传入数据列表 [evr_data, decay_data]，对应颜色 ["black", "#D81525"]
            # stacked=True 实现“加和”效果，柱子会上下堆叠，总高度即为总 counts
            n, bins_out, patches = ax.hist([evr_data, decay_data], bins=bins, 
                                           color=["black", "#D81525"], stacked=True, 
                                           alpha=0.9, edgecolor='white', linewidth=0.5)
        else:
            # 针对第二和第三张图：保持普通的单色直方图绘制
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
    ax.text(0.05, 0.85, label, transform=ax.transAxes, fontsize=FONT_SIZE_MAIN, fontweight=FONT_WEIGHT)
    
    # 2. 添加事件统计 (如 "13 events")
    # 0.05, 0.75: Y位置比上面那个低 0.1，形成两行字的效果
    count_text = f"{len(data)} events"
    ax.text(0.05, 0.75, count_text, transform=ax.transAxes, fontsize=FONT_SIZE_COUNT, fontweight=FONT_WEIGHT)

    # 3. 添加右上角子图编号 (a)(b)(c)
    ax.text(0.95, 0.90, PANEL_LABELS[i], transform=ax.transAxes,
            ha='right', va='top', fontsize=FONT_SIZE_MAIN, fontweight=FONT_WEIGHT)
    
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
    ax.tick_params(direction='in', top=True, right=True, which='both', labelsize=FONT_SIZE_TICK)
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight(FONT_WEIGHT)
    
    # --- E. 隐藏不需要的 X 轴标签 ---
    # 如果不是最后一张图 (i < 2)，就把 X 轴的数字隐藏掉，避免中间重叠
    if i < 2:
        ax.tick_params(labelbottom=False)
    else:
        # 最后一张图，加上 X 轴名称
        ax.set_xlabel(r"$E_{\alpha} $/ MeV", fontsize=FONT_SIZE_AXIS, fontweight=FONT_WEIGHT)

# ==============================================================================
# --- 7. 全局标签与保存 ---
# ==============================================================================
# fig.text: 在整张画布上写字，不属于任何一个子图
# 0.10, 0.5: 【用户配置】X, Y 坐标。
# - 如果觉得字离图太近，就把 0.10 改成 0.05 (往左移)
# - 如果觉得字被切掉了，就结合 subplots_adjust 的 left 参数一起改
fig.text(0.10, 0.5, 'Counts / 20 keV', va='center', rotation='vertical',
         fontsize=FONT_SIZE_AXIS, fontweight=FONT_WEIGHT)

# subplots_adjust: 调整图表边缘留白
# left=0.15: 左边留 15% 空白 (给 Y 轴标题腾位置)
# hspace=0:再次确保垂直间距为 0
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.08)

# savefig: 保存图片
# dpi=300: 分辨率 (300是出版级标准)
# bbox_inches='tight': 自动裁剪掉多余的白边
plt.savefig('HsData_WithEvents.pdf', dpi=300, bbox_inches='tight')
plt.savefig('HsData_WithEvents.png', dpi=600, bbox_inches='tight')
# show: 在屏幕上显示出来
plt.show()
