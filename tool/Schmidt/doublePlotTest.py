import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import norm
from matplotlib.gridspec import GridSpec

# =================================================================
# === USER CONFIGURATION (用户配置区) ===
# =================================================================

# --- 1. 画布设置 ---
FIG_WIDTH = 14      # 画布宽度
FIG_HEIGHT = 7      # 画布高度
TITLE = r"Rf $\alpha$ vs $\beta$" # 图标题
NUCLIDE = r"$^{261}$Rf"    # 核素标注

# --- 2. 能量轴 (左图) 设置 ---
E_XMIN = 8.0        # 能量 X轴 最小值
E_XMAX = 9.2       # 能量 X轴 最大值
E_YMAX = 25.0       # Y轴 最大值 (两图共用)
E_BINS = 50         # 能量轴切分成多少个柱子

# --- 3. 时间轴 (右图) 设置 (对数坐标) ---
T_XMIN = 1e2        # 时间 X轴 最小值 (100 ms)
T_XMAX = 2e5        # 时间 X轴 最大值 (200,000 ms)
T_BINS = 50         # 时间轴切分成多少个柱子

# --- 4. 分组拟合开关 ---
# 列表顺序对应: [Present, Dubna, RIKEN, GSI, Chem]
# True = 显示该组的拟合曲线, False = 不显示
SHOW_FIT_CURVES = [True, True, True, True, False] 

# --- 5. 样式与颜色 ---
# 建议使用透明度 (alpha) 来防止大柱子完全遮挡小柱子
BAR_ALPHA = 0.7     # 柱子透明度 (0.0-1.0)
labels = ['Sghigh', 'Sglow', 'Rfhigh', 'Rflow', 'Chem']
# 颜色配置 (对应上面5组)
colors = ["#A714AC", "#2D10D4", "#DB7611", "#4C10AD", '#E8F1F2']
# 拟合线颜色 (通常深一点好看)
fit_colors = ["#6A0DAD", '#0B3E4D', '#104E63', '#5D8C99', "#9E1F1F"]

# =================================================================

# -----------------------------------------------------------------
# --- 设置全局字体 ---
# -----------------------------------------------------------------
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['STIXGeneral', 'Times New Roman', 'DejaVu Serif', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'

# -----------------------------------------------------------------
# --- 1. 辅助函数 ---
# -----------------------------------------------------------------
def theoretical_curve(t, lambda_est):
    return lambda_est * t * np.exp(-lambda_est * t)

def paste_data(raw_string):
    """将从Excel粘贴的多行文本转换为numpy数组"""
    return np.fromstring(raw_string, sep=' ')

# -----------------------------------------------------------------
# --- 2. 数据定义 (在此处粘贴数据) ---
# -----------------------------------------------------------------

# >>> Group 1: Present <<<
hs_en_present = paste_data("""
9.03
8.96
8.93
8.93
8.93
8.91
8.902
8.9
8.89
8.88
8.88
8.87
8.86
8.86
8.86
8.86
8.86
8.85
8.85
8.85
8.85
8.85
8.85
8.85
8.84
8.84
8.84
8.84
8.83
8.82
8.82
8.81
8.81
8.8
8.79
8.79
8.77
8.77

""")
hs_hist_present = paste_data("""
28.8
12.8
19.1
9.39
7


17.078
8.8
9.79
13.33
17.5
21.7
0.6
1.64
18.51
25.8
17.65
5.6
19.1
17.81
0.22
4.78

8.29
4.51
5.64
12
2.65
6.03
27.3
2.8

18.27
15.1
56.18
18.82
6.8

""") * 1000.0

# >>> Group 2: Dubna <<<
hs_en_dubna = paste_data("""
8.76
8.76
8.76
8.75
8.74
8.74
8.74
8.73
8.72
8.72
8.72
8.72
8.72
8.72
8.71
8.71
8.71
8.7
8.7
8.7
8.7
8.7
8.695
8.69
8.69
8.69
8.69
8.69
8.69
8.69
8.69
8.68
8.68
8.68
8.68
8.68
8.68
8.68
8.68
8.671
8.67
8.67
8.67
8.66
8.66
8.66
8.66
8.66
8.65
8.65
8.65
8.64
8.63
8.63
8.63
8.62
8.62
8.61
8.6
8.6
8.6
8.59
8.54
8.52
8.509
8.397
8.35

""")
hs_hist_dubna = paste_data("""
1.4
0.14
6.6
27.31

3.5
9.04
51.32
6.82
13.62
3.82
23.22
33.14
11.07
23
8.7

79.9
4.8
10.82
5.74
7.37


3.809
7.75
4.409
1.43
33.54
32.5
1.4
7.61
2.48
9.325
85.6
68.32
2.73
3.3
8.41
11.663
43.57
14.77
3.34
25.7
13.8
6.339
84.9
22.08
6.75
6.87
24.886
11.9
3.33
1.2
52
95.33
10.9
3.35
7.7




48.9
58.19
9.076
48.86

""") * 1000.0

# >>> Group 3: RIKEN <<<
hs_en_riken = paste_data("""8.57
8.56
8.53
8.52
8.52
8.52
8.51
8.5
8.5
8.47
""")* 1000.0
hs_hist_riken = paste_data("""
14.561
7.52
1.7
3.14
6.453
2.4
0.4
1.29
4.68
4.56
4.35
2.26
1.06
2.97
0.58
8.3
3.7
0.8
3.98
13.599
0.36
0.875
0.36
2.25
7.09
7.92
4.44
2
2.07
4.38
7.24
1.55
0.199
3.73
5.297
7
3.04
6.69
1.8
1.216
3.86
7.32
6.42
0.191
2.011
0.748
0.215
2.8
0.3923
3.032
14.5
3.124
15.88
31.9
0.8
1.94
4.7
35
4.89
0.846
2.365
2.54
 """)* 1000.0

# >>> Group 4: GSI <<<
hs_en_gsi = paste_data("""8.44
8.43
8.41
8.41
8.39
8.37
8.37
8.37
8.37
8.37
8.36
8.35
8.35
8.35
8.34
8.33
8.33
8.33
8.33
8.33
8.31
8.31
8.31
8.31
8.31
8.3
8.3
8.3
8.29
8.29
8.29
8.29
8.29
8.29
8.29
8.28
8.28
8.28
8.28
8.28
8.27
8.27
8.26
8.26
8.25
8.25
8.25
8.25
8.25
8.25
8.24
8.24
8.24
8.23
8.23
8.23
8.22
8.22
8.21
8.21
8.2
8.2
8.2
8.19
8.19
8.185
8.18
8.18
8.18
8.17
8.16
8.16
8.15
8.14
8.14
8.12
8.04
""")* 1000.0
hs_hist_gsi = paste_data("""7.22




93.65
164.33
18.29
8.35
31


11.63
21.72
64.51
108.89
127.29
40.65
36.71

35.67
13.33
4.78
49.81

49.74
24.98

70.77
16.8
83.53
16.95
104.91
32.1
58.69

12.77
43.53
86.31
67.1
95.7
93.7
5
13.83
28.51
47.56
17.81
70.37
46.2
17.07
17.67
18.76

59.05
20.12


37.56
74.91

48.39
7.5
100.6
85.47


44.48
10.72
18.08
37.86
17.79




32.62
""")* 1000.0

# >>> Group 5: Chem <<<
hs_en_chem = np.array([])
hs_hist_chem = np.array([])

# 打包数据以便循环
all_hist_data_energy = [hs_en_present, hs_en_dubna, hs_en_riken, hs_en_gsi, hs_en_chem]
all_hist_data_time = [hs_hist_present, hs_hist_dubna, hs_hist_riken, hs_hist_gsi, hs_hist_chem]

# -----------------------------------------------------------------
# --- 3. 绘图主逻辑 ---
# -----------------------------------------------------------------
fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
gs = GridSpec(1, 2, figure=fig, wspace=0, width_ratios=[1, 1])

ax1 = fig.add_subplot(gs[0, 0])  # 能量
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)  # 时间

fig.suptitle(TITLE, fontsize=18, y=0.98)

# ==========================================
# 左图：能量 (Energy) - 叠加绘制
# ==========================================
bins_energy = np.linspace(E_XMIN, E_XMAX, E_BINS)

# 循环绘制每一组数据，而不是一次性绘制
# 这样可以控制 "小的在大的前面" (通过 zorder 或 alpha)
# 这里使用 alpha 透明度叠加法，这是最通用的不遮挡方法
for i, data in enumerate(all_hist_data_energy):
    if len(data) == 0: continue
    
    # 绘制直方图
    ax1.hist(data,
             bins=bins_energy,
             stacked=False, # 关闭堆叠，改为覆盖
             color=colors[i],
             alpha=BAR_ALPHA, # 透明度
             edgecolor='black',
             linewidth=0.8,
             label=labels[i])

ax1.set_xlim(E_XMIN, E_XMAX)
ax1.set_ylim(bottom=0, top=E_YMAX)
ax1.set_xlabel("Energy / MeV", fontsize=14)
ax1.set_ylabel("Event Counts", fontsize=14)

ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
ax1.text(0.05, 0.9, NUCLIDE, transform=ax1.transAxes, fontsize=20)

# 强制显示图例 (因为我们是分开画的)
# ax1.legend(loc='upper left', frameon=False) 


# ==========================================
# 右图：时间 (Time) - 叠加绘制 + 独立拟合
# ==========================================
ax2.set_xscale('log')
ax2.set_xlim(T_XMIN, T_XMAX)
bins_time = np.logspace(np.log10(T_XMIN), np.log10(T_XMAX), num=T_BINS)

# 1. 先画所有的直方图
for i, data in enumerate(all_hist_data_time):
    if len(data) == 0: continue
    
    # 绘制直方图
    # zorder 越高越靠前，这里可以让数据量少的稍微靠前一点，或者直接依赖透明度
    n_counts, _, _ = ax2.hist(data,
                              bins=bins_time,
                              stacked=False,
                              color=colors[i],
                              alpha=BAR_ALPHA,
                              edgecolor='black',
                              linewidth=0.8,
                              label=labels[i])

    # 2. 绘制各自的拟合曲线
    if SHOW_FIT_CURVES[i] and len(data) > 1:
        # 计算 Mean Lifetime
        t_bar_ms = np.mean(data)
        lambda_est = 1.0 / t_bar_ms
        print(f"[{labels[i]}] Mean Life: {t_bar_ms:.2f} ms")

        # 生成曲线
        t_axis = np.logspace(np.log10(T_XMIN), np.log10(T_XMAX), 500)
        y_curve_unscaled = theoretical_curve(t_axis, lambda_est)

        # 归一化曲线高度 (匹配当前这一组数据的最高柱子)
        max_count = np.max(n_counts)
        if max_count == 0: max_count = 1
        
        if y_curve_unscaled.max() > 0:
            scale_factor = max_count / y_curve_unscaled.max()
            y_scaled_curve = y_curve_unscaled * scale_factor
            
            # 绘制曲线
            ax2.plot(t_axis, y_scaled_curve, '-', 
                     linewidth=2.5, 
                     color=fit_colors[i], # 使用专用拟合色
                     label=f'{labels[i]} Fit')

ax2.set_xlabel("Time (ms)", fontsize=14)

# 合并图例：因为分开画会导致图例重复或不全，这里统一收集
handles, labs = ax2.get_legend_handles_labels()
# 去重 (如果需要)
by_label = dict(zip(labs, handles))
ax2.legend(by_label.values(), by_label.keys(), frameon=False, fontsize=11, loc='upper right')

# -----------------------------------------------------------------
# --- 4. 边框与美化 ---
# -----------------------------------------------------------------
ax1.spines['top'].set_visible(True)
ax1.spines['right'].set_visible(False)
ax1.tick_params(direction='in', top=True, right=False, which='both', labelsize=12)

ax2.spines['top'].set_visible(True)
ax2.spines['left'].set_linestyle('--')
ax2.spines['left'].set_linewidth(1.5)
ax2.spines['left'].set_color('gray')
ax2.tick_params(direction='in', top=True, right=True, which='both', labelsize=12)
ax2.tick_params(axis='y', labelleft=False, left=False)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()