import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ==========================================
# 1. 全局格式设置 (PRL / Origin 风格)
# ==========================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',     
    'font.size': 12,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 10,
    'xtick.direction': 'in',        
    'ytick.direction': 'in',
    'xtick.top': True,              
    'ytick.right': True,            
    'axes.linewidth': 1.2,          
    'lines.linewidth': 1.5,         
    'errorbar.capsize': 3           
})

# ==========================================
# 2. 读取数据
# ==========================================
try:
    # 🔴【重要修改】：更换了文件名，并使用 \t (Tab) 作为分隔符，防止中间的空缺列导致数据错位
    df_hl = pd.read_csv('halfLifeDataPro.dat', sep='\t')
    df_hl = df_hl.sort_values('N')
except FileNotFoundError:
    print("Error: 找不到 halfLifeDataPro.dat 文件。请确保文件名和路径正确。")
    exit()

try:
    df_qa = pd.read_csv('QaData.dat', sep='\s+')
    df_qa = df_qa.sort_values('N')
except FileNotFoundError:
    print("Error: 找不到 QaData.dat 文件。")
    exit()

# ==========================================
# 3. 创建画布与子图布局
# ==========================================
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7), sharex=True)
fig.subplots_adjust(hspace=0.0) 

# ==========================================
# 4. 绘制上方子图：Q_alpha vs N
# ==========================================
valid_exp = df_qa.dropna(subset=['EXP', 'error'])
exp_N = valid_exp['N']
exp_Q = valid_exp['EXP'] / 1000.0
exp_err = valid_exp['error'] / 1000.0

ax1.errorbar(exp_N, exp_Q, yerr=exp_err, fmt='ro', markersize=4, label='EXP', zorder=3)

# 绘制实心五角星 isomer
ax1.errorbar(163, 11.020, yerr=0.020, fmt='*', color="#FAF6F6", 
             markersize=8, markeredgecolor='black', markeredgewidth=0.8, 
             label=r"EXP ($^{273}\mathrm{Ds}$ i.s.)", zorder=4)

# 2. 理论线 1 (UNEDF)：黑色虚线，与下半图的 URF 呼应
ax1.plot(df_qa['N'], df_qa['UNEDF'], '--', color='#7F7F7F', label='UNEDF1', zorder=2)

# 3. 理论线 2 (WS4+RBF)：深灰色虚线，作为辅助对比
ax1.plot(df_qa['N'], df_qa['WS4+RBF'], '--', color='black', label='WS4+RBF', zorder=2)

ax1.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
ax1.text(163, 9.5, r'$N=162$', fontsize=12, color='black')

ax1.set_ylabel(r'$Q_{\alpha}$ (MeV)')
ax1.legend(loc='upper right', frameon=False) 

# ==========================================
# 5. 绘制下方子图：T_1/2 vs N
# ==========================================
ax2.set_yscale('log')

# 过滤出有实验半衰期的数据点进行绘制，避免 NaN 报错
valid_hl = df_hl.dropna(subset=['T0.5', '+', '-'])
N_hl_exp = valid_hl['N']
T05 = valid_hl['T0.5']
err_plus = valid_hl['+']
err_minus = valid_hl['-']

# 绘制实验半衰期红球
ax2.errorbar(N_hl_exp, T05, yerr=[err_minus, err_plus], fmt='ro', markersize=4, label='EXP', zorder=3)

# 绘制 URF 理论虚线 (贯穿所有数据)
ax2.plot(df_hl['N'], df_hl['URF'], '--', color='black', label='URF', zorder=2)

# 🌟【新增核心逻辑】：筛选 idx == 1 的行，并在理论线上打出黑点
mask_idx1 = df_hl['idx'] == 1
ax2.plot(df_hl.loc[mask_idx1, 'N'], df_hl.loc[mask_idx1, 'URF'], 'ko', markersize=4, zorder=4)

ax2.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
ax2.text(163, 1e3, r'$N=162$', fontsize=12, color='black')

ax2.set_xlabel(r'Neutron number $N$')
ax2.set_ylabel(r'$T_{1/2}$ (s)')
ax2.legend(loc='lower right', frameon=False)

ax2.set_ylim(1e-2, 10**6.3) 
ax2.set_yticks([10**i for i in range(-2, 7)])
ax2.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)))


# ==========================================
# 6. 细节修饰与导出
# ==========================================
min_N = min(df_hl['N'].min(), df_qa['N'].min())
max_N = max(df_hl['N'].max(), df_qa['N'].max())
ax1.set_xlim(min_N - 1, max_N + 1)

for ax in [ax1, ax2]:
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='major', length=6, width=1.2)
    ax.tick_params(which='minor', length=3, width=1.0)
    
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
plt.setp(ax1.get_xticklabels(), visible=False) 

# plt.savefig('Decay_Properties_Pro.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Decay_Properties_Pro.png', dpi=600, bbox_inches='tight')

plt.show()