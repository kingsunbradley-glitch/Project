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
    df_hl = pd.read_csv('halfLifeData.dat', sep='\s+')
    df_hl = df_hl.sort_values('N')
except FileNotFoundError:
    print("Error: 找不到 halfLifeData.dat 文件。")
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
# 紧密贴合：设置垂直间距为 0
fig.subplots_adjust(hspace=0.0) 

# ==========================================
# 4. 绘制上方子图：Q_alpha vs N
# ==========================================
valid_exp = df_qa.dropna(subset=['EXP', 'error'])
exp_N = valid_exp['N']
exp_Q = valid_exp['EXP'] / 1000.0
exp_err = valid_exp['error'] / 1000.0

# 红球大小 markersize=4
ax1.errorbar(exp_N, exp_Q, yerr=exp_err, fmt='ro', markersize=4, label='EXP', zorder=3)
ax1.errorbar(163, 11.020, yerr=0.020, fmt='r*', markersize=6, label='EXP (isomer)', zorder=4,markerfacecolor='none',)

ax1.plot(df_qa['N'], df_qa['UNEDF'], '--', color='#1f77b4', label='UNEDF', zorder=2)
ax1.plot(df_qa['N'], df_qa['WS4+RBF'], '--', color='#2ca02c', label='WS4+RBF', zorder=2)

ax1.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
ax1.text(163, 9.5, r'$N=162$', fontsize=18, color='black')

ax1.set_ylabel(r'$Q_{\alpha}$ (MeV)')
ax1.legend(loc='upper right', frameon=False) 

# ==========================================
# 5. 绘制下方子图：T_1/2 vs N
# ==========================================
ax2.set_yscale('log')

N_hl = df_hl['N']
T05 = df_hl['T0.5']
err_plus = df_hl['+']
err_minus = df_hl['-']

# 下半部分红球大小 markersize=4
ax2.errorbar(N_hl, T05, yerr=[err_minus, err_plus], fmt='ro', markersize=4, label='EXP', zorder=3)
ax2.plot(N_hl, df_hl['URF'], '--', color='black', label='URF', zorder=2)

ax2.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
ax2.text(163, 1e3, r'$N=162$', fontsize=18, color='black')

ax2.set_xlabel(r'Neutron number $N$')
ax2.set_ylabel(r'$T_{1/2}$ (s)')
ax2.legend(loc='lower right', frameon=False)

# 🔴【核心修改 1】：实现 10 的 6.3 次方
ax2.set_ylim(1e-2, 10**6.3) 

# 🔴【核心修改 2】：强制规定 Y 轴显示每一个 10 的整数次幂 (-2 到 6)
# 生成 [10^-2, 10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6]
ax2.set_yticks([10**i for i in range(-2, 7)])

# 🔴 强制补齐对数坐标下的次要刻度 (minor ticks)，让图表看起来更专业
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

# plt.savefig('Decay_Properties_Updated.pdf', dpi=300, bbox_inches='tight')
plt.savefig('Decay_Properties_Updated.png', dpi=600, bbox_inches='tight')

plt.show()