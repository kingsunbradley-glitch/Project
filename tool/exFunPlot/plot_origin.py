#!/usr/bin/env python3
"""Origin 风格：4n、5n 激发函数及实验能量误差。"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator, LogLocator, NullFormatter
from plot_excitation_functions import read_calculation, clip_log_polyline

HERE = Path(__file__).resolve().parent
# 加载真正的 Times New Roman，避免字体回退。
for name in ('times.ttf','timesi.ttf','timesbd.ttf','timesbi.ttf'):
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/'+name)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],
    'mathtext.fontset':'custom','mathtext.rm':'Times New Roman:bold',
    'mathtext.it':'Times New Roman:bold:italic','mathtext.bf':'Times New Roman:bold',
    'font.size':13,'font.weight':'bold','axes.labelsize':16,'axes.labelweight':'bold',
    'axes.linewidth':.8,'xtick.direction':'in','ytick.direction':'in',
    'xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'svg.fonttype':'none'})
data = read_calculation(HERE/'ExFunCalc.csv')
with (HERE/'experimental_data.csv').open(newline='',encoding='utf-8') as f:
    rows = list(csv.DictReader(f))
for row in rows:
    row['E_star_err_minus_MeV'] = '1'
    row['E_star_err_plus_MeV'] = '1.2' if row['kind']=='measurement' else '1'
with (HERE/'experimental_data_with_energy_errors.csv').open('w',newline='',encoding='utf-8') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
fig,ax=plt.subplots(figsize=(190/25.4,120/25.4))
fig.subplots_adjust(left=.12,right=.975,bottom=.17,top=.965)
ax.set(xlim=(30,70),ylim=(1.1e-3,3),yscale='log')
for isotope,color,style,label in [('274Ds','#000000','--',r'$4n$: $^{274}$Ds (calc.)'),('273Ds','#D32F2F','-',r'$5n$: $^{273}$Ds (calc.)')]:
    x,y=clip_log_polyline(data['E*'],data[isotope],ax.get_xlim(),ax.get_ylim())
    ax.plot(x,y,color=color,linestyle=style,lw=1.2,label=label)
# 两个道在 42、45 MeV 的上限及能量误差完全相同：共用黑色符号，能量不偏移。
seen=set()
for row in rows:
    x,y=float(row['E_star_MeV']),float(row['sigma_pb'])
    xe=np.array([[float(row['E_star_err_minus_MeV'])],[float(row['E_star_err_plus_MeV'])]])
    if row['kind']=='upper_limit':
        # 按用户要求，不显示 50 MeV 的 4n 上限；原始实验数据仍保留。
        if row['channel']=='4n' and np.isclose(x,50):
            continue
        key=(x,y)
        if key in seen: continue
        seen.add(key)
        color='#000000'
        ax.errorbar(x,y,xerr=xe,fmt='none',ecolor=color,elinewidth=1.0,capsize=3,zorder=4)
        # 箭长仅用于表示上限，不表示误差。
        end=y/2.1
        ax.annotate('',xy=(x,end),xytext=(x,y),arrowprops=dict(arrowstyle='-|>',color=color,lw=1.0,mutation_scale=9,shrinkA=0,shrinkB=0),zorder=4)
    else:
        ax.errorbar(x,y,xerr=xe,yerr=[[float(row['err_minus_pb'])],[float(row['err_plus_pb'])]],fmt='o',color='#D32F2F',ms=4.5,capsize=3,elinewidth=1.0,zorder=5)
# 用户指定的势垒位置；箭头是能量标记，不代表截面测量。
# VBass = 43.8 MeV，来自实验报告的图 3。
ax.annotate('', xy=(43.8, 1.22e-3), xytext=(43.8, 2.8e-3),
            arrowprops=dict(arrowstyle='-|>', color='black', lw=1.4, mutation_scale=12))
ax.text(44.3, 2.65e-3, r'$V_{\mathrm{BASS}}$', fontsize=16, fontweight='bold', va='center')
ax.text(36.3, 0.18, '4n', fontsize=17, fontweight='bold', color='black')
ax.text(54.2, 0.16, '5n', fontsize=17, fontweight='bold', color='#D32F2F')
ax.set_xlabel(r'Excitation energy $E^{*}$ (MeV)')
ax.set_ylabel(r'Cross section $\sigma$ (pb)')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(LogLocator(base=10,subs=np.arange(2,10)))
ax.yaxis.set_minor_formatter(NullFormatter())
ax.tick_params(which='major',length=5,width=.8,pad=6)
ax.tick_params(which='minor',length=2.5,width=.6)
out=HERE/'figures'; out.mkdir(exist_ok=True)
for ext in ('pdf','svg','png'):
    fig.savefig(out/f'excitation_functions_origin.{ext}',dpi=600)
plt.close(fig)
print(out/'excitation_functions_origin.pdf')
