'''
python OverlayTriplePlot_DsHs_special.py dataDsPro.dat dataHsPro.dat \
  --special-ds specialDataDsPro.dat \
  --special-hs specialDataHsPro.dat
'''
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D


# =========================
# 绘图常量区（建议只改这里）
# =========================
# 字体相关：总字号、坐标轴标题字号、刻度字号、图例字号
BASE_FONT_SIZE = 16
AXES_LABEL_SIZE = 19
XTICK_LABEL_SIZE = 16
YTICK_LABEL_SIZE = 16
LEGEND_FONT_SIZE = 14

# 线条相关：坐标轴线宽、曲线线宽、误差棒帽宽度/线宽
AXES_LINEWIDTH = 1.6
CURVE_LINEWIDTH = 2.0
ERRORBAR_CAPSIZE = 4
ERRORBAR_LINEWIDTH = 1.4

# 画布与子图间距：控制整张图大小与边界留白
FIGSIZE = (14.5, 12.0)
SUBPLOT_LEFT = 0.08
SUBPLOT_RIGHT = 0.92
SUBPLOT_BOTTOM = 0.08
SUBPLOT_TOP = 0.95
SUBPLOT_HSPACE = 0.0
SUBPLOT_WSPACE = 0.0

# 刻度相关：主/副刻度间隔、长度、线宽
X_MAJOR_TICK_STEP = 2
X_MINOR_SUBDIV = 2
MAJOR_TICK_LENGTH = 7
MAJOR_TICK_WIDTH = 1.6
MINOR_TICK_LENGTH = 4
MINOR_TICK_WIDTH = 1.2

# N=162 壳闭合线与标注文本位置
SHELL_X = 162.2
SHELL_LINEWIDTH = 1.3
N162_FONT_SIZE = 17
N162_LINEAR_FRAC = {'down': 0.12, 'middle': 0.35, 'up': 0.82}
N162_LOG_FRAC = {'down': 0.20, 'middle': 0.65, 'up': 0.82}

# 子图角标 (a)-(f) 位置与字号
PANEL_TAG_X_LEFT = 0.04
PANEL_TAG_X_RIGHT = 0.96
PANEL_TAG_Y = 0.95
PANEL_TAG_FONT_SIZE = 24

# 数据点样式：实验点、special 空心点、星标增量等
EXP_MARKER_SIZE = 8.0
SPECIAL_MARKER_SIZE = 10.5
HIGHLIGHT_STAR_DELTA = 1.0
HOLLOW_EDGE_WIDTH = 1.6
HOLLOW_EDGE_WIDTH_LEGEND = 1.4
URF_IDX_MARKER_SIZE = 6.0
LEGEND_URF_MARKER_SIZE = 5.0
LEGEND_EXP_MARKER_SIZE = 6.0
LEGEND_SPECIAL_MARKER_SIZE = 10.5
LEGEND_STAR_EXP_MARKER_SIZE = 10.5
LEGEND_STAR_SPECIAL_MARKER_SIZE = 10.5

# 左右列标题字号与上边距
TITLE_FONT_SIZE = 22
TITLE_PAD = 10

# x 轴左右扩展边距（避免数据贴边）
X_MARGIN = 0.35

# 底部 x 轴标题（英文）
LEFT_COLUMN_XLABEL = r'Neutron number of Ds isotopes '
RIGHT_COLUMN_XLABEL = r'Neutron number of Hs isotopes '

# y 轴上下扩展边距（用于图内“留白”）
QALPHA_YPAD_LOW = 0.08
QALPHA_YPAD_HIGH = 0.22
HALFLIFE_LOG_PAD_LOW = 0.20
HALFLIFE_LOG_PAD_HIGH = 0.30
DELTA2_LOG_PAD_LOW = 0.20
DELTA2_LOG_PAD_HIGH = 0.25

# 图例布局参数
LEGEND_COL_SPACING = 0.9
LEGEND_HANDLE_TEXT_PAD = 0.5

# 左右大面板之间的分割线样式
SEPARATOR_LINEWIDTH = 1.4
SEPARATOR_ALPHA = 0.85


plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': BASE_FONT_SIZE,
    'font.weight': 'bold',
    'axes.labelsize': AXES_LABEL_SIZE,
    'axes.labelweight': 'bold',
    'axes.titleweight': 'bold',
    'xtick.labelsize': XTICK_LABEL_SIZE,
    'ytick.labelsize': YTICK_LABEL_SIZE,
    'legend.fontsize': LEGEND_FONT_SIZE,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.linewidth': AXES_LINEWIDTH,
    'lines.linewidth': CURVE_LINEWIDTH,
    'errorbar.capsize': ERRORBAR_CAPSIZE,
    'axes.unicode_minus': False,
})


# 固定 Ds / Hs 的点型、线型与需要高亮为五角星的 N 值
DS_STYLE = {'marker': 'o', 'linestyle': '--', 'label': r'', 'highlight_n': 163}
HS_STYLE = {'marker': 's', 'linestyle': '-.', 'label': r'', 'highlight_n': 161}


def load_data(filename: str) -> pd.DataFrame:
    """读取并清洗输入数据；统一列名并做必要列校验。"""
    df = pd.read_csv(filename, sep='\t')
    df = df.rename(columns={
        '+': 'T_plus',
        '-': 'T_minus',
        '+.1': 'delta2_plus',
        '-.1': 'delta2_minus',
    })
    need = [
        'N', 'A', 'idx', 'EXP', 'error', 'UNEDF1', 'WS4+RBF',
        'T0.5', 'T_plus', 'T_minus', 'URF',
        'delta2', 'delta2_plus', 'delta2_minus'
    ]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f'{filename} 缺少列: {missing}')
    for c in need:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df.sort_values(['N', 'idx', 'A'], ascending=[True, False, True]).reset_index(drop=True)


def unique_n_for_line(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    """每个 N 只保留一个点用于连线，避免同一 N 多次连接。"""
    sub = df[df[ycol].notna()].copy()
    sub = sub.sort_values(['N', 'idx', 'A'], ascending=[True, False, True])
    sub = sub.drop_duplicates(subset=['N'], keep='first')
    return sub.sort_values('N').reset_index(drop=True)


def set_axis_side(ax, side: str):
    """设置 y 轴在左侧或右侧显示，并隐藏内侧 spine。"""
    if side == 'left':
        ax.yaxis.set_label_position('left')
        ax.yaxis.tick_left()
        ax.tick_params(axis='y', which='both', left=True, labelleft=True, right=False, labelright=False)
        ax.spines['right'].set_visible(False)
    else:
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
        ax.tick_params(axis='y', which='both', left=False, labelleft=False, right=True, labelright=True)
        ax.spines['left'].set_visible(False)


def add_common_style(ax):
    """给每个子图应用统一刻度样式。"""
    ax.xaxis.set_major_locator(ticker.MultipleLocator(X_MAJOR_TICK_STEP))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(X_MINOR_SUBDIV))
    ax.tick_params(which='major', length=MAJOR_TICK_LENGTH, width=MAJOR_TICK_WIDTH)
    ax.tick_params(which='minor', length=MINOR_TICK_LENGTH, width=MINOR_TICK_WIDTH)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight('bold')


def annotate_shell(ax, mode='linear', text_pos='middle'):
    """绘制 N=162 竖线与标注文本，支持上线性/对数坐标。"""
    ax.axvline(x=162, color='gray', linestyle=':', linewidth=SHELL_LINEWIDTH, zorder=1)
    x = SHELL_X
    if mode == 'log':
        y0, y1 = ax.get_ylim()
        frac = N162_LOG_FRAC.get(text_pos, N162_LOG_FRAC['middle'])
        y = 10 ** (np.log10(y0) + frac * (np.log10(y1) - np.log10(y0)))
    else:
        y0, y1 = ax.get_ylim()
        frac = N162_LINEAR_FRAC.get(text_pos, N162_LINEAR_FRAC['middle'])
        y = y0 + frac * (y1 - y0)
    ax.text(x, y, r'$N=162$', fontsize=N162_FONT_SIZE, fontweight='bold', color='black')


def add_panel_tag(ax, tag: str, side='left'):
    """添加子图角标 (a)-(f)，side 控制放左上还是右上。"""
    if side == 'left':
        x, ha = PANEL_TAG_X_LEFT, 'left'
    else:
        x, ha = PANEL_TAG_X_RIGHT, 'right'
    ax.text(
        x, PANEL_TAG_Y, tag,
        transform=ax.transAxes,
        fontsize=PANEL_TAG_FONT_SIZE,
        fontweight='bold',
        va='top',
        ha=ha
    )


def draw_exp(ax, df: pd.DataFrame, ycol: str, ylow: str, yhigh: str, marker: str,
             yscale=1.0, color='red', markersize=EXP_MARKER_SIZE, highlight_n=None):
    """绘制实验点（含误差棒）；若命中 highlight_n 则改为更大的五角星。"""
    exp = df[(df['idx'] == 1) & df[ycol].notna()].copy()
    if exp.empty:
        return exp
    x = exp['N'].to_numpy()
    y = (exp[ycol] * yscale).to_numpy()
    low = (exp[ylow] * yscale).to_numpy() if ylow in exp.columns else None
    high = (exp[yhigh] * yscale).to_numpy() if yhigh in exp.columns else None
    if low is not None and high is not None:
        ax.errorbar(
            x, y, yerr=[low, high], fmt='none', ecolor=color,
            elinewidth=ERRORBAR_LINEWIDTH, capsize=ERRORBAR_CAPSIZE, zorder=3, label='_nolegend_'
        )
    if highlight_n is None:
        ax.plot(
            x, y, linestyle='none', marker=marker, color=color, markersize=markersize,
            zorder=4, label='_nolegend_'
        )
    else:
        # 带红星的点不再画同位置的红圈/红方框
        x_arr = np.asarray(x)
        mask_star = np.isclose(x_arr, float(highlight_n))
        if np.any(~mask_star):
            ax.plot(
                x_arr[~mask_star], y[~mask_star], linestyle='none', marker=marker, color=color,
                markersize=markersize, zorder=4, label='_nolegend_'
            )
        if np.any(mask_star):
            ax.plot(
                x_arr[mask_star], y[mask_star], linestyle='none', marker='*', color=color,
                markersize=markersize + HIGHLIGHT_STAR_DELTA, zorder=5, label='_nolegend_'
            )
    return exp


def draw_special(ax, df: pd.DataFrame | None, ycol: str, ylow: str, yhigh: str, marker: str,
                 yscale=1.0, edgecolor='red', markersize=SPECIAL_MARKER_SIZE, highlight_n=None):
    """绘制 special 点（空心）；若命中 highlight_n 同样改为空心五角星。"""
    if df is None or df.empty or ycol not in df.columns:
        return
    sub = df[df[ycol].notna()].copy()
    if sub.empty:
        return
    x = sub['N'].to_numpy()
    y = (sub[ycol] * yscale).to_numpy()
    low = (sub[ylow] * yscale).to_numpy() if ylow in sub.columns else None
    high = (sub[yhigh] * yscale).to_numpy() if yhigh in sub.columns else None
    if low is not None and high is not None:
        ax.errorbar(
            x, y, yerr=[low, high], fmt='none', ecolor=edgecolor,
            elinewidth=ERRORBAR_LINEWIDTH, capsize=ERRORBAR_CAPSIZE, zorder=5, label='_nolegend_'
        )
    if highlight_n is None:
        ax.plot(
            x, y, linestyle='none', marker=marker, markersize=markersize,
            markerfacecolor='white', markeredgecolor=edgecolor, markeredgewidth=HOLLOW_EDGE_WIDTH,
            zorder=6, label='_nolegend_'
        )
    else:
        # 带空心红星的点不再画同位置空心圈/空心方框
        x_arr = np.asarray(x)
        mask_star = np.isclose(x_arr, float(highlight_n))
        if np.any(~mask_star):
            ax.plot(
                x_arr[~mask_star], y[~mask_star], linestyle='none', marker=marker, markersize=markersize,
                markerfacecolor='white', markeredgecolor=edgecolor, markeredgewidth=HOLLOW_EDGE_WIDTH,
                zorder=6, label='_nolegend_'
            )
        if np.any(mask_star):
            ax.plot(
                x_arr[mask_star], y[mask_star], linestyle='none', marker='*',
                markersize=markersize + HIGHLIGHT_STAR_DELTA,
                markerfacecolor='white', markeredgecolor=edgecolor, markeredgewidth=HOLLOW_EDGE_WIDTH,
                zorder=7, label='_nolegend_'
            )


def qalpha_panel(ax, df, style, special=None, side='left', show_legend=True, ylim=None,
                 legend_loc=None, shell_pos='middle', legend_ncol=2):
    """上排面板：Qα。"""
    unedf = unique_n_for_line(df, 'UNEDF1')
    ws4 = unique_n_for_line(df, 'WS4+RBF')

    ax.plot(unedf['N'], unedf['UNEDF1'], linestyle=style['linestyle'], color='gray', zorder=2)
    ax.plot(ws4['N'], ws4['WS4+RBF'], linestyle=style['linestyle'], color='black', zorder=2)
    draw_exp(
        ax, df, 'EXP', 'error', 'error', style['marker'], yscale=1/1000.0,
        highlight_n=style.get('highlight_n')
    )
    draw_special(
        ax, special, 'EXP', 'error', 'error', style['marker'], yscale=1/1000.0,
        highlight_n=style.get('highlight_n')
    )

    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_ylabel(r'$Q_{\alpha}$ (MeV)', fontweight='bold')
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    set_axis_side(ax, side)
    annotate_shell(ax, mode='linear', text_pos=shell_pos)

    if show_legend:
        # 图例点型与线型单独定义，避免受实际数据中星标高亮影响
        if legend_loc is None:
            legend_loc = 'lower left' if side == 'left' else 'lower right'
        handles = [
            # a/d 图要求：先 EXP，再 This work 与 This work(i.s.)，最后理论线
            Line2D([0], [0], color='red', linestyle='none', marker=style['marker'],
                   markersize=LEGEND_EXP_MARKER_SIZE, label='EXP'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_EXP_MARKER_SIZE,
                   label='This work'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_SPECIAL_MARKER_SIZE,
                   markerfacecolor='white', markeredgecolor='red',
                   markeredgewidth=HOLLOW_EDGE_WIDTH_LEGEND, label='This work(i.s.)'),
            Line2D([0], [0], color='gray', linestyle=style['linestyle'], label='UNEDF1'),
            Line2D([0], [0], color='black', linestyle=style['linestyle'], label='WS4+RBF'),
        ]
        ax.legend(
            handles=handles,
            loc=legend_loc,
            ncol=legend_ncol,
            frameon=False,
            columnspacing=LEGEND_COL_SPACING,
            handletextpad=LEGEND_HANDLE_TEXT_PAD,
        )


def half_life_panel(ax, df, style, special=None, side='left', show_legend=True, ylim=None,
                    legend_loc=None, shell_pos='middle'):
    """中排面板：半衰期 T1/2（对数坐标）。"""
    ax.set_yscale('log')

    urf = unique_n_for_line(df[df['URF'] > 0], 'URF')
    if not urf.empty:
        ax.plot(urf['N'], urf['URF'], linestyle=style['linestyle'], color='black', zorder=2)

    draw_exp(
        ax, df, 'T0.5', 'T_minus', 'T_plus', style['marker'], yscale=1.0,
        highlight_n=style.get('highlight_n')
    )

    urf_idx1 = df[(df['idx'] == 1) & df['URF'].notna() & (df['URF'] > 0)]
    if not urf_idx1.empty:
        ax.plot(
            urf_idx1['N'], urf_idx1['URF'], linestyle='none', marker=style['marker'],
            color='black', markersize=URF_IDX_MARKER_SIZE, zorder=4, label='_nolegend_'
        )

    draw_special(
        ax, special, 'T0.5', 'T_minus', 'T_plus', style['marker'], yscale=1.0,
        highlight_n=style.get('highlight_n')
    )

    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_ylabel(r'$T_{1/2}$ (s)', fontweight='bold')
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    set_axis_side(ax, side)
    annotate_shell(ax, mode='log', text_pos=shell_pos)

    if show_legend:
        if legend_loc is None:
            legend_loc = 'upper left' if side == 'left' else 'upper right'
        handles = [
            Line2D([0], [0], color='black', linestyle=style['linestyle'], marker=style['marker'],
                   markersize=LEGEND_URF_MARKER_SIZE, label='URF'),
            Line2D([0], [0], color='red', linestyle='none', marker=style['marker'],
                   markersize=LEGEND_EXP_MARKER_SIZE, label='EXP'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_EXP_MARKER_SIZE,
                   label='This work'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_SPECIAL_MARKER_SIZE,
                   markerfacecolor='white', markeredgecolor='red',
                   markeredgewidth=HOLLOW_EDGE_WIDTH_LEGEND, label='This work(i.s.)'),
        ]
        ax.legend(
            handles=handles,
            loc=legend_loc,
            ncol=1,
            frameon=False,
            columnspacing=LEGEND_COL_SPACING,
            handletextpad=LEGEND_HANDLE_TEXT_PAD,
        )


def delta2_panel(ax, df, style, special=None, side='left', show_legend=True, ylim=None,
                 legend_loc=None, shell_pos='middle'):
    """下排面板：delta^2。"""
    ax.set_yscale('log')

    # 对数坐标下仅绘制正值
    df_pos = df[df['delta2'].notna() & (df['delta2'] > 0)].copy()
    special_pos = None
    if special is not None and not special.empty:
        special_pos = special[special['delta2'].notna() & (special['delta2'] > 0)].copy()

    line_df = unique_n_for_line(df_pos, 'delta2')
    if not line_df.empty:
        ax.plot(line_df['N'], line_df['delta2'], linestyle=style['linestyle'], color='black', zorder=2)
    draw_exp(
        ax, df_pos, 'delta2', 'delta2_minus', 'delta2_plus', style['marker'], yscale=1.0,
        highlight_n=style.get('highlight_n')
    )
    draw_special(
        ax, special_pos, 'delta2', 'delta2_minus', 'delta2_plus', style['marker'], yscale=1.0,
        highlight_n=style.get('highlight_n')
    )

    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_ylabel(r'$\delta^2$', fontweight='bold')
    ax.set_xlabel(r'Neutron Number $N$', fontweight='bold')
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    set_axis_side(ax, side)
    annotate_shell(ax, mode='log', text_pos=shell_pos)

    if show_legend:
        if legend_loc is None:
            legend_loc = 'upper left' if side == 'left' else 'upper right'
        handles = [
            Line2D([0], [0], color='black', linestyle=style['linestyle'], marker=style['marker'],
                   markersize=LEGEND_URF_MARKER_SIZE, label=r'$\delta^2$'),
            Line2D([0], [0], color='red', linestyle='none', marker=style['marker'],
                   markersize=LEGEND_EXP_MARKER_SIZE, label='EXP'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_EXP_MARKER_SIZE,
                   label='This work'),
            Line2D([0], [0], color='red', linestyle='none', marker='*',
                   markersize=LEGEND_STAR_SPECIAL_MARKER_SIZE,
                   markerfacecolor='white', markeredgecolor='red',
                   markeredgewidth=HOLLOW_EDGE_WIDTH_LEGEND, label='This work(i.s.)'),
        ]
        ax.legend(
            handles=handles,
            loc=legend_loc,
            ncol=1,
            frameon=False,
            columnspacing=LEGEND_COL_SPACING,
            handletextpad=LEGEND_HANDLE_TEXT_PAD,
        )


def infer_output_name(ds_file: str, hs_file: str) -> str:
    """根据输入文件名生成输出文件前缀。"""
    a = os.path.splitext(os.path.basename(ds_file))[0]
    b = os.path.splitext(os.path.basename(hs_file))[0]
    return f'{a}_{b}_final_3x2'


def get_x_limits(df: pd.DataFrame, special: pd.DataFrame | None, margin=X_MARGIN):
    """计算单列 x 轴范围，并在两侧额外留出 margin。"""
    mins = [df['N'].min()]
    maxs = [df['N'].max()]
    if special is not None and not special.empty:
        mins.append(special['N'].min())
        maxs.append(special['N'].max())
    xmin = float(min(mins))
    xmax = float(max(maxs))
    if xmax <= xmin:
        return xmin - 0.5, xmax + 0.5
    return xmin - margin, xmax + margin


def calc_qalpha_ylim(ds_df, hs_df, ds_special=None, hs_special=None):
    """计算 Qα 行共享 y 轴范围（左右两图统一）。"""
    vals = []
    for subdf in [ds_df, hs_df, ds_special, hs_special]:
        if subdf is None:
            continue
        for col, scale in [('UNEDF1', 1.0), ('WS4+RBF', 1.0), ('EXP', 1/1000.0)]:
            if col in subdf.columns:
                s = pd.to_numeric(subdf[col], errors='coerce') * scale
                vals.extend(s.dropna().tolist())
    if not vals:
        return None
    ymin, ymax = min(vals), max(vals)
    span = (ymax - ymin) if ymax > ymin else 1.0
    return ymin - QALPHA_YPAD_LOW * span, ymax + QALPHA_YPAD_HIGH * span


def calc_half_life_ylim(ds_df, hs_df, ds_special=None, hs_special=None):
    """计算半衰期行共享 y 轴范围（对数坐标）。"""
    vals = []
    for subdf in [ds_df, hs_df, ds_special, hs_special]:
        if subdf is None:
            continue
        for col in ['T0.5', 'URF']:
            if col not in subdf.columns:
                continue
            s = pd.to_numeric(subdf[col], errors='coerce').dropna()
            vals.extend(s[s > 0].tolist())
    if not vals:
        return None
    vals = np.array(vals)
    return (
        10 ** (np.log10(vals.min()) - HALFLIFE_LOG_PAD_LOW),
        10 ** (np.log10(vals.max()) + HALFLIFE_LOG_PAD_HIGH)
    )


def calc_delta2_ylim(ds_df, hs_df, ds_special=None, hs_special=None):
    """计算 delta^2 行共享 y 轴范围（左右两图统一）。"""
    vals = []
    for subdf in [ds_df, hs_df, ds_special, hs_special]:
        if subdf is None:
            continue
        d = pd.to_numeric(subdf['delta2'], errors='coerce')
        vals.extend(d[d > 0].tolist())
        if all(c in subdf.columns for c in ['delta2', 'delta2_minus', 'delta2_plus']):
            exp = subdf[subdf['delta2'].notna()].copy()
            if not exp.empty:
                low = pd.to_numeric(exp['delta2'] - exp['delta2_minus'], errors='coerce')
                high = pd.to_numeric(exp['delta2'] + exp['delta2_plus'], errors='coerce')
                vals.extend(low[low > 0].dropna().tolist())
                vals.extend(high[high > 0].dropna().tolist())
    if not vals:
        return None
    vals = np.array(vals)
    return (
        10 ** (np.log10(vals.min()) - DELTA2_LOG_PAD_LOW),
        10 ** (np.log10(vals.max()) + DELTA2_LOG_PAD_HIGH)
    )


def main():
    """主流程：读取数据 -> 计算范围 -> 绘图 -> 导出 PNG/PDF。"""
    parser = argparse.ArgumentParser(description='生成 3x2 对比图：左列 ^{273}Ds，右列 ^{269}Hs（支持 special 空心标记）。')
    parser.add_argument('ds_file', nargs='?', default='dataDsPro.dat')
    parser.add_argument('hs_file', nargs='?', default='dataHsPro.dat')
    parser.add_argument('--special-ds', default='specialDataDsPro.dat')
    parser.add_argument('--special-hs', default='specialDataHsPro.dat')
    parser.add_argument('--no-show', action='store_true')
    args = parser.parse_args()

    # 1) 读取主数据与 special 数据
    ds_df = load_data(args.ds_file)
    hs_df = load_data(args.hs_file)

    ds_special = load_data(args.special_ds) if args.special_ds and os.path.exists(args.special_ds) else None
    hs_special = load_data(args.special_hs) if args.special_hs and os.path.exists(args.special_hs) else None

    # 2) 先计算每一行共享的 y 轴范围（左/右统一）
    q_ylim = calc_qalpha_ylim(ds_df, hs_df, ds_special, hs_special)
    hl_ylim = calc_half_life_ylim(ds_df, hs_df, ds_special, hs_special)
    d2_ylim = calc_delta2_ylim(ds_df, hs_df, ds_special, hs_special)

    fig, axes = plt.subplots(
        3, 2,
        figsize=FIGSIZE,
        sharex='col',
        sharey='row',
        gridspec_kw={'hspace': SUBPLOT_HSPACE, 'wspace': SUBPLOT_WSPACE}
    )

    # 3) 绘制左列（Ds）与右列（Hs）
    qalpha_panel(
        axes[0, 0], ds_df, DS_STYLE, ds_special, side='left', show_legend=True,
        ylim=q_ylim, shell_pos='middle', legend_ncol=1
    )
    half_life_panel(
        axes[1, 0], ds_df, DS_STYLE, ds_special, side='left', show_legend=True, ylim=hl_ylim, shell_pos='middle'
    )
    delta2_panel(
        axes[2, 0], ds_df, DS_STYLE, ds_special, side='left', show_legend=True,
        ylim=d2_ylim, shell_pos='up', legend_loc='lower left'
    )

    qalpha_panel(
        axes[0, 1], hs_df, HS_STYLE, hs_special, side='right', show_legend=True, ylim=q_ylim,
        legend_loc='upper right', shell_pos='up', legend_ncol=1
    )
    half_life_panel(
        axes[1, 1], hs_df, HS_STYLE, hs_special, side='right', show_legend=True, ylim=hl_ylim,
        legend_loc='lower right', shell_pos='down'
    )
    delta2_panel(
        axes[2, 1], hs_df, HS_STYLE, hs_special, side='right', show_legend=True,
        ylim=d2_ylim, shell_pos='up', legend_loc='lower right'
    )

    # 4) 不显示顶部 273Ds / 269Hs 标题（按当前要求）

    # 5) 左右列分别设置 x 轴范围
    ds_xmin, ds_xmax = get_x_limits(ds_df, ds_special, margin=X_MARGIN)
    hs_xmin, hs_xmax = get_x_limits(hs_df, hs_special, margin=X_MARGIN)
    for i in range(3):
        axes[i, 0].set_xlim(ds_xmin, ds_xmax)
        axes[i, 1].set_xlim(hs_xmin, hs_xmax)

    for i in range(3):
        for j in range(2):
            add_common_style(axes[i, j])

    # 底部 x 轴标题使用常量字符串
    axes[2, 0].set_xlabel(LEFT_COLUMN_XLABEL, fontweight='bold')
    axes[2, 1].set_xlabel(RIGHT_COLUMN_XLABEL, fontweight='bold')

    plt.setp(axes[0, 0].get_xticklabels(), visible=False)
    plt.setp(axes[1, 0].get_xticklabels(), visible=False)
    plt.setp(axes[0, 1].get_xticklabels(), visible=False)
    plt.setp(axes[1, 1].get_xticklabels(), visible=False)

    # 6) 调整整图边距
    fig.subplots_adjust(
        left=SUBPLOT_LEFT,
        right=SUBPLOT_RIGHT,
        bottom=SUBPLOT_BOTTOM,
        top=SUBPLOT_TOP,
        hspace=SUBPLOT_HSPACE,
        wspace=SUBPLOT_WSPACE
    )

    # 7) 添加子图角标：(a)(b)(c) 在左列，(d)(e)(f) 在右列
    tags = [('(a)', '(d)'), ('(b)', '(e)'), ('(c)', '(f)')]
    for i in range(3):
        add_panel_tag(axes[i, 0], tags[i][0], side='right')
        add_panel_tag(axes[i, 1], tags[i][1], side='left')

    # 8) 添加左右大面板分割线
    left_bbox = axes[0, 0].get_position()
    right_bbox = axes[0, 1].get_position()
    x_sep = 0.5 * (left_bbox.x1 + right_bbox.x0)
    y_top = max(axes[0, 0].get_position().y1, axes[0, 1].get_position().y1)
    y_bot = min(axes[2, 0].get_position().y0, axes[2, 1].get_position().y0)
    fig.add_artist(
        Line2D(
            [x_sep, x_sep], [y_bot, y_top],
            transform=fig.transFigure,
            color='gray',
            linewidth=SEPARATOR_LINEWIDTH,
            alpha=SEPARATOR_ALPHA,
            zorder=10
        )
    )

    # 9) 导出图片与 PDF
    out_base = infer_output_name(args.ds_file, args.hs_file)
    out_png = out_base + '.png'
    out_pdf = out_base + '.pdf'
    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f'已保存: {out_png}')
    print(f'已保存: {out_pdf}')

    if args.no_show:
        plt.close(fig)
    else:
        plt.show()


if __name__ == '__main__':
    main()
