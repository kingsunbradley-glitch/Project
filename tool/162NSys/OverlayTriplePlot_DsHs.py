import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'STIXGeneral', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'font.size': 11,
    'axes.labelsize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 9,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.linewidth': 1.2,
    'lines.linewidth': 1.5,
    'errorbar.capsize': 3,
    'axes.unicode_minus': False,
})

DS_STYLE = {'marker': 'o', 'linestyle': '--', 'label': 'Ds'}
HS_STYLE = {'marker': 's', 'linestyle': '-.', 'label': 'Hs'}


def load_data(filename: str) -> pd.DataFrame:
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
    sub = df[df[ycol].notna()].copy()
    sub = sub.sort_values(['N', 'idx', 'A'], ascending=[True, False, True])
    sub = sub.drop_duplicates(subset=['N'], keep='first')
    return sub.sort_values('N').reset_index(drop=True)


def add_common_style(ax):
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='major', length=6, width=1.2)
    ax.tick_params(which='minor', length=3, width=1.0)


def annotate_shell(ax, mode='linear'):
    ax.axvline(x=162, color='gray', linestyle=':', linewidth=1.1, zorder=1)
    x = 162.2
    if mode == 'log':
        y0, y1 = ax.get_ylim()
        y = 10 ** (np.log10(y0) + 0.55 * (np.log10(y1) - np.log10(y0)))
    else:
        y0, y1 = ax.get_ylim()
        y = y0 + 0.35 * (y1 - y0)
    ax.text(x, y, r'$N=162$', fontsize=11, color='black')


def draw_exp(ax, df: pd.DataFrame, xcol: str, ycol: str, ylow: str, yhigh: str, marker: str, label: str,
             yscale=1.0):
    exp = df[(df['idx'] == 1) & df[ycol].notna()].copy()
    if exp.empty:
        return exp
    x = exp[xcol]
    y = exp[ycol] * yscale
    low = exp[ylow] * yscale if ylow in exp.columns else None
    high = exp[yhigh] * yscale if yhigh in exp.columns else None
    if low is not None and high is not None:
        ax.errorbar(
            x, y, yerr=[low, high], fmt='none', ecolor='red',
            elinewidth=1.2, capsize=3, zorder=3, label='_nolegend_'
        )
    ax.plot(x, y, linestyle='none', marker=marker, color='red', markersize=4.5,
            zorder=4, label=label)
    return exp


def qalpha_panel(ax, ds_df, hs_df):
    for df, style in [(ds_df, DS_STYLE), (hs_df, HS_STYLE)]:
        unedf = unique_n_for_line(df, 'UNEDF1')
        ws4 = unique_n_for_line(df, 'WS4+RBF')
        ax.plot(unedf['N'], unedf['UNEDF1'], linestyle=style['linestyle'], color='gray', zorder=2)
        ax.plot(ws4['N'], ws4['WS4+RBF'], linestyle=style['linestyle'], color='black', zorder=2)
        draw_exp(ax, df, 'N', 'EXP', 'error', 'error', style['marker'], f"{style['label']} EXP", yscale=1/1000.0)

    vals = []
    for df in [ds_df, hs_df]:
        for col, scale in [('UNEDF1', 1.0), ('WS4+RBF', 1.0), ('EXP', 1/1000.0)]:
            s = pd.to_numeric(df[col], errors='coerce') * scale
            vals.extend(s.dropna().tolist())
    if vals:
        ymin, ymax = min(vals), max(vals)
        pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - pad, ymax + pad)

    ax.set_ylabel(r'$Q_{\alpha}$ (MeV)')
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    annotate_shell(ax, mode='linear')

    handles = [
        Line2D([0], [0], color='red', linestyle='none', marker='o', markersize=5, label='Ds EXP'),
        Line2D([0], [0], color='gray', linestyle='--', marker='o', markersize=4, label='Ds UNEDF1'),
        Line2D([0], [0], color='black', linestyle='--', marker='o', markersize=4, label='Ds WS4+RBF'),
        Line2D([0], [0], color='red', linestyle='none', marker='s', markersize=5, label='Hs EXP'),
        Line2D([0], [0], color='gray', linestyle='-.', marker='s', markersize=4, label='Hs UNEDF1'),
        Line2D([0], [0], color='black', linestyle='-.', marker='s', markersize=4, label='Hs WS4+RBF'),
    ]
    ax.legend(handles=handles, loc='lower left', ncol=2, frameon=False, columnspacing=1.0, handletextpad=0.6)


def half_life_panel(ax, ds_df, hs_df):
    ax.set_yscale('log')
    for df, style in [(ds_df, DS_STYLE), (hs_df, HS_STYLE)]:
        urf = unique_n_for_line(df[df['URF'] > 0], 'URF')
        if not urf.empty:
            ax.plot(urf['N'], urf['URF'], linestyle=style['linestyle'], color='black', zorder=2)
        exp = draw_exp(ax, df, 'N', 'T0.5', 'T_minus', 'T_plus', style['marker'], f"{style['label']} EXP", yscale=1.0)
        urf_idx1 = df[(df['idx'] == 1) & df['URF'].notna() & (df['URF'] > 0)]
        if not urf_idx1.empty:
            ax.plot(urf_idx1['N'], urf_idx1['URF'], linestyle='none', marker=style['marker'],
                    color='black', markersize=4.5, zorder=4, label='_nolegend_')

    vals = []
    for df in [ds_df, hs_df]:
        for col in ['T0.5', 'URF']:
            s = pd.to_numeric(df[col], errors='coerce').dropna()
            vals.extend(s[s > 0].tolist())
    if vals:
        vals = np.array(vals)
        ax.set_ylim(10 ** np.floor(np.log10(vals.min()) - 0.45),
                    10 ** np.ceil(np.log10(vals.max()) + 0.45))

    ax.set_ylabel(r'$T_{1/2}$ (s)')
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    annotate_shell(ax, mode='log')

    handles = [
        Line2D([0], [0], color='black', linestyle='--', marker='o', markersize=4, label='Ds URF'),
        Line2D([0], [0], color='red', linestyle='none', marker='o', markersize=5, label='Ds EXP'),
        Line2D([0], [0], color='black', linestyle='-.', marker='s', markersize=4, label='Hs URF'),
        Line2D([0], [0], color='red', linestyle='none', marker='s', markersize=5, label='Hs EXP'),
    ]
    ax.legend(handles=handles, loc='upper left', ncol=2, frameon=False, columnspacing=1.0, handletextpad=0.6)


def delta2_panel(ax, ds_df, hs_df):
    for df, style in [(ds_df, DS_STYLE), (hs_df, HS_STYLE)]:
        line_df = unique_n_for_line(df, 'delta2')
        if not line_df.empty:
            ax.plot(line_df['N'], line_df['delta2'], linestyle=style['linestyle'], color='black', zorder=2)
        draw_exp(ax, df, 'N', 'delta2', 'delta2_minus', 'delta2_plus', style['marker'], f"{style['label']} idx=1", yscale=1.0)

    vals = []
    for df in [ds_df, hs_df]:
        vals.extend(pd.to_numeric(df['delta2'], errors='coerce').dropna().tolist())
        exp = df[(df['idx'] == 1) & df['delta2'].notna()]
        if not exp.empty:
            vals.extend((exp['delta2'] - exp['delta2_minus']).dropna().tolist())
            vals.extend((exp['delta2'] + exp['delta2_plus']).dropna().tolist())
    if vals:
        ymin, ymax = min(vals), max(vals)
        pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - pad, ymax + pad)

    ax.set_ylabel(r'$\delta^2$')
    ax.set_xlabel(r'Neutron number $N$')
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    annotate_shell(ax, mode='linear')

    handles = [
        Line2D([0], [0], color='black', linestyle='--', marker='o', markersize=4, label=r'Ds $\delta^2$'),
        Line2D([0], [0], color='red', linestyle='none', marker='o', markersize=5, label='Ds EXP'),
        Line2D([0], [0], color='black', linestyle='-.', marker='s', markersize=4, label=r'Hs $\delta^2$'),
        Line2D([0], [0], color='red', linestyle='none', marker='s', markersize=5, label='Hs EXP'),
    ]
    ax.legend(handles=handles, loc='upper left', ncol=2, frameon=False, columnspacing=1.0, handletextpad=0.6)


def infer_output_name(ds_file: str, hs_file: str) -> str:
    a = os.path.splitext(os.path.basename(ds_file))[0]
    b = os.path.splitext(os.path.basename(hs_file))[0]
    return f'{a}_{b}_overlay_triple'


def main():
    parser = argparse.ArgumentParser(description='在一张三联图中叠加 Ds 与 Hs 数据。Ds 用圆点，Hs 用方块。')
    parser.add_argument('ds_file', nargs='?', default='dataDsPro.dat')
    parser.add_argument('hs_file', nargs='?', default='dataHsPro.dat')
    parser.add_argument('--no-show', action='store_true')
    args = parser.parse_args()

    ds_df = load_data(args.ds_file)
    hs_df = load_data(args.hs_file)

    fig, axes = plt.subplots(3, 1, figsize=(8.0, 12.0), sharex=True,
                             gridspec_kw={'hspace': 0.0})

    qalpha_panel(axes[0], ds_df, hs_df)
    half_life_panel(axes[1], ds_df, hs_df)
    delta2_panel(axes[2], ds_df, hs_df)

    xmin = int(min(ds_df['N'].min(), hs_df['N'].min())) - 1
    xmax = int(max(ds_df['N'].max(), hs_df['N'].max())) + 1
    axes[0].set_xlim(xmin, xmax)

    for ax in axes:
        add_common_style(ax)
    plt.setp(axes[0].get_xticklabels(), visible=False)
    plt.setp(axes[1].get_xticklabels(), visible=False)

    fig.align_ylabels(axes)
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
