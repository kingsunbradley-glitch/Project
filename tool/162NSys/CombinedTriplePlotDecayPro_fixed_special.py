'''
python CombinedTriplePlotDecayPro_fixed_special.py dataHsPro.dat dataDsPro.dat \
  --special-left specialDataHsPro.dat \
  --special-right specialDataDsPro.dat
  
'''
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
    'lines.linewidth': 1.4,
    'errorbar.capsize': 3,
    'axes.unicode_minus': False,
})


def load_data(filename: str) -> pd.DataFrame:
    df = pd.read_csv(filename, sep='\t')
    df = df.rename(columns={
        '+': 'T_plus',
        '-': 'T_minus',
        '+.1': 'delta2_plus',
        '-.1': 'delta2_minus',
    })

    required_cols = [
        'N', 'A', 'idx', 'EXP', 'error', 'UNEDF1', 'WS4+RBF',
        'T0.5', 'T_plus', 'T_minus', 'URF',
        'delta2', 'delta2_plus', 'delta2_minus'
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f'{filename} 缺少列: {missing}')

    for col in required_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df.sort_values(['N', 'idx', 'A'], ascending=[True, False, True]).reset_index(drop=True)


def load_optional_data(filename: str | None) -> pd.DataFrame | None:
    if not filename:
        return None
    if not os.path.exists(filename):
        return None
    try:
        df = load_data(filename)
        return df
    except Exception:
        return None


def unique_N_for_line(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
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
        y = 10 ** (np.log10(y0) + 0.14 * (np.log10(y1) - np.log10(y0)))
    else:
        y0, y1 = ax.get_ylim()
        y = y0 + 0.10 * (y1 - y0)
    ax.text(x, y, r'$N=162$', fontsize=11, color='black')


def draw_exp_points_with_errors(ax, x, y, yerr_low=None, yerr_high=None, label='EXP'):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if yerr_low is not None and yerr_high is not None:
        yerr_low = np.asarray(yerr_low, dtype=float)
        yerr_high = np.asarray(yerr_high, dtype=float)
        ax.errorbar(
            x, y, yerr=[yerr_low, yerr_high],
            fmt='none', ecolor='red', elinewidth=1.2, capsize=3,
            zorder=3, label='_nolegend_'
        )
    ax.plot(x, y, 'o', color='red', markersize=4, zorder=4, label=label)


def draw_special_points_with_errors(ax, x, y, yerr_low=None, yerr_high=None):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if yerr_low is not None and yerr_high is not None:
        yerr_low = np.asarray(yerr_low, dtype=float)
        yerr_high = np.asarray(yerr_high, dtype=float)
        ax.errorbar(
            x, y, yerr=[yerr_low, yerr_high],
            fmt='none', ecolor='red', elinewidth=1.2, capsize=3,
            zorder=5, label='_nolegend_'
        )
    ax.plot(
        x, y, linestyle='none', marker='o', markersize=5.2,
        markerfacecolor='white', markeredgecolor='red',
        markeredgewidth=1.3, zorder=6, label='_nolegend_'
    )


def plot_qalpha(ax, df, special_df=None, show_ylabel=True, legend=True):
    unedf = unique_N_for_line(df, 'UNEDF1')
    ws4 = unique_N_for_line(df, 'WS4+RBF')

    ax.plot(unedf['N'], unedf['UNEDF1'], '--', color='#7F7F7F', label='UNEDF1', zorder=2)
    ax.plot(ws4['N'], ws4['WS4+RBF'], '--', color='black', label='WS4+RBF', zorder=2)

    exp = df[(df['idx'] == 1) & df['EXP'].notna()]
    if len(exp):
        draw_exp_points_with_errors(
            ax,
            exp['N'], exp['EXP'] / 1000.0,
            exp['error'] / 1000.0, exp['error'] / 1000.0,
            label='EXP'
        )

    sp = None
    if special_df is not None:
        sp = special_df[(special_df['idx'] == 1) & special_df['EXP'].notna()]
        if len(sp):
            draw_special_points_with_errors(
                ax,
                sp['N'], sp['EXP'] / 1000.0,
                sp['error'] / 1000.0, sp['error'] / 1000.0
            )

    vals = []
    for s in [unedf['UNEDF1'], ws4['WS4+RBF'], exp['EXP'] / 1000.0 if len(exp) else pd.Series(dtype=float)]:
        vals.extend(pd.to_numeric(s, errors='coerce').dropna().tolist())
    if sp is not None and len(sp):
        vals.extend((sp['EXP'] / 1000.0).dropna().tolist())
    if vals:
        ymin, ymax = min(vals), max(vals)
        pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - pad, ymax + pad)

    annotate_shell(ax, mode='linear')
    if show_ylabel:
        ax.set_ylabel(r'$Q_{\alpha}$ (MeV)')
    if legend:
        ax.legend(loc='upper right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def plot_half_life(ax, df, special_df=None, show_ylabel=True, legend=True):
    ax.set_yscale('log')
    exp = df[(df['idx'] == 1) & df['T0.5'].notna() & (df['T0.5'] > 0)]
    if len(exp):
        draw_exp_points_with_errors(
            ax,
            exp['N'], exp['T0.5'],
            exp['T_minus'], exp['T_plus'],
            label='EXP'
        )

    sp = None
    if special_df is not None:
        sp = special_df[(special_df['idx'] == 1) & special_df['T0.5'].notna() & (special_df['T0.5'] > 0)]
        if len(sp):
            draw_special_points_with_errors(
                ax,
                sp['N'], sp['T0.5'],
                sp['T_minus'], sp['T_plus']
            )

    urf = unique_N_for_line(df[df['URF'] > 0], 'URF')
    if len(urf):
        ax.plot(urf['N'], urf['URF'], '--', color='black', label='URF', zorder=2)

    urf_idx1 = df[(df['idx'] == 1) & df['URF'].notna() & (df['URF'] > 0)]
    if len(urf_idx1):
        ax.plot(urf_idx1['N'], urf_idx1['URF'], 'ko', markersize=4, zorder=4, label='_nolegend_')

    vals = pd.concat([
        exp['T0.5'] if len(exp) else pd.Series(dtype=float),
        urf['URF'] if len(urf) else pd.Series(dtype=float),
        sp['T0.5'] if (sp is not None and len(sp)) else pd.Series(dtype=float)
    ]).dropna()
    vals = vals[vals > 0]
    if len(vals):
        ax.set_ylim(10 ** np.floor(np.log10(vals.min()) - 0.45),
                    10 ** np.ceil(np.log10(vals.max()) + 0.45))

    annotate_shell(ax, mode='log')
    if show_ylabel:
        ax.set_ylabel(r'$T_{1/2}$ (s)')
    if legend:
        ax.legend(loc='lower right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))


def plot_delta2(ax, df, special_df=None, show_ylabel=True, legend=True):
    line_df = unique_N_for_line(df, 'delta2')
    if len(line_df):
        ax.plot(line_df['N'], line_df['delta2'], '--', color='black', label=r'$\delta^2$', zorder=2)

    exp = df[(df['idx'] == 1) & df['delta2'].notna()]
    if len(exp):
        draw_exp_points_with_errors(
            ax,
            exp['N'], exp['delta2'],
            exp['delta2_minus'], exp['delta2_plus'],
            label='EXP'
        )

    sp = None
    if special_df is not None:
        sp = special_df[(special_df['idx'] == 1) & special_df['delta2'].notna()]
        if len(sp):
            draw_special_points_with_errors(
                ax,
                sp['N'], sp['delta2'],
                sp['delta2_minus'], sp['delta2_plus']
            )

    vals = pd.concat([
        line_df['delta2'] if len(line_df) else pd.Series(dtype=float),
        (exp['delta2'] - exp['delta2_minus']) if len(exp) else pd.Series(dtype=float),
        (exp['delta2'] + exp['delta2_plus']) if len(exp) else pd.Series(dtype=float),
        (sp['delta2'] - sp['delta2_minus']) if (sp is not None and len(sp)) else pd.Series(dtype=float),
        (sp['delta2'] + sp['delta2_plus']) if (sp is not None and len(sp)) else pd.Series(dtype=float),
    ]).dropna()
    if len(vals):
        ymin, ymax = vals.min(), vals.max()
        pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
        ax.set_ylim(ymin - pad, ymax + pad)

    annotate_shell(ax, mode='linear')
    if show_ylabel:
        ax.set_ylabel(r'$\delta^2$')
    ax.set_xlabel(r'Neutron number $N$')
    if legend:
        ax.legend(loc='upper right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())


def infer_label(path):
    base = os.path.basename(path)
    if 'Hs' in base:
        return 'Hs'
    if 'Ds' in base:
        return 'Ds'
    return os.path.splitext(base)[0]


def infer_special_path(data_path: str) -> str | None:
    directory = os.path.dirname(data_path) or '.'
    label = infer_label(data_path)
    candidate = os.path.join(directory, f'specialData{label}Pro.dat')
    return candidate if os.path.exists(candidate) else None


def make_combined_plot(file_left: str, file_right: str,
                       special_left: str | None = None,
                       special_right: str | None = None,
                       show=True):
    df_left = load_data(file_left)
    df_right = load_data(file_right)
    sp_left = load_optional_data(special_left)
    sp_right = load_optional_data(special_right)

    label_left = infer_label(file_left)
    label_right = infer_label(file_right)

    fig, axes = plt.subplots(
        3, 2, figsize=(10.8, 9.2), sharex='col',
        gridspec_kw={'hspace': 0.0, 'wspace': 0.12}
    )

    plot_qalpha(axes[0, 0], df_left, special_df=sp_left, show_ylabel=True, legend=True)
    plot_half_life(axes[1, 0], df_left, special_df=sp_left, show_ylabel=True, legend=True)
    plot_delta2(axes[2, 0], df_left, special_df=sp_left, show_ylabel=True, legend=True)

    plot_qalpha(axes[0, 1], df_right, special_df=sp_right, show_ylabel=False, legend=True)
    plot_half_life(axes[1, 1], df_right, special_df=sp_right, show_ylabel=False, legend=True)
    plot_delta2(axes[2, 1], df_right, special_df=sp_right, show_ylabel=False, legend=True)

    axes[0, 0].set_title(label_left, fontsize=14, pad=8)
    axes[0, 1].set_title(label_right, fontsize=14, pad=8)

    axes[0, 0].set_xlim(int(df_left['N'].min()) - 1, int(df_left['N'].max()) + 1)
    axes[0, 1].set_xlim(int(df_right['N'].min()) - 1, int(df_right['N'].max()) + 1)

    for i in range(3):
        for j in range(2):
            add_common_style(axes[i, j])

    plt.setp(axes[0, 0].get_xticklabels(), visible=False)
    plt.setp(axes[1, 0].get_xticklabels(), visible=False)
    plt.setp(axes[0, 1].get_xticklabels(), visible=False)
    plt.setp(axes[1, 1].get_xticklabels(), visible=False)

    out_png = f'{label_left}_{label_right}_combined_triple_fixed_special.png'
    out_pdf = f'{label_left}_{label_right}_combined_triple_fixed_special.pdf'
    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f'已保存: {out_png}')
    print(f'已保存: {out_pdf}')

    if show:
        plt.show()
    else:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description='左右拼接 Hs/Ds 三联图，并自动叠加 specialDataHsPro.dat / specialDataDsPro.dat 中的空心特殊点。'
    )
    parser.add_argument('file_left')
    parser.add_argument('file_right')
    parser.add_argument('--special-left', default=None,
                        help='左图特殊数据文件；默认自动寻找 specialDataHsPro.dat 或 specialDataDsPro.dat')
    parser.add_argument('--special-right', default=None,
                        help='右图特殊数据文件；默认自动寻找 specialDataHsPro.dat 或 specialDataDsPro.dat')
    parser.add_argument('--no-show', action='store_true')
    args = parser.parse_args()

    special_left = args.special_left or infer_special_path(args.file_left)
    special_right = args.special_right or infer_special_path(args.file_right)

    make_combined_plot(
        args.file_left, args.file_right,
        special_left=special_left,
        special_right=special_right,
        show=not args.no_show
    )


if __name__ == '__main__':
    main()
