import os
import sys
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
})


def load_data(filename: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(filename, sep='\t')
    except FileNotFoundError:
        print(f"Error: 找不到文件 {filename}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: 读取文件 {filename} 失败: {e}")
        sys.exit(1)

    rename_map = {
        '+': 'T_plus',
        '-': 'T_minus',
        '+.1': 'delta2_plus',
        '-.1': 'delta2_minus',
    }
    df = df.rename(columns=rename_map)

    required_cols = [
        'N', 'A', 'idx', 'EXP', 'error', 'UNEDF1', 'WS4+RBF',
        'T0.5', 'T_plus', 'T_minus', 'URF',
        'delta2', 'delta2_plus', 'delta2_minus'
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f'Error: {filename} 缺少以下列:')
        print(', '.join(missing))
        sys.exit(1)

    for col in required_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df.sort_values(['N', 'idx', 'A'], ascending=[True, False, True]).reset_index(drop=True)


def unique_N_for_line(df: pd.DataFrame, ycol: str) -> pd.DataFrame:
    """
    线图去掉同一 N 的重复点，避免在同一 x 处出现竖直“尾巴”。
    优先保留 idx==1 的记录；若同为 idx，再保留先出现的那一行。
    """
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
    if yerr_low is not None and yerr_high is not None:
        ax.errorbar(
            x, y, yerr=[yerr_low, yerr_high],
            fmt='none', ecolor='red', elinewidth=1.2, capsize=3,
            zorder=3, label='_nolegend_'
        )
    ax.plot(x, y, 'o', color='red', markersize=4, zorder=4, label=label)


def plot_qalpha(ax, df, show_ylabel=True, legend=True):
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

    vals = []
    for s in [unedf['UNEDF1'], ws4['WS4+RBF'], exp['EXP'] / 1000.0 if len(exp) else pd.Series(dtype=float)]:
        vals.extend(pd.to_numeric(s, errors='coerce').dropna().tolist())
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


def plot_half_life(ax, df, show_ylabel=True, legend=True):
    ax.set_yscale('log')
    exp = df[(df['idx'] == 1) & df['T0.5'].notna() & (df['T0.5'] > 0)]
    if len(exp):
        draw_exp_points_with_errors(
            ax,
            exp['N'], exp['T0.5'],
            exp['T_minus'], exp['T_plus'],
            label='EXP'
        )

    urf = unique_N_for_line(df[df['URF'] > 0], 'URF')
    if len(urf):
        ax.plot(urf['N'], urf['URF'], '--', color='black', label='URF', zorder=2)

    urf_idx1 = df[(df['idx'] == 1) & df['URF'].notna() & (df['URF'] > 0)]
    if len(urf_idx1):
        ax.plot(urf_idx1['N'], urf_idx1['URF'], 'ko', markersize=4, zorder=4, label='_nolegend_')

    vals = pd.concat([
        exp['T0.5'] if len(exp) else pd.Series(dtype=float),
        urf['URF'] if len(urf) else pd.Series(dtype=float)
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


def plot_delta2(ax, df, show_ylabel=True, legend=True):
    line_df = unique_N_for_line(df, 'delta2')
    if len(line_df):
        ax.plot(line_df['N'], line_df['delta2'], '--', color='black', label=r'$\delta^2$', zorder=2)

    exp = df[(df['idx'] == 1) & df['delta2'].notna()]
    if len(exp):
        draw_exp_points_with_errors(
            ax,
            exp['N'], exp['delta2'],
            exp['delta2_minus'], exp['delta2_plus'],
            label='idx=1'
        )

    vals = pd.concat([
        line_df['delta2'] if len(line_df) else pd.Series(dtype=float),
        (exp['delta2'] - exp['delta2_minus']) if len(exp) else pd.Series(dtype=float),
        (exp['delta2'] + exp['delta2_plus']) if len(exp) else pd.Series(dtype=float),
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


def make_combined_plot(file_left: str, file_right: str, show=True):
    df_left = load_data(file_left)
    df_right = load_data(file_right)

    label_left = infer_label(file_left)
    label_right = infer_label(file_right)

    fig, axes = plt.subplots(
        3, 2, figsize=(10.8, 9.2), sharex='col',
        gridspec_kw={'hspace': 0.0, 'wspace': 0.12}
    )

    plot_qalpha(axes[0, 0], df_left, show_ylabel=True, legend=True)
    plot_half_life(axes[1, 0], df_left, show_ylabel=True, legend=True)
    plot_delta2(axes[2, 0], df_left, show_ylabel=True, legend=True)

    plot_qalpha(axes[0, 1], df_right, show_ylabel=False, legend=True)
    plot_half_life(axes[1, 1], df_right, show_ylabel=False, legend=True)
    plot_delta2(axes[2, 1], df_right, show_ylabel=False, legend=True)

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

    out_png = f'{label_left}_{label_right}_combined_triple_fixed.png'
    out_pdf = f'{label_left}_{label_right}_combined_triple_fixed.pdf'
    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
    print(f'已保存: {out_png}')
    print(f'已保存: {out_pdf}')

    if show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('用法: python CombinedTriplePlotDecayPro_fixed.py dataHsPro.dat dataDsPro.dat')
        print('或  : python CombinedTriplePlotDecayPro_fixed.py dataHsPro.dat dataDsPro.dat --no-show')
        sys.exit(1)

    file_left = sys.argv[1]
    file_right = sys.argv[2]
    show_flag = '--no-show' not in sys.argv
    make_combined_plot(file_left, file_right, show=show_flag)
