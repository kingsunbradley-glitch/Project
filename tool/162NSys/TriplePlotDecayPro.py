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

    # 统一列名，避免重复的 + / - 列混淆
    rename_map = {
        '+': 'T_plus',
        '-': 'T_minus',
        '+.1': 'delta2_plus',
        '-.1': 'delta2_minus',
    }
    df = df.rename(columns=rename_map)

    # 检查必要列
    required_cols = [
        'N', 'idx', 'EXP', 'error', 'UNEDF1', 'WS4+RBF',
        'T0.5', 'T_plus', 'T_minus', 'URF',
        'delta2', 'delta2_plus', 'delta2_minus'
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print('Error: 数据文件缺少以下列:')
        print(', '.join(missing))
        sys.exit(1)

    for col in required_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.sort_values('N').reset_index(drop=True)
    return df



def add_common_style(ax):
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='major', length=6, width=1.2)
    ax.tick_params(which='minor', length=3, width=1.0)



def annotate_shell(ax, y_pos=None):
    ax.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
    if y_pos is not None:
        ax.text(162.25, y_pos, r'$N=162$', fontsize=12, color='black')



def plot_qalpha(ax, df: pd.DataFrame):
    # 理论线
    ax.plot(df['N'], df['UNEDF1'], '--', color='#7F7F7F', label='UNEDF1', zorder=2)
    ax.plot(df['N'], df['WS4+RBF'], '--', color='black', label='WS4+RBF', zorder=2)

    # idx == 1 的实验点
    exp = df[(df['idx'] == 1) & df['EXP'].notna()]
    q_exp = exp['EXP'] / 1000.0
    q_err = exp['error'] / 1000.0
    ax.errorbar(
        exp['N'], q_exp, yerr=q_err,
        fmt='o', color='red', markersize=4,
        label='EXP', zorder=3
    )

    ymin = min(df['UNEDF1'].min(), df['WS4+RBF'].min(), np.nanmin(q_exp) if len(q_exp) else np.inf)
    ymax = max(df['UNEDF1'].max(), df['WS4+RBF'].max(), np.nanmax(q_exp) if len(q_exp) else -np.inf)
    pad = 0.08 * (ymax - ymin if ymax > ymin else 1.0)
    annotate_shell(ax, y_pos=ymin + 0.10 * (ymax - ymin if ymax > ymin else 1.0))
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_ylabel(r'$Q_{\alpha}$ (MeV)')
    ax.legend(loc='upper right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())



def plot_half_life(ax, df: pd.DataFrame):
    ax.set_yscale('log')

    # 实验点：idx == 1
    exp = df[(df['idx'] == 1) & df['T0.5'].notna()]
    exp = exp[exp['T0.5'] > 0]
    ax.errorbar(
        exp['N'], exp['T0.5'], yerr=[exp['T_minus'], exp['T_plus']],
        fmt='o', color='red', markersize=4,
        label='EXP', zorder=3
    )

    # URF 理论线
    urf = df[df['URF'].notna() & (df['URF'] > 0)]
    ax.plot(urf['N'], urf['URF'], '--', color='black', label='URF', zorder=2)

    # 仅在 idx == 1 的位置给理论值打黑点
    urf_idx1 = df[(df['idx'] == 1) & df['URF'].notna() & (df['URF'] > 0)]
    ax.plot(urf_idx1['N'], urf_idx1['URF'], 'ko', markersize=4, zorder=4)

    y_all = pd.concat([
        exp['T0.5'],
        urf['URF']
    ]).dropna()
    if len(y_all):
        y_min = y_all[y_all > 0].min()
        y_max = y_all.max()
        ax.set_ylim(10 ** np.floor(np.log10(y_min) - 0.5),
                    10 ** np.ceil(np.log10(y_max) + 0.5))

    annotate_shell(ax, y_pos=(10 ** (np.log10(ax.get_ylim()[0]) + 0.15 * (np.log10(ax.get_ylim()[1]) - np.log10(ax.get_ylim()[0])))))
    ax.set_ylabel(r'$T_{1/2}$ (s)')
    ax.legend(loc='lower right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))



def plot_delta2(ax, df: pd.DataFrame):
    # 所有中心值都用黑色虚线连接；idx==0 只连线不打点
    all_valid = df[df['delta2'].notna()]
    ax.plot(all_valid['N'], all_valid['delta2'], '--', color='black', label=r'$\delta^2$', zorder=2)

    # idx == 1 的实验/有效点加误差棒
    exp = df[(df['idx'] == 1) & df['delta2'].notna()]
    ax.errorbar(
        exp['N'], exp['delta2'],
        yerr=[exp['delta2_minus'], exp['delta2_plus']],
        fmt='o', color='red', markersize=4,
        label='idx=1', zorder=3
    )

    y_all = pd.concat([
        all_valid['delta2'],
        (exp['delta2'] - exp['delta2_minus']).dropna(),
        (exp['delta2'] + exp['delta2_plus']).dropna(),
    ]).dropna()
    if len(y_all):
        y_min = y_all.min()
        y_max = y_all.max()
        pad = 0.08 * (y_max - y_min if y_max > y_min else 1.0)
        ax.set_ylim(y_min - pad, y_max + pad)

    annotate_shell(ax, y_pos=ax.get_ylim()[0] + 0.10 * (ax.get_ylim()[1] - ax.get_ylim()[0]))
    ax.set_ylabel(r'$\delta^2$')
    ax.set_xlabel(r'Neutron number $N$')
    ax.legend(loc='upper right', frameon=False)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())



def make_plot(filename: str, show: bool = True):
    df = load_data(filename)

    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(6.2, 9.4), sharex=True,
        gridspec_kw={'hspace': 0.0}
    )

    plot_qalpha(ax1, df)
    plot_half_life(ax2, df)
    plot_delta2(ax3, df)

    nmin = int(df['N'].min())
    nmax = int(df['N'].max())
    ax1.set_xlim(nmin - 1, nmax + 1)

    for ax in (ax1, ax2, ax3):
        add_common_style(ax)

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    base = os.path.splitext(os.path.basename(filename))[0]
    out_png = f'{base}_triple.png'
    out_pdf = f'{base}_triple.pdf'

    plt.savefig(out_png, dpi=600, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=300, bbox_inches='tight')

    print(f'已保存: {out_png}')
    print(f'已保存: {out_pdf}')

    if show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('用法: python TriplePlotDecayPro.py dataHsPro.dat')
        print('或  : python TriplePlotDecayPro.py dataDsPro.dat --no-show')
        sys.exit(1)

    input_file = sys.argv[1]
    show_flag = '--no-show' not in sys.argv
    make_plot(input_file, show=show_flag)
