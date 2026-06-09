'''
py -m http.server 8000
'''


import argparse
import json
import os
import re
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/matplotlib-cache')
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
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
    'axes.unicode_minus': False,
})


UNIT_TO_MS = {
    '': 1.0,
    'ms': 1.0,
    's': 1.0e3,
    'sec': 1.0e3,
    'min': 60.0e3,
}

ELEMENT_STYLES = {
    'Rg': {'color': '#D62728', 'marker': 'o', 'plotly_symbol': 'circle'},
    'Ds': {'color': '#1F77B4', 'marker': 's', 'plotly_symbol': 'square'},
    'Ds^m': {'color': '#17BECF', 'marker': '*', 'plotly_symbol': 'star'},
    'Mt': {'color': '#2CA02C', 'marker': '^', 'plotly_symbol': 'triangle-up'},
    'Hs': {'color': '#9467BD', 'marker': 'D', 'plotly_symbol': 'diamond'},
    'Bh': {'color': '#FF7F0E', 'marker': 'v', 'plotly_symbol': 'triangle-down'},
    'Sg': {'color': '#8C564B', 'marker': 'P', 'plotly_symbol': 'cross'},
}

SPECIAL_POINTS = [
    {
        'A': 273,
        'Z': 110,
        'nucl': 'Ds^m',
        'N': 163,
        'element': 'Ds^m',
        'T_half': 17.34,
        'T_plus': 31.55,
        'T_minus': 6.83,
    },
]


INTERACTIVE_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Half-life by nuclide</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {
    margin: 0;
    background: #f3f4f6;
    color: #111827;
    font-family: "Times New Roman", "STIXGeneral", "DejaVu Serif", serif;
  }
  .layout {
    display: grid;
    grid-template-columns: minmax(720px, 1fr) 280px;
    gap: 16px;
    min-height: 100vh;
    padding: 16px;
    box-sizing: border-box;
  }
  .plot-wrap,
  .controller {
    background: #ffffff;
    border: 1px solid #d1d5db;
    border-radius: 6px;
  }
  .plot-wrap {
    min-width: 0;
  }
  #plot {
    width: 100%;
    height: calc(100vh - 34px);
    min-height: 620px;
  }
  .controller {
    max-height: calc(100vh - 34px);
    overflow: auto;
    padding: 12px;
    box-sizing: border-box;
  }
  .controller-title {
    margin: 0 0 10px;
    font-size: 18px;
    font-weight: 700;
  }
  .actions {
    display: flex;
    gap: 8px;
    margin-bottom: 12px;
  }
  button {
    font: inherit;
    cursor: pointer;
  }
  .action {
    border: 1px solid #9ca3af;
    border-radius: 4px;
    background: #ffffff;
    padding: 4px 9px;
  }
  .element-button {
    width: 100%;
    display: flex;
    align-items: center;
    justify-content: space-between;
    margin: 6px 0;
    border: 1px solid #d1d5db;
    border-left-width: 5px;
    border-radius: 4px;
    background: #ffffff;
    padding: 7px 9px;
    text-align: left;
  }
  .element-button.active {
    background: #f3f4f6;
    border-color: #6b7280;
    font-weight: 700;
  }
  .element-meta {
    color: #6b7280;
    font-size: 12px;
  }
  @media (max-width: 980px) {
    .layout {
      grid-template-columns: 1fr;
    }
    #plot {
      height: 620px;
    }
  }
</style>
</head>
<body>
<main class="layout">
  <section class="plot-wrap">
    <div id="plot"></div>
  </section>
  <aside class="controller">
    <p class="controller-title">Element visibility</p>
    <div class="actions">
      <button class="action" type="button" id="hide-all">Hide all</button>
      <button class="action" type="button" id="show-all">Show all</button>
    </div>
    <div id="element-list"></div>
  </aside>
</main>
<script>
const markerTraces = __MARKER_TRACES__;
const elementCounts = __ELEMENT_COUNTS__;
const styles = __STYLES__;
const elementOrder = __ELEMENT_ORDER__;
const layout = __LAYOUT__;
const config = {responsive: true, displaylogo: false};
const activeElements = new Set(elementOrder);
const buttonsByElement = new Map();
const traceByElement = new Map(markerTraces.map((trace, index) => [trace.name, index]));

Plotly.newPlot("plot", markerTraces, layout, config).then(() => {
  buildController();
  document.getElementById("plot").on("plotly_legendclick", event => {
    const trace = markerTraces[event.curveNumber];
    if (!trace) return false;
    toggleElement(trace.name);
    return false;
  });
});

function setElement(element, visible) {
  const traceIndex = traceByElement.get(element);
  if (traceIndex === undefined) return;
  Plotly.restyle("plot", {visible: visible}, [traceIndex]);
  const button = buttonsByElement.get(element);
  if (visible) {
    activeElements.add(element);
    if (button) button.classList.add("active");
  } else {
    activeElements.delete(element);
    if (button) button.classList.remove("active");
  }
}

function toggleElement(element) {
  setElement(element, !activeElements.has(element));
}

function buildController() {
  const host = document.getElementById("element-list");
  elementOrder.forEach(element => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = "element-button active";
    button.style.borderLeftColor = styles[element].color;
    button.innerHTML = `<span>${element}</span><span class="element-meta">${elementCounts[element] || 0} points</span>`;
    button.addEventListener("click", () => toggleElement(element));
    buttonsByElement.set(element, button);
    host.appendChild(button);
  });

  document.getElementById("hide-all").addEventListener("click", () => {
    elementOrder.forEach(element => setElement(element, false));
  });
  document.getElementById("show-all").addEventListener("click", () => {
    elementOrder.forEach(element => setElement(element, true));
  });
}
</script>
</body>
</html>
"""


def parse_value(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if not text:
        return np.nan

    match = re.fullmatch(r'([+-]?\d+(?:\.\d+)?)(ms|min|sec|s)?', text)
    if not match:
        return np.nan

    number = float(match.group(1))
    unit = match.group(2) or ''
    return number * UNIT_TO_MS[unit]


def parse_half_life(value):
    if pd.isna(value):
        return np.nan, np.nan, np.nan
    text = str(value).strip()
    if not text:
        return np.nan, np.nan, np.nan

    match = re.fullmatch(
        r'([+-]?\d+(?:\.\d+)?)\+([+-]?\d+(?:\.\d+)?)_-([+-]?\d+(?:\.\d+)?)(ms|min|sec|s)?',
        text
    )
    if match:
        unit = match.group(4) or ''
        scale = UNIT_TO_MS[unit]
        return (
            float(match.group(1)) * scale,
            float(match.group(2)) * scale,
            float(match.group(3)) * scale,
        )

    return parse_value(text), np.nan, np.nan


def element_symbol(nuclide):
    match = re.match(r'([A-Za-z]+)', str(nuclide))
    return match.group(1) if match else str(nuclide)


def add_special_points(df):
    if not SPECIAL_POINTS:
        return df

    special = pd.DataFrame(SPECIAL_POINTS)
    for col in df.columns:
        if col not in special.columns:
            special[col] = np.nan
    special = special[df.columns]
    return pd.concat([df, special], ignore_index=True)


def load_data(filename):
    try:
        df = pd.read_csv(filename, sep=r'\s+', engine='python')
    except FileNotFoundError:
        print(f'Error: 找不到文件 {filename}')
        sys.exit(1)
    except Exception as exc:
        print(f'Error: 读取文件 {filename} 失败: {exc}')
        sys.exit(1)

    required = ['A', 'Z', 'nucl', 'N', 'T0.5', 'T+', 'T-']
    missing = [col for col in required if col not in df.columns]
    if missing:
        print('Error: 数据文件缺少以下列:')
        print(', '.join(missing))
        sys.exit(1)

    df['element'] = df['nucl'].map(element_symbol)
    for col in ['A', 'Z', 'N']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    parsed = df['T0.5'].map(parse_half_life)
    df['T_half'] = [item[0] for item in parsed]
    df['T_plus_embedded'] = [item[1] for item in parsed]
    df['T_minus_embedded'] = [item[2] for item in parsed]
    df['T_plus'] = df['T+'].map(parse_value)
    df['T_minus'] = df['T-'].map(parse_value)

    missing_plus = df['T_plus'].isna() & df['T_plus_embedded'].notna()
    missing_minus = df['T_minus'].isna() & df['T_minus_embedded'].notna()
    df.loc[missing_plus, 'T_plus'] = df.loc[missing_plus, 'T_plus_embedded']
    df.loc[missing_minus, 'T_minus'] = df.loc[missing_minus, 'T_minus_embedded']

    df = df[df['N'].notna() & df['T_half'].notna() & (df['T_half'] > 0)].copy()
    df = add_special_points(df)
    df = df.sort_values(['Z', 'N']).reset_index(drop=True)
    return df


def add_common_style(ax):
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='major', length=6, width=1.2)
    ax.tick_params(which='minor', length=3, width=1.0)


def annotate_shell(ax):
    ax.axvline(x=162, color='gray', linestyle=':', linewidth=1.2, zorder=1)
    y0, y1 = ax.get_ylim()
    y = 10 ** (np.log10(y0) + 0.15 * (np.log10(y1) - np.log10(y0)))
    ax.text(162.25, y, r'$N=162$', fontsize=12, color='black')


def plot_half_life(df, output_prefix):
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.set_yscale('log')

    for element, group in df.groupby('element', sort=False):
        group = group.sort_values('N')
        style = ELEMENT_STYLES.get(element, {'color': 'black', 'marker': 'o'})
        x = group['N'].astype(float).to_numpy()
        y = group['T_half'].astype(float).to_numpy()
        yerr_low = group['T_minus'].astype(float).to_numpy()
        yerr_high = group['T_plus'].astype(float).to_numpy()
        ax.errorbar(
            x,
            y,
            yerr=[yerr_low, yerr_high],
            fmt=style['marker'],
            color=style['color'],
            ecolor=style['color'],
            markersize=5,
            elinewidth=1.2,
            capsize=3,
            linestyle='none',
            label=rf'${element}$',
            zorder=3,
        )

    y_all = df['T_half'].dropna()
    if len(y_all):
        ax.set_ylim(
            10 ** np.floor(np.log10(y_all.min()) - 0.5),
            10 ** np.ceil(np.log10(y_all.max()) + 0.5),
        )

    ax.set_xlim(int(df['N'].min()) - 1, int(df['N'].max()) + 2)
    annotate_shell(ax)
    add_common_style(ax)
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1))
    ax.set_xlabel(r'Neutron number $N$')
    ax.set_ylabel(r'$T_{1/2}$ (ms)')
    ax.legend(loc='lower right', frameon=False, ncol=2)

    out_png = f'{output_prefix}.png'
    out_pdf = f'{output_prefix}.pdf'
    fig.savefig(out_png, dpi=600, bbox_inches='tight')
    fig.savefig(out_pdf, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存: {out_png}')
    print(f'已保存: {out_pdf}')


def finite_or_none(value):
    if value is None or pd.isna(value):
        return None
    return float(value)


def element_order(df):
    present = list(dict.fromkeys(df['element'].tolist()))
    preferred = [el for el in ELEMENT_STYLES if el in present]
    extras = [el for el in present if el not in ELEMENT_STYLES]
    return preferred + extras


def make_marker_trace(element, group):
    style = ELEMENT_STYLES.get(element, {
        'color': '#000000',
        'plotly_symbol': 'circle',
    })
    group = group.sort_values('N')
    return {
        'type': 'scatter',
        'mode': 'markers',
        'name': element,
        'x': group['N'].astype(float).tolist(),
        'y': group['T_half'].astype(float).tolist(),
        'customdata': [[row['id'], row['label']] for _, row in group.iterrows()],
        'hovertemplate': (
            '%{customdata[1]}<br>'
            'N=%{x}<br>'
            'T<sub>1/2</sub>=%{y:.4g} ms'
            '<extra></extra>'
        ),
        'marker': {
            'color': style['color'],
            'symbol': style['plotly_symbol'],
            'size': 10,
            'line': {'color': style['color'], 'width': 1.2},
        },
        'error_y': {
            'type': 'data',
            'symmetric': False,
            'array': [finite_or_none(v) for v in group['T_plus']],
            'arrayminus': [finite_or_none(v) for v in group['T_minus']],
            'color': style['color'],
            'thickness': 1.4,
            'width': 4,
            'visible': True,
        },
    }


def make_interactive(df, output_prefix):
    df = df.copy()
    df['label'] = df['A'].astype(int).astype(str) + df['element']
    df['id'] = df['label'] + '_' + df['N'].astype(int).astype(str)

    orders = element_order(df)
    marker_traces = [
        make_marker_trace(element, df[df['element'] == element])
        for element in orders
    ]
    element_counts = {
        element: int((df['element'] == element).sum())
        for element in orders
    }

    y_min = float(df['T_half'].min())
    y_max = float(df['T_half'].max())
    y_range = [
        np.floor(np.log10(y_min) - 0.5),
        np.ceil(np.log10(y_max) + 0.5),
    ]
    x_min = int(df['N'].min()) - 1
    x_max = int(df['N'].max()) + 2

    layout = {
        'template': 'plotly_white',
        'margin': {'l': 88, 'r': 28, 't': 26, 'b': 78},
        'font': {'family': 'Times New Roman, STIXGeneral, DejaVu Serif, serif', 'size': 16},
        'xaxis': {
            'title': {'text': 'Neutron number <i>N</i>', 'font': {'size': 24}},
            'range': [x_min, x_max],
            'dtick': 2,
            'ticks': 'inside',
            'mirror': 'ticks',
            'showline': True,
            'linecolor': '#000000',
            'linewidth': 1.6,
            'showgrid': False,
            'zeroline': False,
        },
        'yaxis': {
            'title': {'text': '<i>T</i><sub>1/2</sub> (ms)', 'font': {'size': 24}},
            'type': 'log',
            'range': y_range,
            'ticks': 'inside',
            'mirror': 'ticks',
            'showline': True,
            'linecolor': '#000000',
            'linewidth': 1.6,
            'showgrid': False,
            'zeroline': False,
        },
        'legend': {
            'x': 0.82,
            'y': 0.08,
            'xanchor': 'left',
            'yanchor': 'bottom',
            'bgcolor': 'rgba(255,255,255,0)',
            'borderwidth': 0,
            'orientation': 'v',
            'font': {'size': 16},
        },
        'shapes': [{
            'type': 'line',
            'xref': 'x',
            'yref': 'paper',
            'x0': 162,
            'x1': 162,
            'y0': 0,
            'y1': 1,
            'line': {'color': '#7f7f7f', 'width': 1.5, 'dash': 'dot'},
        }],
        'annotations': [{
            'x': 162.25,
            'y': 10 ** (y_range[0] + 0.15 * (y_range[1] - y_range[0])),
            'xref': 'x',
            'yref': 'y',
            'text': '<i>N</i>=162',
            'showarrow': False,
            'font': {'size': 18, 'color': '#000000'},
            'xanchor': 'left',
        }],
    }

    style_payload = {
        element: {'color': ELEMENT_STYLES[element]['color']}
        for element in orders
        if element in ELEMENT_STYLES
    }
    for element in orders:
        style_payload.setdefault(element, {'color': '#000000'})

    html = INTERACTIVE_TEMPLATE
    replacements = {
        '__MARKER_TRACES__': marker_traces,
        '__ELEMENT_COUNTS__': element_counts,
        '__STYLES__': style_payload,
        '__ELEMENT_ORDER__': orders,
        '__LAYOUT__': layout,
    }
    for token, payload in replacements.items():
        html = html.replace(token, json.dumps(payload, ensure_ascii=False))

    out_html = f'{output_prefix}.html'
    with open(out_html, 'w', encoding='utf-8') as handle:
        handle.write(html)
    print(f'已保存: {out_html}')


def main():
    parser = argparse.ArgumentParser(
        description='按元素绘制 162Sys/data.dat 的半衰期图，并生成可交互 HTML。'
    )
    parser.add_argument('input', nargs='?', default='data.dat')
    parser.add_argument('-o', '--output-prefix')
    parser.add_argument('--no-static', action='store_true', help='不生成 PNG/PDF')
    parser.add_argument('--no-html', action='store_true', help='不生成交互 HTML')
    args = parser.parse_args()

    input_path = args.input
    if not os.path.isabs(input_path):
        input_path = os.path.abspath(input_path)

    output_prefix = args.output_prefix
    if output_prefix is None:
        base = os.path.splitext(os.path.basename(input_path))[0]
        output_prefix = os.path.join(os.path.dirname(input_path), f'{base}_halflife_by_element')

    df = load_data(input_path)
    if not args.no_static:
        plot_half_life(df, output_prefix)
    if not args.no_html:
        make_interactive(df, output_prefix)


if __name__ == '__main__':
    main()
