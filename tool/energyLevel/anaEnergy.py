import csv
import os
import re
import pandas as pd


def clean_text(x):
    """清理文本中的 BOM / 不可见空格 / 首尾空白"""
    if x is None:
        return ""
    s = str(x)
    s = s.replace('\ufeff', '').replace('\xa0', ' ')
    s = s.strip()
    return s


def canonical_col(col):
    """把列名标准化，解决 \xa0E / 空格 / 奇怪编码问题"""
    c = clean_text(col)
    c_no_space = re.sub(r'\s+', '', c)

    mapping = {
        'K': 'K',
        '宇称': '宇称',
        '量子数[N,nz,Λ]': '量子数',
        '量子数[N,nz,Λ': '量子数',
        '量子数': '量子数',
        'configuration': 'configuration',
        'β2': 'β2',
        'γ': 'γ',
        'β4': 'β4',
        'E': 'E',
        'Ex': 'Ex',
    }
    return mapping.get(c_no_space, c_no_space)


def load_data(filename):
    """
    新版稳健读取：
    1. 按 tab 读取，保留空列（避免 configuration 空列导致错位）
    2. 清理 NBSP / BOM
    3. 统一列名
    4. 数值列强制转数值
    """
    rows = []
    with open(filename, 'r', encoding='utf-8-sig', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            # 保留空列，但要清理内容
            cleaned = [clean_text(x) for x in row]
            # 去掉行尾多余空列
            while cleaned and cleaned[-1] == '':
                cleaned.pop()
            if cleaned:
                rows.append(cleaned)

    if not rows:
        raise ValueError(f'{filename} 为空或无法读取')

    header = [canonical_col(x) for x in rows[0]]
    data_rows = rows[1:]

    # 对齐列数：短的补空，长的截断（一般不会长）
    ncol = len(header)
    aligned_rows = []
    for r in data_rows:
        if len(r) < ncol:
            r = r + [''] * (ncol - len(r))
        elif len(r) > ncol:
            r = r[:ncol]
        aligned_rows.append(r)

    df = pd.DataFrame(aligned_rows, columns=header)

    # 如果有空的 configuration 列就保留；不影响后续
    # 统一数值列
    for col in ['K', '量子数', 'β2', 'γ', 'β4', 'E', 'Ex']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    if '宇称' in df.columns:
        df['宇称'] = df['宇称'].astype(str).str.strip()

    required = ['K', '宇称', '量子数', 'E']
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'{filename} 缺少必要列: {missing}; 当前列为: {list(df.columns)}')

    return df


def calculate_internal_transitions(df, nucleus_name, output_filename, only_positive=False):
    """计算同一核内能级差。优先用 Ex；若 Ex 缺失则退回 E。"""
    transitions = []

    use_ex = 'Ex' in df.columns and df['Ex'].notna().any()

    for _, initial_row in df.iterrows():
        for _, final_row in df.iterrows():
            if use_ex and pd.notna(initial_row['Ex']) and pd.notna(final_row['Ex']):
                e_diff = initial_row['Ex'] - final_row['Ex']
            else:
                e_diff = initial_row['E'] - final_row['E']

            if pd.isna(e_diff):
                continue
            if only_positive and e_diff <= 0:
                continue

            transitions.append({
                f'{nucleus_name}_初态_K': initial_row['K'],
                f'{nucleus_name}_初态_宇称': initial_row['宇称'],
                f'{nucleus_name}_初态_量子数': int(initial_row['量子数']) if pd.notna(initial_row['量子数']) else None,
                f'{nucleus_name}_初态_Ex(MeV)': initial_row['Ex'] if 'Ex' in df.columns else None,
                f'{nucleus_name}_末态_K': final_row['K'],
                f'{nucleus_name}_末态_宇称': final_row['宇称'],
                f'{nucleus_name}_末态_量子数': int(final_row['量子数']) if pd.notna(final_row['量子数']) else None,
                f'{nucleus_name}_末态_Ex(MeV)': final_row['Ex'] if 'Ex' in df.columns else None,
                '能级差_E(MeV)': round(float(e_diff), 4),
            })

    out = pd.DataFrame(transitions)
    if not out.empty:
        out = out.sort_values('能级差_E(MeV)', ascending=False)
    out.to_csv(output_filename, index=False, encoding='utf-8-sig')
    print(f'已生成能级差数据表: {output_filename} (共 {len(out)} 条)')
    return out


def calculate_alpha_decay(ds_file='mm273Ds.txt', hs_file='mm269Hs.txt'):
    print('正在加载数据...')
    ds_df = load_data(ds_file)
    hs_df = load_data(hs_file)

    print('273Ds 列名:', list(ds_df.columns))
    print('269Hs 列名:', list(hs_df.columns))

    # 定义质量数，用于计算反冲动能
    A_parent = 273.0
    A_daughter = 269.0

    # 273Ds -> 269Hs alpha 衰变
    alpha_results = []
    for _, ds_row in ds_df.iterrows():
        for _, hs_row in hs_df.iterrows():
            if pd.isna(ds_row['E']) or pd.isna(hs_row['E']):
                continue

            q_value = ds_row['E'] - hs_row['E'] + 28.3
            e_alpha = q_value * (A_daughter / A_parent)

            alpha_results.append({
                '273Ds_K': ds_row['K'],
                '273Ds_宇称': ds_row['宇称'],
                '273Ds_量子数': int(ds_row['量子数']) if pd.notna(ds_row['量子数']) else None,
                '273Ds_Ex(MeV)': ds_row['Ex'] if 'Ex' in ds_df.columns else None,
                '269Hs_K': hs_row['K'],
                '269Hs_宇称': hs_row['宇称'],
                '269Hs_量子数': int(hs_row['量子数']) if pd.notna(hs_row['量子数']) else None,
                '269Hs_Ex(MeV)': hs_row['Ex'] if 'Ex' in hs_df.columns else None,
                'Q值(MeV)': round(float(q_value), 4),
                'Alpha动能(MeV)': round(float(e_alpha), 4),
            })

    df_alpha = pd.DataFrame(alpha_results)
    df_alpha = df_alpha.sort_values('Alpha动能(MeV)', ascending=False)
    df_alpha.to_csv('Alpha_Decay_Energies_Corrected.csv', index=False, encoding='utf-8-sig')
    print(f'已生成修正后的 Alpha 衰变数据表: Alpha_Decay_Energies_Corrected.csv (共 {len(df_alpha)} 条数据)')

    # 同核内部能级差
    calculate_internal_transitions(ds_df, '273Ds', '273Ds_Internal_Transitions.csv', only_positive=False)
    calculate_internal_transitions(hs_df, '269Hs', '269Hs_Internal_Transitions.csv', only_positive=False)

    # 额外给一个只保留正能级差的版本，更像“初态->低末态跃迁”
    calculate_internal_transitions(ds_df, '273Ds', '273Ds_Internal_Transitions_PositiveOnly.csv', only_positive=True)
    calculate_internal_transitions(hs_df, '269Hs', '269Hs_Internal_Transitions_PositiveOnly.csv', only_positive=True)


if __name__ == '__main__':
    calculate_alpha_decay()
