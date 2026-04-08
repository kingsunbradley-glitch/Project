import pandas as pd

def load_data(filename):
    # 读取制表符分割的txt文件
    df = pd.read_csv(filename, sep='\t')
    df.columns = df.columns.str.strip()
    return df

def calculate_alpha_decay():
    print("正在加载数据...")
    try:
        ds_df = load_data('mm273Ds.txt')
        hs_df = load_data('mm269Hs.txt')
    except FileNotFoundError as e:
        print(f"找不到文件，请确保文本文件与脚本在同一目录: {e}")
        return

    # 找到具体的量子数列名
    ds_qnum_col = [col for col in ds_df.columns if '量子数' in col][0]
    hs_qnum_col = [col for col in hs_df.columns if '量子数' in col][0]

    # 定义质量数，用于计算反冲动能
    A_parent = 273.0
    A_daughter = 269.0

    # ==========================================
    # 任务 1: 计算 273Ds 到 269Hs 的 Q值 与 alpha 动能
    # ==========================================
    alpha_results = []
    
    for _, ds_row in ds_df.iterrows():
        for _, hs_row in hs_df.iterrows():
            # 1. 计算总衰变能 (Q值)
            q_value = ds_row['E'] - hs_row['E'] + 28.3
            
            # 2. 计算分配给 Alpha 粒子的实际动能 (排除子核反冲动能)
            e_alpha = q_value * (A_daughter / A_parent)
            
            alpha_results.append({
                '273Ds_K': ds_row['K'],
                '273Ds_宇称': ds_row['宇称'],
                '273Ds_量子数': ds_row[ds_qnum_col],
                '269Hs_K': hs_row['K'],
                '269Hs_宇称': hs_row['宇称'],
                '269Hs_量子数': hs_row[hs_qnum_col],
                'Q值(MeV)': round(q_value, 4),
                'Alpha动能(MeV)': round(e_alpha, 4) # 这是探测器实际测到的能量
            })
            
    df_alpha = pd.DataFrame(alpha_results)
    df_alpha.to_csv('Alpha_Decay_Energies_Corrected.csv', index=False, encoding='utf-8-sig')
    print(f"已生成修正后的 Alpha 衰变数据表: Alpha_Decay_Energies_Corrected.csv (共 {len(df_alpha)} 条数据)")

    # ==========================================
    # 任务 2: 计算 273Ds 和 269Hs 内部的能级差 (此部分与之前一致，不涉及反冲)
    # ==========================================
    def calculate_internal_transitions(df, qnum_col, output_filename):
        transitions = []
        for _, initial_row in df.iterrows():
            for _, final_row in df.iterrows():
                # 能级跃迁只看能量差
                e_diff = initial_row['E'] - final_row['E']
                
                transitions.append({
                    '初态_K': initial_row['K'],
                    '初态_宇称': initial_row['宇称'],
                    '初态_量子数': initial_row[qnum_col],
                    '末态_K': final_row['K'],
                    '末态_宇称': final_row['宇称'],
                    '末态_量子数': final_row[qnum_col],
                    '能级差_E(MeV)': round(e_diff, 4)
                })
        
        df_transitions = pd.DataFrame(transitions)
        df_transitions.to_csv(output_filename, index=False, encoding='utf-8-sig')
        print(f"已生成能级差数据表: {output_filename}")

    # 生成内部能级差表
    calculate_internal_transitions(ds_df, ds_qnum_col, '273Ds_Internal_Transitions.csv')
    calculate_internal_transitions(hs_df, hs_qnum_col, '269Hs_Internal_Transitions.csv')

if __name__ == "__main__":
    calculate_alpha_decay()