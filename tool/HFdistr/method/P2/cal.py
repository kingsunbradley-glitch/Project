import pandas as pd
import numpy as np

def get_magic_interval(value, magic_list):
    """
    获取 N 和 Z 所在的 magic-plus-one 区间 (Ni, Ni+1] 或 (Zi, Zi+1]
    """
    for i in range(len(magic_list) - 1):
        if magic_list[i] < value <= magic_list[i+1]:
            return magic_list[i], magic_list[i+1]
    
    return magic_list[-2], magic_list[-1]

def calculate_log10_T(row):
    """
    根据 Poenaru 等人 (1980) 的半经验公式计算半衰期的对数 log10(T)
    """
    Z = row['Z']
    A = row['A']
    E_alpha = row['Ealpha_MeV'] 
    
    A_d = A - 4
    Z_d = Z - 2
    N = A - Z
    
    # 根据公式 Q = Ea * A / Ad 将 E_alpha 转换为 Q 值
    Q = E_alpha * A / A_d
    
    # 论文中给出的 magic-plus-one 序列
    N_i_list = [3, 9, 21, 29, 51, 83, 127, 185, 250]
    Z_i_list = [3, 9, 21, 29, 51, 83, 115, 121, 150]
    
    # 获取区间
    n_i, n_i1 = get_magic_interval(N, N_i_list)
    z_i, z_i1 = get_magic_interval(Z, Z_i_list)
    
    # 计算无量纲变量 y 和 z
    y = (N - n_i) / (n_i1 - n_i)
    z_val = (Z - z_i) / (z_i1 - z_i)
    
    # 计算 x
    x = 0.4253 * Q * (1.5874 + A_d**(1/3)) / Z_d
    
    # 物理意义上 x 必须属于 (0, 1]
    if x <= 0 or x > 1:
        return np.nan
    
    # 计算 Ks
    term1 = 2.52956 * Z_d * np.sqrt(A_d / (A * Q))
    term2 = np.arccos(np.sqrt(x)) - np.sqrt(x * (1 - x))
    K_s = term1 * term2
    
    # 论文给出的拟合参数 B1~B6
    B1 = 0.988662
    B2 = 0.016314
    B3 = 0.020433
    B4 = 0.027896
    B5 = -0.003033
    B6 = -0.16820
    
    # 多项式求和
    poly = B1 + B2*y + B3*z_val + B4*(y**2) + B5*y*z_val + B6*(z_val**2)
    
    # 计算 log10(T) (注意除以 ln(10) 以转换为底数为 10 的对数)
    log10_T = (poly * K_s) / np.log(10) - 20.446
    
    return log10_T

def main():
    # 1. 读取您提供的 input.csv 文件
    df = pd.read_csv('input.csv')
        
    # 2. 应用公式进行计算，将结果放入新列 'log10T_calc'
    df['log10T_calc'] = df.apply(calculate_log10_T, axis=1)
        
    # 3. 将对数半衰期转换为实际的秒数 T_calc_s
    df['T_calc_s'] = 10 ** df['log10T_calc']
    
    # 4. 计算阻碍因子 HF = T_exp / T_calc (利用 input.csv 中的 Texp_s 列)
    df['HF'] = df['Texp_s'] / df['T_calc_s']
        
    # 5. 将包含原数据和新计算结果的 DataFrame 保存为 CSV
    df.to_csv('output_Poenaru_HF.csv', index=False)
    
    print("计算完毕！前几行计算结果如下：")
    print(df[['Name', 'Texp_s', 'T_calc_s', 'HF']].head())

if __name__ == "__main__":
    main()