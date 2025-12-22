import pandas as pd
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

# 读取 Excel 文件 (请修改文件名为您的实际文件名)
file_path = 'data.xlsx' 
# 如果没有文件，这里创建一个模拟数据演示
df = pd.read_excel(file_path) 

def calculate_single_width(row):
    # 1. 提取参数
    Am = row['A']
    Zm = row['Z']
    Q_input = row['Q']
    T_ms = row['T_ms']
    BR = row['BR']
    
    # 2. 基础物理常数
    A_dau = Am - 4.0
    Z_dau = Zm - 2.0
    Hplk = 4.13566727e-21  # MeV s
    Hbarc = 197.32696      # MeV fm
    eSqua = 1.4399644      # MeV fm
    Rmass = 4.0 * A_dau / (4.0 + A_dau) * 931.494013 # 折合质量
    
    # 3. 能量计算 (加上电子屏蔽 Escr)
    # Escr 公式: (65.3*Z^(7/5) - 80*Z^(2/5)) * 1e-6 MeV
    Escr = (65.3 * np.power(Zm, 7./5.) - 80. * np.power(Zm, 2./5.)) * 1.0E-6
    Etot = Q_input + Escr # 隧穿总能量
    
    # 4. 定义势能函数
    def Vtot(r):
        # Nuclear Potential (Igo)
        Vn = -1100.0 * np.exp(-((r - 1.17 * np.power(A_dau, 1./3.)) / 0.574))
        # Coulomb Potential
        Vc = 2.0 * Z_dau * eSqua / r
        # Centrifugal (L=0)
        Vcent = 0
        return Vn + Vc + Vcent

    # 5. 寻找转折点 (Vtot = Etot)
    # 辅助函数: V - E
    f = lambda r: Vtot(r) - Etot
    
    # 外转折点 Rout (库仑近似作为初猜)
    Rout_guess = 2.0 * Z_dau * eSqua / Etot
    try:
        Rout = brentq(f, 10.0, Rout_guess + 20.0)
    except:
        Rout = Rout_guess # 如果数值求解失败，回退到近似值
        
    # 内转折点 Rin (在核表面附近搜索)
    try:
        Rin = brentq(f, 2.0, Rout - 1.0)
    except:
        Rin = 1.2 * np.power(A_dau, 1./3.) # 失败时的兜底
        
    # 6. WKB 积分
    # 积分 sqrt(2*mu*(V-E)) dr
    def integrand(r):
        val = Vtot(r) - Etot
        return np.sqrt(2.0 * Rmass * val) if val > 0 else 0

    G_integral, _ = quad(integrand, Rin, Rout)
    
    # 7. 穿透因子 P
    # P = exp( -2/hbar * Integral )
    # 注意：integrand 已经是动量 p，除以 hbar 得到波数积分
    Gamow_factor = 2.0 * G_integral / Hbarc
    P = np.exp(-Gamow_factor)
    
    # 8. 计算宽度 Delta^2 (keV)
    # Delta^2 = h * ln2 * BR / (T * P)
    # T 需要转为秒: T_ms * 1e-3
    T_sec = T_ms * 1e-3
    Width_keV = Hplk * (np.log(2) / T_sec * BR) / P * 1000.0
    
    return pd.Series({
        'Escr': Escr,
        'Etot': Etot,
        'Rin': Rin,
        'Rout': Rout,
        'Penetration': P,
        'Width_keV': Width_keV
    })

# 批量应用计算
print("正在计算...")
results = df.apply(calculate_single_width, axis=1)
df_final = pd.concat([df, results], axis=1)

# 计算误差界限
# 这里的逻辑是：半衰期越短，宽度越大。所以 T-err 对应 Width+err
def calc_error(row):
    T = row['T_ms']
    W = row['Width_keV']
    T_min = T - row['Err_minus']
    T_max = T + row['Err_plus']
    
    # 简单比例缩放误差 (Width ~ 1/T)
    # W_max = W * (T / T_min)
    # W_min = W * (T / T_max)
    
    W_upper = W * (T / T_min) if T_min > 0 else np.nan
    W_lower = W * (T / T_max)
    
    return pd.Series({
        'Err_W_Plus': W_upper - W,
        'Err_W_Minus': W - W_lower
    })

errors = df_final.apply(calc_error, axis=1)
df_final = pd.concat([df_final, errors], axis=1)

print(df_final[['A', 'Z', 'Q', 'T_ms', 'Width_keV', 'Err_W_Plus', 'Err_W_Minus']])

# 保存结果
df_final.to_excel("result_width.xlsx", index=False)