#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
专门绘制 Z=110 (Ds) 的 Qalpha 对比图
1. 包含所有理论模型
2. 包含现有的实验数据
3. 特别添加 N=163 的新实验点 (红点)，且不显示文字标注
4. 图例字号调大
"""

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# ==========================================
# 1. 数据读取辅助函数 (之前缺失的部分)
# ==========================================

def _try_parse_numeric(parts):
    try:
        return [float(x) for x in parts]
    except Exception:
        return None

def read_qalpha_table(path: str, model_name: str) -> pd.DataFrame:
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                # 跳过空行、注释或非数字开头的行
                if not s or s.startswith("#") or s.startswith("[") or not s[0].isdigit():
                    continue
                parts = s.split()
                nums = _try_parse_numeric(parts)
                if nums is None:
                    continue
                
                # 兼容不同的数据格式
                if len(nums) >= 4: 
                    # 格式: Z N A Q ...
                    Z = int(round(nums[0]))
                    N = int(round(nums[1]))
                    A = int(round(nums[2]))
                    Q = float(nums[3])
                    rows.append((Z, N, A, Q))
                elif len(nums) == 3: 
                    # 格式: A Z Q (WS4等可能格式)
                    A = int(round(nums[0]))
                    Z = int(round(nums[1]))
                    Q = float(nums[2])
                    N = A - Z
                    rows.append((Z, N, A, Q))
    except FileNotFoundError:
        print(f"[提示] 未找到文件: {path} (跳过模型 {model_name})")
        return pd.DataFrame()
    
    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha"])
    df["model"] = model_name
    return df

def read_exp_table(path: str) -> pd.DataFrame:
    rows = []
    try:
        with open(path, "r", errors="replace") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#") or s.startswith("[") or not s[0].isdigit():
                    continue
                parts = s.split()
                nums = _try_parse_numeric(parts)
                if nums is None or len(nums) < 4:
                    continue
                # 假设 EXP 格式: A Z N Q(keV) Err(keV)
                if len(nums) >= 5:
                    A = int(round(nums[0]))
                    Z = int(round(nums[1]))
                    N = int(round(nums[2]))
                    Q_kev = float(nums[3])
                    Err_kev = float(nums[4])
                    # 转换为 MeV
                    rows.append((Z, N, A, Q_kev / 1000.0, Err_kev / 1000.0))
    except FileNotFoundError:
        print(f"[提示] 未找到实验数据文件: {path}")
        return pd.DataFrame()
    
    df = pd.DataFrame(rows, columns=["Z", "N", "A", "Qalpha", "Error"])
    df["model"] = "EXP"
    return df

# ==========================================
# 2. 绘图主逻辑
# ==========================================

def plot_single_z_special(target_Z, element_sym, nmin=145, nmax=180):
    # --- 配置要显示的模型 ---
    # 如果不想看某个模型，在行首加 # 注释掉即可
    model_files = {
        #"FRDM(2012)": "Q@a_FRDM.dat",
        #"HFB-32":     "Q@a_HFB.dat",
        "WS4+RBF":    "Q@a_WS4+RBF.dat",
        "UNEDF1":     "Q@a_UNEDF1.dat",
        #"SKMS":       "Q@a_SKMS.dat",
        #"SLY4":       "Q@a_SLY4.dat",
        "SV-MIN":     "Q@a_SV-MIN.dat"
    }
    
    # 定义颜色和标记风格
    model_styles = {
        "FRDM(2012)": {"c": "blue", "m": "o"},
        "HFB-32":     {"c": "green", "m": "v"},
        "WS4+RBF":    {"c": "red", "m": "^"},
        "UNEDF1":     {"c": "purple", "m": "s"},
        "SKMS":       {"c": "cyan", "m": "p"},
        "SLY4":       {"c": "magenta", "m": "*"},
        "SV-MIN":     {"c": "orange", "m": "D"},
    }
    # 指定图例排序
    model_order = ["FRDM(2012)", "HFB-32", "WS4+RBF", "UNEDF1", "SKMS", "SLY4", "SV-MIN"]

    # --- 加载数据 ---
    dfs = []
    print("正在读取模型数据...")
    for name, path in model_files.items():
        df = read_qalpha_table(path, name)
        if not df.empty:
            dfs.append(df)
    
    df_models = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    df_exp = read_exp_table("Q@a_EXP.dat")

    # --- 开始绘图 ---
    print(f"正在绘制 {element_sym} (Z={target_Z})...")
    fig, ax = plt.subplots(figsize=(10, 7))

    # 1. 绘制理论模型曲线
    if not df_models.empty:
        d = df_models[(df_models["Z"] == target_Z) & (df_models["N"] >= nmin) & (df_models["N"] <= nmax)]
        for model_name in model_order:
            if model_name not in d["model"].unique(): continue
            dm = d[d["model"] == model_name].sort_values("N")
            if dm.empty: continue
            
            st = model_styles.get(model_name, {"c": "gray", "m": "."})
            ax.plot(dm["N"], dm["Qalpha"], marker=st["m"], color=st["c"],
                    markersize=6, linewidth=1.5, label=model_name, alpha=0.7)

    # 2. 绘制已有的实验数据 (空心黑方块)
    if not df_exp.empty:
        de = df_exp[(df_exp["Z"] == target_Z) & (df_exp["N"] >= nmin) & (df_exp["N"] <= nmax)].sort_values("N")
        if not de.empty:
            ax.errorbar(de["N"], de["Qalpha"], yerr=de["Error"], 
                        fmt='s', color='black', label='Existing EXP', 
                        capsize=4, markersize=7, zorder=10, 
                        markerfacecolor='none', markeredgewidth=2)

    # =================================================================
    # 3. --- 添加新的实验点 (红点) ---
    # 数据: Z=110, N=163, Q=11270 keV, Error=70 keV
    # =================================================================
    new_N = 163
    new_Q_MeV = 11270 / 1000.0  # 11.270 MeV
    new_Err_MeV = 70 / 1000.0   # 0.070 MeV
    
    ax.errorbar(new_N, new_Q_MeV, yerr=new_Err_MeV,
                fmt='o', color='red', label='high energy EXP (N=163)',
                capsize=5, markersize=12, zorder=20, # zorder=20 保证在最上层
                markerfacecolor='red', markeredgecolor='darkred', markeredgewidth=2)
    
    # 【已根据要求注释掉】 不显示文字标注
    # ax.annotate(f'New Point\nN={new_N}\nQ={new_Q_MeV:.3f} MeV', 
    #             xy=(new_N, new_Q_MeV), 
    #             xytext=(new_N - 2, new_Q_MeV - 0.3),
    #             arrowprops=dict(facecolor='red', arrowstyle='->', linewidth=1.5),
    #             fontsize=11, color='red', fontweight='bold')
    # =================================================================

    # --- 格式化设置 ---
    ax.set_title(f"$Q_\\alpha$ Values for {element_sym} (Z={target_Z}) ", fontsize=16)
    ax.set_xlabel("Neutron number N", fontsize=14)
    ax.set_ylabel(r"$Q_\alpha$ (MeV)", fontsize=14)
    
    # 设置刻度字号
    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # 强制 X 轴只显示整数刻度
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    # 显示网格
    ax.grid(True, linestyle=":", linewidth=1, alpha=0.6)
    
    # 设置图例大小 (已设置为 large)
    ax.legend(fontsize='large', loc='best', framealpha=0.9)

    plt.tight_layout()
    output_file = f"Qalpha_Z{target_Z}_{element_sym}_special.png"
    plt.savefig(output_file, dpi=300)
    #print(f"绘图完成！图片已保存至: {output_file}")
    plt.show()

if __name__ == "__main__":
    # 执行 Z=110 (Ds) 的绘图，聚焦范围 N=155-175
    plot_single_z_special(target_Z=110, element_sym="Ds", nmin=155, nmax=175)