import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
import scipy.stats as stats
from matplotlib.offsetbox import AnchoredText
import warnings
import re

# --- 用户配置区 ---

# 在此处输入您的寿命数据
# 示例数据为 269Hs (单位: s) 
lifetimes_data = np.array([5.005,9.363, 4.978, 14.2, 0.27, 36, 22, 19.7]) #269Hs 
#lifetimes_data = np.array([0.107,0.184,0.52,0.0399,0.373,0.31,0.11]) #237Ds
# 在此处输入数据的时间单位 (例如: 'ms', 's', 'μs')
# 这个单位将用于所有输出和图表
time_unit = 'ms'

# 在此处输入核素的 LaTeX 格式名称
nuclide_name = r'$^{269}\mathrm{Hs}^{*}$'
#nuclide_name = r'$^{273}\mathrm{Ds}$'

# --- (已移除) ---
# 文件名逻辑已移至函数内部，此处不再需要

# --- 脚本核心功能 ---

def log_likelihood_T(T_val, N_events, sum_t):
    """
    计算以半衰期 T 为参数的对数似然函数
    """
    if T_val <= 0:
        return -np.inf
    lambda_val = np.log(2) / T_val
    if lambda_val <= 0:
        return -np.inf
    return N_events * np.log(lambda_val) - lambda_val * sum_t

def get_robust_brackets(N, S, t_half_mle):
    """
    使用 chi^2 分布为非对称区间估计一个稳健的搜索范围
    """
    df = 2 * N
    alpha = 1.0 - 0.9973 # 3-sigma
    alpha_low = alpha / 2.0
    alpha_high = 1.0 - (alpha / 2.0)
    
    try:
        chi2_lower_3s = stats.chi2.ppf(alpha_low, df)
        chi2_upper_3s = stats.chi2.ppf(alpha_high, df)
        
        tau_L_3s_est = (2 * S) / chi2_upper_3s
        tau_U_3s_est = (2 * S) / chi2_lower_3s
        
        t_half_L_3s_est = tau_L_3s_est * np.log(2)
        t_half_U_3s_est = tau_U_3s_est * np.log(2)
        
        return t_half_L_3s_est * 0.5, t_half_U_3s_est * 1.5
        
    except Exception as e:
        print(f"警告: 计算 chi^2 边界失败 ({e}). 使用基于 MLE 的回退方案。")
        return t_half_mle / 10.0, t_half_mle * 10.0

def find_intercepts(level, N_events, sum_t, t_half_mle, bracket_low, bracket_high):
    """
    寻找 log_likelihood_T(T) = level 的两个根
    """
    func_to_solve = lambda T: log_likelihood_T(T, N_events, sum_t) - level
    
    T_L, T_U = np.nan, np.nan
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            safe_bracket_low = max(bracket_low, 1e-9)
            if safe_bracket_low >= t_half_mle:
                safe_bracket_low = t_half_mle / 10.0
            sol_lower = root_scalar(func_to_solve, bracket=[safe_bracket_low, t_half_mle])
            T_L = sol_lower.root
        except ValueError:
            pass # T_L 保持 nan

        try:
            safe_bracket_high = max(bracket_high, t_half_mle * 1.1)
            if safe_bracket_high <= t_half_mle:
                 safe_bracket_high = t_half_mle * 10.0
            sol_upper = root_scalar(func_to_solve, bracket=[t_half_mle, safe_bracket_high])
            T_U = sol_upper.root
        except ValueError:
            pass # T_U 保持 nan
            
    return T_L, T_U

def calculate_mle_and_plot(data, unit, nuclide_name_latex):
    """
    主函数：计算 MLE, 置信区间, 并绘制似然函数图
    """
    
    # 1. 基本统计
    N = len(data)
    if N == 0:
        print("错误: 寿命数组为空。")
        return
    S = np.sum(data)
    
    # 2. 计算 MLE
    tau_mle = S / N
    t_half_mle = tau_mle * np.log(2)
    
    # 3. 计算最大对数似然值
    log_L_max = log_likelihood_T(t_half_mle, N, S)
    
    print("--- MLE 计算结果 ---")
    print(f"观测事件数 (N): {N}")
    print(f"平均寿命 (τ_MLE): {tau_mle:.3f} {unit}")
    print(f"半衰期 (T_1/2_MLE): {t_half_mle:.3f} {unit}")
    print(f"最大对数似然 (ln(L_max)): {log_L_max:.3f}")
    
    # 4. 定义置信区间水平
    level_1s = log_L_max - 0.5
    level_2s = log_L_max - 2.0
    level_3s = log_L_max - 4.5
    
    # 5. 寻找置信区间 (根)
    bracket_L, bracket_U = get_robust_brackets(N, S, t_half_mle)
    
    T_L_1s, T_U_1s = find_intercepts(level_1s, N, S, t_half_mle, bracket_L, bracket_U)
    T_L_2s, T_U_2s = find_intercepts(level_2s, N, S, t_half_mle, bracket_L, bracket_U)
    T_L_3s, T_U_3s = find_intercepts(level_3s, N, S, t_half_mle, bracket_L, bracket_U)
    
    # 6. 打印置信区间结果
    print("\n--- 置信区间 (C.I.) ---")
    print(f"MLE 半衰期 = {t_half_mle:.3f} {unit}")

    # 定义误差变量，默认为 nan
    err_low_1s, err_high_1s = np.nan, np.nan
    err_low_2s, err_high_2s = np.nan, np.nan
    err_low_3s, err_high_3s = np.nan, np.nan

    if not (np.isnan(T_L_1s) or np.isnan(T_U_1s)):
        err_low_1s = t_half_mle - T_L_1s
        err_high_1s = T_U_1s - t_half_mle
        print(f"1-sigma (68.3%) C.I.: [{T_L_1s:.3f}, {T_U_1s:.3f}] {unit} (即 +{err_high_1s:.3f} / -{err_low_1s:.3f} {unit})")
    
    if not (np.isnan(T_L_2s) or np.isnan(T_U_2s)):
        err_low_2s = t_half_mle - T_L_2s
        err_high_2s = T_U_2s - t_half_mle
        print(f"2-sigma (95.4%) C.I.: [{T_L_2s:.3f}, {T_U_2s:.3f}] {unit} (即 +{err_high_2s:.3f} / -{err_low_2s:.3f} {unit})")

    if not (np.isnan(T_L_3s) or np.isnan(T_U_3s)):
        err_low_3s = t_half_mle - T_L_3s
        err_high_3s = T_U_3s - t_half_mle
        print(f"3-sigma (99.7%) C.I.: [{T_L_3s:.3f}, {T_U_3s:.3f}] {unit} (即 +{err_high_3s:.3f} / -{err_low_3s:.3f} {unit})")

    # 7. 绘制图像
    print("\n正在生成图像...")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 设置绘图范围
    plot_min = bracket_L
    plot_max = bracket_U
    if np.isnan(T_L_3s) or np.isnan(T_U_3s):
        plot_min = t_half_mle / 5.0
        plot_max = t_half_mle * 5.0
        
    T_values = np.linspace(plot_min, plot_max, 500)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        log_L_values = np.array([log_likelihood_T(T, N, S) for T in T_values])
        
    norm_log_L = log_L_values - log_L_max
    
    # 绘制似然曲线 (移除 label)
    ax.plot(T_values, norm_log_L, color='black', linewidth=2)
    
    # 绘制置信区间水平线 (移除 label)
    ax.axhline(-0.5, linestyle='--', color='#d62728')
    ax.axhline(-2.0, linestyle=':', color='#2ca02c')
    ax.axhline(-4.5, linestyle='-.', color='#1f77b4')
    
    # 绘制 MLE 垂线 (移除 label)
    ax.axvline(t_half_mle, linestyle='-', color='gray', alpha=0.8)
    
    # 设置坐标轴 (使用 f-string 和 unit 变量)
    ax.set_xlabel(rf'Half-life $T_{{1/2}}$ ({unit})', fontsize=14)
    ax.set_ylabel(r'Normalized Log-Likelihood $\ln(L) - \ln(L_{\mathrm{max}})$', fontsize=14)
    
    # 设置 y 轴范围
    ax.set_ylim(-5.0, 0.5)
    # 设置 x 轴范围
    if not np.isnan(T_L_3s) and not np.isnan(T_U_3s):
         ax.set_xlim(T_L_3s * 0.8, T_U_3s * 1.2)
    else:
         ax.set_xlim(plot_min, plot_max)

    # 放置核素框
    at = AnchoredText(nuclide_name_latex, 
                      loc='upper right', 
                      prop=dict(size=14), 
                      frameon=True, 
                      bbox_to_anchor=(0.98, 0.98), 
                      bbox_transform=ax.transAxes)
    at.patch.set(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9)
    ax.add_artist(at)
    
    # --- 添加新的文本标注 ---
    
    # 获取 x 轴 55% 的位置 (使用您修改后的 0.55)
    plot_xlim = ax.get_xlim()
    x_pos_annot = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * 0.55

    # 1-sigma text
    if not np.isnan(err_low_1s):
        text_1s = rf'$1\sigma: T_{{1/2}} = {t_half_mle:.2f}^{{+{err_high_1s:.2f}}}_{{-{err_low_1s:.2f}}}$ {unit}'
        ax.text(x_pos_annot, -0.5 + 0.1, text_1s, ha='left', va='bottom', fontsize=11, color='#d62728')

    # 2-sigma text
    if not np.isnan(err_low_2s):
        text_2s = rf'$2\sigma: T_{{1/2}} = {t_half_mle:.2f}^{{+{err_high_2s:.2f}}}_{{-{err_low_2s:.2f}}}$ {unit}'
        ax.text(x_pos_annot, -2.0 + 0.1, text_2s, ha='left', va='bottom', fontsize=11, color='#2ca02c')

    # 3-sigma text
    if not np.isnan(err_low_3s):
        text_3s = rf'$3\sigma: T_{{1/2}} = {t_half_mle:.2f}^{{+{err_high_3s:.2f}}}_{{-{err_low_3s:.2f}}}$ {unit}'
        ax.text(x_pos_annot, -4.5 + 0.1, text_3s, ha='left', va='bottom', fontsize=11, color='#1f77b4')

    # 科学规范格式
    ax.tick_params(axis='both', which='major', direction='in', labelsize=12, top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', top=True, right=True)
    ax.minorticks_on()
    ax.grid(True, which='major', linestyle=':', linewidth=0.5, alpha=0.7)
    
    fig.tight_layout()
    
    # --- (已修复) 保存图像 ---
    
    # 1. 创建一个安全的文件名
    #    nuclide_name_latex 是传入的 LaTeX 字符串，例如 r'$^{269}\mathrm{Hs}$'
    #    我们移除所有 LaTeX 特殊字符 (\, $, ^, {, }) 和 \mathrm
    safe_filename_base = re.sub(r'[\$\^\{\}\\]|mathrm', '', nuclide_name_latex)
    #    结果 (示例): '269Hs'
    
    # 2. 组合文件名
    output_filename = f'{safe_filename_base}_mle_likelihood_plot.png'
    #    结果 (示例): '269Hs_mle_likelihood_plot.png'
    
    # 3. 保存图像
    fig.savefig(output_filename, dpi=300)
    plt.close(fig)
    
    print(f"图像已保存到: {output_filename}") # <-- 现在会打印正确的文件名

# --- 运行脚本 ---
if __name__ == "__main__":
    # 我们将用户配置区的变量传递给函数
    calculate_mle_and_plot(lifetimes_data, time_unit, nuclide_name)