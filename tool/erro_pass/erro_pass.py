import numpy as np

# ==================================================
# 地区         用户配置区 (USER CONFIG AREA)
# ==================================================
#
# 请在此处以 (E_i, sigma_i) 的格式输入您的数据
# E_i 是能量值, sigma_i 是对应的误差
#data = [(10857, 14),(10929, 22),] # Ds273a

data = [(11650, 14),(11017, 22),(11140,35),(11150,50),(11030,50),(11200,35),(11080,35),] # Ds273b

# 设置 Sigma 筛选的阈值
# 1.0 代表 1-sigma。
# 2.0 代表 2-sigma。
# 设为 None 或 0 将禁用筛选。
SIGMA_CLIPPING_THRESHOLD = 3.0
# SIGMA_CLIPPING_THRESHOLD = None # 如需禁用筛选，请取消本行注释

# ==================================================
# 地区            程序主逻辑 (DO NOT EDIT)
# ==================================================

def perform_sigma_clipping(data_points, threshold):
    """
    执行迭代的 Sigma Clipping 来剔除异常值。
    """
    if threshold is None or threshold <= 0:
        print("Sigma 筛选被禁用。\n")
        return data_points

    print(f"--- 正在执行 {threshold}σ 筛选 ---")
    
    current_data = np.array(data_points)
    iteration = 0
    
    while True:
        iteration += 1
        e_values = current_data[:, 0]
        
        if len(e_values) < 2:
            print("数据点过少，无法进行 Sigma 筛选。")
            break
            
        # 1. 计算 E_i 的简单平均值和样本标准差
        mean_e = np.mean(e_values)
        std_e = np.std(e_values, ddof=1) # ddof=1 使用样本标准差
        
        if std_e == 0:
            print("所有 E_i 值均相同，无法进行 Sigma 筛选。")
            break

        # 2. 定义筛选范围
        lower_bound = mean_e - threshold * std_e
        upper_bound = mean_e + threshold * std_e
        
        # 3. 筛选，找出在范围内的点 (mask)
        mask = (e_values >= lower_bound) & (e_values <= upper_bound)
        
        num_removed = len(current_data) - np.sum(mask)
        
        print(f"第 {iteration} 轮 (μ={mean_e:.3f}, S={std_e:.3f}): 范围 [{lower_bound:.3f}, {upper_bound:.3f}]")
        
        if num_removed == 0:
            print(f"筛选完成，没有点被剔除。剩余 {len(current_data)} 个数据点。\n")
            break
        else:
            # 打印被剔除的点
            removed_points = current_data[~mask]
            for point in removed_points:
                print(f"  [剔除] E={point[0]:<8} (σ={point[1]})")
            
            # 4. 保留“好”数据，进入下一轮
            current_data = current_data[mask]

    return list(current_data)

def calculate_weighted_mean(data_points):
    """
    计算加权平均值和其误差
    """
    if not data_points:
        print("错误：没有有效数据点用于计算加权平均。")
        return None, None, 0
    
    # 转换为 numpy 数组以便于计算
    data_array = np.array(data_points)
    e_values = data_array[:, 0]
    sigma_values = data_array[:, 1]

    # 1. 计算权重 w_i = 1 / sigma_i^2
    weights = 1.0 / (sigma_values**2)
    
    # 2. 计算加权平均 E
    # E = (Σ w_i * E_i) / (Σ w_i)
    E_weighted = np.sum(weights * e_values) / np.sum(weights)
    
    # 3. 计算加权平均误差 σ
    # σ^2 = 1 / (Σ w_i)
    sigma_weighted = np.sqrt(1.0 / np.sum(weights))
    
    return E_weighted, sigma_weighted, len(e_values)

# --- 主程序执行 ---

# 1. (可选) 执行 Sigma 筛选
cleaned_data = perform_sigma_clipping(data, SIGMA_CLIPPING_THRESHOLD)

# 2. 用筛选后的数据计算加权平均
E, sigma, n_points = calculate_weighted_mean(cleaned_data)

# --- 输出最终结果 ---
if E is not None:
    print("--- 最终计算结果 (加权平均法) ---")
    print(f"使用的数据点 (n): {n_points}")
    print(f"加权平均能量 (E): {E:.4f}")
    print(f"加权平均误差 (σ): {sigma:.4f}")