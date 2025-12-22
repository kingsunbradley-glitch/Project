import matplotlib.pyplot as plt

# --- 1. 定义能级坐标和标签 ---
# (x_start, y)
levels = {
    'Rf_b': {'pos': (0, 0), 'label': '$^{261}$Rf', 'side_label': 'b'},
    'Sg_b': {'pos': (3.5, 3), 'label': '$^{265}$Sg', 'side_label': 'b'},
    'Sg_a': {'pos': (3.5, 4.5), 'label': '', 'side_label': 'a'},
    'Hs_g': {'pos': (7, 7), 'label': '$^{269}$Hs', 'side_label': ''},
    'Ds_a': {'pos': (10.5, 10), 'label': '$^{273}$Ds', 'side_label': 'a'},
    'Ds_b': {'pos': (10.5, 11.5), 'label': '', 'side_label': 'b'},
}

# 每条能级线的水平长度
level_len = 2.0

# --- 2. 创建画布 ---
fig, ax = plt.subplots(figsize=(10, 8))

# --- 3. 绘制能级线和标签 ---
for key, data in levels.items():
    x, y = data['pos']
    label = data['label']
    side_label = data['side_label']
    
    # 绘制能级线
    ax.plot([x, x + level_len], [y, y], color='black', linewidth=3)
    
    # 绘制核素标签 (在右侧)
    if label:
        ax.text(x + level_len + 0.2, y, label, 
                fontsize=14, color='darkblue', va='center')
        
    # 绘制 'a'/'b' 标签 (在左侧)
    if side_label:
        ax.text(x - 0.2, y, side_label, 
                fontsize=14, color='darkblue', ha='right', va='center')

# --- 4. 绘制衰变箭头 ---
# 我们使用 ax.annotate 来绘制箭头

# 定义箭头连接点 (x, y)
points = {
    'Rf_b_in': (levels['Rf_b']['pos'][0] + 1.8, levels['Rf_b']['pos'][1]),
    'Rf_b_out': (levels['Rf_b']['pos'][0] + 0.2, levels['Rf_b']['pos'][1]),
    'Sg_b_in': (levels['Sg_b']['pos'][0] + 1.8, levels['Sg_b']['pos'][1]),
    'Sg_b_out': (levels['Sg_b']['pos'][0] + 0.2, levels['Sg_b']['pos'][1]),
    'Sg_a_in': (levels['Sg_a']['pos'][0] + 1.8, levels['Sg_a']['pos'][1]),
    'Sg_a_mid': (levels['Sg_a']['pos'][0] + 0.8, levels['Sg_a']['pos'][1]),
    'Sg_b_mid': (levels['Sg_b']['pos'][0] + 0.8, levels['Sg_b']['pos'][1]),
    'Hs_g_in1': (levels['Hs_g']['pos'][0] + 1.8, levels['Hs_g']['pos'][1]),
    'Hs_g_in2': (levels['Hs_g']['pos'][0] + 1.0, levels['Hs_g']['pos'][1]),
    'Hs_g_out': (levels['Hs_g']['pos'][0] + 0.2, levels['Hs_g']['pos'][1]),
    'Ds_a_out': (levels['Ds_a']['pos'][0] + 0.2, levels['Ds_a']['pos'][1]),
    'Ds_b_out': (levels['Ds_b']['pos'][0] + 0.2, levels['Ds_b']['pos'][1]),
}

# 箭头样式
arrow_solid_red = dict(arrowstyle="->", color='red', linewidth=2)
arrow_dashed_red = dict(arrowstyle="->", color='red', linewidth=2, linestyle='--')
arrow_dashed_blue = dict(arrowstyle="->", color='darkblue', linewidth=2, linestyle='--')
arrow_solid_sf = dict(arrowstyle="->", color='darkcyan', linewidth=2)

# Ds(a) -> Hs(g)
ax.annotate("", xy=points['Hs_g_in1'], xytext=points['Ds_a_out'], arrowprops=arrow_solid_red)

# Ds(b) -> Hs(g)
ax.annotate("", xy=points['Hs_g_in2'], xytext=points['Ds_b_out'], arrowprops=arrow_dashed_red)

# Hs(g) -> Sg(a)
ax.annotate("", xy=points['Sg_a_in'], xytext=points['Hs_g_out'], arrowprops=arrow_solid_red)

# Sg(b) -> Rf(b)
ax.annotate("", xy=points['Rf_b_in'], xytext=points['Sg_b_out'], arrowprops=arrow_solid_red)

# Sg(a) -> Sg(b) [Internal]
ax.annotate("", xy=points['Sg_b_mid'], xytext=points['Sg_a_mid'], arrowprops=arrow_dashed_blue)

# Rf(b) -> SF
sf_target = (points['Rf_b_out'][0] - 3, points['Rf_b_out'][1] + 3)
ax.annotate("", xy=sf_target, xytext=points['Rf_b_out'], arrowprops=arrow_solid_sf)
ax.text(sf_target[0] - 0.2, sf_target[1], 'SF', 
        fontsize=14, color='darkcyan', ha='right', weight='bold')


# --- 5. 清理和显示 ---
ax.set_ylim(-1, 13)
ax.set_xlim(-2, 14)
ax.axis('off')  # 关闭坐标轴
plt.show()