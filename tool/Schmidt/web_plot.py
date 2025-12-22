import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.ticker as ticker
import io
import json
import os

# ==========================================
# --- 1. 配置保存/加载系统 (核心升级) ---
# ==========================================
CONFIG_FILE = "plot_config.json"
default_groups = ["Present", "Dubna", "RIKEN", "GSI", "Chem"]

def load_config():
    """启动时尝试加载配置文件到 session_state"""
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                saved_config = json.load(f)
                # 将保存的配置更新到 session_state
                for key, value in saved_config.items():
                    st.session_state[key] = value
            # 标记已加载，防止重复提示
            if "config_loaded" not in st.session_state:
                st.session_state["config_loaded"] = True
        except Exception as e:
            st.error(f"加载配置文件失败: {e}")

def save_current_config():
    """将当前所有控件的值保存到本地文件"""
    config_data = {}
    # 1. 收集全局设置
    global_keys = [
        "fig_width", "fig_height", "preview_dpi", 
        "e_min", "e_max", "e_bins", "e_ymax",
        "t_min_exp", "t_max_exp", "t_bins"
    ]
    for k in global_keys:
        if k in st.session_state:
            config_data[k] = st.session_state[k]

    # 2. 收集每一组的数据和设置
    for i in range(len(default_groups)):
        group_keys = [
            f"show_{i}", f"lbl_{i}", f"c_bar_{i}", f"c_fit_{i}", f"fit_{i}",
            f"data_e_{i}", f"data_t_{i}", f"u_{i}"
        ]
        for k in group_keys:
            if k in st.session_state:
                config_data[k] = st.session_state[k]
    
    try:
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump(config_data, f, indent=4, ensure_ascii=False)
        st.toast("✅ 设置已保存！下次打开将自动加载。", icon="💾")
    except Exception as e:
        st.error(f"保存失败: {e}")

# --- 在页面渲染前加载配置 ---
load_config()

# ==========================================
# --- 2. 网页基础配置 ---
# ==========================================
st.set_page_config(page_title="Nuclear Decay Plotter V5", layout="wide")

st.title("Nuclear Decay Data Visualizer (V5 - AutoSave)")
st.markdown("支持 **配置保存**，关闭网页后数据不会丢失。")

# ==========================================
# --- 3. 侧边栏：控制与数据输入 ---
# ==========================================
st.sidebar.header("1. 全局控制")

# --- 保存/加载按钮区 ---
col_save, col_info = st.sidebar.columns([1, 1])
if col_save.button("💾 保存当前设置", type="primary", use_container_width=True):
    save_current_config()
    
# 显示加载状态提示
if st.session_state.get("config_loaded"):
    st.sidebar.success("已自动加载上次的配置")
    st.session_state["config_loaded"] = False # 只显示一次

st.sidebar.markdown("---")

# --- 状态初始化 ---
# 初始化 checkboxes 的默认状态 (如果配置文件没覆盖它们)
for i in range(len(default_groups)):
    key_name = f"show_{i}"
    if key_name not in st.session_state:
        st.session_state[key_name] = (i < 2)

# --- 全选 / 全不选 ---
col_sel1, col_sel2 = st.sidebar.columns(2)
def set_all_state(state):
    for i in range(len(default_groups)):
        st.session_state[f"show_{i}"] = state
if col_sel1.button("✅ 全选显示", use_container_width=True):
    set_all_state(True)
if col_sel2.button("⬜ 全部隐藏", use_container_width=True):
    set_all_state(False)

# --- 画布设置 ---
# 注意：所有 widget 必须添加 key 参数，以便 session_state 能够捕获并保存
with st.sidebar.expander("画布与清晰度设置", expanded=False):
    c1, c2 = st.columns(2)
    # 使用 session_state.get 设置默认值，确保加载配置后生效
    fig_width = c1.number_input("宽度 (inch)", value=st.session_state.get("fig_width", 14.0), step=0.5, key="fig_width")
    fig_height = c2.number_input("高度 (inch)", value=st.session_state.get("fig_height", 7.0), step=0.5, key="fig_height")
    preview_dpi = st.slider("预览清晰度 (DPI)", 100, 400, st.session_state.get("preview_dpi", 250), key="preview_dpi")

# --- 坐标轴设置 ---
with st.sidebar.expander("坐标轴 (Axis) 设置", expanded=False):
    st.markdown("**Energy Axis**")
    ec1, ec2 = st.columns(2)
    e_min = ec1.number_input("E Min", value=st.session_state.get("e_min", 8.0), step=0.1, key="e_min")
    e_max = ec2.number_input("E Max", value=st.session_state.get("e_max", 9.0), step=0.1, key="e_max")
    e_bins = st.slider("Energy Bins", 10, 100, st.session_state.get("e_bins", 40), key="e_bins")
    e_ymax = st.number_input("E Y-Max (0=Auto)", value=st.session_state.get("e_ymax", 25.0), step=1.0, key="e_ymax")
    
    st.markdown("---")
    st.markdown("**Time Axis (Log)**")
    tc1, tc2 = st.columns(2)
    t_min_exp = tc1.slider("10^x Min", 1, 5, st.session_state.get("t_min_exp", 2), key="t_min_exp") 
    t_max_exp = tc2.slider("10^x Max", 3, 8, st.session_state.get("t_max_exp", 5), key="t_max_exp") 
    t_bins = st.slider("Time Bins", 10, 100, st.session_state.get("t_bins", 30), key="t_bins")

# --- 数据输入区域 ---
st.sidebar.header("2. 数据输入")

def parse_input(text_input):
    if not text_input or not text_input.strip(): return np.array([])
    try:
        clean_text = text_input.replace(",", "\n").replace(";", "\n")
        return np.fromstring(clean_text, sep=' ')
    except: return np.array([])

default_colors_list = ["#A714AC", '#126782', '#2A9BC4', '#A8D5E2', '#E8F1F2']
default_fit_colors_list = ["#6A0DAD", '#0B3E4D', '#104E63', '#5D8C99', '#999999']

data_store = []

st.sidebar.caption("👇 调整好后记得点击顶部的 '保存当前设置'")

# 创建循环输入框
for i, default_name in enumerate(default_groups):
    c_check, c_exp = st.sidebar.columns([0.15, 0.85])
    
    # 1. 外部勾选框 (Key 已绑定 session_state)
    is_show = c_check.checkbox("", key=f"show_{i}")
    
    # 获取当前的标签作为 Expander 标题
    current_label = st.session_state.get(f"lbl_{i}", default_name)
    if not current_label: current_label = default_name
    
    # 2. 折叠面板
    with c_exp.expander(f"{current_label}", expanded=False):
        
        # 图例名称
        custom_label = st.text_input("图例名称", value=default_name, key=f"lbl_{i}")
        
        # 颜色与拟合
        cc1, cc2, cc3 = st.columns([1, 1, 1])
        # 使用 defaults 列表防止 crash，优先读取 session state
        def_bar = default_colors_list[i] if i < len(default_colors_list) else "#000000"
        def_fit = default_fit_colors_list[i] if i < len(default_fit_colors_list) else "#000000"
        
        bar_col = cc1.color_picker("柱色", value=st.session_state.get(f"c_bar_{i}", def_bar), key=f"c_bar_{i}")
        fit_col = cc2.color_picker("线色", value=st.session_state.get(f"c_fit_{i}", def_fit), key=f"c_fit_{i}")
        fit_chk = cc3.checkbox("拟合", value=st.session_state.get(f"fit_{i}", (i==0)), key=f"fit_{i}")

        # 数据输入 (默认值仅在第一次且无配置时显示)
        def_e_val = "8.57 8.56 8.53 8.52 8.52 8.52 8.51 8.5 8.5 8.47 8.44 8.43 8.41 8.41" if i==0 else ""
        def_t_val = "14.561 7.52 1.7 3.14 6.453 2.4 0.4 1.29 4.68 4.56 4.35 2.26 1.06 2.97 0.58 8.3 3.7 0.8 3.98" if i==0 else ""
        
        # 从 session_state 获取，如果没有则用默认值
        val_e = st.session_state.get(f"data_e_{i}", def_e_val)
        val_t = st.session_state.get(f"data_t_{i}", def_t_val)
        
        raw_e = st.text_area("Energy", value=val_e, height=60, key=f"data_e_{i}")
        raw_t = st.text_area("Time", value=val_t, height=60, key=f"data_t_{i}")
        # --- 修复后的时间单位选择逻辑 ---
        # 1. 定义选项
        u_options = [1.0, 1000.0]
        
        # 2. 获取当前状态 (可能是从json加载的浮点数 1000.0，也可能是默认索引 1)
        val_in_state = st.session_state.get(f"u_{i}", 1000.0)
        
        # 3. 智能计算 index (确保它是整数 0 或 1)
        if val_in_state in u_options:
            # 如果存的是数值 (1.0 或 1000.0)，找到它在列表里的位置
            current_index = u_options.index(val_in_state)
        elif isinstance(val_in_state, int) and val_in_state in [0, 1]:
            # 如果存的是旧版本的索引，直接用
            current_index = val_in_state
        else:
            # 默认选第2个 (s / 1000.0)
            current_index = 1

        # 4. 生成组件 (注意 index 现在一定是整数了)
        unit_mult = st.radio(
            "单位", 
            u_options, 
            index=current_index, 
            format_func=lambda x: "ms" if x==1 else "s", 
            key=f"u_{i}", 
            horizontal=True
        )
        
        data_store.append({
            "label": custom_label,
            "show": is_show,
            "fit": fit_chk,
            "bar_color": bar_col,
            "fit_color": fit_col,
            "energy": parse_input(raw_e),
            "time": parse_input(raw_t) * unit_mult
        })

# ==========================================
# --- 4. 绘图逻辑 ---
# ==========================================

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif', 'serif']
plt.rcParams['mathtext.fontset'] = 'stix'

def theoretical_curve(t, lambda_est):
    return lambda_est * t * np.exp(-lambda_est * t)

fig = plt.figure(figsize=(fig_width, fig_height), dpi=preview_dpi)
gs = GridSpec(1, 2, figure=fig, wspace=0, width_ratios=[1, 1])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)

# --- 绘图 (Energy) ---
bins_e = np.linspace(e_min, e_max, e_bins)
for d in data_store:
    if d["show"] and len(d["energy"]) > 0:
        ax1.hist(d["energy"], bins=bins_e, color=d["bar_color"], 
                 alpha=0.7, label=d["label"], edgecolor='black', histtype='stepfilled', linewidth=0.5)

ax1.set_xlim(e_min, e_max)
if e_ymax > 0: ax1.set_ylim(0, e_ymax)
ax1.set_xlabel("Energy / MeV", fontsize=14); ax1.set_ylabel("Counts", fontsize=14)
ax1.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
if any(d["show"] and len(d["energy"]) > 0 for d in data_store):
    ax1.legend(loc='upper left', frameon=False, fontsize=10)

# --- 绘图 (Time) ---
t_xmin, t_xmax = 10**t_min_exp, 10**t_max_exp
ax2.set_xscale('log'); ax2.set_xlim(t_xmin, t_xmax)
bins_t = np.logspace(t_min_exp, t_max_exp, t_bins)

for d in data_store:
    if d["show"] and len(d["time"]) > 0:
        n, _, _ = ax2.hist(d["time"], bins=bins_t, color=d["bar_color"], 
                           alpha=0.7, label=d["label"], edgecolor='black', histtype='stepfilled', linewidth=0.5)
        if d["fit"] and len(d["time"]) > 1:
            t_bar = np.mean(d["time"])
            t_axis = np.logspace(t_min_exp, t_max_exp, 500)
            y_curve = theoretical_curve(t_axis, 1.0/t_bar)
            scale = (np.max(n) / y_curve.max()) if y_curve.max() > 0 else 1
            ax2.plot(t_axis, y_curve*scale, color=d["fit_color"], linewidth=2, label=f"{d['label']} Fit")

ax2.set_xlabel("Time (ms)", fontsize=14)
if any(d["show"] and len(d["time"]) > 0 for d in data_store):
    handles, labs = ax2.get_legend_handles_labels()
    by_label = dict(zip(labs, handles))
    ax2.legend(by_label.values(), by_label.keys(), loc='upper right', frameon=False, fontsize=10)

for ax in [ax1, ax2]:
    ax.spines['top'].set_visible(True)
    ax.tick_params(direction='in', top=True, which='both')
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.tick_params(labelleft=False)
ax2.axvline(x=t_xmin, color='gray', linestyle='--', linewidth=0.8)

plt.tight_layout()

# ==========================================
# --- 5. 显示与下载 ---
# ==========================================
st.pyplot(fig)

st.write("### 📥 下载区域")
c_dl1, c_dl2, c_dl3 = st.columns(3)
img_png = io.BytesIO()
plt.savefig(img_png, format='png', dpi=300, bbox_inches='tight')
c_dl1.download_button("下载 PNG (300 DPI)", data=img_png, file_name="plot_300dpi.png", mime="image/png")
img_pdf = io.BytesIO()
plt.savefig(img_pdf, format='pdf', bbox_inches='tight')
c_dl2.download_button("下载 PDF (矢量)", data=img_pdf, file_name="plot_vector.pdf", mime="application/pdf")
img_svg = io.BytesIO()
plt.savefig(img_svg, format='svg', bbox_inches='tight')
c_dl3.download_button("下载 SVG (矢量)", data=img_svg, file_name="plot_vector.svg", mime="image/svg+xml")