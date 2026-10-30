# Ds 激发函数与实验数据

## 绘图约定

使用 Python / matplotlib，按 nature-figure 的流程生成 183 × 92 mm 的三面板图。
图的科学问题是：各蒸发道的激发函数如何随激发能变化，4n、5n 计算值是否与给定实验约束相容？
面板 a 展示全部四条计算曲线的完整能区及所有正截面；面板 b、c 分别对照 4n、5n 实验数据，
采用相同的局部坐标范围，避免同一能量处不同道的符号相互遮挡。
这种分工属于定量比较图，不能仅凭相容性推断唯一核素或核态指认。

- 输入计算数据：`ExFunCalc.csv`，保持原样。用户确认截面单位为 pb。
- 实验数据：从提供的表格逐项录入 `experimental_data.csv`，共 6 条道/能量记录。
- 用户确认复合核为 278Ds；275Ds、274Ds、273Ds、272Ds 分别对应 3n、4n、5n、6n。
- 横轴直接使用两套数据各自给出的 E*，单位 MeV，不转换或平移能量。
  计算 CSV 中 E*=50 MeV 对应 Elab=214.227 MeV，而实验表中对应 225 MeV；
  在两者能量定义核实前，不绘制共享 Elab 副轴，也不据此修改输入数据。
- 计算点依次用线段连接，不拟合、不平滑、不添加理论不确定度。
- 对数轴中零截面被掩蔽，保留断点，不替换成任意小正数；完整统计记录在 QA 数据报告。
- 实验上限以横帽和向下箭头表示，箭头长度只用于显示，不代表测量误差。
- 5n 在 E*=50 MeV 的测量值为 0.14 (+0.18/-0.08) pb，误差棒端点为 0.06、0.32 pb。
  图表未说明误差或上限的置信水平，因此不将其标成 1σ 或 95% CI；没有给出横向误差。
- 同一实验点的两组靶厚和束流剂量用 `|` 保留，不当作两个独立截面测量。

## 运行环境

优先复用已有 Conda 环境 `py_env`：

```bash
cd /home/evalie/Project/Function/Project/tool/exFunPlot
/home/evalie/miniconda3/envs/py_env/bin/python plot_excitation_functions.py
```

绘图需要 numpy、matplotlib。nature-figure 的 PDF 碰撞检查另需 PyMuPDF。
当前检查库单独位于临时目录，因此本次完整复现命令为：

```bash
PYTHONPATH=/tmp/exfunplot-qa-deps /home/evalie/miniconda3/envs/py_env/bin/python plot_excitation_functions.py
```

临时目录清理后，可在自行管理的依赖位置安装 `PyMuPDF>=1.24,<2`；不需要重装绘图库。
脚本通过 `CODEX_HOME` 或默认 Codex 目录定位 nature-figure，也可用 `--skill-dir` 指定。

## 输出和检查

`figures/excitation_functions.{pdf,svg,png}` 为最终图；PDF/SVG 保留文本，PNG 为 600 dpi。
`figures/qa/` 保存源代码检查、对齐、PDF 字号、碰撞和输入数据核对结果。
面板 a 全范围为 30–80 MeV，面板 b、c 放大至 30–65 MeV；局部坐标裁切不改变输入数据。
后两幅图的对数纵轴下界为 0.001 pb，低于该显示范围的理论尾部仍可在面板 a 查阅。

## 图注

**Excitation functions for the 278Ds compound system.**
**a**, Calculated cross sections for the 3n–6n evaporation channels, connected at the supplied energy grid.
**b,c**, Calculated 4n and 5n cross sections compared with the experimental values from the supplied table.
Capped downward arrows denote upper limits; the filled circle and asymmetric error bar show
the measured 5n cross section at E*=50 MeV. The confidence levels of the limits and errors
are not specified in the supplied table. Each dataset uses its supplied excitation energies;
no common laboratory-energy conversion is applied. Zero calculated cross sections are not
shown on the logarithmic axes. Panels b and c use matching zoomed axes.

技能来源：[Yuan1z0825/nature-skills](https://github.com/Yuan1z0825/nature-skills)，
使用 `nature-figure` 及其支持文件；本图使用其排版和检查流程，不宣称已满足某刊的全部投稿要求。
