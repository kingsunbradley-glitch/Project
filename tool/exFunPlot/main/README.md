# 各模型蒸发道激发函数（Origin 风格）

`plot_all_models.py` 统一读取以下三套现有计算数据，并为每个模型分别绘制一张包含全部可用蒸发道的图：

- Havip：3n–6n；
- DNS_Li：2n–5n；
- WangECC3：2n–5n。

脚本不修改原始数据。Havip、DNS_Li 直接读取以 pb 为单位的激发能数据；WangECC3 使用用户给定的 `Q = -132.481 MeV`，按 `E* = Ecm + Q` 计算激发能，并使用 `1 mb = 10^9 pb` 换算截面。零值和空值在对数坐标中保留为曲线断点。

## 用户配置区

编辑 `plot_all_models.py` 顶部“用户配置区”即可：

- `X_LIMITS_MEV`、`Y_LIMITS_PB`：坐标范围；当前横轴为 `30`–`62` MeV，纵轴为 `1e-2`–`3` pb。
- `WANG_ECC3_Q_MEV`：WangECC3 的 Q 值，当前为 `-132.481 MeV`。
- `EXPERIMENTAL_POINTS`：红、蓝点坐标，以及左、右、下、上四个独立误差。
- `UPPER_LIMITS`：截面上限坐标及左、右能量误差；向下箭头只表示上限。
- `BARRIER_ENERGY_MEV`：默认不画势垒箭头；填入数值后启用。
- `FIGURE_SIZE_MM`、`OUTPUT_FORMATS`、`PNG_DPI`：图幅和输出设置。

由于纵轴下限为 `1e-2` pb，WangECC3 的 2n、3n 全部计算值位于显示范围以下。脚本仍读取并绘制这些数据，图例也保留相应蒸发道，但坐标框内看不到曲线。

## 运行

```bash
cd /home/evalie/Project/Function/Project/tool/exFunPlot/main
/home/evalie/miniconda3/envs/py_env/bin/python plot_all_models.py
```

只绘制一个模型：

```bash
/home/evalie/miniconda3/envs/py_env/bin/python plot_all_models.py --model wang_ecc3
```

结果保存到 `figures/`，每张图同时输出 PDF、SVG 和 600 dpi PNG。
