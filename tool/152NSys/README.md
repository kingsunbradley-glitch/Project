# Z=100--104, N=148--156 的 Qα 与 α 部分半衰期系统学图

这个小工具用于从 IAEA LiveChart 的 `ground_states` CSV API 下载评价数据，筛选：

$$
Z=100\text{--}104,\qquad N=148\text{--}156,
$$

并绘制：

1. $Q_\alpha$ vs. $N$；
2. 考虑 α 分支比后的部分半衰期

$$
T_{1/2,\alpha}=\frac{T_{1/2,\mathrm{tot}}}{b_\alpha},
$$

其中 $b_\alpha$ 使用小数形式，例如 30% 写成 `0.30`。

## 安装依赖

```bash
conda activate py_env
pip install pandas numpy matplotlib
```

## 直接自动下载并绘图

```bash
python make_alpha_systematics.py --outdir out
```

输出：

```text
out/raw_livechart_ground_states.csv
out/alpha_systematics_Z100_104_N148_156.csv
out/alpha_systematics_Z100_104_N148_156.pdf
out/alpha_systematics_Z100_104_N148_156.png
```

## 推荐的严谨流程

先运行一次：

```bash
python make_alpha_systematics.py --outdir out
```

然后检查：

```text
out/alpha_systematics_Z100_104_N148_156.csv
```

重点核查这些列：

```text
nuclide, Qalpha_keV, T12_total_s, alpha_branch, T12_alpha_s, decay_summary, note
```

如果你想用 ENSDF/NuDat/NUBASE 的人工核查值覆盖某些点，先生成模板：

```bash
python make_alpha_systematics.py --outdir out --write-template
```

编辑：

```text
out/overrides_template.csv
```

然后重画：

```bash
python make_alpha_systematics.py --outdir out --overrides out/overrides_template.csv
```

## 覆盖文件格式

`alpha_branch_override` 是小数，不是百分数：

```csv
Z,N,A,symbol,nuclide,qa_keV_override,T12_total_s_override,alpha_branch_override,alpha_branch_flag_override,note_override
102,152,254,No,254No,8100,51.0,0.90,exact,manual checked with ENSDF/NuDat
```

`alpha_branch_flag_override` 可用：

```text
exact   精确/正常评价值
approx  近似值
lower   b_alpha > value，因此 T_alpha 是上限
upper   b_alpha < value，因此 T_alpha 是下限
```

## 物理上需要注意

不要把总半衰期直接当成 α 半衰期。这个区域 SF/EC/α 竞争明显，应该使用：

$$
T_{1/2,\alpha}=T_{1/2,\mathrm{tot}}/b_\alpha.
$$

主图默认只使用 LiveChart `ground_states` 表。若要讨论异构态，建议另做一张图，或者把异构态用空心符号单独标出，不要和 ground-state 系统学混画。
