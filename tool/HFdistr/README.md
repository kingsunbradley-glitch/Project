# HF 核素图与交互网页

本目录以 `ackermann_ground_states.csv` 为数据库，调用
`method/all/run_all.py` 所登记的 7 种 α 衰变半衰期方法，统一计算阻碍因子
（hindrance factor, HF），并生成可离线打开的交互式核素图网页。

## HF 定义

计算严格采用数据库给出的实验 α **部分半衰期**：

\[
HF = \frac{T_{\alpha,\mathrm{partial}}^{\mathrm{exp}}}
          {T_{\alpha}^{\mathrm{calc}}}
\]

- 分子直接读取 `T_alpha_half_s`。
- `T_half_s`（总半衰期）和 `b_alpha_percent` 只用于数据质量核对，绝不用于
  重建或替换 `T_alpha_half_s`。
- 理论计算直接使用数据库的 `Qalpha_MeV`，不做
  \(E_\alpha\rightarrow Q_\alpha\) 二次反冲修正。
- 对只提供实验 `Ealpha_MeV` 的新增状态，写入数据库前按本项目方法约定
  \(Q_\alpha=E_\alpha A/(A-4)\) 计算并保存 `Qalpha_MeV`；运行计算时仍只读取
  已保存的 `Qalpha_MeV`。

## 生成网页与结果

从仓库根目录运行：

```bash
python3 tool/HFdistr/build_hf_dashboard.py
```

主要输出位于 `tool/HFdistr/results/`：

- `hf_nuclide_chart.html`：自包含交互网页，可直接双击打开；无需网络或 CDN。
- `hf_results_long.csv`：每条核素记录 × 每种方法的长格式结果。
- `hf_results_wide.csv`：便于横向比较 7 种方法的宽格式结果。
- `hf_data_quality.csv`：全部数据库记录及质量标记。
- `hf_excluded_records.csv`：因缺少正的 `Qalpha_MeV` 或
  `T_alpha_half_s` 而无法定义 HF 的记录。
- `hf_method_summary.csv`：各方法 HF 范围和统计量。
- `hf_consistency_audit.csv`：总/部分半衰期、分支比、Eα/Qα、单位与误差字段的
  全库一致性审计结果。
- `hf_lt_0p1_source_audit.csv`：Viola–Seaborg 与 E. Rurarz 数值结果中
  `HF<0.1` 候选的逐核素来源核对、修订动作和修订后结果。
- `hf_extreme_source_audit.csv`：Viola–Seaborg 或 E. Rurarz 中初始满足
  `HF<1` 或 `HF>100` 的全部 50 条核素/状态记录及修订后结果。
- `261Sg_source_audit.csv`：`261Sg` 半衰期单位、α 跃迁能量、部分半衰期及
  Viola–Seaborg/E. Rurarz 重算结果的专项来源核对。
- `paper_lookup_raw.json`：本轮原始论文检索使用的 OpenAlex/Crossref 原始响应。
- `261Sg_paper_lookup_raw.json`：`261Sg` 原始衰变研究论文的 Crossref/OpenAlex
  检索响应。
- `hf_extreme_paper_lookup_raw.json`：极端 HF 全表核对使用的 NUBASE2020 与
  关键新实验论文的 Crossref/OpenAlex 原始响应。
- `hf_nuclide_chart_overview.png/.pdf`：7 种方法的静态总览。
- `run_metadata.json`：输入、HF 定义、方法、假设和色标元数据。

如只需网页和 CSV、跳过静态图：

```bash
python3 tool/HFdistr/build_hf_dashboard.py --no-static-overview
```

## 网页说明

- 横轴为中子数 \(N\)，纵轴为质子数 \(Z\)。
- 颜色直接按 HF 离散分档：`HF<0.1`、`0.1≤HF<4`、`4≤HF<10`、
  `10≤HF≤100` 和 `HF>100`；同一档内颜色完全相同。
- 7 种方法共享同一套固定颜色，格内文字和悬停框继续显示未截断的原始 HF。
- 悬停、键盘聚焦或点击核素格可查看 `Qalpha_MeV`、实验 α 部分半衰期、
  理论半衰期、HF、质量标记及模型状态。
- 灰色“—”表示数据库中存在该核素，但没有足够数据定义 HF。
- 虚线边框和 `*` 表示使用估算的 Qα；HF 的 `≈`、`<`、`>`、`≤`、`≥`
  直接继承 α 部分半衰期限定符。
- `261Rf^a/b`、`265Sg^a/b`、`273Ds`/`273Ds^m`，以及 `262Bh` 的两条记录，
  均在同一 N–Z 格内拆分显示，不做平均。
- 状态分辨数据额外保留 `state_label`、`Ealpha_MeV`、`alpha_daughter` 和
  `alpha_daughter_population_raw`，网页悬停框会显示这些衰变链信息。

## 方法与假设

默认方法集为：Poenaru、DZR、Ismail 2022 Formula E、Qi 2009 UDL、
Viola–Seaborg、Xu 2022 Unified 和 E. Rurarz/Poenaru target。

数据库没有 α 跃迁角动量 `l`。DZR、Ismail 和 Xu 默认采用 `l=0`，并在网页
与结果元数据中明确标记为假设。Xu 方法中当前数据库的全部核素均位于本地
脚本 `C(Z,N)` 的定义区之外，因此结果沿用脚本规则 `C=0`，属于域外外推。

## 当前数据库覆盖

默认数据库共 210 条记录。176 条记录具有正的 `Qalpha_MeV` 和
`T_alpha_half_s`，对应 172 个唯一 N–Z 格点；其余 34 条仍显示在网页中，但
HF 标为不可计算。7 种方法合计生成 1232 条 HF 结果。

当前 `273Ds^a` 使用论文表中的 `T_alpha_half_s=17.34 ms`、
`Ealpha_MeV=10.89 MeV`，对应 `Qalpha_MeV=11.051933... MeV`。原始字符串读取
失败的 `264Lr` 也已修正为 `T_half_s=T_SF_half_s=17280 s`。

针对 Viola–Seaborg 与 E. Rurarz 的 `HF<0.1` 候选，已根据 ENSDF、
NUBASE2020 和原始论文修正 `251Es`、`255Es`、`245Fm`、`256Fm`、`258Rf`、
`260Rf`、`290Mc`。其中包括衰变分支次序、秒/毫秒单位和未观测 α 分支的
评估标记。修订后仍小于 0.1 的数值均保留 `<`/`>` 限定符，不能当成无条件的
精确 HF。

`261Sg` 也已专项复核。Ackermann 汇总表中的 `183(5) s` 是单位误植，采用
ENSDF 的 `T1/2=184(5) ms`、`BRalpha=98.1(4)%`，得到 α 部分半衰期
`0.1875637 s`。HF 对应实际布居 `257Rf` 基态的 α 跃迁：由
`Ealpha(lab)=9.410(10) MeV` 反冲修正得到 `Qtransition=9.556459 MeV`；质量表
基态差 `Qmass=9.714(15) MeV` 仅保留在原始说明中，不用于该跃迁的 HF。

对 Viola–Seaborg 或 E. Rurarz 中初始 `HF<1` 或 `HF>100` 的 50 条记录也已
完成逐项来源核对。确认并修正：`248Fm` 的 α 分支小数/百分数误读，`252Fm`
的 Qα 数字次序及主 SF 分支选择，`286Fl` 的秒/毫秒单位，`246Fm` 的 SF/EC
分支次序，以及 `254Fm` 的主 SF 分支选择。修订后剩余 48 条极端记录，其中
10 条使用估算 Qα、15 条的 α 部分半衰期带上下限或近似限定。`281Rg` 的
14/88 中心值合计 102%，因来源不足以唯一裁决而保留为显式警告。
