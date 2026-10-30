# 偶然符合概率：Schmidt (1984)

依据 K.-H. Schmidt, C.-C. Sahm, K. Pielenz, H.-G. Clerc,
*Some Remarks on the Error Analysis in the Case of Poor Statistics*,
Z. Phys. A **316**, 19–26 (1984)，第 2、3、8 节。
[论文 DOI / 出版社页面](https://doi.org/10.1007/BF01415656)。
公式已按用户提供的 PDF 第 20–21 页逐项核对。

本工具从背景计数率、测量时长、关联时间窗和观测链数计算随机符合概率；
不读取 ROOT 事件树，也不自动筛选衰变链。实现式 (2)、(3)、(5)–(8)，
不包含论文后半部分的能量/寿命相容性检验。

## 直接使用

依赖 Python 3 和 SciPy，现有 `py_env` 已满足。终端中运行：

```bash
cd /home/evalie/Project/Function/Project/tool/accidental_coincidence
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py
```

默认读取脚本同目录的 `input.json`。**当前文件装的是论文 Th-214、5 秒窗示例，
不是用户的实验参数**；实际使用时修改该文件，或指定另一个 JSON 文件：

```bash
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py --config input.json
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py --config input.json --json
```

脚本只向终端输出，不自动覆盖任何数据或结果文件。

## 参数

全部时间用秒，计数率用 Hz（s⁻¹）。例如 5 ms 写成 `0.005`，250 h 写成 `900000`。

| 参数 | 含义 |
| --- | --- |
| `mode` | `fixed`：固定顺序，式 (5)；`partial`：仅首事件固定，式 (7)。默认 `fixed`。 |
| `measurement_time_s` | 本次搜索随机链的有效测量时长 $T>0$，不是背景样本时长。 |
| `rates_hz` | 各类事件的背景率数组 $[R_1,\ldots,R_K]$，与 `background_counts` 二选一。 |
| `background_counts` | 各类事件的非负整数计数，利用下面的背景采样时长换算成率。 |
| `background_time_s` | 背景样本时长：一个数（各类共用），或长度为 $K$ 的数组。省略时用 $T$；只与 `background_counts` 配合使用。 |
| `windows_s` | 长度为 $K-1$ 的非负时间窗数组，定义见下文。 |
| `observed_chains` | 满足同一组筛选条件的实际完整链数 $n_m$，默认 1。 |
| `position_bins` | 等效、均匀、独立的位置分区数 $M$，默认 1；见位置符合说明。 |
| `sequence` | 可选，链步骤对应计数率数组的下标（从 0 编号），例如 `[0,1,2,2]`。默认各类别顺序各出现一次；重复类别仅支持 `fixed`。 |

计数率应已应用相应的能量、粒子类型等筛选条件。这里的 λ 是背景事件率，
不是核素的衰变常数 `ln(2)/T_half`。默认不传播背景率估计的不确定度。

例如，搜索 `ER → α → SF`，依次输入三个互不重叠事件类别的背景率：

```json
{
  "mode": "fixed",
  "measurement_time_s": 900000,
  "rates_hz": [0.00032, 0.00003472222222, 0.01333333333333],
  "windows_s": [0.005, 12.9],
  "observed_chains": 1,
  "position_bins": 50
}
```

以上仅是格式演示，请替换成实验参数。

### 时间窗与事件顺序

- `fixed`：事件必须为 $E_1\to E_2\to\cdots\to E_K$，时间窗为
  $[\Delta t_{12},\Delta t_{23},\ldots,\Delta t_{K-1,K}]$。原文式 (4) 要求相邻链成员之间
  没有这 $K$ 类中的其他事件；允许插入其他候选事件的“所有组合”配对不属于这个模型。
- `partial`：固定 $E_1$ 为首事件，后续类别顺序可变；时间窗为
  $[\Delta t_{12},\Delta t_{13},\ldots,\Delta t_{1K}]$，**全部相对于首事件**。
  按原文式 (7) 的乘积形式实现，每个因子含 $(\lambda_1+\lambda_i)$。
  这不是任意宽窗、多次重复事件、共享起点条件下的通用精确联合概率；
  应按论文所讨论的稀有链情形使用。

例如相邻实测时间间隔为 `[0.005, 0.0273, 12.9]`，用于 `partial` 时对应首事件的
时间上限应是 `[0.005, 0.0323, 12.9323]`；不能原样照搬。
内置 109 号元素示例直接采用论文式 (21) 列出的三个首事件时间窗。
窗口的下限均为零；有死时间/非零下限的窗口需要另行处理。

同一种事件类别重复出现（例如共用同一能窗的两个 α）时，计数率数组只列出这个类别一次，
通过 `sequence` 重复引用。例如 ER-alpha1-alpha23-alpha23 使用三个背景计数率和
`sequence: [0,1,2,2]`，三个相邻时间窗。每步所需的事件率可以重复使用，
但指数项中的总率只对独立背景类别求和一次。这是在原文相邻事件、无插入其他类别事件模型下
利用泊松过程无记忆性对重复类别的推广，链成员仍是两个不同的事件。
不同能窗互不重叠、且可视为独立背景过程时，才可分别作为不同类别输入。

## 计算公式与输出

先考虑一个位置区，背景率为 $\lambda_i$，总率 $\Lambda=\sum_i\lambda_i$。

固定顺序，式 (5)：

\[
n_b=T\lambda_1\prod_{i=1}^{K-1}
\frac{\lambda_{i+1}}{\Lambda}
\left(1-e^{-\Lambda\Delta t_{i,i+1}}\right).
\]

部分自由顺序，式 (7)：

\[
n_b=T\lambda_1\prod_{i=2}^{K}
\frac{\lambda_i}{\lambda_1+\lambda_i}
\left(1-e^{-(\lambda_1+\lambda_i)\Delta t_{1,i}}\right).
\]

低计数率展开，式 (6)/(8)：

\[
n_b\simeq T\left(\prod_{i=1}^{K}\lambda_i\right)
\left(\prod \Delta t\right).
\]

最终按式 (2) 计算整个测量期间出现至少 $n_m$ 条随机链的概率：

\[
P_{\rm err}=P(N\ge n_m\mid N\sim\operatorname{Poisson}(n_b))
=\sum_{n=n_m}^{\infty}\frac{n_b^n e^{-n_b}}{n!}.
\]

这是给定随机背景模型的尾概率，不是“这条链为假”的后验概率，也不是单个起始事件的
条件符合概率。$n_b$ 是期望链数，可以大于 1，不能一般性地将它当作概率。
当 $n_m=1$ 时，$P_{\rm err}=1-e^{-n_b}$；只有 $n_b\ll1$ 才有 $P_{\rm err}\simeq n_b$。
式 (3) 的一般近似为 $P_{\rm err}\simeq n_b^{n_m}/n_m!$，同样要求 $n_b\ll1$。

| 输出 | 含义 |
| --- | --- |
| `n_b` | 原文式 (5)/(7) 的随机链期望数（保留指数项）。 |
| `P_err` | 上述 $n_b$ 代入式 (2) 的泊松上尾概率。 |
| `log10_P_err` | 概率的十进制对数，保留极小概率的信息。 |
| `n_b_low_rate` | 式 (6)/(8) 的低计数率近似期望数。 |
| `P_err_low_rate` | 低计数率近似 $n_b$ 代入式 (2)，不是直接使用式 (3)。 |
| `log10_P_eq3` | 保留指数项的 $n_b$ 代入式 (3) 后的对数，仅供小 $n_b$ 时参考。 |
| `occupancy` | 各时间窗的 $\Lambda\Delta t$ 或 $(\lambda_1+\lambda_i)\Delta t$。 |

通过 `expm1` 和泊松 survival function 避免极小概率被 `1-CDF` 的数值消减抹成零。
浮点下溢时仍计算 `log10_P_err`。严格为零的量，其对数在 JSON 中用字符串 `"-inf"` 表示。
`max(λΔt)>=0.1`、`n_b>=0.1` 的提示阈值只是数值使用提醒，不是论文规定的显著性阈值。

### 位置符合

`rates_hz` 和 `background_counts` 始终覆盖所分析的整个区域。
若将该区域视为 $M$ 个均匀、独立的有效位置区，并且要求链中所有事件在同一区，
脚本先令每区 $\lambda_i=R_i/M$，再对各区求和：

\[
n_b^{\rm total}=M\,n_b^{\rm one\ bin},\qquad
n_b^{\rm total,low}\simeq
\frac{T\prod_i R_i\prod\Delta t}{M^{K-1}}.
\]

这里不会再把起始事件总数除以 $M$。已输入单条带/单像素背景率且只分析该位置时，
设 `position_bins=1`。非均匀探测器应逐区计算 $n_{b,j}$，然后将其求和，再对总期望数
求泊松尾概率；不要对各区 p 值求和。允许相邻条带的重叠位置窗也不能简单当作独立分区。

### 适用范围

采用原文的平稳、独立泊松背景假设，以及以泊松分布描述随机链数的处理。
保留指数项只是相对于式 (6)/(8) 不做低率展开，并不使任意事件关联算法的链计数严格服从泊松分布。
共享事件、大量重叠窗口或高占有率时，需要针对实际筛选算法做时间打乱或蒙特卡洛检验。
本程序未修正测量起止边界、束流开关、采集间断和死时间，默认窗口远小于有效测量区段。
背景观测为零时输出的是零率代入结果，不能据此声称真实偶然概率为零。
扫描多个能窗、时间窗、位置条件后挑最小概率，还需考虑搜索次数的影响；此处输出的是给定筛选的局部概率。

## 原文示例复算

```bash
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py --example th214-5s
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py --example th214-40s
/home/evalie/miniconda3/envs/py_env/bin/python calculate.py --example element109
```

| 论文例子 | 本程序 $n_b$ | 本程序 $P_{\rm err}$ | 低率 $n_b$ 代入式 (2) |
| --- | --- | --- | --- |
| Th-214，$T=6080$ s，$n_1=28,n_2=38$，5 s 窗、12 条链 | 0.8516779582 | $1.3879261702\times10^{-10}$ | $1.8786903400\times10^{-10}$ |
| 同上，40 s 窗、18 条链 | 5.6782799581 | $2.8451964606\times10^{-5}$ | $3.6178431660\times10^{-4}$ |
| 109 号元素，250 h，50 个位置区、1 条链 | $2.0868390726\times10^{-18}$ | $2.0868390726\times10^{-18}$ | $2.0869333333\times10^{-18}$ |

前两行满足原文第 24 页的 $P<10^{-9}$、$P<10^{-3}$；原文这里使用式 (6) 的近似。
第三行复现第 25 页约 $2\times10^{-18}$ 的估计，其中后续三类事件总区域率分别是
`1/(8 h)`、`1/(75 s)`、`1/(250 h)`，首类是 `288/(250 h)`。
原文额外要求裂变为最后事件，而式 (7) 只固定首事件，所以论文将此估计作为上限。

## 验证

```bash
/home/evalie/miniconda3/envs/py_env/bin/python -m unittest discover -s . -p 'test_*.py' -v
```

检查论文例子、两类事件的独立积分、位置分区求和一致性、计数/率输入等价、极小尾概率、零值和无效输入。
