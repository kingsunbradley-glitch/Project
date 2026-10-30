# Berkeley 方法

默认参数已改为用户使用的**逐级窗口内至少一次**模型。原有总窗计数保留为另一种明确命名的模型，
不能拿它替代实际的三个 `search_win`。

## 逐级窗口：本次实际使用

设置 `model="sequential"`。每个窗口均相对于前一个选定事件重新开始，
每步使用自己的能量筛选、时间窗与位置窗内背景率：

$$
\mu_j=\frac{R_j}{M}t_j,\qquad
p_j=1-e^{-\mu_j},\qquad
P_{chain\mid ER}=\prod_j p_j,\qquad n_b=N_{ER}P_{chain\mid ER}.
$$

| 参数 | 含义 |
| --- | --- |
| `step_background_counts` | 各步骤背景计数，例如 `[65110,154000,154000]`。同筛选类别在不同步骤使用相同计数，不相加或平分。 |
| `background_time_s` | 上述背景计数的共同有效采集时长，单位秒。 |
| `step_rates_hz` | 也可直接提供各步骤全区域背景率；与计数输入二选一，使用率时不设置 `background_time_s`。 |
| `step_windows_s` | 各步骤窗口上限，如 `[0.1,100,100]`，不是相对于首个 ER 的累计窗口。 |
| `position_bins` | 有效等价位置单元数，默认 1。输入单位置率时设为 1。 |
| `implantations` | 满足入口条件的 ER 数量。 |
| `observed_chains` | 用于附加泊松尾概率的观测链数，默认 1。 |

这正是用户后续截图的乘积公式。alpha2 与 alpha3 同筛选时，
分别计算它们在各自重新开始的 100 s 内至少一次的概率，得到相同因子的乘积；
不额外除以 2!。阶乘适用于**同一个固定总窗中**的计数质量函数。
逐级模型假设平稳泊松背景并在选定事件后重新计时；按首个合格事件推进时，
可利用无记忆性逐步相乘。复杂回溯配对、时间/位置非均匀背景应另与实际搜索算法核对。

运行默认配置：`python calculate.py`；也可 `python calculate.py --config input_sequential.json`。
本次 T=655.2 h、M=6144、ER=13250000，输出 n_b=6.715226456e−6。
`P_per_step` 列出每步概率；`P_per_implantation` 是它们的乘积。
`P_at_least_observed_chains` 另按整次随机链数的泊松模型计算，默认约 6.715203909e−6。

## 固定总窗计数：另一种定义

设置 `model="total_window"`，或沿用只含原有总窗字段的旧配置。
先计算单次注入后、同一位置区和固定总窗内的 alpha 本底期望数，再计算泊松概率：

$$
\mu=\frac{R_\alpha}{N_{pixel}}\Delta t_{max},\qquad
P_{exact}=\frac{\mu^{n_\alpha}}{n_\alpha!}e^{-\mu},\qquad
n_b=N_{ER}P_{exact}.
$$

截图式 (4.29) 分母印为 n_alpha；本程序使用 **n_alpha!**。
该形式也见 Berkeley 原始资料附录 A 的式 (A2)/(A3)：
[LBNL 原始论文档案](https://escholarship.org/uc/item/1w75714h)。

| 参数 | 含义 |
| --- | --- |
| `alpha_rate_hz` | 所选能量范围内，整个分析区域的 alpha 总背景率 R_alpha。 |
| `position_bins` | 区域中等效均匀位置单元数；已输入单个位置的率并只分析该位置时设 1。 |
| `window_s` | 相对注入时刻的整个观察窗，非每步的相邻时间间隔。 |
| `n_alpha` | 需要的 alpha 数，正整数。 |
| `implantations` | 满足入口筛选的注入事件总数。 |
| `count_mode` | 默认 `exact`，复现上述公式；`at_least` 则对泊松分布求 N≥n_alpha 的上尾。 |
| `observed_chains` | 可选，默认 1；用于由 n_b 估计整个搜索的泊松尾概率。 |

主要输出为 `P_per_implantation`（一次注入的概率）、`n_b`（期望随机链数）。
`P_at_least_observed_chains` 是附加的整个搜索概率估计，另假设随机链数服从泊松分布。
重叠的注入窗会影响独立性，不能把这一步视为无条件精确结论。

`input_total_window.json` 示例：R_alpha=0.01 Hz，384 个位置区，总窗 60 s，3 个 alpha，10000 次注入，
得到 mu=0.0015625、单次概率约 6.35e−10、n_b 约 6.35e−6。
这些是格式演示参数。

总时间窗可按核素的预计寿命与搜索设计指定；截图提到最长半衰期的 5–6 倍，
程序不自动替用户选择这个倍数。此计数模型不额外要求每步能量顺序、逐步寿命或 SF 条件；
如果实际筛选包含这些条件，应使用相应细分模型或 Monte Carlo。

运行旧总窗示例：本目录 `python calculate.py --config input_total_window.json`，支持 `--json`。
推荐解释器为 `/home/evalie/miniconda3/envs/py_env/bin/python`。
