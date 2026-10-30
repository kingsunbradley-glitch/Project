# Dubna / Zlokazov 方法

依据截图式 (4.32)–(4.41) 及 [Zlokazov, Statistical analysis of rare events—synthesis of the element 114, EPJ A 8, 81–86 (2000)](https://doi.org/10.1007/s100500070121)。
已核对用户本地同名 PDF。配置中的 `counts`、`time`、`decay`、`position` 四块可独立使用，至少保留一块。
输入时间全部为秒，率全部为 Hz。不同检验针对不同假设，程序分别输出，不自动相乘或给链贴真假标签。

## 分类计数 counts

固定窗 `window_s`；`categories` 中每类填写 `name`、`rate_hz`、`count`：

$$
Q_{K_i}=\frac{(l_i t)^{K_i}}{K_i!}e^{-l_i t},\qquad
P_{configuration}=\prod_i Q_{K_i}.
$$

率是**该链位置门内**的率；程序不再除以像素数。各类要求独立、互不重叠。
输出为恰好出现这些计数的概率质量，不是“至少这么多”的尾概率，不能单凭概率质量很小就赋予显著性。
若以已有 SF 为锚，按其前方窗口检验 alpha 和注入本底，则条件化在 SF 上，不再乘 SF 计数概率。
若要无条件检验全部类别，需自行把 SF 也列入类别并定义曝光。

模板按原文：34 min 窗，alpha 率 1/h、恰好 3 个，注入率 1.3/h、恰好 1 个，
输出 `P_exact_configuration=0.0060683733693`，复现截图约 0.00607。

## 时间 time

填写 `background_rate_hz` 与 `event_times_s`（递增时间）或 `gaps_s`（相邻间隔），二选一。
令 m 为**间隔数**，对 m 个独立指数样本，输出

$$
E[T_{min}]=\frac{1}{ml},\qquad
E[T_{max}]=\frac{H_m}{l},\qquad H_m=\sum_{j=1}^m\frac1j.
$$

五个成员对应四个间隔。模板复现原文：本底平均间隔 26.09 min、链均值 8.5 min，
本底最小/最大间隔期望约 6.52/54.35 min，链拟合值约 2.125/17.708 min。
这些值不是置信区间；程序不会把二者是否重叠当作一个严格显著性检验。
另输出 `P_mean_gap_le_observed_if_iid_exponential`，用 Gamma(m, rate=l) 的间隔总和分布计算。
该尾概率要求未被时间窗/终点条件截断的固定 m 个独立指数间隔；预先筛选过的链一般不满足此条件。

## 串联衰变 decay

填写 `time_s`，以及 `half_lives_s` 或 `half_life_bounds_s` 二选一。
按独立串联衰变等待时间和的 CDF 计算截止 t 已完成所有步骤的概率。
采用吸收态连续时间马尔可夫链的矩阵指数，避免衰变常数相等或接近时解析分母的奇点。
单步时退化为 1−exp(−lambda*t)，相等速率时退化为 Erlang CDF。

`half_life_bounds_s` 每项为 `[正下限,正上限]`。CDF 随每个半衰期增大而减小，
故所有上限组成的角点给最小值，所有下限角点给最大值，无需随机遍历参数范围。
`detection_efficiency` 是可选的**整链总效率**，默认 1；结果另给乘以它的 `detected_probabilities`。
不自动补分支比、漏失衰变、单步效率或核素产额等额外权重。

模板的半衰期区间来自原文，但只计算清楚定义的串联完成 CDF。
原文还包含特定探测效率与不同解释的附加因素；本程序不声称直接复现文中最终的 0.0083–0.3364。
截图等速率的卷积式有指数符号问题，程序使用正确的负指数形式。

## 位置 position

输入 `standardized_differences`，即各位置差除以其**位置差标准差**；
例如 `[ (p_ER-p_alpha1)/sigma_ERalpha, ... ]`。这些除法需先算成 JSON 数值。
输入对应的 `correlation_matrix`；其对角必须为 1，对称且半正定。
`samples` 默认 100000、`seed` 默认固定、`confidence` 默认 0.67。

$$
S=\sum_i x_i^2,\quad E[S]=n,\quad Var(S)=2\operatorname{tr}(C^2).
$$

相关标准高斯差分下，S 与 sum(eigenvalue(C)_i × Z_i²) 同分布。
数值模拟这一分布，输出等尾中心参考区间、S 的上尾概率及模拟误差区间。
不能直接使用独立 chi-square(n)。邻接差分相关 −0.5、非邻接为 0 时，方差为 3n−1，
模板 n=4 得到 11。

相关矩阵必须匹配实际选取的位置差。若位置测量相互独立、误差方差为 v_j，
差分矩阵为 A，则差分协方差为 A diag(v_j) A^T，再按各差分标准差归一化得到相关矩阵。
截图式 (4.41) 的最后一项是 ER−SF，因此不能直接套用全相邻差分的三对角相关矩阵。
模板的位置差和三对角矩阵仅作格式演示，不是截图的 S=9.56 数据。

原文给出的 67% 区间 (1.8,10.8) 不是这里的等尾中心分位数区间；其表格本身对应的
等尾 67% 区间约为 (1.2,6.6)。程序采用明确的分位数定义，不把历史经验界限当作精确分位数。
位置相容性以同一位置的真实衰变链假设为条件，并非随机背景形成假链的概率。

运行：本目录 `python calculate.py`，支持 `--config`、`--json`。
推荐解释器为 `/home/evalie/miniconda3/envs/py_env/bin/python`。
