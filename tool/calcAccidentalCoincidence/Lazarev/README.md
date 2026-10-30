# Lazarev 方法

实现截图式 (4.28)：对 `ER → alpha1 → … → alpha(K-1)`，

$$
n_b=C_{ER,\alpha_1}\prod_{i=1}^{K-2}\frac{C_{\alpha_i,\alpha_{i+1}}}{N_{\alpha_i}}.
$$

| 参数 | 含义 |
| --- | --- |
| `pair_counts` | 按链顺序排列的相邻类别配对计数，共 K−1 个；第一个是 ER-alpha1。允许已按曝光归一化的非整数背景估计。 |
| `intermediate_counts` | 中间事件类别的总计数，共 K−2 个；依次为 N_alpha1、N_alpha2 等，必须为正。 |
| `observed_chains` | 实际完整链数，默认 1。 |

两类事件的链只需一个 `pair_counts`，`intermediate_counts` 用空数组或省略。
模板四类事件的例子：`12 × (5/1000) × (2/800) = 0.00015`。

配对计数必须按相同时间、能量、位置条件以及相容的曝光得到，并代表**偶然配对本底**。
输入包含真实强关联的配对计数会污染本底估计。程序不自动选择旁带，也不再除以位置分区数；
位置门已经包含在配对统计中。中间计数为零时，方法不可定义；配对计数为零时只得到零值代入估计。
允许多个候选后继导致 C/N 大于 1 时，该因子是平均配对数，不应解读为一次 Bernoulli 概率。

`n_b` 得到后，按截图式 (4.23) 对泊松分布求上尾，输出 `P_at_least_observed_chains`。
这一步沿用随机链数的泊松假设，不传播配对计数的统计误差。

运行：在本目录执行 `python calculate.py`，也可 `python calculate.py --config other.json --json`。
推荐解释器为 `/home/evalie/miniconda3/envs/py_env/bin/python`。

公式直接依据用户提供的截图。相关背景文献：
Lazarev et al., *Discovery of Enhanced Nuclear Stability near the Deformed Shells N=162 and Z=108*,
[Phys. Rev. Lett. 73, 624 (1994)](https://doi.org/10.1103/PhysRevLett.73.624)；
Stoyer et al. 对不同方法的讨论见 [NIM A 455, 433–441 (2000)](https://doi.org/10.1016/S0168-9002(00)00503-9)。
截图未附参考文献编号对应表，因此不把其中的 [176]、[177] 擅自对应到具体条目。
