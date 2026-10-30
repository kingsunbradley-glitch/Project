# Monte Carlo 随机插入方法

实现截图式 (4.31)：将一个符合指定链步骤入口条件的标记事件随机放入实际数据，
搜索包含这个标记的完整链；重复后计算

$$
P_{insertion}=N_{accepted}/N_{inserted}.
$$

方法来源：[Stoyer et al., Random probability analysis of heavy-element data, NIM A 455, 433–441 (2000)](https://doi.org/10.1016/S0168-9002(00)00503-9)。
作者的公开说明见 [Random Probability Analysis of Recent 48Ca Experiments](https://www.osti.gov/servlets/purl/877859)。
这里实现明确列出的时间/能量/位置门，没有实现文献中的全部核素系统学判据。

## 两种用法

`mode="simulate"`：读取 `events_csv`，真正运行随机插入和链搜索。`input.json` 和
`example_events.csv` 是人工演示数据；默认 10000 次插入，固定随机种子。

`mode="counts"`：已有自己的关联搜索程序时，输入它给出的 `inserted_events` 和
`correlated_events`。后者应是**至少找到一条完整链的插入事件数**，一个插入事件最多计一次。
示例：`python calculate.py --config input_counts.json`。

两者都输出 Clopper–Pearson 二项置信区间。零命中仍给出非零上限，不能证明偶然概率为零。

## CSV 与参数

CSV 列：`time_s,position_bin,event_type,energy_mev`。前三列必需；使用能量门时最后一列必需。
`position_bin` 从 0 编号，时间单位为秒；能量列及门均为 MeV。允许未排序输入。

| 参数 | 含义 |
| --- | --- |
| `events_csv` | 实际列表数据路径，相对配置文件目录解析。 |
| `live_intervals_s` | 排好序且不重叠的连续采集区段，例如 `[[0,1000],[1100,2000]]`。不跨采集间断组链。 |
| `position_bins` | 有效位置单元个数；实际 CSV 和注入事件均用相同编号。 |
| `stages` | 固定顺序链的步骤。每步 `event_type` 与 CSV 精确匹配；可加 `energy_mev:[下限,上限]`。 |
| `windows_s` | K−1 个相邻步骤时间窗，每个用 `[非负下限,上限]` 表示。 |
| `injection_index` | 标记事件所在步骤，从 0 编号。默认 0 注入 ER；也支持在链中间/末端插入并搜索两侧。 |
| `trials` | 插入次数。 |
| `seed` | 随机种子；相同配置、输入和 NumPy 版本下可复现。 |
| `injection_interval_weights` | 可选，各采集区段的注入概率质量，程序归一化。默认按区段长度抽样；区段内时间均匀。 |
| `position_weights` | 可选，各位置的注入概率质量。默认均匀，可用真实入口分布指定非均匀权重。 |
| `confidence` | 二项区间的置信水平，默认 0.95。 |
| `reference_events` | 可选，实验中对应插入步骤的真实入口事件总数，用于估计 n_b。 |
| `observed_chains` | 与 `reference_events` 配套，默认 1。 |

每次仅插入一个标记事件，背景列表保持固定，标记事件之间不互相组链。
时间严格递增，要求同一位置单元；允许链成员之间有其他事件；满足任一候选组合即判为命中。
同一真实事件可供不同试验使用。与实际实验算法比较时，应使这些约定和实际搜索一致。
目前位置采用离散单元；连续位置差门、邻条带门及束流触发规则应在实际搜索程序中实现，
再通过 `counts` 模式汇总结果，或扩展这里的 `has_chain`。

实际数据的聚集、时间变化和空间不均匀会保留。区段内均匀插入是可配置抽样设计的一部分，
应按要检验的事件分布选取权重；时间分布与位置有关时，可对不同区域分别运行再按入口概率加权。
入口标记默认已经满足该步骤的类别和能量门；没有模拟该能量与后续寿命的额外耦合。
固定数据中的真实关联也会参与匹配，因此本方法估计标记与观测数据的偶然组链率。

若提供 `reference_events`，输出 `poisson_plugin_estimate`：先用 `N_reference × P_insertion`
估计随机链期望数，再按泊松假设给出整个搜索的尾概率。
这是把模拟估计代入解析模型，非逐次完整实验模拟；概率区间只包括有限插入次数的抽样误差，
不包括真实背景数据量有限、门选择或抽样分布估计的误差。

在演示配置中，搜索概率可独立求得：两个位置上可接受的 ER 插入区间总长为 15 s，
总时空曝光为 100×2，因此精确插入概率为 0.075。模拟值应在其统计误差范围内。

运行：本目录 `python calculate.py`，支持 `--config`、`--json`。
推荐解释器为 `/home/evalie/miniconda3/envs/py_env/bin/python`。
