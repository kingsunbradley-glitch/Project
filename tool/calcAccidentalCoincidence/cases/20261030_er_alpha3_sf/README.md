# ER–alpha1–alpha2–alpha3–SF 的逐级偶然符合

沿用 Berkeley 逐级窗口至少一次的公式，按用户确认在 alpha3 后加一个 100 s 的 SF 窗口。

共用输入：有效曝光 655.2 h = 2358720 s；位置单元 128×48=6144；
ER=13250000，alpha1=65110，alpha2、alpha3 各使用同一筛选类别的 154000；
四个相邻时间窗依次为 0.1、100、100、100 s。
SF 总计数也视为同一曝光、全探测器范围内的计数。
背景按时间、位置均匀模型归一化，每步从前一个选定事件重新计时。

原先不含 SF 的整次随机链期望数为 n_base=6.715226456307335e−6。
追加 SF 后

$$
\mu_{SF}=\frac{N_{SF}}{6144\times2358720}\times100,\qquad
p_{SF}=1-e^{-\mu_{SF}},\qquad
n_b=n_{base}\,p_{SF}.
$$

| SF-like 全实验总数 | mu_SF | 最后一个 100 s 窗中 SF 至少一次的概率 | 单个 ER 的完整链概率 | 整次随机链期望数 n_b | 整次至少一条概率（泊松模型） |
| --- | ---: | ---: | ---: | ---: | ---: |
| 188600 | 0.001301409857 | 0.001300563391 | 6.591379388e−16 | 8.733577689e−9 | 8.733577651e−9 |
| 194 | 1.338671857e−6 | 1.338670961e−6 | 6.784512193e−19 | 8.989478656e−12 | 8.989478656e−12 |

SF 计数改为 194 后，期望随机链数降低约 971.53 倍。
整次至少一条的概率按 1−exp(−n_b) 计算；两种情况下 n_b 很小，与该概率数值近乎相同。

`SF_188600.json`、`SF_194.json` 为两组独立配置，`results.json` 保存完整结果。
没有改变原先三步 alpha 例子的默认输入。

从仓库根目录复算：

```bash
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/Berkeley/calculate.py --config tool/calcAccidentalCoincidence/cases/20261030_er_alpha3_sf/SF_188600.json
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/Berkeley/calculate.py --config tool/calcAccidentalCoincidence/cases/20261030_er_alpha3_sf/SF_194.json
```

本次结果已与直接计算各步 1−exp(−mu) 的乘积作独立数值校验。

## 其他方法比较

使用相同曝光、位置数和四个相邻窗口，结果如下。表内为 n_b；整次至少一条的泊松概率
为 1−exp(−n_b)，在当前极小 n_b 下与表内数字在所示精度上相同。

| 方法及具体模型 | SF-like=188600 | SF-like=194 |
| --- | ---: | ---: |
| GSI：固定顺序、相邻成员间无其他候选类别事件 | 7.603332645e−9 | 7.836064943e−12 |
| GSI 低率展开 / Lazarev 率模型配对估计 | 8.748554820e−9 | 8.999043664e−12 |
| Berkeley：每步自己的窗口内至少一次 | 8.733577689e−9 | 8.989478656e−12 |
| Dubna 计数式加全部四步时间门（见下述额外条件） | 8.675002898e−9 | 8.958268516e−12 |
| 实际列表 Monte Carlo 插入法 | 缺输入数据 | 缺输入数据 |

GSI 的总率对 ER、alpha1、alpha23、SF 四类各计一次，链顺序为 `[0,1,2,2,3]`。
加入 SF 后，SF 率也进入每一步的总竞争率，因此需重算全部指数项，不能只在先前 GSI
结果后乘一个 SF 因子。这里假设四类互不重叠、背景均匀。

Lazarev 尚缺实测的相邻配对数。本次按率模型估计四个配对期望数：

| 配对 | SF-like=188600 | SF-like=194 |
| --- | ---: | ---: |
| ER-alpha1 | 5.953001296 | 5.953001296 |
| alpha1-alpha23 | 69.18959996 | 69.18959996 |
| alpha23-alpha23 | 163.64918437 | 163.64918437 |
| alpha23-SF | 200.41711801 | 0.2061554660 |

各配对均满足对应相邻时间门和位置门，分母依次使用 65110、154000、154000。
这些是模型期望值，不是声称实际观察到了小数条配对。代入后恰好退化为 GSI 低率式。

### Dubna 如何施加相同时间门

以 ER 为给定锚点，最大总跨度 W=0.1+100+100+100=300.1 s。
固定窗计数条件为：恰好一个 alpha1、两个 alpha23、一个 SF。
两个 alpha23 是同一类别的 K=2 计数，包含 2! 分母。

计数概率本身不含各步时间门，不能直接作为本次逐级搜索结果。
在均匀泊松背景且给定这些计数后，alpha1 和 SF 的时间分别均匀分布，
两个有序 alpha23 的联合时间密度为 2/W²。四个有序时间的合格体积为四个窗口长度的乘积，
故需再乘条件接受度

$$
A_t=\frac{2(0.1)(100)(100)(100)}{300.1^4}
=2.465846363\times10^{-5}.
$$

最终

$$
n_b^{Dubna+time}=N_{ER}Q_1(r_1W)Q_2(r_BW)Q_1(r_{SF}W)A_t
=n_b^{low}\exp[-(r_1+r_B+r_{SF})W].
$$

这是对原始固定窗计数式显式添加时间选择的扩展；它另外要求固定总窗内没有额外
alpha1、alpha23 或 SF，但允许额外 ER。它与 Berkeley 的逐级至少一次、GSI 的相邻事件
无其他候选类别插入各有不同事件选择条件，不能把差值解释为同一概率的纯数值误差。
原始固定窗计数输出保留在结果文件的 reference 字段，仅供核对，未代替逐级门结果。

实际 MCRP 仍需事件时间、位置、类别列表，或已有随机插入的总数/命中数。
Dubna 的时间、寿命及位置相容性分支也需要实际链数据；搜索窗上限不能代替观测值。

`compare_methods.py` 复算并打印所有可计算分支；`comparison_results.json` 保存各方法输入、
完整输出与模型条件，原来的 `results.json` 仍为 Berkeley 两种 SF 计数的结果。

```bash
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/cases/20261030_er_alpha3_sf/compare_methods.py
```

复算时额外核对了 GSI 下一事件积分的直接乘积、Lazarev 与低率式的一致性，以及
Dubna 加时间门的两种等价表达式。
