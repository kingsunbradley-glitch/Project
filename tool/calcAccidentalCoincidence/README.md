# 偶然符合与关联链统计

按用户提供的第 4.4.1–4.4.5 节截图整理。各方法独立放置，入口均为本目录下各子目录的 `calculate.py`，默认读取该脚本同目录的 `input.json`。

| 目录 | 主要输入 | 主要输出 |
| --- | --- | --- |
| [GSI](GSI/README.md) | 各类背景率、相邻/首事件时间窗、位置分区 | 随机链期望数、整个搜索的泊松尾概率 |
| [Lazarev](Lazarev/README.md) | 相邻类别的偶然配对计数、中间事件计数 | 由配对计数估计的随机完整链期望数、泊松尾概率 |
| [Berkeley](Berkeley/README.md) | 各步骤背景计数/率、各自时间窗、像素数、注入数 | 逐级至少一次的乘积、期望随机链数；另保留旧总窗计数模型 |
| [MonteCarlo](MonteCarlo/README.md) | 实际事件 CSV、搜索条件、随机插入分布，或已有模拟命中计数 | 每个插入事件形成完整链的频率及其二项置信区间 |
| [Dubna](Dubna/README.md) | 分类背景率与计数、时间间隔、半衰期、位置差相关矩阵 | 分类泊松概率，以及时间、衰变与位置相容性统计 |

新增方法使用 `common.py` 提供通用统计与参数检查，不依赖 GSI 的实现。依赖现有 `py_env` 的 NumPy、SciPy。
所有时间用 **s**，计数率用 **Hz**；若原数据是每分钟/每小时计数，分别除以 60/3600。
模板中的数字已注明是文献例子还是人工示例。

本次用户参数的具体复算见 [128×48 位置单元、655.2 h、ER–alpha1–alpha2–alpha3](cases/20261030_er_alpha3/README.md)。

运行示例（在仓库根目录）：

```bash
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/Lazarev/calculate.py
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/Berkeley/calculate.py
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/MonteCarlo/calculate.py
/home/evalie/miniconda3/envs/py_env/bin/python tool/calcAccidentalCoincidence/Dubna/calculate.py
```

均支持 `--config /path/to/parameters.json` 和 `--json`。配置中的 CSV 相对路径以配置文件所在目录解析。
程序只读取输入、向终端输出；不写入原始事件数据。

这些输出的统计问题各有定义。Berkeley 当前默认是逐级窗口“至少一次”的乘积；
其旧总窗模型与 Dubna 的 `exact` 是“恰好出现指定个数”，
GSI/Lazarev 的尾概率是“出现至少观测链数”，Monte Carlo 是“一个随机插入事件被搜索接受”的频率。
比较方法时，应先统一搜索条件、条件事件、位置范围和输出定义，不能仅比较名称相近的概率。
任何单项输出都不自动给出一条链为真的后验概率。

截图中有两处会影响实现的公式排版问题：Berkeley 式 (4.29) 的泊松分母应为阶乘；
Dubna 的等衰变常数卷积式中指数应为负号。实现采用归一化的泊松分布和指数等待时间卷积。
Dubna 的最小/最大时间是次序统计量的**期望值**，并非统计置信界。
各目录 README 记录了实现的具体约定与可复核来源。

新增方法的验证：

```bash
/home/evalie/miniconda3/envs/py_env/bin/python -m unittest discover -s tool/calcAccidentalCoincidence -p test_methods.py -v
```
