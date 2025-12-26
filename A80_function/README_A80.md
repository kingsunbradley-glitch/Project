# A80 / X80 / Y80 批处理 QA 管线（ROOT）

本文件夹提供一套基于 ROOT 的小型 C++ 程序 + 并行运行脚本，用于批量处理实验数据 `run%05d_map.root`，并输出：

- **每个 run 一份 QA PDF（多页）**：`run%05d_QA_A80.pdf`
- 一份合并后的 **数值汇总 ROOT**：`summary.root`（每个 run 一行）
- 一份合并后的 **QA 对象 ROOT**：`qa_plots.root`（保存谱图/热力图/投影/Canvas 等）

核心目的：在三个能量峰附近定义 ROI，构建 XY 热力图，并用 **A80 / X80 / Y80** 量化“聚焦/弥散”程度，便于你画 run-by-run 的漂移/色散趋势（`mux/muy` 等）与对比不同设置。

---

## 1. 推荐目录结构

建议把代码放到实验数据目录下的 `function/`，数据放在上一级，避免污染原始数据：

```
<DATA_DIR>/
  run00001_map.root
  run00002_map.root
  ...
  function/
    *.cpp *.h
    Makefile
    mulicore_run.sh
    README.md
    out/                (自动创建；所有输出都放这里)
```

默认所有结果写入 `function/out/`，不会把 `summary.root / qa_plots.root / pdf/` 等文件堆到实验数据同名目录下。

---

## 2. 环境依赖

- Linux / WSL / macOS shell
- 已正确安装 ROOT（要求命令可用：`root-config`, `root`, `hadd`）
- 支持 C++17 的编译器（推荐 `g++`）

快速检查：

```bash
root-config --version
which root
which hadd
g++ --version
```

---

## 3. 编译

### 3.1 使用 Makefile（推荐）

在 `function/` 目录执行：

```bash
make -j5
```

生成可执行文件：

- `./process_runs_A80`

清理编译产物：

```bash
make clean
```

### 3.2 不用 Makefile（脚本自动编译）

如果你的 `mulicore_run.sh` 已写好“运行前编译”，那么直接跑脚本也可以（见下一节）。

---

## 4. 运行（处理 run 区间）

### 4.1 单进程运行（用于调试 / 验证）

程序接口：

```bash
./process_runs_A80 <INDIR> <RUN_FIRST> <RUN_LAST> <outSummary.root> <outQA.root> <pdfDir>
```

示例：数据在上一级目录（`..`），只跑 run 16：

```bash
./process_runs_A80 .. 16 16 out_summary.root out_qa.root out_pdf
```

单跑一个 run 的意义：
- 验证能打开 `run%05d_map.root`
- 验证 `tr_map` / branch 名称匹配
- 验证 PDF 能输出、图正常

### 4.2 并行批处理（推荐）

脚本接口（以本项目脚本为准；常见形式如下）：

```bash
./mulicore_run.sh NPROCS INDIR RUN_FIRST RUN_LAST [OUTDIR] [--no-clean]
```

典型用法：你在 `function/` 目录，数据在上一级：

```bash
./mulicore_run.sh 8 ../ 1 385
```

默认输出目录为：

```
function/out/
```

#### 清理策略（强烈建议开启）

默认会在运行前清理 `OUTDIR` 的历史输出（避免混淆旧结果）。  
如果你明确想保留旧输出，使用：

```bash
./mulicore_run.sh 8 ../ 1 385 ./out --no-clean
```

> 重要提示：如果脚本里会 `cd jobs/job_xxx` 运行子任务，**INDIR 必须是绝对路径**或脚本要做绝对化，否则相对路径（比如 `.`、`..`）在子目录里会失效。  
> 推荐脚本内做一次：`INDIR="$(cd "$INDIR" && pwd)"`

---

# 5. 输入数据说明（程序读什么）

## 5.1 输入文件命名规则

对 run 编号 `r`，程序会尝试打开：

```
<INDIR>/run%05d_map.root
```

例如 run 16：

```
<INDIR>/run00016_map.root
```

若文件缺失/无法打开，则该 run **跳过**（不会生成该 run 的 PDF，也不会写该 run 的 QA 对象）。

## 5.2 输入 TTree

每个输入 ROOT 文件中必须存在：

- `TTree *tr_map`

若不存在 `tr_map`，该 run 跳过。

## 5.3 使用到的 Branch（分支）及意义

本流程只读取以下分支（名字必须严格一致）：

| Branch 名称 | 读取类型 | 在本流程中的用途 |
|---|---|---|
| `DSSDX_mul` | `UShort_t` | X 侧 multiplicity（多重性）；要求 **== 1** |
| `DSSDY_mul` | `UShort_t` | Y 侧 multiplicity；要求 **== 1** |
| `DSSDX_Ch[256]` | `Double_t[256]` | X 命中条/像素 index；用 **`DSSDX_Ch[0]`** 作为 X |
| `DSSDY_Ch[256]` | `Double_t[256]` | Y 命中条/像素 index；用 **`DSSDY_Ch[0]`** 作为 Y |
| `DSSDX_E[256]` | `Double_t[256]` | 能量（keV）；用 **`DSSDX_E[0]`** 作为能量 |

> 说明：本流程专门针对**DSSD单触发**事件，先用 `DSSDX_mul==1 && DSSDY_mul==1` 筛选，  
> 然后只取数组第 0 个元素（`[0]`）作为该事件的 X、Y、E。

---

# 6. 分析逻辑

## 6.1 全局门（所有后续都必须满足）

后续能谱、峰搜索、XY 热力图、A80/X80/Y80 都只使用满足以下条件的事件：

- `DSSDX_mul == 1`
- `DSSDY_mul == 1`
- `DSSDX_E[0] ∈ [5000, 6000] keV`

该门用于压制噪声与多重性事件（避免噪声进入 ROI）。

## 6.2 三个峰的搜索窗口（固定范围）

在全局门下的能谱中，分别在三个固定窗口独立找峰：

- Peak1：**5100–5200 keV**（$^{239}Pu:5156.6~\text{keV}$）
- Peak2：**5435–5535 keV**（$^{241}Am:5485.6~\text{keV}$）
- Peak3：**5754–5854 keV**（$^{244}Cm:5804.8~\text{keV}$）

### 峰寻找与拟合步骤

对每个峰窗口：
1. 在窗口内找最大 bin 作为 `mu_seed`
2. 以 `mu_seed` 为中心，用局部拟合窗口 **±80 keV** 做高斯拟合（范围裁剪到该峰窗口内）
3. 用第一次拟合得到的 `mu, sigma`，再在更紧的 **`[mu ± 2σ]`** 里二次拟合（同样裁剪到窗口内）
4. 定义用于填 XY 的 ROI（你最关心的）：

```
ROI = [mu - 3σ, mu + 3σ]
```

若某个 run 里某个峰不存在（例如磁场设置导致只剩一个峰），该峰会被标记为“未找到”，并且不会参与 XY 填充与指标计算；程序仍可正常跑完该 run。

## 6.3 XY 热力图与投影

对通过全局门的事件，如果其能量落在某个峰的 ROI 内，则填入对应 XY 图：

- 命中 Peak1 ROI → 填 `hXY_p1`
- 否则命中 Peak2 ROI → 填 `hXY_p2`
- 否则命中 Peak3 ROI → 填 `hXY_p3`
- 三者任意命中 → 同时也会填入 ROI 并集图 `hXY_all`

坐标定义：

- `x = DSSDX_Ch[0]`
- `y = DSSDY_Ch[0]`

XY binning（像素化）：

- X：128 bins，范围 `[-0.5, 127.5]`（对应 0–127）
- Y：48 bins，范围 `[-0.5, 47.5]`（对应 0–47）

并对每张 XY 图做投影（用于 X80/Y80）：
- `hX_*`：对 XY 做 X 投影
- `hY_*`：对 XY 做 Y 投影

---

# 7. A80 / X80 / Y80 指标定义

指标思想：**覆盖 80% 计数需要多少个 bin/pixel？**

## 7.1 A80（二维）

对某个 2D 直方图（例如 `hXY_all`）：

1. 取所有非零 bin 的计数集合 `{c_i}`
2. 按计数从大到小排序
3. 累加，直到达到总计数的 80%

输出三类量（常用的是 `A80_eq`）：

- `A80_pix`：达到 ≥80% 所需的**整 bin 数**
- `A80_eq`：最后一个 bin 做线性插值后的“等效 bin 数”（可为小数）
- `A80_cov`：使用 `A80_pix` 个 bin 时实际覆盖率（≥0.8，通常略大于 0.8）

**解释**：
- `A80_eq` 越小 → 命中更集中（更“聚焦”）
- `A80_eq` 越大 → 分布更散（更“弥散/漂移”）

## 7.2 X80 / Y80（一维）

对 `hX_*` 或 `hY_*`（投影直方图）做同样的排序累计，得到：
- `X80_eq / X80_pix / X80_cov`
- `Y80_eq / Y80_pix / Y80_cov`

**解释**：
- X80 反映 X 向展开程度
- Y80 反映 Y 向展开程度

---

# 8. 输出文件说明（看完就能用输出，不用看代码）

默认输出目录：

```
function/out/
```

其中包含：

- `summary.root`：数值汇总（每个 run 一行，最适合画趋势）
- `qa_plots.root`：QA 对象库（每个 run 一个目录，存直方图/Canvas）
- `pdf/`：每个 run 的 QA PDF
- `jobs/`：并行任务日志
- `hadd_summary.log`, `hadd_qa.log`：合并日志

下面详细解释每个输出文件“干什么、里面有什么”。

---

## 8.1 out/summary.root（最重要：趋势分析用）

### 文件用途
这是“只存数字”的输出，便于：
- 画 `A80_all vs run`
- 画 `mux/muy vs run`（漂移趋势）
- 画峰位 `muE vs run`（能量标定/漂移）
- 快速筛选异常 run（比如 `nPeaks<3` 或 `N_all==0`）

### 内部对象
包含一个树：

- `TTree summary`：**每个 run 一条 entry**

### summary 树的 Branch（及用途）

> 注：具体分支名以你的版本为准，但本项目的设计意图如下（你画趋势主要用这些）。

```text
runnum (Int)
  run编号，例如 16 对应 run00016_map.root

nPeaks (Int)
  成功找到的峰数（0..3）

muE[3] (Double[3])
sigE[3] (Double[3])
  三个峰拟合得到的 mu 和 sigma（keV）
  如果该峰未找到，通常会填 -999 或保持无效值

N_all (Long64)
  ROI 并集（hXY_all）里的总计数（entries）
  很小或为 0 往往意味着该run在[5000,6000]门内统计量少或峰不在窗口里

A80_eq_all / A80_pix_all / A80_cov_all
X80_eq_all / X80_pix_all / X80_cov_all
Y80_eq_all / Y80_pix_all / Y80_cov_all
  ROI 并集（all）对应的 A80 / X80 / Y80

Npk[3] (Long64[3])
  三个峰各自 ROI 的计数

mux[3], muy[3] (Double[3])
  三个峰各自 XY 的均值（用于漂移/色散趋势）
  注意：这是“在 ROI 事件上做位置均值”，不额外加假设

A80eq[3], X80eq[3], Y80eq[3]（以及对应 pix/cov）
  三个峰各自 ROI 的 A80 / X80 / Y80（用于分峰趋势）
```

### 如何查看 summary.root（ROOT 命令示例）

打开：

```bash
root -l out/summary.root
```

在 ROOT 里：

```cpp
TTree* t = (TTree*)gFile->Get("summary");
t->Print(); // 打印所有 branch

// 看前 30 行：run、峰数、总计数、A80/X80/Y80
t->Scan("run:nPeaks:N_all:A80_eq_all:X80_eq_all:Y80_eq_all", "", "", 30);

// 找异常：缺峰 或 没统计量
t->Scan("run:nPeaks:N_all", "nPeaks<3 || N_all==0");

// 画趋势：A80 vs run
t->Draw("A80_eq_all:run", "", "P");

// 画漂移：Peak1 的 mux/muy vs run（条件：该峰有统计量）
t->Draw("mux[0]:run", "Npk[0]>0", "P");
t->Draw("muy[0]:run", "Npk[0]>0", "P");

// 画峰位漂移：muE vs run
t->Draw("muE[0]:run", "muE[0]>0", "P");
```

> 实用建议：先用 `summary.root` 找异常 run 编号，再去 `qa_plots.root` 或 PDF 精查该 run 的谱与 XY。

---

## 8.2 out/qa_plots.root（QA 对象库：直方图/画布都在这里）

### 文件用途
这是“所有 QA 图的仓库”，用于：
- 不看 PDF、直接程序化读取直方图并重画
- 批量导出某些 run 的 XY 热力图/投影
- 做更细致的二次分析（比如你要对某些 run 重新算别的指标）

### 内部结构
通常为：

- 顶层：每个 run 一个目录

```
run00016/
run00017/
...
```

进入某个 run 目录后，你会看到：

- `hE_5000_6000`：能谱（应用全局门后：mx==1,my==1 且 E∈[5000,6000]）
- `hXY_all`：三个峰 ROI 并集的 XY 热力图
- `hXY_p1 / hXY_p2 / hXY_p3`：分峰 ROI 的 XY 热力图
- `hX_* / hY_*`：对应的投影
- `c_*`：用于生成 PDF 的画布（Canvas），例如能谱页/并集页/分峰页等

### 如何查看 qa_plots.root（ROOT 命令示例）

打开：

```bash
root -l out/qa_plots.root
```

在 ROOT 里：

```cpp
gDirectory->ls();            // 列出 run 目录
gDirectory->cd("run00016");  // 进入 run00016
gDirectory->ls();            // 列出该 run 的对象

// 画 XY 并集热力图
TH2D* h = (TH2D*)gDirectory->Get("hXY_all");
h->Draw("COLZ");

// 画能谱
TH1D* e = (TH1D*)gDirectory->Get("hE_5000_6000");
e->Draw();
```

如果 ROOT 支持 GUI 浏览器：

```bash
rootbrowse out/qa_plots.root
```

---

## 8.3 out/pdf/run%05d_QA_A80.pdf（每个 run 的 QA PDF）

### 文件用途
这是“人工快速 QA”的结果：你不需要打开 ROOT，也能直接翻页确认：
- 能谱里三个峰是否存在、拟合是否合理
- ROI 是否框到正确区域
- XY 热力图是否集中/漂移/扩展
- A80/X80/Y80 是否异常

### PDF 常见页面（典型）
1. 能谱页：
   - 全局门下的能谱（5000–6000）
   - 标出三个搜索窗口
   - 标出拟合得到的 ROI 边界（若该峰存在）
   - 显示 mu/sigma/ROI 数值

2. “并集页”（All peaks union）：
   - `hXY_all` + X/Y 投影
   - A80/X80/Y80（并集）

3–5. 分峰页（Peak1/2/3）：
   - `hXY_p?` + 投影
   - 计数 N、`mux/muy`
   - 分峰 A80/X80/Y80

若某个峰未找到，对应分峰页会提示 “Peak NOT FOUND”。

---

## 8.4 out/jobs/ 与 out/hadd_*.log（排错用）

- `out/jobs/job_*/run.log`
  - 每个并行 worker 的完整输出（最常用排错入口）
- `out/hadd_summary.log`, `out/hadd_qa.log`
  - 合并时 `hadd` 的日志

常见排错：
- “某些 run 的 PDF 缺失” → 检查该 run 的输入文件是否存在，或者在 log 里看是否被跳过
- “看起来卡在第 2 步不动” → 脚本在 `wait` 等后台进程，进度都在 `run.log` 里

---

# 9. 常见问题（FAQ）

## Q1：为什么某些 run 没有 PDF？
最常见原因：
1) 该 run 的 `run%05d_map.root` 缺失/打不开 → 该 run 直接跳过  
2) 该 run 的 `tr_map` 不存在 → 跳过  
3) 该 run 在全局门（E∈[5000,6000] 且 mx=my=1）下几乎无统计量 → 可能 PDF 仍生成，但 `nPeaks=0`、`N_all=0`

建议：
- 先查 `out/pdf/` 是否缺该 run
- 再查对应 worker 的 `out/jobs/.../run.log`

## Q2：为什么 `nPeaks=0` 或 `N_all=0`？
典型原因：
- 在 `[5000,6000]` keV + 单击门下统计量极低（或峰不在窗口）
- 磁场/设置导致峰漂移出固定搜索窗口（5100–5200等）

解决思路：
- 先用 `summary.root` 找出这些 run 编号
- 打开对应 PDF 看能谱（是否峰跑飞）
- 需要时调整峰窗口或能量门后重新跑

## Q3：跑的时候内存占用很大怎么办？
ROOT 多进程 + 画图 + 写 PDF 容易吃内存。建议：
- 把并行数 `NPROCS` 从 8 降到 6/4
- 确保没有内存泄漏（本项目应避免常见的 Projection 泄漏）

---

## 10. 输出文件“怎么用”的典型工作流（推荐）

1) 先从 `out/summary.root` 做“全局趋势”：
   - `A80_eq_all vs run`
   - `mux/muy vs run`（按峰分别画）
   - `muE vs run`（峰位漂移）
2) 用 `summary.root` 过滤出异常 run（缺峰、A80很大、mux/muy突变）
3) 再打开对应 run 的：
   - PDF：快速人工确认
   - 或 `qa_plots.root`：做更细的直方图检查/导出

---

如需我给你配一份“自动画趋势并导出图片”的 ROOT 宏（读 `summary.root`，输出一组 png/pdf），告诉我你想画哪些量（A80_all、三峰 mux/muy、三峰 muE 等），我可以直接给你一键脚本。
