# ROOT `TCanvas` 波形读取与短寿命 pile-up 分析程序 —— Codex 实现任务说明

## 1. 任务目标

编写一个 ROOT/C++ 程序，用于读取用户输入的 ROOT 文件，例如：

```text
wave_242Fm.root
```

该 ROOT 文件的顶层对象格式类似：

```text
TFile**         wave_242Fm.root
 TFile*         wave_242Fm.root
  KEY: TCanvas  ChainNum_1;1    ChainNum_1: run=6, map_entry=43971794
  KEY: TCanvas  ChainNum_2;1    ChainNum_2: run=9, map_entry=42197111
  KEY: TCanvas  ChainNum_3;1    ChainNum_3: run=..., map_entry=...
  ...
```

也就是说：

- 一个 ROOT 文件中包含多个 `TCanvas`；
- Canvas 名称一般为：

```text
ChainNum_1
ChainNum_2
ChainNum_3
...
```

- Canvas title 中可能包含事件信息，例如：

```text
ChainNum_1: run=6, map_entry=43971794
```

程序需要自动遍历整个 ROOT 文件，读取每个 `TCanvas` 内部真正保存的波形对象，并将波形数据提取出来，供后续 pile-up 分析使用。

---

# 2. 第一阶段：先实现“通用 ROOT Canvas 波形读取器”

不要预先假设 Canvas 内部一定保存的是某一种固定对象。

程序必须递归检查：

```text
TCanvas
  └── TPad
       ├── TGraph
       ├── TGraphErrors
       ├── TGraphAsymmErrors
       ├── TH1 / TH1D / TH1F
       ├── TF1
       ├── TLine
       ├── TLatex
       └── 其他 primitive
```

真正需要提取波形的对象优先级为：

1. `TGraph`
2. `TGraphErrors`
3. `TGraphAsymmErrors`
4. `TH1`

程序必须能够递归进入所有 `TPad`，不能只检查 Canvas 最外层的 primitive list。

建议实现函数：

```cpp
void ScanPad(TVirtualPad* pad,
             int chainIndex,
             const std::string& canvasName,
             const std::string& canvasTitle);
```

递归逻辑：

```cpp
TIter next(pad->GetListOfPrimitives());
TObject* obj = nullptr;

while ((obj = next())) {

    if (obj->InheritsFrom(TPad::Class())) {
        ScanPad((TVirtualPad*)obj, ...);
    }

    else if (obj->InheritsFrom(TGraph::Class())) {
        // 提取波形
    }

    else if (obj->InheritsFrom(TH1::Class())) {
        // 提取波形
    }
}
```

---

# 3. 输入方式

程序必须支持命令行输入 ROOT 文件名。

例如编译后：

```bash
./wave_reader wave_242Fm.root
```

或者如果实现成 ROOT macro：

```bash
root -l
.x wave_reader.C("wave_242Fm.root")
```

更推荐最终写成一个普通 C++ 可执行程序：

```bash
g++ wave_reader.cpp $(root-config --cflags --libs) -O2 -o wave_reader
```

使用：

```bash
./wave_reader wave_242Fm.root
```

---

# 4. ROOT 文件遍历

必须使用 ROOT 的 key 系统自动遍历文件，而不是硬编码 `ChainNum_1`、`ChainNum_2`。

推荐：

```cpp
TFile* fin = TFile::Open(filename.c_str(), "READ");

TIter nextKey(fin->GetListOfKeys());
TKey* key = nullptr;

while ((key = (TKey*)nextKey())) {

    std::string className = key->GetClassName();

    if (className == "TCanvas") {
        TCanvas* c = dynamic_cast<TCanvas*>(key->ReadObj());
        ...
    }
}
```

程序必须输出类似：

```text
============================================================
Canvas: ChainNum_1
Title : ChainNum_1: run=6, map_entry=43971794
============================================================

Primitive found:
  class = TGraph
  name  = Graph
  title = DSSD waveform
  N     = 4096

Primitive found:
  class = TGraph
  name  = Graph1
  title = SSD waveform
  N     = 4096
```

这样第一次运行时就可以知道真正的数据结构。

---

# 5. 从 Canvas title 解析事件信息

Canvas title 可能类似：

```text
ChainNum_1: run=6, map_entry=43971794
```

程序尝试自动解析：

```cpp
struct EventInfo {
    int chainNum = -1;
    int run = -1;
    long long mapEntry = -1;
};
```

需要支持：

```text
ChainNum_1
```

解析：

```text
chainNum = 1
```

以及：

```text
run=6
map_entry=43971794
```

建议使用 `std::regex`。

例如：

```cpp
std::regex reChain("ChainNum_([0-9]+)");
std::regex reRun("run=([0-9]+)");
std::regex reEntry("map_entry=([0-9]+)");
```

如果解析失败，不要退出程序，只需要把对应值设置为 `-1`。

---

# 6. 波形统一数据结构

不管原始对象是 `TGraph` 还是 `TH1`，最后都转换成统一的数据结构：

```cpp
struct Waveform {
    int chainNum;
    int run;
    long long mapEntry;

    std::string canvasName;
    std::string canvasTitle;

    std::string padName;
    std::string objectName;
    std::string objectTitle;
    std::string objectClass;

    std::vector<double> x;
    std::vector<double> y;
};
```

其中：

```text
x = 时间轴 / sample index

y = ADC amplitude
```

如果对象是 `TGraph`：

```cpp
for (int i = 0; i < graph->GetN(); ++i) {
    double x, y;
    graph->GetPoint(i, x, y);
}
```

如果对象是 `TH1`：

```cpp
for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double x = h->GetBinCenter(i);
    double y = h->GetBinContent(i);
}
```

---

# 7. 第一次运行必须提供 `--inspect` 模式

这一点非常重要。

因为当前只知道 ROOT 文件顶层保存的是 `TCanvas`，并不知道 Canvas 里面的具体对象名字与结构。

所以程序首先必须支持：

```bash
./wave_reader wave_242Fm.root --inspect
```

输出整个 ROOT 文件的结构，例如：

```text
[Canvas]
name  = ChainNum_1
title = ChainNum_1: run=6, map_entry=43971794

  [TPad]
  name = pad1

      [TGraph]
      name  = Graph0
      title = DSSD waveform
      N     = 4096

      [TGraph]
      name  = Graph1
      title = SSD waveform
      N     = 4096

  [TPad]
  name = pad2

      [TH1D]
      name = hWave
      bins = 4096
```

要求：

- 显示 primitive class；
- 显示 object name；
- 显示 title；
- 若为 `TGraph`，显示 `N`；
- 若为 `TH1`，显示 `NbinsX`；
- 显示 primitive 所属的 pad；
- 显示 Canvas title。

这是程序开发的第一优先级。

---

# 8. 不要把 Canvas 当图片读取

注意：

这个程序不是从截图中 OCR 波形。

必须读取 ROOT Canvas 内部保存的真正对象：

```text
TGraph
TH1
...
```

然后直接获得原始数值点。

禁止：

```text
canvas -> png -> image processing
```

正确流程：

```text
TFile
  -> TCanvas
      -> TPad
          -> TGraph / TH1
              -> x[i], y[i]
```

---

# 9. 第二阶段：绘制提取出来的 waveform

读取成功之后，程序提供：

```bash
./wave_reader wave_242Fm.root --draw 1
```

表示显示：

```text
ChainNum_1
```

或者：

```bash
./wave_reader wave_242Fm.root --draw-all
```

程序重新绘制提取的 waveform，而不是直接 Draw 原来的 Canvas。

例如：

```cpp
TGraph* g = new TGraph(wave.x.size(),
                       wave.x.data(),
                       wave.y.data());
```

这样可以验证提取的数据与原 Canvas 一致。

---

# 10. 第三阶段：自动寻找主要 pulse

在成功读取 waveform 后，实现第一版非常简单的 pulse finder。

先不要立即做复杂 fit。

处理流程：

```text
raw waveform
    ↓
baseline subtraction
    ↓
optional smoothing
    ↓
derivative / matched filter
    ↓
find pulse candidates
```

---

# 11. Baseline 估计

默认使用波形前面一段数据估计 baseline。

例如：

```text
前 10% samples
```

或者：

```cpp
Nbaseline = min(200, N/10);
```

baseline：

$$
B = \frac{1}{N_B}\sum_{i=1}^{N_B} y_i
$$

噪声：

$$
\sigma_B =
\sqrt{
\frac{1}{N_B-1}
\sum_i(y_i-B)^2
}
$$

得到：

```cpp
ycorr[i] = y[i] - baseline;
```

输出：

```text
Baseline = xxxx ADC
Noise RMS = xx ADC
```

最好同时实现 robust baseline，例如 median + MAD，作为可选模式。

---

# 12. Pulse polarity 必须自动判断

不要假设 pulse 一定向上。

程序比较：

```cpp
max(ycorr)
abs(min(ycorr))
```

如果：

```text
abs(min) > max
```

则认为是负脉冲，并自动乘以 `-1`。

最终分析内部统一使用“正脉冲”。

---

# 13. 初始 pulse finder

第一版可以使用导数：

$$
d_i = y_{i+1}-y_i
$$

寻找明显的 rising edge。

阈值可取：

$$
d_i > k\sigma_d
$$

例如：

```text
k = 5
```

同时加入 minimum separation，防止一个 pulse 被重复找到。

程序输出：

```text
Pulse candidate 1:
    sample = 1523
    x      = 7.615 us
    height = 10342 ADC

Pulse candidate 2:
    sample = 1578
    x      = 7.890 us
    height = 8321 ADC

Delta_t = 0.275 us
```

这一阶段只是 pulse candidate finder，不作为最终精确结果。

---

# 14. 第四阶段：构造单脉冲模板

为了后面做 pile-up 分析，需要支持从“干净单脉冲”生成 template。

建议允许用户指定一个或多个 Canvas，例如：

```bash
./wave_reader wave_242Fm.root --make-template 3 5 8 12
```

或者读取全部满足“只有一个 pulse”的 waveform。

处理步骤：

```text
baseline subtraction
↓
pulse polarity normalization
↓
time alignment
↓
amplitude normalization
↓
average
```

得到模板：

```cpp
std::vector<double> pulseTemplate;
```

将 template 保存为：

```text
pulse_template.root
```

ROOT 文件中保存：

```text
TGraph pulse_template
```

同时输出：

```text
pulse_template.csv
```

便于调试。

---

# 15. 第五阶段：pile-up 检测思路

后续不要直接使用一个高维 MINUIT 同时拟合所有参数。

采用：

```text
matched filter
+
variable projection
```

模型：

$$
y(t)=A_1h(t-t_1)+A_2h(t-t_2)+b_0+b_1t
$$

其中：

- $h(t)$：单脉冲模板；
- $A_1,A_2$：两个 pulse 的幅度；
- $t_1,t_2$：两个 pulse 的时间；
- $b_0,b_1$：残余 baseline。

定义：

$$
\Delta t=t_2-t_1
$$

---

# 16. Variable Projection

关键要求：

不要把：

```text
A1
A2
b0
b1
```

作为 MINUIT 非线性参数。

对于给定：

```text
t1
t2
```

模型可以写为：

$$
\mathbf y=X(t_1,t_2)\boldsymbol\beta
$$

其中：

$$
\boldsymbol\beta=
\begin{pmatrix}
A_1\\
A_2\\
b_0\\
b_1
\end{pmatrix}
$$

而：

$$
X=
\begin{pmatrix}
h(t_1) & h(t_2) & 1 & t
\end{pmatrix}
$$

因此：

$$
\hat\beta=(X^TX)^{-1}X^T y
$$

或者使用更加稳定的：

```text
QR decomposition
```

或：

```text
SVD
```

禁止直接显式计算矩阵逆，如果 ROOT 自带线性代数支持可以使用：

```cpp
TDecompQRH
TDecompSVD
TMatrixD
TVectorD
```

优先采用 QR 或 SVD。

这样非线性搜索只剩：

```text
t1, t2
```

如果第一个 pulse 已经由 trigger 很好地确定，则进一步变成：

```text
Delta_t
```

一个参数。

---

# 17. 一维 `Delta_t` 扫描模式

优先实现一个非常稳定的 brute-force scan：

```bash
./wave_reader wave_242Fm.root --fit 1
```

对某个 waveform 扫：

```text
Delta_t = 0 -> 2000 ns
```

步长首先使用一个 ADC sample。

例如：

```text
Delta_t       chi2
0 ns          182.4
10 ns         160.2
20 ns         131.5
30 ns          96.2
40 ns          71.4
50 ns          60.3
60 ns          57.1   <-- minimum
70 ns          59.0
80 ns          66.1
```

找到最小值后，在 minimum 附近进行 sub-sample refinement。

可以采用：

1. parabolic interpolation；
2. cubic spline；
3. Brent minimization。

首选 Brent method，因为只有一维。

---

# 18. 单 pulse vs 双 pulse

程序必须同时计算：

## 单脉冲模型

$$
y(t)=Ah(t-t_1)+b_0+b_1t
$$

得到：

$$
\chi^2_1
$$

## 双脉冲模型

$$
y(t)=A_1h(t-t_1)+A_2h(t-t_2)+b_0+b_1t
$$

得到：

$$
\chi^2_2
$$

输出：

$$
\Delta\chi^2=\chi^2_1-\chi^2_2
$$

同时计算：

```text
AIC
BIC
```

例如：

```text
Single pulse:
    chi2 = 241.3
    NDF  = 995
    BIC  = 268.9

Double pulse:
    chi2 = 172.1
    NDF  = 993
    BIC  = 213.5

Delta BIC = 55.4

Result: strong evidence for pile-up
```

不要一开始把某条波形强行假设为双 pulse。

---

# 19. 极短时间间隔：加入模板导数分析

对于：

$$
\Delta t\rightarrow0
$$

有：

$$
h(t-\Delta t)
\approx
h(t)-\Delta t h'(t)
$$

因此双脉冲近似：

$$
A_1h(t)+A_2h(t-\Delta t)
\approx
(A_1+A_2)h(t)-A_2\Delta t h'(t)
$$

所以程序后期可以实现：

```text
pulse template h(t)
derivative template h'(t)
```

计算投影：

$$
R_0=\langle y,h\rangle
$$

$$
R_1=\langle y,h'\rangle
$$

`R1` 可作为 unresolved pile-up 的快速判据。

这一部分作为第二阶段算法优化，不要求第一版立即完成。

---

# 20. 输出结果 ROOT 文件

程序最终输出：

```text
wave_242Fm_analysis.root
```

至少保存一个 `TTree`：

```text
WaveAnalysis
```

branch 建议：

```cpp
int chainNum;
int run;
long long mapEntry;

char canvasName[128];
char objectName[128];

int nSamples;

double baseline;
double noiseRMS;

int nPulseCandidate;

double t1;
double t2;
double deltaT;

double A1;
double A2;

double chi2_1pulse;
double chi2_2pulse;

double deltaChi2;
double AIC1;
double AIC2;
double BIC1;
double BIC2;
```

额外保存：

```text
TGraph chi2_vs_deltaT
```

对于每条分析过的 chain，可保存在目录：

```text
ChainNum_1/
ChainNum_2/
...
```

例如：

```text
ChainNum_1/raw_waveform
ChainNum_1/baseline_subtracted
ChainNum_1/best_single_fit
ChainNum_1/best_double_fit
ChainNum_1/residual_single
ChainNum_1/residual_double
ChainNum_1/chi2_vs_deltaT
```

---

# 21. CSV 输出

同时生成：

```text
wave_242Fm_analysis.csv
```

例如：

```text
chainNum,run,mapEntry,objectName,nSamples,baseline,noiseRMS,nPulse,t1,t2,deltaT,A1,A2,chi2_1,chi2_2,deltaBIC
1,6,43971794,Graph,4096,8192.4,3.2,2,7.132,7.401,0.269,10234,8421,251.4,174.2,58.1
```

---

# 22. 程序命令建议

要求最终至少支持：

```bash
./wave_reader FILE.root --inspect
```

检查 ROOT 文件结构。

```bash
./wave_reader FILE.root --list
```

列出所有 Canvas。

```bash
./wave_reader FILE.root --draw 3
```

提取并绘制 `ChainNum_3`。

```bash
./wave_reader FILE.root --analyze 3
```

分析 `ChainNum_3`。

```bash
./wave_reader FILE.root --analyze-all
```

分析全部 Canvas。

```bash
./wave_reader FILE.root --make-template 1 3 5 7
```

使用指定 waveform 建立模板。

```bash
./wave_reader FILE.root --template pulse_template.root --analyze-all
```

使用已有 template 做分析。

---

# 23. 特别注意：一个 Canvas 内可能有多个 waveform

不能认为：

```text
一个 Canvas = 一个 TGraph
```

一个 Canvas 可能包含：

```text
DSSD waveform
SSD waveform
MWPC waveform
fitted curve
marker
legend
```

所以第一次 `--inspect` 必须把所有 primitive 打出来。

之后程序允许通过：

```text
object title
object name
pad name
```

选择真正需要分析的 waveform。

建议命令：

```bash
./wave_reader FILE.root --object-pattern "DSSD"
```

或者：

```bash
./wave_reader FILE.root --graph-index 0
```

如果无法自动判断波形对象，不要静默选择。

程序必须把候选对象全部列出。

---

# 24. 对未知 ROOT 文件结构要尽量鲁棒

原则：

```text
先 introspection
后 analysis
```

不要在代码最开始写死：

```cpp
TGraph* g = (TGraph*)c->GetPrimitive("Graph");
```

因为不同 Canvas 的 primitive 名字可能变化。

应该递归自动寻找所有可能 waveform 对象。

---

# 25. 第一版程序的最低完成标准

第一阶段暂时不要急着做 pile-up fit。

第一版必须先可靠完成以下功能：

1. 从命令行读取 ROOT 文件名；
2. 打开 `TFile`；
3. 遍历所有 `TCanvas`；
4. 递归遍历所有 `TPad`；
5. 找出所有 `TGraph / TH1`；
6. 提取 `(x,y)` 数值；
7. 解析 `ChainNum / run / map_entry`；
8. 输出 ROOT 文件内部结构；
9. 能重新绘制提取的 waveform；
10. 能把 waveform 保存为 CSV。

例如：

```bash
./wave_reader wave_242Fm.root --export-waveforms
```

生成：

```text
waveforms/
  ChainNum_1_Graph0.csv
  ChainNum_1_Graph1.csv
  ChainNum_2_Graph0.csv
  ...
```

每个 CSV：

```text
sample,x,y
0,0.000,8193
1,0.005,8191
2,0.010,8192
...
```

---

# 26. 推荐代码结构

建议拆成：

```text
wave_reader/
│
├── CMakeLists.txt
│
├── include/
│   ├── Waveform.h
│   ├── RootCanvasReader.h
│   ├── Baseline.h
│   ├── PulseFinder.h
│   ├── PulseTemplate.h
│   └── PileupFitter.h
│
├── src/
│   ├── main.cpp
│   ├── RootCanvasReader.cpp
│   ├── Baseline.cpp
│   ├── PulseFinder.cpp
│   ├── PulseTemplate.cpp
│   └── PileupFitter.cpp
│
└── README.md
```

如果希望先快速验证，也可以先写一个单文件版本：

```text
wave_reader.cpp
```

先保证功能正确，再模块化。

---

# 27. ROOT 依赖

建议至少使用：

```cpp
#include <TFile.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TList.h>
#include <TClass.h>
#include <TROOT.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
```

编译：

```bash
g++ -std=c++17 -O2 wave_reader.cpp \
    $(root-config --cflags --libs) \
    -o wave_reader
```

---

# 28. 调试输出

程序增加：

```bash
--verbose
```

例如：

```text
[INFO] Open file: wave_242Fm.root
[INFO] Found canvas: ChainNum_1
[INFO] Parsed run = 6
[INFO] Parsed map_entry = 43971794
[INFO] Enter pad: pad1
[INFO] Found TGraph: Graph, N=4096
[INFO] Extract waveform successfully.
```

错误使用：

```text
[WARNING]
[ERROR]
```

例如：

```text
[WARNING] ChainNum_7 contains no TGraph or TH1 waveform.
```

---

# 29. 绝对不要做的事情

不要：

```text
1. 把 ChainNum 数量写死；
2. 把 Graph 名字写死；
3. 假定 Canvas 只有一个 pad；
4. 假定 waveform 一定是 TH1；
5. 假定 waveform 一定是 TGraph；
6. 假定 pulse polarity；
7. 一开始就用复杂神经网络；
8. 一开始就同时拟合十几个 nonlinear 参数；
9. 只保存图片而不保存数值结果；
10. 修改原始 ROOT 文件。
```

原始文件始终只读：

```cpp
TFile::Open(filename, "READ")
```

---

# 30. 开发顺序

严格按照以下顺序开发。

## Step 1

只做：

```text
ROOT file introspection
```

确认 `wave_242Fm.root` 中：

```text
TCanvas
TPad
TGraph / TH1
```

的真实结构。

## Step 2

提取所有 waveform 数值并导出 CSV。

## Step 3

重新绘制 waveform，确认与原 Canvas 完全一致。

## Step 4

baseline subtraction + noise estimation。

## Step 5

pulse candidate finder。

## Step 6

建立 clean single-pulse template。

## Step 7

single-pulse fit。

## Step 8

double-pulse variable-projection fit。

## Step 9

扫描 $\Delta t$，输出：

```text
chi2(Delta_t)
```

## Step 10

single vs double pulse model comparison。

---

# 31. 最终物理目标

最终程序不是简单“画波形”，而是为超短寿命核素的 pulse pile-up 提供：

```text
waveform extraction
↓
pulse identification
↓
sub-sample timing
↓
energy reconstruction
↓
single/double pulse hypothesis test
↓
pile-up resolving limit
```

最终每个事件希望得到：

$$
E_1,
\quad
E_2,
\quad
\Delta t,
\quad
\sigma_{E_1},
\quad
\sigma_{E_2},
\quad
\sigma_{\Delta t},
$$

以及：

$$
\Delta\chi^2,
\quad
\Delta\mathrm{BIC}.
$$

这里最重要的设计原则是：

> **先可靠读取 ROOT Canvas 中的真实 waveform，再逐步增加分析算法。不要在不知道 Canvas 内部真实 primitive 结构之前写死读取逻辑。**

