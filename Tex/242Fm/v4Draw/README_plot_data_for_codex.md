# `plot_data.root` 读取说明

本文档用于指导 Codex 读取 `plot_data.root` 并绘制后续图像。

## 1. 文件用途

`plot_data.root` 已经完成原始 ROOT 数据的筛选与拼接。

**不要重新读取原始 `runXXXXX_Rechain.root` 文件。**

后续绘图只需要读取：

```text
plot_data.root
```

即可。

---

## 2. 数据来源规则

在生成 `plot_data.root` 时，已经按照能量区域自动选择不同的数据源：

```text
SumE < 50 MeV
    -> ../../dataFinal/runXXXXX_Rechain.root

SumE >= 50 MeV
    -> ../../data0/runXXXXX_Rechain.root
```

即以

```text
50 MeV = 50000 keV
```

为分界。

因此，后续绘图程序 **不要再次根据 50 MeV 切换输入文件**。

`tree_full_green` 和 `tree_full_red` 已经是正确拼接完成的数据。

---

## 3. Run 分组

实验数据分为两组：

```text
GREEN:
Run 4-11

RED:
Run 12-40
```

在绘图中：

```text
green -> Run 4-11
red   -> Run 12-40
```

---

## 4. ROOT 文件中的主要对象

`plot_data.root` 中包含以下四个主要 `TTree`：

```text
tree_low_green
tree_low_red
tree_full_green
tree_full_red
```

以及对应的 `TGraph`：

```text
g_low_green
g_low_red
g_full_green
g_full_red
```

推荐绘图时优先读取 `TTree`，因为 `TTree` 保留了完整变量。

---

# 5. 三张图应该读取什么

## Figure 1：RED，7.0-9.5 MeV

读取：

```text
tree_low_red
```

绘图范围：

```text
X: 7000-9500 keV
```

推荐：

```text
X = SumE
Y = log10_DeltaT_ns
```

对应：

```text
Run 12-40
```

---

## Figure 2：GREEN，7.0-9.5 MeV

读取：

```text
tree_low_green
```

绘图范围：

```text
X: 7000-9500 keV
```

推荐：

```text
X = SumE
Y = log10_DeltaT_ns
```

对应：

```text
Run 4-11
```

---

## Figure 3：RED + GREEN，0-250 MeV

读取：

```text
tree_full_green
tree_full_red
```

绘图范围：

```text
X: 0-250000 keV
```

两组数据画在同一个坐标系中：

```text
tree_full_green -> GREEN
tree_full_red   -> RED
```

推荐：

```text
X = SumE
Y = log10_DeltaT_ns
```

注意：

`tree_full_green` 和 `tree_full_red` 已经完成：

```text
0-50 MeV       -> dataFinal
50-250 MeV     -> data0
```

的拼接。

**禁止重新从原始文件读取或重新拼接。**

---

# 6. TTree 中的 Branch

四个 tree 都具有相同的 branch。

## `SumE`

```text
SumE
```

类型：

```text
Double_t
```

单位：

```text
keV
```

因此：

```text
7000 keV   = 7 MeV
9500 keV   = 9.5 MeV
50000 keV  = 50 MeV
250000 keV = 250 MeV
```

---

## `DeltaT_ns`

```text
DeltaT_ns
```

原始时间差。

单位：

```text
ns
```

---

## `DeltaT_s`

```text
DeltaT_s
```

单位：

```text
s
```

定义：

```text
DeltaT_s = DeltaT_ns / 1e9
```

---

## `log10_DeltaT_ns`

```text
log10_DeltaT_ns
```

定义：

```text
log10_DeltaT_ns = log10(DeltaT_ns)
```

这是目前推荐用于散点图 Y 坐标的变量。

例如：

```text
DeltaT = 10 ns      -> log10_DeltaT_ns = 1
DeltaT = 1 us       -> log10_DeltaT_ns = 3
DeltaT = 1 ms       -> log10_DeltaT_ns = 6
DeltaT = 1 s        -> log10_DeltaT_ns = 9
DeltaT = 10 s       -> log10_DeltaT_ns = 10
```

---

## `log10_DeltaT_s`

```text
log10_DeltaT_s
```

定义：

```text
log10_DeltaT_s = log10(DeltaT_s)
```

例如：

```text
1 ns   -> -9
1 us   -> -6
1 ms   -> -3
1 s    ->  0
10 s   ->  1
```

如果最终图的纵轴希望直接表示秒的数量级，可以使用这个 branch。

---

## `pass_10s`

```text
pass_10s
```

类型：

```text
Bool_t
```

含义：

```text
true  -> DeltaT < 10 s
false -> DeltaT >= 10 s
```

如果某张图只需要 10 s 内的事件，可以进一步筛选：

```cpp
pass_10s == true
```

---

## `source`

```text
source
```

类型：

```text
Int_t
```

含义：

```text
source = 0
    -> ../../dataFinal/

source = 1
    -> ../../data0/
```

这个 branch 主要用于检查数据来源。

正常绘图不需要根据 `source` 再做筛选。

---

# 7. ROOT C++ 读取示例

```cpp
TFile *f = TFile::Open("plot_data.root");

TTree *red =
    (TTree *)f->Get("tree_low_red");

TTree *green =
    (TTree *)f->Get("tree_low_green");

TTree *full_red =
    (TTree *)f->Get("tree_full_red");

TTree *full_green =
    (TTree *)f->Get("tree_full_green");
```

例如查看 RED 的低能散点：

```cpp
red->Draw(
    "log10_DeltaT_ns:SumE",
    "",
    ""
);
```

GREEN：

```cpp
green->Draw(
    "log10_DeltaT_ns:SumE",
    "",
    ""
);
```

---

# 8. ROOT RDataFrame 读取示例

也可以使用：

```cpp
ROOT::RDataFrame df_red(
    "tree_low_red",
    "plot_data.root"
);

ROOT::RDataFrame df_green(
    "tree_low_green",
    "plot_data.root"
);
```

读取完整能区：

```cpp
ROOT::RDataFrame df_full_red(
    "tree_full_red",
    "plot_data.root"
);

ROOT::RDataFrame df_full_green(
    "tree_full_green",
    "plot_data.root"
);
```

---

# 9. Python / uproot 读取示例

如果 Codex 使用 Python 绘图，推荐使用 `uproot`：

```python
import uproot

f = uproot.open("plot_data.root")

red = f["tree_low_red"]
green = f["tree_low_green"]

full_red = f["tree_full_red"]
full_green = f["tree_full_green"]
```

读取数组：

```python
red_data = red.arrays(
    [
        "SumE",
        "DeltaT_ns",
        "DeltaT_s",
        "log10_DeltaT_ns",
        "log10_DeltaT_s",
        "pass_10s",
        "source",
    ],
    library="np",
)
```

例如：

```python
x = red_data["SumE"]
y = red_data["log10_DeltaT_ns"]
```

然后直接：

```python
ax.scatter(x, y)
```

---

# 10. 推荐的绘图读取逻辑

Codex 应按照以下逻辑绘图。

```text
Figure 1
--------
input:
    tree_low_red

X:
    SumE

X range:
    7000-9500 keV

Y:
    log10_DeltaT_ns

group:
    RED / Run 12-40
```

```text
Figure 2
--------
input:
    tree_low_green

X:
    SumE

X range:
    7000-9500 keV

Y:
    log10_DeltaT_ns

group:
    GREEN / Run 4-11
```

```text
Figure 3
--------
input:
    tree_full_red
    tree_full_green

X:
    SumE

X range:
    0-250000 keV

Y:
    log10_DeltaT_ns

groups:
    GREEN = Run 4-11
    RED   = Run 12-40
```

---

# 11. 关于 Y 轴的重要说明

当前文件同时保存了：

```text
DeltaT_ns
DeltaT_s
log10_DeltaT_ns
log10_DeltaT_s
```

因此有两种合理的画法。

### 方法 A：直接绘制 log10 后的时间

推荐用于当前三张散点图：

```text
Y = log10_DeltaT_ns
```

此时 Y 轴本身是普通线性轴，但变量已经做了 `log10`。

例如纵轴可以标为：

```text
log10(Delta t / ns)
```

---

### 方法 B：绘制真实时间，然后设置 log Y

例如：

```text
Y = DeltaT_s
```

然后：

```cpp
gPad->SetLogy();
```

或 Python：

```python
ax.set_yscale("log")
```

此时纵轴可以直接标：

```text
Delta t (s)
```

**不要同时对 `log10_DeltaT_ns` 再设置 log Y。**

也就是说：

```text
正确：
DeltaT_s + log-scale axis

或者：

log10_DeltaT_ns + linear axis
```

不要：

```text
log10_DeltaT_ns + log-scale axis
```

---

# 12. Codex 必须遵守的规则

1. 只读取：

```text
plot_data.root
```

2. 不重新读取：

```text
../../dataFinal/runXXXXX_Rechain.root
../../data0/runXXXXX_Rechain.root
```

3. 不重新进行 50 MeV 数据拼接。

4. 低能图直接读取：

```text
tree_low_red
tree_low_green
```

5. 全能区图直接读取：

```text
tree_full_red
tree_full_green
```

6. `SumE` 单位始终为：

```text
keV
```

7. GREEN 与 RED 的定义始终为：

```text
GREEN = Run 4-11
RED   = Run 12-40
```

8. 如果使用：

```text
log10_DeltaT_ns
```

则 Y 轴使用普通线性坐标。

9. 如果需要真正的 logarithmic Y axis，则使用：

```text
DeltaT_s
```

并设置：

```text
log Y
```

---

# 13. 最简 Codex 提示

如果只需要快速告诉 Codex 如何读取，可以直接使用下面这段：

```text
Read only plot_data.root.

The file already contains all selected and merged data.

Trees:

tree_low_red:
    Run 12-40
    SumE = 7.0-9.5 MeV

tree_low_green:
    Run 4-11
    SumE = 7.0-9.5 MeV

tree_full_red:
    Run 12-40
    SumE = 0-250 MeV

tree_full_green:
    Run 4-11
    SumE = 0-250 MeV

For tree_full_red/tree_full_green the data source has already been merged:
    SumE < 50 MeV  -> dataFinal
    SumE >= 50 MeV -> data0

Do not read the original Rechain ROOT files again.

Available branches:
    SumE                 [keV]
    DeltaT_ns            [ns]
    DeltaT_s             [s]
    log10_DeltaT_ns
    log10_DeltaT_s
    pass_10s
    source

source:
    0 = dataFinal
    1 = data0

For scatter plots either use:

    X = SumE
    Y = log10_DeltaT_ns

with a linear Y axis,

or use:

    X = SumE
    Y = DeltaT_s

with a logarithmic Y axis.

Do NOT apply a logarithmic axis to log10_DeltaT_ns.
```
