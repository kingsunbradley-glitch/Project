# SSD 内置α源刻度程序 (V2.0 - 健壮版)

这是一个基于 ROOT 和 TSpectrum 的 C++ 程序工作流，用于对 SSD (Silicon Strip Detector) 探测器进行能量刻度。该版本在专业版的基础上进行了全面的代码审查和重构，显著提升了程序的健壮性、安全性和用户友好性。

## 核心组件

1. `config.h`: **唯一的配置中心**。所有文件路径、物理参数、程序行为都由此文件控制。
2. `SSD_calibration.cpp`: (第1步) 分析原始数据，执行背景扣除、自动寻峰、线性拟合，最终生成刻度参数文件和一系列诊断图。
3. `apply_calibrationv1.cpp`: (第2步) 读取参数文件，将刻度应用于原始数据，并生成一个全新的、已刻度的 ROOT 文件。
4. `Makefile`: 自动化编译和执行脚本。

## 主要特性

### 集中化配置
用户仅需修改 `config.h` 文件即可控制整个工作流，无需接触 `.cpp` 源码。

### 高健壮性
- 全面的错误处理机制，能捕获文件、TTree、TBranch的丢失等常见问题。
- 通过智能指针管理内存，杜绝内存泄漏风险。
- 增加了数组边界检查，防止因数据异常导致程序崩溃。
- 启动时进行参数验证，确保配置合理。

### 清晰的输出
增加了详细的运行日志和一份最终的刻度效果总结报告。

### 自动化流程
使用 Makefile 实现一键编译、运行和清理。

---

## 使用方法

### 1. 环境准备
确保您的系统中已正确安装 ROOT，并且 `root-config` 命令在终端中可用。

### 2. **参数配置 (关键步骤)**
打开 `config.h` 文件，根据您的实验数据和文件存放路径，仔细修改以下几个部分：
- **文件路径配置**: `input_root_file`, `calibration_param_file` 等。
- **TTree 名称**: `tree_name`。
- **刻度峰配置**: `peaks_of_interest` 向量，定义参考峰的能量和大致道址范围。
- **寻峰和物理参数**: 根据需要微调 `peak_search_sigma`, `peak_search_threshold` 等。

### 3. 编译
将所有文件 (`config.h`, `SSD_calibration.cpp`, `apply_calibrationv1.cpp`, `Makefile`) 放在同一目录下。打开终端，执行：
```bash
make
```
或
```bash
make all
```
该命令会自动编译生成 `SSD_calibration` 和 `apply_calibration` 两个可执行文件。

### 4. 运行
推荐使用 `make` 提供的 `run` 命令，它会自动按顺序执行整个工作流：
```bash
make run
```
程序将首先运行 `./SSD_calibration` 生成参数，随后运行 `./apply_calibration` 应用参数。

### 5. 清理
当您需要重新编译或清理生成的文件时，执行：
```bash
make clean
```
该命令会删除所有生成的可执行文件和中间文件。

---

## 输出文件详解

### `calibration_plots_pro/` 目录:
- `pos_*_1_background.png`: 显示原始能谱和 `TSpectrum` 估算的本底曲线。
- `pos_*_2_peak_search.png`: 显示扣除本底后的能谱，并标记找到的峰（红色用于刻度，绿色为其他）。
- `pos_*_3_linear_fit.png`: 显示刻度点、线性拟合结果及拟合优度。

### `..._params_pro.txt` 文件:
- 包含四列的刻度参数：`位置(Pos)` `斜率(k)` `截距(b)` `置信度(confidence)`。
- `confidence` 代表用于拟合的峰的个数。**2或更大**表示拟合成功；**1**表示只找到一个峰无法拟合；**0**表示未找到可用峰或拟合失败；**-1**表示该位置的数据量不足。

### `..._SSDcal_pro.root` 文件:
- **最终的数据产品**。该文件是原始 ROOT 文件的完整克隆，但 `SSD_E` 分支的内容已被**替换为刻度后**的能量值（单位：keV）。

---

## 故障排查与参数微调

如果刻度效果不佳，请首先检查 `calibration_plots_pro/` 中的诊断图，然后参考以下建议调整 `config.h`：

### 问题：没有找到足够的刻度峰 (used_markers 数量 < 2)
- 检查 `pos_*_2_peak_search.png` 图，确认能谱中是否存在清晰的峰。
- 放宽 `peaks_of_interest` 中的 `win_min` 和 `win_max` 范围。
- 适当降低 `peak_search_threshold` (例如从 0.05 降到 0.03) 来提高寻峰灵敏度。
- 检查 `MIN_HIST_ENTRIES`，确保统计量足够。

### 问题：找到了太多无关的杂峰
- 适当提高 `peak_search_threshold` (例如从 0.05 提高到 0.08)。
- 适当增大 `peak_search_sigma` (单位: bin)，如果峰形本身较宽的话。

### 问题：背景扣除效果不佳
- 在 `pos_*_1_background.png` 中，如果背景曲线与实际谱的连续部分偏离过大，可以尝试增减 `background_iterations` 的值。

---

## 后期优化思路

### 利用时间窗去除干扰峰

在处理计数率较高的实验数据时，能谱中可能混入大量由随机事件构成的本底，这会干扰甚至淹没真实的α刻度峰。通过对事件的时间戳施加筛选条件（即"卡时间窗"），可以极大地提高信噪比。

例如，如果我们知道刻度用的α粒子来自于某个特定的衰变链，我们可以只选择那些在注入事件之后特定时间范围内发生的事件。

#### 实现方法

**1. 修改 `config.h` 文件**

在文件末尾（ErrorCode 枚举之前）添加一个新的配置项，用于定义全局的时间削减条件。

```cpp
// --- 粘贴到 config.h ---

// --- 7. 附加全局削减 ---
// 应用于TChain::Draw的全局事件削减条件 (ROOT TCut format)
// 例如，只选择第一个alpha衰变时间小于5秒的事件: "Delta_Ts[1] > 0 && Delta_Ts[1] < 5000e6"
// 留空字符串 "" 则不应用任何额外削减。
const char* global_event_cut = "Delta_Ts[1] > 0 && Delta_Ts[1] < 5000e6";
```

**2. 修改 `SSD_calibration.cpp` 文件**

在 `calibrate_ssd_professional` 函数的 for 循环内部，找到绘制直方图的部分，将其修改为如下代码，以合并基础的位置削减和全局的时间削减。

```cpp
// --- 找到并替换 SSD_calibration.cpp 中的对应部分 ---

// ... 在 for (int pos = 0; pos < NUM_SSD_POS; ++pos) 循环内部 ...

// --- 原有代码 ---
// auto h_sum_E = std::make_unique<TH1D>(...);
// TCut cut = TCut(TString::Format("SSD_E > 0 && DSSD_E > 0 && TMath::Nint(SSD_Pos) == %d", pos));
// chain->Draw(TString::Format("DSSD_E + SSD_E >> %s", h_sum_E->GetName()), cut, "gooff");
// -----------------

// --- 替换为以下代码 ---
auto h_sum_E = std::make_unique<TH1D>(TString::Format("h_sum_E_pos%d", pos), TString::Format("Energy Spectrum for Pos %d;Energy Sum (channels);Counts", pos), 700, 3000, 10000);

// 基础的位置削减
TCut pos_cut = TCut(TString::Format("SSD_E > 0 && DSSD_E > 0 && TMath::Nint(SSD_Pos) == %d", pos));

// 合并全局事件削减
TCut final_cut = pos_cut;
if (std::string(global_event_cut).length() > 0) {
    final_cut = final_cut && TCut(global_event_cut);
    if (enable_debug_output) {
        std::cout << "    Applying combined cut for Pos " << pos << ": " << final_cut.GetTitle() << std::endl;
    }
}

// 使用最终的削减条件绘制能谱
chain->Draw(TString::Format("DSSD_E + SSD_E >> %s", h_sum_E->GetName()), final_cut, "gooff");
// -------------------------

// ... 后续代码不变 ...
```

完成以上两步修改并重新用 `make` 编译后，`SSD_calibration` 程序在生成能谱时就会自动应用您在 `config.h` 中定义的时间窗条件了。