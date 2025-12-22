# DSSD 归一化与刻度程序使用说明（README）

> 本说明文档介绍整套 DSSD 处理流程：**Step 0 预筛选 → Step 1 归一化 → Step 2 绝对刻度（内源刻度）与诊断**，以及相关配置、输入输出文件和常见问题。  
> 可以直接把本段内容作为工程的 README 或任务说明书使用。

---

## 1. 程序整体结构

本套程序主要由以下源文件组成：

- Step 0：预筛选  
  - `Preselect_Main.cpp`  
    从原始 map ROOT 文件中筛出 *implant-only* 与 *decay-only* 的小 ROOT 文件，减小后续处理的数据量。

- Step 1：归一化（Normalization）  
  - `Normalize_Main.cpp`  
    利用 implant+decay 的符合事件，对 X/Y 面各条进行条间归一化，输出归一化系数表和诊断 ROOT 文件。

- Step 2：绝对刻度与诊断（Calibration）  
  - `Calibrator.h` / `Calibrator.cpp`  
    - 读入归一化系数  
    - 生成 X / Y / YH 总能谱  
    - 对平面做整体绝对刻度  
    - 对每条做 FWHM 诊断  
    - 输出最终的能量刻度表 `ener_cal_Recal.dat`
  - `Calibrator_Main.cpp`  
    - 程序入口，创建 `Calibrator` 对象并执行 `Run()`。

- 统一配置  
  - `Config.h`  
    所有运行区间、文件路径、树名、归一化与刻度参数、参考峰和诊断峰等都集中在这里配置。

- 辅助  
  - `makefile`（如果存在）  
    提供统一编译入口，可以通过 `make` 命令一键编译各个可执行程序。

---

## 2. 依赖环境

1. 操作系统  
   推荐 Linux 环境。

2. 编译器  
   - `g++`，需支持 C++20（例如编译选项使用 `-std=c++20` 或 `-std=gnu++20`）。

3. ROOT 环境  
   需要以下组件：
   - 基本 I/O：`TFile`, `TTree`, `TChain`
   - 直方图和拟合：`TH1D`, `TH2D`, `TProfile`, `TF1`, `TFitResult`
   - 框架：`TApplication`, `ROOT::RDataFrame`
   - 峰搜索与本底扣除：`TSpectrum`

   运行标定程序时，会在 `Calibrator_Main.cpp` 中初始化 `TApplication`。

4. 多线程  
   - 使用 `ROOT::EnableImplicitMT(NUM_THREADS)` 启动多线程。  
   - 线程数在 `Config.h` 中通过 `NUM_THREADS` 设置。

---

## 3. 配置文件 `Config.h`

所有关键参数都集中在 `Config.h` 中，需要根据实际实验环境进行一次性修改。

### 3.1 线程与输入文件

示意：
``` C++
    inline constexpr int NUM_THREADS = 10;

    inline constexpr char INPUT_DIR_PATTERN[] =
        "/path/to/data/SS032%05d_map.root";
    inline constexpr int RUN_START = 2;
    inline constexpr int RUN_END   = 178;

    inline constexpr char TREE_NAME[] = "tr_map";
```
含义：

- `NUM_THREADS`：RDataFrame 隐式多线程线程数，-1的话默认最大，谨慎使用-1。
- `INPUT_DIR_PATTERN`：原始 map 文件路径模式（使用 `sprintf` 用 run 号填充）。
- `RUN_START` / `RUN_END`：run 编号范围（闭区间）。
- `TREE_NAME`：原始 ROOT 文件中的 TTree 名称（例如 `"tr_map"`）。

### 3.2 中间与输出文件名

示意：
``` C++
    inline constexpr char DECAY_ONLY_FILE[]   = "SS032_decay_only.root";
    inline constexpr char IMPLANT_ONLY_FILE[] = "SS032_implant_only.root";

    inline constexpr char NORM_PARAM_FILE[]   = "./SS032_Normalize_Params.txt";
    inline constexpr char NORM_DIAG_ROOT[]    = "Diagnose_Normalize.root";
    inline constexpr char CALIB_DIAG_ROOT[]   = "Diagnose_Calibration.root";

    inline constexpr char OUTPUT_DAT_FILE[]   = "ener_cal_Recal.dat";
```
说明：

- `DECAY_ONLY_FILE` / `IMPLANT_ONLY_FILE`：预筛选输出的小 ROOT 文件。
- `NORM_PARAM_FILE`：归一化系数文本表。
- `NORM_DIAG_ROOT`：归一化诊断 ROOT 文件。
- `CALIB_DIAG_ROOT`：最终刻度诊断 ROOT 文件。
- `OUTPUT_DAT_FILE`：最终能量刻度参数与逐条 FWHM 结果。

可根据工程目录结构进行适当修改。

### 3.3 通道定义与归一化参数

示意：
``` C++
    #define TOTAL_CH 224

    inline constexpr int NUM_DSSDX_POS = 128; // X 面条数：0–127
    inline constexpr int NUM_DSSDY_POS = 48;  // Y 面条数：128–175
    // YH 面：176–223
```
归一化时使用的能量范围及条件：
``` C++
    inline constexpr double NORM_MIN_E = 4000.0;
    inline constexpr double NORM_MAX_E = 14000.0;
    inline constexpr double NORM_E_DIFF = 200.0;

    inline constexpr int NORM_MIN_ENTRIES = 500;
    inline constexpr int NORM_REF_MIN_COUNTS = 2000;

    // 可选：指定参考条
    // inline constexpr int REF_STRIP_X = 39;
    // inline constexpr int REF_STRIP_Y = 25;
```
说明：

- `NORM_MIN_E` / `NORM_MAX_E`：归一化选用的能量范围。
- `NORM_E_DIFF`：归一化时要求 X/Y 能量差的最大允许值。
- `NORM_MIN_ENTRIES`：每条进行归一化拟合所需的最小统计量。
- `NORM_REF_MIN_COUNTS`：作为参考条所需的最小 hit 数。
- 若未显式指定 `REF_STRIP_X` / `REF_STRIP_Y`，则自动选取统计量最大的条作为参考条。
- 建议X面不选择39 89附近的，在Veto峰中有没有筛选干净的信号
- XY尽量选择FWHM比较好的条，可以尝试运行make Cali查看每一条的FWHM。

### 3.4 平面绝对刻度峰与诊断峰

示意定义：
``` C++
    struct PeakDef {
        double energy;   // keV
        double window;   // keV
    };

    inline const std::vector<PeakDef> REF_PEAKS = {
        {6113, 50.0},//242Cm
        {7686, 50.0},//214Po
        {8784, 50.0}//212Po
    };

    inline const std::vector<PeakDef> DIAG_PEAKS = {
        {6113, 40.0},//242Cm
        {7686, 40.0},//214Po
        {5307, 40.0}//212Po
    };
```
说明：

- `REF_PEAKS`：用于 X/Y/YH 三个平面的整体绝对刻度。
- `DIAG_PEAKS`：用于刻度后诊断，统计峰位偏差 RMS 和最大 FWHM。

### 3.5 总谱直方图参数与本底扣除

示意：
``` C++
    inline constexpr int HIST_BINS = 2000;
    inline constexpr double HIST_MIN = 4000.0;
    inline constexpr double HIST_MAX = 14000.0;

    inline constexpr bool ENABLE_BKG_SUB = true;
    inline constexpr int BKG_ITERATIONS = 20;
``` 
说明：

- `HIST_BINS` / `HIST_MIN` / `HIST_MAX`：平面总谱及条谱直方图的 bin 设置。
- `ENABLE_BKG_SUB`：是否启用本底扣除。
- `BKG_ITERATIONS`：`TSpectrum::Background` 的迭代次数。

---

## 4. 编译方式

### 4.1 使用命令行直接编译（示例）

在已加载 ROOT 环境的前提下，可以使用如下命令：
``` shell
    # Step 0：预筛选
    g++ -O2 -std=gnu++20 -o Preselect Preselect_Main.cpp `root-config --cflags --glibs`

    # Step 1：归一化
    g++ -O2 -std=gnu++20 -o Normalize Normalize_Main.cpp `root-config --cflags --glibs`

    # Step 2：绝对刻度与诊断
    g++ -O2 -std=gnu++20 -o DSSD_Calib \
        Calibrator_Main.cpp Calibrator.cpp \
        `root-config --cflags --glibs` -lSpectrum
``` 
说明：

- `root-config --cflags --glibs`：自动添加 ROOT 所需的编译与链接选项。
- `-lSpectrum`：用于链接 `TSpectrum` 所在的库。

### 4.2 使用 makefile 编译（推荐）

如果目录中提供了 `makefile`，则可以直接在工程目录下执行：
``` shell
    make
``` 
常见的目标：
``` shell
    make Preselect
    make Normalize
    make DSSD_Calib
    make all
    make clean
``` 
这样可以简化日常编译过程。

---

## 5. 运行流程总览

标准处理流程分为三步，建议严格按顺序执行：
``` shell
     make Pres      # Step 0：生成 implant-only / decay-only ROOT 文件
     make Norm      # Step 1：根据 implant+decay 符合事件做条间归一化
     make Cali      # Step 2：平面绝对刻度 + 逐条 FWHM 诊断 + 输出最终刻度表
``` 
在每个步骤之间建议检查输出文件和诊断 ROOT 文件，以便及时修正问题。

---

## 6. Step 0：预筛选原始数据（Preselect）

对应程序：`Preselect_Main.cpp`

### 6.1 主要作用

- 读取原始 map 文件（由 `INPUT_DIR_PATTERN` + `RUN_START/RUN_END` 决定）。
- 按事件类型和多重性条件，筛选出：
  - implant-only 事件 → `IMPLANT_ONLY_FILE`
  - decay-only 事件 → `DECAY_ONLY_FILE`

通过过滤无需参与后续分析的事件，显著减小数据量。

### 6.2 典型筛选条件（示意）

- Implant（注入）事件：

  - `MWPC_mul > 0`
  - `DSSDX_mul == 1 && DSSDY_mul == 1`
  - `Veto_mul == 0 && SSD_mul == 0`

- Decay（衰变）事件：

  - `MWPC_mul == 0`
  - `DSSDX_mul == 1 && DSSDY_mul == 1`
  - `Veto_mul == 0 && SSD_mul == 0`

并只保留与 DSSD 能量与通道相关的分支，以压缩文件体积。

### 6.3 运行方式
``` shell
    make Pres
``` 
运行结束后，目录中应出现：

- `SS032_implant_only.root`（或 `IMPLANT_ONLY_FILE` 指定的文件名）
- `SS032_decay_only.root`（或 `DECAY_ONLY_FILE` 指定的文件名）

建议检查一下文件大小是否与预期相符，以确认筛选是否成功。

---

## 7. Step 1：条间归一化（Normalize）

对应程序：`Normalize_Main.cpp`

### 7.1 主要作用

- 读取 `IMPLANT_ONLY_FILE` 和 `DECAY_ONLY_FILE`。
- 利用 X/Y 面的符合事件，在给定能量范围和能量差条件下，构建 `(ChX, EX)` 与 `(ChY, EY)` 的二维分布。
- 对每条 X/Y 条，相对于参考条做 profile 拟合，得到线性归一化系数 `(k_norm, b_norm)`。
- 对 Y 面使用桥接方法（bridge correction），将 X/Y 统一到同一能量标度。
- 输出：
  - 文本归一化参数表 `NORM_PARAM_FILE`
  - 诊断 ROOT 文件 `NORM_DIAG_ROOT`

### 7.2 核心选择条件（示意）

对符合事件施加以下限制：

- DSSD 多重性：

  - `DSSDX_mul == 1 && DSSDY_mul == 1`
  - `Veto_mul == 0 && SSD_mul == 0`

- 能量范围：

  - `EX ∈ [NORM_MIN_E, NORM_MAX_E]`
  - `EY ∈ [NORM_MIN_E, NORM_MAX_E]`

- X/Y 能量接近：

  - `|EX - EY| < NORM_E_DIFF`

在这些条件下，构建每条与参考条之间的 profile，用于线性拟合。

### 7.3 参考条选择与桥接

1. 自动或手动选择参考条

   - 程序会统计每条的 hits 数目，自动选取统计量最大的条作为 `refX` / `refY`。
   - 若在 `Config.h` 中定义了 `REF_STRIP_X` / `REF_STRIP_Y`，则优先使用用户指定的参考条。

2. X 面归一化

   - 对每条 X 条，使用与 `refY` 的关系（或与 `refX` 的关系）建立线性拟合，得到 `(k_norm_x, b_norm_x)`。

3. Y 面归一化与桥接（bridge）

   - 首先相对于 `refX` 建立线性关系。
   - 再利用 `refX` 与 `refY` 的交叉关系，对 Y 面系数做桥接变换，使其转换到与 `refY` 相同的标度。

### 7.4 输出文件格式

- 文本归一化参数：`NORM_PARAM_FILE`，示意形式：

      # StripID   k_norm        b_norm
      0          k0            b0
      1          k1            b1
      ...
      127        k127          b127
      128        k128          b128
      ...
      175        k175          b175

  文件头或末尾可能附带“坏条”标记（统计量不够或拟合失败的条）。

- 诊断 ROOT 文件：`NORM_DIAG_ROOT`
  - X/Y 各条的 hits 统计直方图
  - 各条 2D 分布 (EX vs EY) 图
  - Profile 曲线与拟合结果
  - 参考条桥接关系诊断图

### 7.5 运行方式

    make Norm

运行结束后，建议：

- 打开 `SS032_Normalize_Params.txt` 检查坏条数量是否合理。
- 用 ROOT 打开 `Diagnose_Normalize.root`，确认：
  - 大部分条的 profile 拟合较为线性；
  - 桥接后 X/Y 的能量关系合理。

---

## 8. Step 2：绝对刻度与 FWHM 诊断（Calibrator）

对应：`Calibrator.h` / `Calibrator.cpp` + `Calibrator_Main.cpp`

### 8.1 总体流程（Calibrator::Run 概要）

1. 读取归一化参数

   - 从 `NORM_PARAM_FILE`（如 `SS032_Normalize_Params.txt`）中读取每条 `(k_norm, b_norm)`。

2. 加载 decay-only 数据

   - 用 `TChain` 读取 `DECAY_ONLY_FILE`（如 `SS032_decay_only.root`）。

3. 生成归一化后的总谱

   - 对每个 decay 事件：
     - 根据通道号映射到 strip ID：
       - X：`DSSDX_Ch ∈ [0, 127]` → strip 0–127
       - Y：`DSSDY_Ch ∈ [0, 47]` → strip 128–175
       - YH：`DSSDYH_Ch ∈ [0, 47]` → strip 176–223
     - 使用归一化系数 `(k_norm, b_norm)` 对对应条的能量做一次线性变换。
     - 灵活构建 X / Y / YH 各自的总能谱直方图。

4. 平面整体绝对刻度

   - 对 X / Y / YH 三个平面总谱：
     - 在 `REF_PEAKS` 指定的能量附近搜索峰（如 6113 / 7686 / 8784 keV）。
     - 对每个峰在指定窗口内做本底扣除 + 高斯+多项式拟合。
     - 用三个峰的峰位拟合得到平面刻度 `(K_plane, B_plane)`。

5. 诊断峰 RMS 与 FWHM

   - 使用 `DIAG_PEAKS`（如 6113 / 7686 / 5307 keV）：
     - 计算刻度后峰位与理论值的差值，得到各平面的峰位 RMS。
     - 从拟合结果中记录每个平面的最大 FWHM（作为整体分辨率指标）。

6. 逐条 FWHM 计算

   - 为每条 strip（尤其是 X/Y 条）在刻度后能量上建立条谱。
   - 在以诊断峰为中心的窗口内进行本底扣除与高斯拟合，得到条级 FWHM。
   - 将条谱和 FWHM 分布写入 `CALIB_DIAG_ROOT` 中。

7. 输出最终刻度表

   - 将平面刻度 `(K_plane, B_plane)` 与归一化 `(k_norm, b_norm)` 合并：
     - `k_total = K_plane * k_norm`
     - `b_total = K_plane * b_norm + B_plane`
   - 把 `(b_total, k_total, FWHM, PlaneRMS)` 等写入 `OUTPUT_DAT_FILE`（如 `ener_cal_Recal.dat`）。

### 8.2 通道与 strip 映射规则

- X 面：
  - `DSSDX_Ch` 在 [0, 127] → strip 0–127
  - 能量来源：`DSSDX_E[0]`

- Y 面：
  - `DSSDY_Ch` 在 [0, 47] → strip 128–175
  - 能量来源：`DSSDY_E[0]`

- YH 面：
  - `DSSDYH_Ch` 在 [0, 47] → strip 176–223
  - 能量来源：`DSSDYH_E[0]`
  - 通常还要求 `DSSDYH_mul == 1` 等条件

### 8.3 峰搜索与本底扣除

- 峰搜索：
  - 使用 `TSpectrum::Search` 在给定能量窗口内寻找峰位。
- 拟合模型：
  - 使用形如 `gaus(0) + pol1(3)` 的函数，对包含峰与本底的区间进行拟合。
- 本底扣除：
  - 若 `ENABLE_BKG_SUB == true`，则使用 `TSpectrum::Background` 对直方图进行本底估计，并从总谱中减去本底。

### 8.4 输出 `ener_cal_Recal.dat` 格式示意

    # Ch     b_total       k_total       FWHM       PlaneRMS
    0        b0            k0            fwhm0      rms_X
    1        b1            k1            fwhm1      rms_X
    ...
    127      b127          k127          fwhm127    rms_X
    128      b128          k128          fwhm128    rms_Y
    ...
    175      b175          k175          fwhm175    rms_Y
    176      b176          k176          fwhm176    rms_YH
    ...
    223      b223          k223          fwhm223    rms_YH

说明：

- `b_total`, `k_total`：该条最终能量刻度系数，可用于后续分析中将 ADC 转换为 keV：
  - `E_cal = k_total * ADC + b_total`
- `FWHM`：该条诊断峰（通常为 6113 keV）拟合得到的 FWHM。
- `PlaneRMS`：该条所属平面诊断峰峰位偏差的 RMS，平面内各条相同。

### 8.5 运行方式

    make Cali

运行结束后生成：

- `Diagnose_Calibration.root`：
  - 平面总谱（刻度前/刻度后）
  - 平面刻度拟合图（TGraph）
  - 各条刻度后能谱
  - FWHM 分布直方图等
- `ener_cal_Recal.dat`：
  - 最终能量刻度系数与逐条 FWHM，用于后续物理分析。

---

## 9. 常见问题与排查建议

1. Step 0 报错：找不到输入文件

   - 检查 `INPUT_DIR_PATTERN` 是否与实际文件名模式匹配。
   - 检查 `RUN_START` / `RUN_END` 是否包含正确的 run 编号。
   - 确认文件路径是否存在，权限是否正确。

2. Step 1 运行后坏条过多

   - 适当降低 `NORM_MIN_ENTRIES`。
   - 放宽 `NORM_E_DIFF` 或扩大 `NORM_MIN_E–NORM_MAX_E` 范围。
   - 检查是否存在物理上统计量极低的条（例如边缘条）。

3. Y 面归一化或桥接不稳定

   - 检查 Y 面参考条统计量是否足够。
   - 尝试在 `Config.h` 中显式设置 `REF_STRIP_X` / `REF_STRIP_Y` 为经验上较好的条。
   - 检查桥接拟合曲线是否近似线性。

4. Step 2 中部分条 FWHM = 0 或缺失

   - 说明该条诊断峰统计量不足或拟合失败。
   - 可在后处理时对这些条进行单独检查、屏蔽或使用邻道插值。

5. 平面刻度后参考峰偏移较大

   - 检查 `REF_PEAKS` 的能量是否与实际标准源一致。
   - 检查拟合窗口设置是否合适（窗口过窄可能截断峰，过宽则可能引入多余背景或其他峰）。
   - 若本底结构复杂，可考虑调节 `BKG_ITERATIONS` 或关闭本底扣除测试差异。

---

## 10. 建议的使用节奏与检查清单

1. 首次配置

   - 修改 `Config.h` 中：
     - 文件路径、run 范围、TTree 名称；
     - 参考峰与诊断峰能量；
     - 归一化能量窗口与选取条件等。
   - 使用 `make` 或命令行完成编译。

2. Step 0：Preselect

   - 执行 `./Preselect`。
   - 检查 `IMPLANT_ONLY_FILE` 与 `DECAY_ONLY_FILE` 的大小是否合理。

3. Step 1：Normalize

   - 执行 `./Normalize`。
   - 检查 `NORM_PARAM_FILE` 中的坏条数量。
   - 用 ROOT 打开 `NORM_DIAG_ROOT` 检查：
     - X/Y hits 分布是否正常；
     - 每条的 profile 拟合是否合理；
     - 参考条桥接关系是否线性。

4. Step 2：Calib

   - 执行 `./DSSD_Calib`。
   - 用 ROOT 打开 `CALIB_DIAG_ROOT` 检查：
     - 平面总谱的刻度是否正确（参考峰位置是否在指定能量附近）；
     - FWHM 分布是否在预期范围内（例如大部分条 FWHM 相近，没有异常坏条）。
   - 使用 `ener_cal_Recal.dat` 在后续分析中统一进行能量刻度。

---

## 11. 后续扩展建议


- 若有多个不同实验周期或 run 集合，可以为每一套单独生成刻度文件，再在分析阶段有选择地使用。
- 对 `Diagnose_Normalize.root` 与 `Diagnose_Calibration.root` 做进一步自动化分析，例如：
  - 生成 FWHM 热力图；
  - 自动列出 FWHM 超出阈值的“坏条”；
  - 自动比较不同 run 的刻度稳定性。



