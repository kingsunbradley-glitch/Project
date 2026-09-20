# ROOT Canvas 波形读取与 pile-up 分析

按 `waveform_reader_codex_spec.md` 实现的独立 ROOT/C++ 程序。直接读取 Canvas 中的数值对象，原始 ROOT 文件始终以 `READ` 打开。支持递归目录和 pad、TGraph 及派生类、1D TH1、波形导出与重绘、基线和脉冲候选、模板以及单／双脉冲比较。

## 编译

在本目录运行：

```bash
make
./wave_reader --help
```

使用当前环境的 `root-config` 获取编译参数，包括 ROOT 要求的 C++ 标准。代码需要 C++17 或更新版本；本机 ROOT 6.36.06 使用 C++20。

也支持：

```bash
cmake -S . -B build
cmake --build build -j2
./build/wave_reader --help
```

## 对 test 数据进行检查和导出

```bash
./wave_reader test/wave_242Fm.root --inspect
./wave_reader test/wave_246Fm.root --list
./wave_reader test/wave_242Fm.root --export-waveforms --output-dir output/242Fm
./wave_reader test/wave_242Fm.root --draw 1 --output-dir output/242Fm
./wave_reader test/wave_242Fm.root --analyze-all --baseline-samples 100 --output-dir output/242Fm
```

`--draw N` 重建所选 ChainNum 的全部候选波形，每条独立输出 PNG、PDF 和 ROOT；`--draw-all` 处理全部 Canvas。程序默认 batch 模式，不依赖图形桌面。输出 ROOT 中的 `extracted_waveform` 是从数值重新构造的 TGraph。

`test/wave_242Fm.root` 有 22 个 Canvas、83 条波形；`test/wave_246Fm.root` 有 8 个 Canvas、26 条波形。真实对象包含 XWave、YWave、YHWave、SSDWave，有的 pad 无数据，有的 pad 有多条波形。`hframe` 是 ROOT 生成的空坐标框，默认在 inspect 中显示但不提取；`--include-frames` 可包含它。其他 1D TH1 保留为候选，TH2/TH3 仅显示结构。每个 key 只读取最新 cycle。

## 对象筛选

```bash
./wave_reader test/wave_242Fm.root --analyze 1 --object-pattern '^YHWave_' --baseline-samples 100
./wave_reader test/wave_242Fm.root --draw 1 --graph-index 0
```

`--object-pattern` 是区分大小写的正则表达式，匹配 object name、title、pad name 或 pad 完整路径。`--graph-index` 是每个 Canvas 内递归遇到的 TGraph 的零起始索引，不含 TH1；两个筛选器同时指定时取交集。命令输出会列出所有候选及是否选中。`--analyze N` 和 `--analyze-all` 都逐条处理选中的对象，不会合并多个探测器的数据。

## 基线与时间单位

默认基线窗口为 `min(200,max(2,N/10))`，使用均值和样本标准差；`--robust-baseline` 使用 median 和 `1.4826*MAD`。实际测试数据中不少脉冲约在 sample 180–190 起跳，默认 200 点可能包含信号，建议先查看波形并设置 `--baseline-samples 100`。检测到起跳与基线窗口重叠时会输出警告。

通过正、负极值自动判断极性。`baseline_subtracted` 保存的是减基线且已统一正极性的波形，原始 y 保存在 `raw_waveform` 和树的 y branch。`--smooth 5` 等奇数移动平均窗口仅影响候选检测，不用于拟合原始数据。阈值默认 5 倍基线段导数噪声，`--min-separation 20` 为候选间最小样本距离。连续超阈值的上升沿合为一个候选；候选 x 对导数峰作抛物线插值。候选 height 是直到下一个候选之前的最大基线修正高度，不等同于分解后的第二脉冲幅度。

默认保持原始 x，单位标签为 `original_x`；本次两份测试数据的 x 是 sample index。知道 ADC 采样周期后才能换算物理时间，例如 **仅在实际采样周期为 10 ns 时**：

```bash
./wave_reader test/wave_242Fm.root --analyze 1 --sample-period-ns 10 --baseline-samples 100
```

转换定义为 `x_ns = 原始 x × 采样周期`，不会重新编号。原始 x 已经是时间时，应使用 `--x-unit ns` 或 `--x-unit us`，不进行第二次换算。`--t1` 和 `--max-delta` 都使用当前 x 单位；`--max-delta` 默认 2000，程序不猜测 ADC 周期。模板也保存并检查单位标签；没有单位元数据的外部模板会发出警告并按当前单位解释。

## 构造模板与拟合

下面示例选择同一类 YHWave；使用前仍需检查各通道的形状是否适合共享模板。

```bash
./wave_reader test/wave_246Fm.root \
  --object-pattern '^YHWave_' --baseline-samples 100 \
  --make-template 1 2 3 4 5 6 7 8 --output-dir output/template_YH

./wave_reader test/wave_242Fm.root \
  --object-pattern '^YHWave_' --baseline-samples 100 \
  --template output/template_YH/pulse_template.root \
  --fit 1 --max-delta 1500 --output-dir output/fit_242Fm
```

`--make-template` 不带编号时，在选中对象中自动寻找只有一个候选的波形。存在多对象 Canvas 时必须指定筛选器，避免不经选择就平均不同探测器。脉冲候选数为 1 只是初筛，不保证信号没有饱和或未分辨的 pile-up。

模板先减基线、统一极性、按导数上升沿位置对齐，逐条归一化幅度，再在共同覆盖的时间区间内插值平均，最后将模板峰值归一到 1。输出 `pulse_template.root` 中的 `TGraph pulse_template` 和同名 CSV；模板 t=0 对应对齐后的上升沿参考。输入模板源要求等间距、相同采样间隔。插值使用分段线性函数，左侧超出范围返回 0，右侧保持最后测得值，以兼容仍在平台上的电荷敏感前放波形。模板尾部不足时需要检查此假设。

`--analyze`/`--analyze-all` 在提供模板后也会拟合；不提供模板时只做基线和候选分析，拟合字段为 NaN，状态为 `not_requested`。`--fit` 必须提供模板，或与 `--make-template` 同时使用。

本版实现固定首脉冲参考的 **一维 Δt 扫描**：默认 t1 来自首个导数候选，也可用 `--t1 VALUE` 提供外部 trigger。扫描步长为一个样本的中位 x 间隔，在最佳点邻域用有界 Brent 细化。每个时间假设下用列缩放的 Householder QR 求解幅度与常数／线性基线；脉冲幅度约束为非负，负幅度解通过边界模型重新求解。Δt=0 单独退化为单脉冲列，避免奇异矩阵。扫描范围不会超过波形末端，最大 100000 步。

`matched_filter_power` 保存相对单脉冲模型的加权残差平方和下降量：对第二模板消去首脉冲与基线分量后，它对应广义匹配滤波的功率响应，幅度边界处按非负约束处理。它与 `chi2_vs_deltaT` 使用同一组扫描时间。

权重使用原始基线噪声；零噪声时明确警告，并以 1 ADC 为诊断权重。计算 χ²、NDF、AIC=χ²+2k、BIC=χ²+k ln N；外部固定 t1 时单／双脉冲分别计 3／5 个参数，数据估计 t1 时分别计 4／6 个参数。ΔBIC=BIC1−BIC2，正值倾向双脉冲。超过 10 且两幅度非零时标为 `pileup_candidate`，扫描最优值落在右边界时标为 `scan_boundary`。

这是对 **固定 t1、指定模板和噪声模型** 的条件比较，尚未实现 t1、t2 联合优化。模板不匹配、ADC 饱和、基线污染、相关噪声或 t1 偏差可能产生很大的 ΔBIC；程序同时输出 reduced χ²，较大时提示检查残差。模型边界附近的 AIC/BIC 是近似判据，不能直接当成发现显著性或最终寿命测量。

## 输出

默认目录 `output/<输入文件名去掉扩展名>/`，可用 `--output-dir` 更改。写入前检查目标文件，已有同名输出需显式 `--overwrite`；即使指定此参数也不允许覆盖输入 ROOT 或输入模板。

- `waveforms/<canvas路径>_<objectName>_w<对象索引>.csv`：`sample,x,y`，17 位精度。
- `drawings/`：重绘 PNG、PDF、ROOT。
- `<输入名>_analysis.csv`：事件、对象、基线、极性、候选数和拟合摘要。
- `<输入名>_analysis.root`：`WaveAnalysis` TTree，包含事件与完整对象路径、原始 x/y、候选样本／时间／高度向量、拟合状态与模型比较。
- 每个波形独立的 ROOT 子目录：`raw_waveform`、`baseline_subtracted`、`best_single_fit`、`best_double_fit`、两种 residual、`chi2_vs_deltaT`、`matched_filter_power`。尚未拟合时不生成拟合图。

对象索引和 Canvas 路径用于区分同名 Graph；若清理文件名字符后仍发生冲突则报错，不覆盖。个别非法波形的分析失败会记录 `analysis_error`，保留其他波形输出并返回退出码 2；命令／输入错误返回 1，成功返回 0。非有限点或非递增 x 可原样导出，但不能分析。

## 验证与当前边界

```bash
make test
```

测试需要当前环境可导入 PyROOT。独立测试脚本递归读取两份原始数据，逐点比较 CSV、ROOT 图和树，并校验输入 SHA-256 未改变。还检查嵌套 pad／目录、TGraphErrors／TGraphAsymmErrors、TH1、同名对象、重复 key cycle、参数错误、模板生成、单位一致性，以及已知正负双脉冲的 Δt 和幅度恢复。每次输出保存在新的 `output/tests/validation_*/`，含命令日志与 `validation.json`。

本版已实现需求中的第一版读取标准，以及后续模板和一维 pile-up 分析。模板导数 h′ 的未分辨脉冲投影、能量标定、参数置信区间和分辨极限扫描属于后续物理分析工作，目前不输出未经标定的 E 或 σE、σΔt。实际验证结果见 [VALIDATION.md](VALIDATION.md)。
