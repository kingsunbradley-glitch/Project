# waveFit 效果检测记录

检测日期：2026-09-09。输入为 `test/wave_242Fm.root` 和 `test/wave_246Fm.root`。本目录 `test` 中提供的是 ROOT 数据，没有额外的文字验收说明，因此使用需求 Markdown 定义的读取、导出和分析行为进行检测。

## 编译与自动检查

`make` 编译成功，使用当前环境 ROOT 6.36.06 的 C++20 参数，无编译警告。自动检查 `python3 tests/validate.py` 全部通过。提供了 CMake 配置，但本环境没有 `cmake` 命令，未验证 CMake 构建。

最新完整记录：[validation.json](output/tests/validation_wq4ichzv/validation.json)。该目录同时保存每条测试命令的 stdout/stderr。

| 测试文件 | Canvas 数 | 真实波形数 | 逐点核对数量 | 结果 |
|---|---:|---:|---:|---|
| wave_242Fm.root | 22 | 83 | 226000 | 全部相等 |
| wave_246Fm.root | 8 | 26 | 70000 | 全部相等 |

独立 PyROOT 读取原始 TGraph，逐点比较导出的 CSV、分析 TTree 的 x/y 和 `raw_waveform`。三个输出与输入的双精度值全部相等；读取前后 SHA-256 相同。Canvas 事件信息和分析条目数正确。坐标框 hframe 没有被当成波形。

补充结构与命令测试通过：深层 TPad、子目录、TGraphErrors、TGraphAsymmErrors、TH1、TH2 排除、TLine 展示、空 Canvas、重复 key cycle、同名 Graph 的独立文件、缺失事件元数据、空格分隔的 run/map_entry、64 位 map_entry、所有波形重绘，以及已有输出、错误正则、错误参数、缺失文件／chain／模板的报错。

## 已知真值的数值检测

采用独立构造的高斯单脉冲模板，叠加确定随机种子的 0.5 ADC 噪声。首沿真值 t1=80.25 samples，双脉冲 Δt=35.4 samples，幅度 1000 与 650 ADC。通过固定首沿的一维扫描、QR 和 Brent 进行恢复。

| 项目 | 真值 | 程序结果 |
|---|---:|---:|
| 正双脉冲 Δt / samples | 35.4 | 35.39995776 |
| 负双脉冲 Δt / samples | 35.4 | 35.40240555 |
| 正双脉冲 A1 / ADC | 1000 | 1000.11927 |
| 正双脉冲 A2 / ADC | 650 | 649.85837 |
| 从含噪单脉冲生成模板后自动检测 Δt / samples | 35.4 | 35.39279877 |

单脉冲对照没有被标为强 pile-up 候选。负脉冲极性正确恢复。Δt=0 扫描点有限，没有重复模板列导致的奇异结果。检查了重建曲线与残差点数、匹配滤波功率与 χ² 扫描的一致性、模板峰值归一化、时间单位转换与不一致单位的拒绝，以及 `--t1` 放在 `--fit` 前后的参数行为。

这些是指定模板、噪声与首沿条件下的算法检验，不是实际探测系统的分辨极限或误差标定。

## 真实波形的模板拟合演示

使用 `wave_246Fm.root` 中 8 个 Canvas 的 YHWave、前 100 点基线构建模板，共同覆盖区间为 2498 点；拟合 `wave_242Fm.root` 的 `ChainNum_1/YHWave_1`。

- 两个导数候选的时间差：约 636.04 samples。
- 模板拟合 Δt：635.90217 samples。
- 双脉冲 χ²：62551.71，NDF=2494，χ²/NDF≈25.08。
- 单脉冲明显不能描述两次阶跃，但双脉冲残差仍远大于基线噪声。

因此该实例证明读取→模板→拟合→输出链条可以运行，不能将 ΔBIC 或拟合间隔直接当成最终物理结果。下一步需检查同类通道模板一致性、首沿偏差、信号相关噪声与可能的饱和，再评估参数误差。

结果文件：

- [拟合 CSV](output/real_fit_final/wave_242Fm_analysis.csv)
- [拟合 ROOT](output/real_fit_final/wave_242Fm_analysis.root)
- [波形、单／双脉冲模型及残差图](output/real_fit_final/fit_validation.png)
- [重绘的原始 YHWave](output/real_fit_final/drawings/ChainNum_1_YHWave_1_w2.png)

尚未校准 sample→ns、ADC→能量，也没有给出 σE、σΔt 或探测器实际 pile-up 分辨极限。
