# Geant4 超导螺线管蒸发余核输运

这是一个独立的 Geant4 11.4 工程，用于模拟轴上点靶产生的完全熔合蒸发余核，从靶点经过漂移距离 `d`、进入有限长超导螺线管并到达出口的过程。默认反应为

\[
^{40}\mathrm{Ar}+^{169}\mathrm{Tm}\rightarrow{}^{209}\mathrm{Fr}^{*}
\rightarrow{}^{205}\mathrm{Fr}+4n / {}^{204}\mathrm{Fr}+5n.
\]

第一版只接受完全熔合后的 `xn` 道；没有模拟真实靶膜厚度、靶内能损、`pxn` 或 `αxn`。入口半径是输运结果，不是发生器输入或预筛选条件。

## 环境与编译

项目使用现有 `root_env`，其中 ROOT 保持为 6.36.06，并增加 Geant4 11.4.2、Qt6/OpenGL 及官方核数据集。

```bash
conda activate root_env
cmake -S . -B build -G Ninja
cmake --build build -j
ctest --test-dir build --output-on-failure
```

若不使用 Ninja，去掉 `-G Ninja` 即可。CMake 会把 `macros/` 复制到构建目录。

## 运行

```bash
# Qt/OpenGL，默认显示 30 条代表轨迹
./build/solenoid_transport macros/vis.mac

# 运行 20 个事件并自动导出 output/open_gl_tracks_0000.pdf
./build/solenoid_transport macros/vis_export.mac

# 10^5 个物理分布事件，不保存逐步轨迹
./build/solenoid_transport macros/batch_Ar40_Tm169.mac

# 60 个角度箱，每箱 1000 个真实余核样本
./build/solenoid_transport macros/angle_scan.mac

# 用 ROOT 绘图，不依赖 Matplotlib
python scripts/analyze_results.py output/angle_scan.root
```

`vis.mac` 使用 Qt/OpenGL 驱动 `OGL`；`tools_sg_vis.mac` 使用 `TSG_QT`。后者要求当前 Geant4 构建实际提供 ToolsSG Qt 驱动。

## 宏参数

- 几何：`/sim/geometry/solenoidRadius`、`solenoidLength`、`targetDistance`、`maxStep`
- 磁场：`/sim/field/centerField`
- 气体：`/sim/gas/enabled`、`pressure`、`temperature`
- 反应：`/sim/reaction/projectileZA Z A`、`targetZA Z A`、`beamEnergyLab`、`excitationEnergy`
- 截面：`/sim/reaction/sigma4n`、`sigma5n`，单位为 μb
- 扫描：`/sim/scan/thetaMin`、`thetaMax`、`bins`、`eventsPerBin`
- 输出：`/sim/output/trackSampleCount`、`postExitTrackLength`、`fileName`、`randomSeed`

束流能量和 `E*` 可以任选其一。若在同一宏内显式设置两者，程序会用 Geant4 质量表检查一致性；反推出的 `E*` 相差超过 0.5 MeV 时直接报错。4n/5n 按相同事件数条件采样；两条截面都大于零时才使用绝对截面权重并报告 `σ_out`，否则按等权归一化只报告效率。

角扫描保留退激模型产生的余核种类、动能和电荷，只把方向旋转到指定角箱。几何参考角为 `atan(R/d)`，不作为硬截断。

## 几何、磁场与终态

靶点为 `(0,0,0)`，螺线管入口和出口分别为 `z=d` 与 `z=d+l`。磁场由有限长圆电流片的完整轴对称场数值积分生成，包含两端边缘场；中心归一为 `B0`，Geant4 内部以笛卡尔坐标用 `G4FieldBuilder` 和 Dormand–Prince 积分器输运。

唯一终态编码为：`1 transmitted`、`2 aperture`、`3 not_entered`、`4 stopped`、`5 backward`、`6 world_exit`。只有穿过出口平面且 `ρ<R` 的事件才是 transmitted。

He 打开时按宏给定的压力和温度由理想气体状态构造密度，参考物理表 `FTFP_BERT` 提供重离子电离、能损和多重散射；动态电荷取自每一步的 `G4DynamicParticle`。关闭气体后输运区为 `G4_Galactic`。

## 输出

每次运行直接生成：

- `*.root`：`events`、`tracks`、`summary` 三棵 ntuple；
- `*_events.csv`：事件级初态、入口/出口量和终态；
- `*_tracks.csv`：抽样轨迹的 `(ρ,φ,z)`、`(vρ,vφ,vz)`、能量、电荷、时间；
- `*_summary.csv`：角箱生成/进入/通过数、效率和 95% Wilson 区间；
- `*_summary.json`：几何角、最大观测通过角、`T≥1%` 外侧区间、通道效率和可选 `σ_out`。

`scripts/analyze_results.py` 额外输出 ROOT 绘制的角分布、接受度、`ρ-z` 轨迹和入口/出口分布，格式为 PNG、PDF 和 ROOT。
