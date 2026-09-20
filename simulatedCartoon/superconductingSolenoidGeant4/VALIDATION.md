# 验收记录

验收环境：`root_env`，Geant4 11.4.2，ROOT 6.36.06，Release 构建，串行固定种子。

| 检查 | 结果 |
|---|---|
| CMake 编译 | 通过 |
| CTest 单元测试 | 通过：质量表、E_lab/E* 互逆、A/Z 与四动量守恒、磁场中心值、轴对称和边缘连续性 |
| CTest 输运测试 | 通过：固定种子复现、无场 `rho_entry=d tan(theta)`、真空动能守恒、He 能损与 `0<=q<=Z` |
| 物理模式 smoke test | 10/10 事件完成，ROOT/CSV/JSON 输出有效 |
| 接受度视觉样本 | 1000 事件、20 角箱；200 transmitted、550 aperture、250 not_entered |
| ROOT 图 | PNG/PDF/ROOT 均成功生成并视觉检查 |
| Qt/OpenGL | 20 事件运行成功；当前 WSL 从 OGLSQt 回退 OGLSX，导出的轨迹 PDF 视觉检查通过 |
| `10^4` 集成运行 | 10000/10000 事件完成；ROOT、CSV 和终态计数均为 10000；344.97 s（约 29.0 event/s） |

默认 `10^5` 宏为 `macros/batch_Ar40_Tm169.mac`。按 `10^4` 的相同 2 mm 步长、He 与磁场配置线性估计约 57.5 分钟；本次未重复等待完整 `10^5` 串行运行。运行时设置 `trackSampleCount 0`，因此大样本不会保存逐步轨迹。
