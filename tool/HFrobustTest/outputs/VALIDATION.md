# 验证记录

验证完成日期：2026-10-26（研究数据截止日期仍为2026-10-22）。

- 15项物理/存储测试通过（新增239Pu归一化/核态回归与近基态异能态匹配歧义测试）：spin_flip 双带头限制与 parity_change 转动成员保留、带头/转动成员分类、两类标签互斥与来源对应、单位、非对称误差、候选Jπ、无允许l、Nilsson K与J区别、部分半衰期归一化、255No显式NB、原式数值回归、HF角动量恒等式、样本去重和SQLite完整性。
- 270组冻结范围与页面动态范围的成员、计数、状态、HF分位数完全一致。
- 页面逻辑验证覆盖84个渲染组合；母核 N 上下限及 same_band 最低允许 l 筛选、重置测试通过。九标签图（含并排图）的 HF 中位数标注与筛选数据一致；same_band 的 l 子集不复用全区 bootstrap，其他标签保持原范围，候选不重复入样，文本转义与CSV转义检查通过。
- 2026-10-23 的本机浏览器HTTP预览验证：初始独立层3条跃迁、DZR范围/角动量切换、全部记录层、质量筛选得到新257Sg单事件、完整详情及来源证据。
- 浏览器实际导出PNG成功，并用Pillow验证图片格式；保留`figures/browser_angular_export.png`作为导出验证样品（导出时点的局部视图，不替代最终冻结全区图）。
- 12幅最终静态范围图由Matplotlib生成SVG/PNG；核查图形排版及数据对应。
- 页面没有外部脚本、样式、字体或运行时fetch，数据内嵌。自动化浏览器策略仅允许HTTP/HTTPS，拒绝file://协议，因此没有通过该自动化工具直接测试双击本地文件；没有绕过限制。双击离线入口仍为标准自包含HTML，附带文件仅用于下载链接。
- 2026-10-25 按用户口径重建 same_intrinsic：相同组态且初末态均为带头；same_band 至少一端非带头，两类互斥。全部 15,655 条公式计算结果与修订前逐项一致；仅结构标签和相关统计变化。重建 270 组范围及12幅静态图，检查图中带头分类说明和中位数标注。
- 2026-10-26 修订：spin_flip 仅限双带头且 Λ−K 改变，149 条候选中保留34条；parity_change 的388条记录、各组成员与点统计保持一致。全部15,655条计算结果保持一致。所有组按固定种子重新执行bootstrap；抽样区间可能随前序组的随机数消耗变化，不将其视为宇称数据变化。
- 原始formula目录及HFdistr只读；交付另存。

详见`validation_tests.txt`、`formula_regression.csv`、`branch_normalization_audit.csv`、`ame_q_crosscheck.csv`。科研内容的覆盖限制与未完成证据复核见README，测试通过不表示研究材料已穷尽。

- 2026-10-26 239Pu核查：仅 ensdf-412f45083d69-a1 的末态轨道、带头与结构标签改变（unknown → same_intrinsic）；所有15,655条结果逐字段一致。NR/NT列显式保存且不进入α计算，IA=0.7077；0.0765与0.0760 keV按3σ合成误差唯一匹配，不合并为0-keV基态。same_intrinsic由27条变为28条，其他已赋结构标签成员保持不变。页面详情专项检查覆盖两份能量、70.77原始百分数和γ归一化说明。全区范围及静态图已重建，新的bootstrap区间随重新采样更新。

- 2026-10-26 nilsson_changed 角动量筛选：验证非零 l=2 子集的成员、HF中位数、表格、显式 l 修正模式与CSV/SVG筛选信息；同same_band控件独立生效。未知和多候选不进入指定l子集，子集不复用全区bootstrap，重置清空两个控件。270组默认范围和84个渲染组合继续通过；本次仅更新页面筛选，数据库数值及冻结统计不变。

- 2026-10-26 三式共同越界审计：固定修订前评估层56个标签成员项/46条跃迁；原始A/L/P/N、单位、分支、核态来源及独立重写公式共598项检查通过。15,655条结果的log10T、T、log10HF、HF及误差、角动量修正数值逐项不变。修正S标记/合峰/能量不一致的统计资格，独立主层3→2条。94条S记录退出单线范围，233U合峰退出单线范围；242Cm的14条推算能量保留来源标记，仅303.42和763.22-keV两条因能量关系不一致退出主/评估层；其他自洽推算能量不被自动排除。270组冻结/动态范围一致，84个渲染组合通过，离线审计页筛选另行验证。

2026-10-26 final delivery checks: all 15,655 formula results unchanged except sample_status; SQLite integrity_check=ok; outlier_audit=46 rows. Real IAB rendering verified: full table 46, search 239Pu 3, combined review status 2, cleared search review status 4. Screenshot figures/outlier_audit_preview.png. Temporary local server and tab closed.
