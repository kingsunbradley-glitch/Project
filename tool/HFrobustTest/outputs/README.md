# α 逐跃迁 HF 数据库：使用说明与研究结论

本交付是可离线使用、可复算的研究数据库。它完整枚举 **ENSDF 2026-10-01 快照内本解析器支持的 Z=89–118 α 衰变记录**，补充 2025 年新报告的 257Sg，提供七套本地公式、九标签和逐条来源。**它尚不是截至 2026-10-22 所有实验论文和全部独立谱学证据的穷尽综述。**

当前有 **1,687 条 α 记录、312 种母核、345 个衰变数据集/实验条目**；其中 **1,494 条具备数值 HF**，形成 **15,655 条公式/角动量/候选结果**。其中包含暂定值和上下限，不能把这个数量等同于精确已测、已分辨跃迁数。

## 直接使用

双击 `index.html`，无需联网、服务器、数据库软件或外部字体。四个视图分别为核素/逐跃迁表、HF—标签范围、角动量对照、方法与覆盖。默认显示独立证据主层；切换“评估指认”或“全部记录”可查看更广的敏感性数据。

可以按公式、角动量、证据、标签、核区、母核 Z 与 N 范围、Z/N 奇偶、初末态类别、数据质量及文字筛选；同一转动带（same_band）与 Nilsson 轨道改变（nilsson_changed）各有独立的最低允许 α 角动量 lₘᵢₙ 筛选，lₘᵢₙ 由 Jπ 的三角条件及宇称选择定则得到。每个控件仅筛选对应标签行；选择该结构标签后，核素图、跃迁表和角动量对照也使用同一子集。多候选及未知 l 不进入指定 l 子集；筛选后重新计算中位数与分位数，不沿用全区 bootstrap。当前筛选值随 CSV 与范围图 SVG/PNG 导出；九标签图直接标注 HF 中位数；点击范围行追溯样本，点击跃迁查看原始量、初末组态和各公式结果。页面可导出筛选 CSV、范围 CSV、SVG、PNG；完整 SQLite 和 JSON 下载链接要求保留同目录文件。

`HF_RANGE_SUMMARY.md` 是七式、各自 native/l0/lmin、九标签、三证据层的完整数值矩阵；`hf_ranges.csv` 包括 270 组统计、样本ID和全区链bootstrap区间。

## 239Pu 逐线核查

见 `PU239_AUDIT.md`。70.77% → 0.7077，γ 归一化 NR/NT=10⁻⁶ 未乘入 α 强度。0.0765 keV 衰变末态与 0.0760 keV Adopted 异能态按测量误差唯一关联后恢复双带头 same_intrinsic 标签；保留两个原始能量。该线七式 HF 为0.569–2.546，计算值未变，结构证据仍在评估层。逐线详情现在同时显示 IA 原始百分数、γ 归一化与 Adopted 能量。

## 目前能成立的结论

独立证据人工复核覆盖 **17 个核态**。在同时满足逐跃迁输入条件后，主层目前只有 **2 条不同跃迁**：229Th→225Ra 基态和241Am→237Np 基态。233U→229Th 的82.7%线经注释核查为基态与低能异能态的合峰，已退出逐线标定。它们提供具体例证，但没有任何类别达到10个不同母核，**不能把主层跨度当作稳定的 HF 结构判别阈值**。

评估/敏感性层有更丰富的同轨道、轨道改变、宇称改变和 Λ−K 改变分布，数值范围存在重叠；其中部分指认依赖原研究中的 α 阻碍或系统学。它们可用于探索，不能反过来声称独立验证了 HF—结构关系。

`shape_changed` 与 `multi_qp_changed` 暂无经本轮复核建立的可标定样本。251Fm 的电磁强度支持粒子—八极声子混合解释，但相应暂定核态等条件未满足主层。空范围的原因见 `missing_range_reasons.csv`；空白不表示物理上不存在该机制，也不表示所有相关文献已查完。

## 输入、公式与角动量

- 总半衰期单位为秒；总 α 分支 `b_alpha` 与单条 α 线内部强度 `I_alpha` 均为0–1比例。`I_alpha` 已包括有明确给出的α强度归一化因子NB；`br=b_alpha*I_alpha`，`T_partial_s=T_total_s/br`。缺失值不补1。
- 旧 HFdistr 的 `T_alpha_half_s` 已是 α 部分半衰期，不再除总α分支。逐线分支仍须单独分辨。
- 主计算采用每条 `Ealpha` 的本地反冲约定 `Q_if=Ealpha*A/(A-4)`；缺少Eα才使用 `Q_gs+Ei−Ef`。`ame_q_crosscheck.csv` 独立保存AME2020质量差比较，未静默混用两种Q。
- 初态与末态匹配采用评估能级能量和偏移符号；已采用较新Adopted半衰期的记录保留旧Decay P值和更新审计。候选Jπ、括号、上下限、未放置α线和未知激发能均保留。
- `l_min` 是满足三角条件、宇称和整数角动量耦合的最小l。0+→0−没有允许解；半整数J与整数J不允许通过自旋0的α相连。多候选分别计算，不当作重复实验样本。
- DZR、Ismail2022 Formula E、Xu2022分别计算l=0和l_min；l=0禁戒时标为形式参考。Royer、Qi2009 UDL、Viola–Seaborg、E. Rurarz没有显式l项，仅输出native，l_min不适用。
- `PoenaruMethod.py`未纳入。E. Rurarz虽含Poenaru型文件名，按要求保留；Royer取自本地DZR脚本的独立Royer函数。原六个脚本逐字复制，SHA256可核查。
- 没有发现充分原文依据支持系数修改，因此没有人为建立“核正版”。Xu区域外C=0保留为本地约定；VS采用偶偶核拟合基准，奇核为外推；Rurarz的A/Ad记号仍待原始扫描进一步核验。详见`formula_audit.csv`。

## 标签与统计

九标签为 same_intrinsic、same_band、nilsson_changed、multi_qp_changed、shape_changed、configuration_mixing、parity_change、spin_flip、unknown。初末轨道标注来自评估带头或明确轨道指认；人工复核的独立证据记录在可复现包的`data/state_evidence_overrides.json`。

**2026-10-25 标签口径修订：** `same_intrinsic` 仅指相同内禀准粒子组态的 **带头 → 带头**，初态与末态必须均有带头指认；任何一端是非带头转动成员，都不归入此类。`same_band` 用于同组态对应转动带、且至少一端为非带头成员的跃迁，这两个标签互斥，其他结构标签仍可叠加。

自动带角色判定要求核态有明确的评估带关联及唯一 Jπ：在名义单准粒子带中，J=K 标为带头，J>K 标为非带头转动成员，且宇称需一致。K=1/2 带的 J>K 基态不因能量最低自动视为 J=K 带头。没有带关联、多候选或指认矛盾时标记“未确定”，不凭相同轨道、基态、相同 Jπ、l=0 或 HF 大小补成带头。暂定指认保留在敏感性层。带头身份不等于基态，激发带头可以纳入。

逐跃迁保存 `initial_band_role`、`final_band_role` 及判定依据、证据等级；核态详情同时给出评估原文和来源。`bandhead_classification_audit.csv` 与 SQLite 同名表可追溯相同轨道候选的分类与排除理由。

**2026-10-26 口径修订：** `spin_flip` 仅限 **带头 → 带头**，且有可识别单准粒子 Nilsson 轨道、初末 Λ−K 不同；任一端为非带头成员或带头身份未确定时均不纳入。`parity_change` 只要求初末态宇称改变，不受带头限制，可包含涉及转动带成员的 α 跃迁。候选保留与排除依据见 `spin_flip_classification_audit.csv`。K为轨道投影，不是转动带成员J；其符号与常用Σ=K−Λ相反。多个可能轨道不强制唯一选取。独立意味着不依赖α-HF，仍可能包含磁矩、转移反应与核模型的解释，并非模型无关测量。

ENSDF A记录第80列S（预期/推断α支路）和已核实未分辨合峰，保留在全部记录表并展示参考计算，但不进入任何逐跃迁HF范围。明确以Q和能级推算的Eα保留derived_alpha_energy来源标记；推算本身不自动排除。能量关系显著不一致者标记energy_balance_inconsistent，退出主/评估层，留在全部记录层供核实。

每组统计按跃迁ID去重；除 same_intrinsic 与 same_band 互斥外，其他标签之间允许重叠。HF空间线性插值求P10、P50、P90；图使用log10HF坐标。小于10个不同母核仅作为小样本跨度；unknown只显示观测分布。上下限和重复备选在全部记录表保留，但不进入点样本范围；l_min多候选也不入范围。

全区范围使用固定种子20261022，按`(N−Z, Z mod 2)`组成α链族，分组重采样1,000次，至少10母核且3链时报告95%区间。这是样本聚类不确定度；实验测量误差按一阶非对称传播单独保存，不含公式系统误差，也未传播未知协方差。筛选子集重新计算点分位数，**不会挪用全区bootstrap区间；子集区间留空**。

HF是实验/指定公式的半衰期比，不等同于唯一的初末波函数重叠或某个形变量。不同公式的归一化与屏障项会改变范围，不能把各式或l0/l_min混合使用。

## 覆盖审计与尚未解决的部分

NUBASE2020中的386个Z89–118含α分支核态里，303个核态能量匹配本库、48个有母核但核态未唯一对应、35个没有找到该ENSDF快照的α衰变数据集。该清单包括估计分支和上限，见`nubase_alpha_inventory.csv`。这些记录没有被伪装成已知末态的α跃迁。

旧库210条核素/状态记录已交叉核查；176条旧计算结果与六套共有公式均通过数值回归（最大log10半衰期差约1.1×10^-11，来自CSV舍入）。Royer另外作独立解析式回归。旧库的核素级记录与本库多条谱线不必一一对应，差异见`legacy_crosscheck.csv`。

2025年257Sg只报告一条α事件，b≈0.1无正式分支误差，末态能量未知；其暂定9/2−[734]直接α喂入态不能与随后裂变的7/2+[624]态混同。253No→249Fm基态带总馈入0.7%未被拆给三个未知强度分支。依赖HF的247Md/257Rf组态指认留在敏感性层。

255No的原始IA以100为最强线相对强度；按ENSDF N记录NB=0.297归一化，十条线强度合计0.995247，未强行补齐到1。对应HF较忽略NB的计算增大1/0.297倍。[NNDC/ICTP归一化讲义第6页](https://indico.ictp.it/event/7641/session/33/contribution/73/material/1/0.pdf)明确给出α强度的NB×BR转换。

原始实验专题复核共11条来源（其中1条仅为理论背景索引），并非所有原论文逐表已读。还需扩大逐核态独立证据复核、高K多准粒子支路及独立形变证据，才能回答全区稳定的经验阈值问题。

## 文件与复算

- `hf_database.sqlite`：transitions、states、sources、results、label_evidence、hf_ranges、methods、coverage、audit、bandhead_classification_audit、spin_flip_classification_audit、formula_audit、literature_evidence、AME/NUBASE/旧库交叉核查及metadata。
- `database.json`与`index.html`：完全对应的内嵌可视化数据。
- 同名CSV：平面导出；JSON数组/嵌套量以JSON文本存储。
- `sources/`：实际使用的ENSDF原始Adopted/Decay文本；`provenance_manifest.json`记录原始缓存和公式SHA256。
- `figures/`：六组范围图，每组SVG和PNG；native/l0与native/l_min、三证据层分别输出。
- `HFrobustTest_reproducible.zip`：源代码、原样公式、冻结原始数据及配置。解压后安装依赖，运行：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python src/run.py
node src/test_viewer.cjs
```

安装依赖后，复算不需要联网。`fetch_data.py`用于另行更新数据缓存；当前`build.py`明确固定`ensdf_261001.zip`，升级快照必须人工修改并重新核查差异，不能悄然覆盖已发布结论。原始计算目录和已有HFdistr数据库保持不变。

## 主要来源

- [ENSDF 数据归档](https://www.nndc.bnl.gov/ensdfarchivals/)及[本次2026-10-01快照](https://www.nndc.bnl.gov/ensdfarchivals/distributions/dist26/ensdf_261001.zip)，评估由NSDD网络提供；NNDC的HF字段仅为文献参考，不用于拟合本项目标签。
- [AME2020：M. Wang等，Chinese Phys. C45,030003](https://doi.org/10.1088/1674-1137/abddaf)；[AMDC官方数据页](https://amdc.impcas.ac.cn/web/masseval.html)。
- [NUBASE2020：F. G. Kondev等，Chinese Phys. C45,030001](https://doi.org/10.1088/1674-1137/abddae)；[官方数据页](https://amdc.impcas.ac.cn/web/nubase_en.html)。
- [257Sg新实验，PRL134,232501(2025)](https://doi.org/10.1103/s7hr-y7zq)。其余实验定位、证据边界及七公式DOI见`literature_evidence.csv`、`formula_audit.csv`和页面“方法与覆盖”。

## 三式共同越界核查（2026-10-26）

`OUTLIER_AUDIT.html` 为可搜索的离线核查表；`OUTLIER_AUDIT.md` 给出结论，`outlier_audit.csv` 给出46条候选的原始输入、各式HF、检查项和证据，`outlier_label_memberships.csv` 保存56个标签成员项的修订前P10–P90阈值。选择固定为评估指认层、全核区、无l筛选、VS/Rurarz/Royer三式共同越界。确认4条候选需改变统计资格：233U合峰1条、239Pu的S支路2条、242Cm能量不一致1条。全库继承同类来源规则，94条S记录和14条242Cm推算能量记录获得相应来源标记；后者仅303.42与763.22-keV两条因能量关系不一致退出主/评估层。所有15,655条公式结果的数值保持不变。

核查覆盖冻结ENSDF转录及本地算术，不等于46条原始实验的重新测量或全部论文的独立复核。Rurarz原论文A/Ad记号仍待核准；VS奇核外推与缺少显式l的影响继续保留。
