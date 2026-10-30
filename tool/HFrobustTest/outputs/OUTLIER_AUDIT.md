# VS / E. Rurarz / Royer：共同越界数据核查

评估指认层、全核区、无l筛选；按修订前每个结构标签各自的P10–P90取三式共同越界的交集。冻结候选为46条跃迁、56个跃迁—标签项、23种母核。P10–P90是中央80%经验范围，不是中心值的置信区间；不以越界作为删除理由。

## 已确认的记录状态问题

- 233U→229Th：0-keV条目的82.7%包括基态与8.355734-eV异能态，不能作为单独基态分支。保留数值为合峰参考，退出逐跃迁范围。
- 239Pu→235U：338.52、821.25-keV条目为ENSDF S标记的预期/推断α支路，不能当作已直接观测、已分辨α线进入标定；IA不是γ归一化误乘。
- 242Cm→238Pu：763.22-keV线的Eα来自能量关系推算，当前Eα与P.Q−Ef相差约4.07 keV。标记energy_balance_inconsistent，保留原值于敏感性层；计算能量本身不自动判错。详表给出改用能级Q的影响，仅作敏感性比较。
- 同类规则覆盖全库：94条S记录获得明确状态；242Cm的14条推算能量记录继承来源标记，其中303.42和763.22-keV两条能量关系明显不一致，退出主/评估层；其余不因推算来源而被自动排除。数值未因越界而调整。独立主标定从3条变为2条（233U合峰退出）。

## 数值与解释

逐项核对原始A/L/P/N记录、百分数、总α分支、半衰期单位、末态自旋来源、反冲Q与部分半衰期；独立重写三式数值表达式复算。失败项：0。这仅证明与冻结源和本地约定一致，不保证实验值、组态纯度或原论文记号完全正确。

- 251Es基态为same_intrinsic唯一共同高端点，bα相对误差40%；计入已知实验误差后与部分公式经验边界重叠，不能仅据点估计越界称为错误。
- 229Th基态IA文献存在0.05%与0.01%分歧，HF会差5倍；原始采用值未改动。
- 大HF常伴弱α支路及较大l；低HF也可能伴Coriolis混合或集体八极带。三式无显式l项，HF不是剥离离心阻碍后的纯内禀量。DZR修正倍数仅供另一模型下的比较，不能直接套入三式。
- configuration_mixing仅4个母核，为小样本；逐线表保留其尾部，不声称是稳定异常阈值。
- VS本地四系数源于偶偶核拟合，奇核使用属于外推；Royer四奇偶参数复算一致；Rurarz本地算式复算一致，但原论文A与Ad的记号仍待核准，未声称原式已完全验证。
- 249Bk低端支路的原注释轨道记号与Adopted带描述不一致，已单列待回原论文核验，未按HF强行改标签。

## 全部候选（修订前固定选择）

| 母核→子核 | Ef keV | 标签 | lmin | HF VS | HF Rurarz | HF Royer | 核查结果 |
|---|---:|---|---:|---:|---:|---:|---|
| 227Ac→223Fr | 515.2 | nilsson_changed, spin_flip | 0 | 2.4905 | 1.5682 | 0.48393 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-3f88e62ed135.txt) |
| 226Th→222Ra | 242.157 | parity_change | 1 | 3.7256 | 4.7575 | 2.5436 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-7f73d69bd986.txt) |
| 228Th→224Ra | 215.985 | parity_change | 1 | 9.7604 | 9.5614 | 5.67 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-4a99f439e3ac.txt) |
| 229Th→225Ra | 0 | nilsson_changed | 2 | 53416 | 49604 | 12281 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-10449e0a1703.txt) |
| 229Th→225Ra | 179.75 | nilsson_changed | 0 | 20.911 | 17.679 | 4.5365 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-10449e0a1703.txt) |
| 229Th→225Ra | 243.56 | nilsson_changed | 2 | 16.782 | 13.687 | 3.5636 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-10449e0a1703.txt) |
| 231Pa→227Ac | 109.989 | nilsson_changed, parity_change | 3 | 22.791 | 19.801 | 6.1044 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ef76c43a3003.txt) |
| 231Pa→227Ac | 330.037 | same_band | 0 | 2.3841 | 1.8235 | 0.58693 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ef76c43a3003.txt) |
| 230U→226Th | 230.37 | parity_change | 1 | 13.393 | 15.436 | 8.8879 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-4f3b1613cbb1.txt) |
| 233U→229Th | 0 | configuration_mixing | 0 | 1.2767 | 1.0232 | 0.28813 | 需降为敏感性/非单线记录；[原始A记录](sources/ensdf-ee6ecdf657a4.txt) |
| 233U→229Th | 42.4333 | same_band, configuration_mixing | 2 | 3.7119 | 2.8992 | 0.82718 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ee6ecdf657a4.txt) |
| 233U→229Th | 195.707 | configuration_mixing | 4 | 517.14 | 366.54 | 110.02 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ee6ecdf657a4.txt) |
| 233U→229Th | 261.967 | configuration_mixing | 2 | 589.64 | 398.26 | 122.68 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ee6ecdf657a4.txt) |
| 233U→229Th | 327.86 | configuration_mixing | 6 | 578.77 | 373.7 | 117.97 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ee6ecdf657a4.txt) |
| 235U→231Th | 236.893 | nilsson_changed | 2 | 11.878 | 6.8723 | 2.513 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-ad0f705d1e04.txt) |
| 235Np→231Pa | 84.2152 | nilsson_changed | 0 | 3.0312 | 2.6973 | 0.90596 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-0ae5bf56c5af.txt) |
| 235Np→231Pa | 101.409 | nilsson_changed | 2 | 5.2129 | 4.5953 | 1.5493 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-0ae5bf56c5af.txt) |
| 237Np→233Pa | 86.482 | nilsson_changed | 0 | 4.6593 | 3.6064 | 1.4031 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-757946240509.txt) |
| 237Np→233Pa | 103.651 | nilsson_changed | 2 | 7.3583 | 5.6338 | 2.2027 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-757946240509.txt) |
| 237Np→233Pa | 109.05 | nilsson_changed | 2 | 16.983 | 12.961 | 5.075 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-757946240509.txt) |
| 237Np→233Pa | 237.904 | configuration_mixing | 0 | 3.1618 | 2.2124 | 0.90216 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-757946240509.txt) |
| 239Pu→235U | 197.119 | same_band | 6 | 1427 | 1203.9 | 385.09 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-412f45083d69.txt) |
| 239Pu→235U | 338.52 | parity_change | 9 | 55404 | 42930 | 14447 | 需降为敏感性/非单线记录；[原始A记录](sources/ensdf-412f45083d69.txt) |
| 239Pu→235U | 821.25 | same_band | 2 | 1226.3 | 665.12 | 280.11 | 需降为敏感性/非单线记录；[原始A记录](sources/ensdf-412f45083d69.txt) |
| 241Pu→237U | 204.3 | same_band | 2 | 3.8081 | 2.9847 | 1.0878 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-6a3f19306351.txt) |
| 241Am→237Np | 267.566 | nilsson_changed | 2 | 10868 | 10677 | 4113.2 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-6713f67fe36e.txt) |
| 241Am→237Np | 305.05 | same_band | 4 | 1385.4 | 1336.9 | 519.36 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-6713f67fe36e.txt) |
| 241Am→237Np | 332.382 | nilsson_changed, parity_change, spin_flip | 3 | 2.2685e+05 | 2.1585e+05 | 84415 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-6713f67fe36e.txt) |
| 241Am→237Np | 395.57 | same_band | 6 | 1331.8 | 1225.8 | 487.12 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-6713f67fe36e.txt) |
| 241Cm→237Pu | 0 | nilsson_changed | 3 | 6080.7 | 7973.7 | 1929.5 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-021dbc42ddbe.txt) |
| 241Cm→237Pu | 369.6 | spin_flip | 2 | 149.88 | 174.09 | 44.907 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-021dbc42ddbe.txt) |
| 242Cm→238Pu | 763.22 | parity_change | 5 | 26067 | 26061 | 23007 | 需降为敏感性/非单线记录；[原始A记录](sources/ensdf-0cf17355062d.txt) |
| 245Cm→241Pu | 95.786 | nilsson_changed | 2 | 6927.3 | 7125.7 | 2337.5 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-01ee4c6bf708.txt) |
| 245Cm→241Pu | 161.73 | nilsson_changed | 4 | 5705.2 | 5682.1 | 1902.6 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-01ee4c6bf708.txt) |
| 247Cm→243Pu | 124.8 | nilsson_changed, parity_change | 1 | 7200.8 | 6258.6 | 2501.3 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-d66fd70886d4.txt) |
| 247Bk→243Am | 84 | parity_change | 1 | 30.641 | 34.548 | 14.792 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-96fa71b358ad.txt) |
| 247Bk→243Am | 108 | parity_change | 3 | 30.266 | 33.808 | 14.544 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-96fa71b358ad.txt) |
| 249Bk→245Am | 19.3 | nilsson_changed | 0 | 6.9783 | 6.7682 | 3.4474 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-de78882abd05.txt) |
| 249Cf→245Cm | 0 | nilsson_changed, parity_change | 1 | 6021.4 | 7841 | 2296.2 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-c6cd481d89aa.txt) |
| 249Cf→245Cm | 54.784 | nilsson_changed, parity_change | 1 | 5998.9 | 7677.2 | 2273.1 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-c6cd481d89aa.txt) |
| 249Cf→245Cm | 121.574 | nilsson_changed, parity_change | 1 | 10843 | 13573 | 4076.3 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-c6cd481d89aa.txt) |
| 249Cf→245Cm | 199.8 | nilsson_changed, parity_change | 3 | 37025 | 45108 | 13789 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-c6cd481d89aa.txt) |
| 249Cf→245Cm | 674.2 | same_band | 4 | 1788.8 | 1790 | 626.71 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-c6cd481d89aa.txt) |
| 251Cf→247Cm | 227.379 | nilsson_changed | 2 | 24.2 | 27.677 | 9.6173 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-5658c4646cf5.txt) |
| 251Es→247Bk | 0 | same_intrinsic | 0 | 3.7667 | 5.1715 | 2.1433 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-fc6312798ee9.txt) |
| 255Fm→251Cf | 319.643 | nilsson_changed | 2 | 4848.3 | 6925.5 | 2100.1 | 未发现转录或换算错误（冻结ENSDF级）；[原始A记录](sources/ensdf-4a74247305f3.txt) |

详见交互表 OUTLIER_AUDIT.html、逐线CSV outlier_audit.csv、逐标签阈值表 outlier_label_memberships.csv。审计保留修订前阈值，不通过反复删除尾部点制造稳定分布。

来源：[ENSDF手册第12节](https://www.nndc.bnl.gov/ensdf/ensdf-manual.pdf)、[Royer与Zhang原文](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevC.77.037602/fulltext)、[VS参数出处](https://periodic-table.com/Preprints/2004/160(E7-2004-160).pdf)。各跃迁的冻结源与在线评估表链接保存在CSV中；本次未重读46条跃迁的全部原始实验论文。
