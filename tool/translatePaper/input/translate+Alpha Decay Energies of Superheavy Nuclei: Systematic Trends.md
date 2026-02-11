# 翻译文稿: Alpha Decay Energies of Superheavy Nuclei: Systematic Trends.pdf

*翻译时间: 2026-02-04 17:44:18*
*模型: gemini-2.5-flash | temperature=0.2*
*输入方式: inline*

---

# 超重核 $\alpha$ 衰变能：系统趋势

E. Olsen¹ 和 W. Nazarewicz¹
¹ 密歇根州立大学物理与天文学系及 FRIB 实验室，美国密歇根州东兰辛 48824
(日期：22019 年 1 月 16 日)

**摘要**

**背景**：新的超重核通常通过其特征性的 $\alpha$ 衰变能来识别，这需要对 $Q_\alpha$ 值进行精确计算。尽管有许多 $Q_\alpha$ 预测可用，但对其不确定性知之甚少，这使得将结果外推到未知系统变得困难。

**目的**：本工作旨在分析几种模型，将其预测与现有实验数据进行比较，并研究它们对当前实验关注的未观测到的 $^{296}120$ 和 $^{298}120$ $\alpha$ 衰变链的性能。我们量化的结果也将作为未来更复杂的统计研究的基准。

**方法**：我们使用核超流密度泛函理论 (DFT) 和几种 Skyrme 能量密度泛函 (EDFs)。为了估计系统模型不确定性，我们采用均匀模型平均。

**结果**：我们评估了从 Fm 到 Z = 120 的偶偶核的 $Q_\alpha$ 值。对于 Fm 到 Ds 之间形变良好的核，我们发现不同模型预测之间具有极好的一致性，并且与实验结果吻合良好。对于 Ds 之外的过渡核，模型间差异增大，导致系统误差显著。特别是，我们的模型低估了最重核 $^{294}Og$ 的 $Q_\alpha$ 值。

**结论**：形变良好的超重核的 DFT 预测的稳健性支持了将实验 $Q_\alpha$ 值与理论预测结合作为合理的 (Z, A) 指示器的想法。不幸的是，这种识别方法预计在接近 N = 184 的形变到球形形状转变区域效果不佳。在开发新的谱学质量 EDF 和更复杂的不确定性量化统计技术方面取得进展，将极大地有益于利用 $Q_\alpha$ 值识别新的超重核。

---

## I. 引言

质子数 Z > 104 的超重核占据了核素图的右上角 [1, 2]。对这些大质量体系的研究，是出于回答核物理、原子物理和化学领域许多基本问题的愿望 [3, 4]。

特别是，几十年来，人们一直在积极寻找自然界中长寿命的超重核。早期的理论计算预测超重幻数（所谓的“稳定岛”[5]）位于 Z = 114 和 N = 184 [6-8]。随着时间的推移和模型的改进，超重幻数被建议为质子数 114、120、124 或 126，中子数 172 或 184 [9-14]。然而，与传统幻数不同，这些对超重核的预测更可能对应于延长的半衰期而非稳定系统 [15]；这归因于大的库仑斥力以及单粒子能级的高密度 [16, 17]，导致弥散的壳结构 [16-18]。

通过冷、热重离子熔合 [19, 20] 的实验技术，在过去十年中发现了 Z = 114 (Fl) 和 Z = 118 (Og) 之间的新元素同位素，并将其添加到核素图中 [21-23]。目前，识别 Og 之外和更富中子体系的核的努力尚未成功 [24-26]。

已知的超重核主要通过 $\alpha$ 衰变和自发裂变 [27-29] 衰变。因此，新同位素通常通过观察其特征性的 $\alpha$ 衰变链 [30] 来识别，这基于实验数据和理论预测。因此，对具有量化不确定性的 $Q_\alpha$ 值进行计算对于未来的超重核搜索非常有用。

有许多关于超重核 $Q_\alpha$ 值的计算，例如参考文献 [16, 31-43]。其中一些研究还包括使用经验公式 [5, 44-51] 计算 $\alpha$ 衰变半衰期，其中半衰期表示为 $Q_\alpha$ 的函数。在这方面，$Q_\alpha$ 值和半衰期包含相同的信息内容。除了最近的综述 [16, 36] 之外，理论研究的重点是特定模型的性能。本文的目的是采取另一种方法：分析和比较几种 Skyrme-DFT 模型预测的 $Q_\alpha$ 值。通过这种方式，可以通过直接分析不同参数化和模型混合，更彻底地估计它们的系统不确定性和稳健性。

本文结构如下。第二节讨论所使用的理论方法。第三节展示了我们对已发现的已知超重核以及未知的 $^{296}120$ 和 $^{298}120$ 的结果。它还展示了模型混合分析的结果。最后，第四节讨论了结论和展望。

## II. 理论方法

我们所有的计算都是在核密度泛函理论 (DFT) [52] 框架内进行的，其中体系的总能量表示为能量密度，能量密度是单体局域密度和流的泛函。核 DFT 是对复杂重核进行全局预测的首选工具。正如参考文献 [53] 所强调的，这种方法足够通用，可以应用于核素图的任何位置；它可以通过内禀对称性破缺来纳入核形变；并且可以为各种可观测量提供量化预测。核 DFT 的主要组成部分是能量密度泛函 (EDF)，它代表了一种有效的介质中核相互作用。EDF 包含许多耦合常数，这些常数根据选定的实验数据和理论伪数据 [52, 54, 55] 进行调整；根据所选择的优化方法和策略，这些低能耦合常数会发生变化，并开发出新的 EDF 参数化。

在这项工作中，我们选择了七种有效的 Skyrme EDF [56, 57]，它们在粒子-空穴通道中通过密度依赖的混合型零程配对项 [58] 增强：SkM* [59]、SLy4 [60]、SV-min [54]、UNEDF0 [55]、UNEDF1 [61]、UNEDF2 [62] 和 UNEDF1SO [63]。SkM* 泛函（侧重于表面能，以使用半经典方法正确解释 $^{240}Pu$ 的裂变位垒）和 SLy4（侧重于富中子核）因其作为传统 Skyrme EDF 的价值而被纳入，并作为新参数化性能的基准。EDF SV-min 通过半幻核的结合能、电荷半径和衍射半径以及表面厚度进行参数化。UNEDF0 参数化针对球形核和形变核的结合能、电荷半径和奇偶结合能差进行了优化。EDF UNEDF1 是为裂变研究而开发的，并通过包含新的质量和裂变同质异能态的激发能扩展了 UNEDF0 的数据集。UNEDF2 泛函考虑了以前 UNEDF 参数化中忽略的张量项；它是为壳结构研究而开发的，并通过双幻核的单粒子能级劈裂扩展了 UNEDF1 的数据集。最后，UNEDF1SO 是一个在超铀区局部优化的 EDF，其自旋-轨道和配对参数经过微调，以在 $^{249}Bk$ 和 $^{251}Cf$ 的激发谱和奇偶质量差方面取得更好的一致性（所有其他参数与 UNEDF1 相同）。选择几种基于不同优化方法的 EDF，可以估计系统误差。

我们进行计算的程序与我们之前关于核滴线的工作 [64] 中使用的程序相同。对于给定的核，我们求解核 DFT 的 Hartree-Fock-Bogoliubov (HFB) 方程 [65] 以找到其基态结合能和其他全局核性质。鉴于形状形变对核结合能的影响，有必要对几种不同的核构型求解 HFB 方程；然后记录每个核对应于最小结合能的核形变（以及其他全局核性质），并用于后续计算。由于需要对许多不同的核进行数千次计算，我们利用高性能计算来加快这一过程。

由于我们的工作重点是超重体系，我们仅限于质子数 98 < Z < 120 的核。此外，我们仅限于偶数质子和中子的核，以避免与奇A和奇奇体系相关的复杂性 [66-68]。为了进行计算，我们使用了 DFT 代码 HFBTHOv300 [69]，它通过在形变谐振子基中直接对角化来求解 HFB 方程。我们对四极形变 $\beta_2$ 施加了约束，以解释长椭球形、扁椭球形和球形形变。为了加快计算速度，我们施加了轴向和反射对称性。尽管该区域存在三轴形变已得到充分证实 [70]，但它们对基态结合能的影响预计很小 [71]，因此这是一个合理的近似。为了近似恢复 HFB 方法中被破坏的粒子数对称性，我们使用了参考文献 [72] 中概述的 Lipkin-Nogami 方法。

虽然对于几种 EDF，参考文献 [64] 中计算的质量表已存储在 MassExplorer 数据库 [73] 中，但为了确保四极形变不会突然跳到极端值，我们重新计算了所有质量表并更新了 MassExplorer。从计算出的结合能中，我们提取了 $Q_\alpha$ 值：

$$Q_\alpha = 28.3 \text{ MeV} – BE(Z, N) + BE(Z – 2, N – 2).$$ (1)

为了评估我们的 $Q_\alpha$ 值的质量，我们采用了两种方法。第一种是直接分析所有 7 个 EDF 的结果，并将它们相互比较以及与现有实验值进行比较。第二种方法用于估计系统不确定性，即混合几种模型。

## III. 结果

我们首先评估我们对选定核的 $\alpha$ 衰变链的计算。在图 1 中，我们将 $^{270}Ds$ $\alpha$ 衰变链的 $Q_\alpha$ 值的核 DFT 结果与实验数据进行比较；选择该核是因为其 $\alpha$ 衰变链中每个核都有实验数据。我们在这里观察到的第一件事是 Skyrme EDF 结果的总体一致性：除了 SkM* 对 $^{266}Hs$ 的结果外，每个预测的 $Q_\alpha$ 值都遵循实验数据的模式，并随着质量数的增加而增加。我们还注意到，当排除 SkM* 时，计算出的 $Q_\alpha$ 值的离散度小于 1 MeV，并且每个数据点都位于此范围内（$^{266}Hs$ 是一个例外，它比 SV-min 最接近的结果高 25 keV）。虽然它几乎总是高估实验数据，但 SkM* 在此处的性能并不令人惊讶，因为后来的 EDF 得到了改进。

**[图 1]**：使用几种 Skyrme EDF 计算的 $^{270}Ds$ $\alpha$ 衰变链中核的 $Q_\alpha$ 值。参考文献 [74] 的实验值用星号表示。
*(图1展示了$^{270}Ds$衰变链中核的$Q_\alpha$值随质量数的变化。不同Skyrme EDF模型的结果与实验值（星号）进行了比较，显示出总体一致性，且$Q_\alpha$值随质量数增加。)*

在图 2 中，我们展示了 $^{292}Lv$ $\alpha$ 衰变链的 $Q_\alpha$ 结果，其中大量外推值突显了该区域实验数据的稀缺性。直到 Ds，该模式与图 1 中给定外推值的模式非常相似。然而，在 $^{284}Cn$ 处，外推的 $Q_\alpha$ 值略有下降，随后在 $^{288}Fl$ 处实验 $Q_\alpha$ 值有所增加。这归因于该区域 [36, 70] 的三轴软度导致在 N = 174 附近从长椭球形到扁椭球形形变的突然转变。结果，EDF 结果的一致性受到影响，最明显的是 UNEDF0，其结果从 Fl 到 Lv 降低，而实验数据增加。形状转变对 UNEDF1 和 UNEDF2 影响的减小似乎突显了在全局 EDF 优化中包含裂变同质异能态和单粒子能级数据的必要性。

**[图 2]**：与图 1 类似，但针对 $^{292}Lv$。参考文献 [74] (Syst) 中系统趋势外推的推荐值也用星号表示。
*(图2展示了$^{292}Lv$衰变链中核的$Q_\alpha$值随质量数的变化。与图1类似，不同Skyrme EDF模型的结果与实验值（星号）进行了比较，但在$^{284}Cn$附近出现形状转变，导致模型一致性下降。)*

我们还希望将分析扩展到尚未观测到的核。为此，我们在图 3 中展示了 $^{296}120$ 和 $^{298}120$ $\alpha$ 衰变链的 $Q_\alpha$ 值结果。就像图 1 和图 2 中一样，从 Fm 到 Ds 的模式对于所有 EDF 和实验都重复出现。在达到长椭球形到扁椭球形转变时，我们看到了与图 2 结果相似的行为：对于 $^{296}120$，我们注意到从系统趋势获得的 $^{280}Cn$ 和 $^{284}Fl$ 值低于 Skyrme EDF 结果的分布，尽管 SLy4、UNEDF1、UNEDF1SO 和 UNEDF2 对 $^{284}Fl$ 的结果相当接近。对于 $^{298}120$，我们再次看到 $^{282}Cn$ 的推荐值低于计算结果的范围，而 $^{286}Fl$ 和 $^{290}Lv$ 的实验数据在其范围内（$^{286}Fl$ 之前缺乏实验数据是由于 $^{282}Cn$ 中裂变占主导地位 [75]，这阻止了该区域 $Q_\alpha$ 值的测量）。这里最引人注目的特征是实验数据点与 $^{294}Og$ 的计算结果之间的比较，其中每个 EDF 都低估了实验值。UNEDF1SO 最接近该值的事实很可能归因于它预测从扁椭球形形变到球形形状转变始于 N = 178 附近，这可能表明自旋-轨道力对该核的重要性。同样，UNEDF0 对实验值的显著低估，与 UNEDF1、UNEDF1SO 和 UNEDF2 相比，说明了在未来 EDF 中纳入更多数据的必要性。我们还注意到，每个计算出的 $Q_\alpha$ 值对于 $^{296}120$ 和 $^{298}120$ 都有增加的趋势；鉴于从 Fl 到 Og 实验数据中看到的模式，这种行为是有希望的。

**[图 3]**：与图 1 和图 2 类似，但针对 $^{296}120$ (a) 和 $^{298}120$ (b)。
*(图3展示了$^{296}120$ (a) 和 $^{298}120$ (b) 衰变链中核的$Q_\alpha$值随质量数的变化。与前两图类似，模型结果与实验/推荐值进行了比较，但在重核区域，特别是$^{294}Og$，模型低估了实验值。)*

图 4 显示了从 Fm 到 Z = 120 的核同位素链的 $Q_\alpha$ 值分析。标记了长椭球形、扁椭球形形变和球形形状区域之间的边界。（有关其他模型中形变预测的讨论，请参见参考文献 [36]。）在 N = 164 附近观察到的不规则性，特别是对于 SLy4、SV-min 和 UNEDF1，是由于长椭球形形变的子壳闭合 [76, 77]。再次，我们观察到每个 EDF 的总体一致性，在每个同位素链的理论计算中出现了相似的模式，即使在形状转变发生的区域也是如此。我们的理论结果与实验值的接近程度，通过均方根 (rms) 偏差 $\delta(Q_\alpha)$ 表示，是相当合理的。如前所述，与实验的最大偏差出现在最重元素 Lv 和 Og 处，它们被预测位于长椭球形到扁椭球形形状转变区域。在检查单个 rms 偏差时，表现最好的模型是 UNEDF1，其 $\delta(Q_\alpha) = 0.31$ MeV，而最早的 Skyrme EDF SkM* 产生了 $\delta(Q_\alpha) = 0.81$ MeV。总的来说，这里的 rms 偏差范围与其他 DFT 模型 [36] 的结果一致。例如，如果考虑参考文献 [16] 的相对论 EDF，$\delta(Q_\alpha)$ 值范围从 PC-PK1 的 0.32 MeV 到 NL3* 的 0.68 MeV。

**[图 4]**：本研究中使用的每个全局 Skyrme EDF 计算的从 Fm 到 Z = 120 的核同位素链的 $Q_\alpha$ 值。标记了长椭球形、扁椭球形和球形形状区域。参考文献 [74] 的实验数据用圆圈表示，并与代表理论计算的相应线条颜色匹配（我们注意到目前没有 Cn 偶偶同位素的实验数据）。每个模型的均方根 (rms) 偏差（以 MeV 为单位）已标出。对于局部 EDF UNEDF1SO（未显示），rms 偏差为 0.46 MeV。
*(图4展示了从Fm到Z=120的核同位素链的$Q_\alpha$值随中子数的变化，并标记了核的形状区域。不同Skyrme EDF模型的结果与实验数据（圆圈）进行了比较，并给出了每个模型的均方根偏差。)*

评估模型在缺乏实验数据区域的预测不确定性是现代核理论的核心问题。到目前为止，我们通过使用许多不同的 EDF 并比较它们各自的性能与实验来估计不确定性。根据参考文献 [64]，我们现在计算几种模型结果的均匀平均值以及相应的标准偏差，以确定系统不确定性。为此，我们选择了 SV-min、UNEDF0、UNEDF1 和 UNEDF2，因为它们是本研究中使用的最新开发的全局 EDF。虽然这个过程可能看起来很幼稚，但在没有额外信息或昂贵的统计计算的情况下，均匀权重的选择本质上是最佳的 [78]。此外，通过在平均中给予每个模型相同的权重，我们可以了解更复杂的模型混合可能如何表现。

在图 5 中，我们展示了 $^{296}120$ 和 $^{298}120$ $\alpha$ 衰变链的模型平均结果。对于这两条链，在 Fm 到 Ds 之间，我们看到了计算值与实验/推荐值之间极好的一致性。然而，从 Cn 及以后，形状转变的影响破坏了这种一致性。如前所述，这种差异在 $^{298}120$ $\alpha$ 衰变链中的 $^{294}Og$ 处尤为明显；这很可能归因于 UNEDF0 预测的低 $Q_\alpha$ 值（参见图 3）。然而，即使将 UNEDF0 从平均中排除（并包含 $Q_\alpha$ 值大得多的 UNEDF1SO），差异仍然很大。在图 6 中，我们展示了存在实验数据的同位素链的 $Q_\alpha$ 值的模型平均结果。从 Fm 到 Fl，与实验数据的接近程度非常好，误差带相对较小。然而，对于 Lv 和 Og，我们看到了与图 5 中类似的行为。

**[图 5]**：使用三种 UNEDF 模型和 SV-min 计算的 $^{296}120$ (a) 和 $^{298}120$ (b) $\alpha$ 衰变链中核的模型平均 $Q_\alpha$ 值。假设均匀模型权重。理论预测的误差棒代表标准偏差。参考文献 [74] 的实验值和推荐值用星号表示，并带有相应的误差棒。
*(图5展示了$^{296}120$ (a) 和 $^{298}120$ (b) 衰变链中核的模型平均$Q_\alpha$值随质量数的变化。模型平均结果与实验/推荐值（星号）进行了比较，在较轻核区域一致性良好，但在重核区域，特别是$^{294}Og$，差异显著。)*

**[图 6]**：使用三种 UNEDF 模型和 SV-min 计算的从 Fm 到 Z = 120 的偶偶核同位素链（不包括 Cn）的模型平均 $Q_\alpha$ 值。假设均匀模型权重。误差棒和误差带代表标准偏差。参考文献 [74] 的实验数据用圆圈表示，并与代表理论预测的线条颜色匹配。模型平均结果与实验数据的均方根偏差为 0.35 MeV。
*(图6展示了从Fm到Z=120的偶偶核同位素链的模型平均$Q_\alpha$值随中子数的变化。模型平均结果与实验数据（圆圈）进行了比较，在较轻核区域吻合良好，但在重核区域，特别是Lv和Og，与实验值存在偏差。)*

## IV. 结论

在本文中，我们研究了在核 DFT 框架内，使用几种不同的 Skyrme EDF 计算的从 Fm 到 Z = 120 的偶偶超重核的 $Q_\alpha$ 值。为了估计系统误差，我们通过比较单个模型和模型平均来分析 $\alpha$ 衰变链的理论预测。在形变良好的超重核区域，理论预测是稳健的，每个 EDF 都给出了相对一致的结果。这种稳健性对于形状过渡核有所降低。总的来说，观察到的与实验数据的一致性是相当合理的。

单个泛函，特别是 UNEDF 家族的行为，也证明了其启发性。在所使用的模型中，表现最好的是 UNEDF1，其均方根偏差 $\delta(Q_\alpha) = 0.31$ MeV。UNEDF1 和 UNEDF2 在形状转变区域的结果优于 UNEDF0，这表明裂变同质异能态和单准粒子态数据在 EDF 优化中的重要性。我们还分析了 UNEDF1SO 泛函的性能，该泛函在超铀同位素 Bk 和 Cf 区域进行了局部优化。鉴于其微调，其 $Q_\alpha$ 值性能与其他 UNEDF 参数化相比相似或略差，这很有趣。

总的来说，通过 $Q_\alpha$ 值识别核素的方法预计在形变到球形形状转变区域效果不佳。在这种背景下，理论将受益于开发新的谱学质量全局 EDF 和更复杂的统计不确定性量化技术方面的进展。实验上，在不使用 $Q_\alpha$ 值的情况下，识别上超重（热熔合）区域新超重核的工作已经在进行中 [79-81]。

随着对 Og 之外元素 [82-85] 的搜索继续进行，对 $Q_\alpha$ 值进行精确计算将变得更加有益。我们模型混合结果在评估和降低不确定性方面的性能似乎很有希望。此外，通过更复杂的模型混合技术，利用贝叶斯模型平均 [78, 86, 87]，其中简单平均通过整合参数空间上的相应似然函数计算出的模型后验概率进行重新加权，预计可预测性将进一步提高。然而，在不久的将来，为了对 $Q_\alpha$ 值进行更可靠的外推，我们打算使用最近参考文献 [88] 中描述的贝叶斯机器学习技术。

## V. 致谢

感谢 Nicolas Schunck 的讨论。劳伦斯利弗莫尔国家实验室的 Quartz 计算设施对这项工作至关重要。这项工作得到了美国能源部以下奖项号的支持：DOE-DE-NA0002847 (NNSA, Stewardship Science Academic Alliances program)、de-sc0013365 (Office of Science) 和 de-sc0018083 (Office of Science, NUCLEI SciDAC-4 collaboration)。

[1] D. C. Hoffman, A. Ghiorso, and G. T. Seaborg, The Transuranium People: The Inside Story (Imperial College Press, 2000).
[2] G. T. Seaborg and W. D. Loveland, The Elements Beyond Uranium (Wiley-Interscience, 1990).
[3] S. Giuliani, Z. Matheson, W. Nazarewicz, E. Olsen, P.-G. Reinhard, J. Sadhukhan, B. Schuetrumpf, N. Schunck, and P. Schwerdtfeger, Rev. Mod. Phys. (2018).
[4] W. Nazarewicz, Nature Physics 14, 537 (2018).
[5] V. Viola and G. Seaborg, J. Inorg. Nucl. Chem. 28, 741 (1966).
[6] W. D. Myers and W. J. Swiatecki, Nucl. Phys. 81, 1 (1966).
[7] S. G. Nilsson, C. F. Tsang, A. Sobiczewski, Z. Szymański, S. Wycech, C. Gustafson, I.-L. Lamm, P. Möller, and B. Nilsson, Nucl. Phys. A 131, 1 (1969).
[8] A. Sobiczewski, F. Gareev, and B. Kalinkin, Phys. Lett. 22, 500 (1966).
[9] A. V. Afanasjev, T. L. Khoo, S. Frauendorf, G. A. Lalazissis, and I. Ahmad, Phys. Rev. C 67, 024309 (2003).
[10] A. V. Afanasjev, Phys. Scr. 2006, 62 (2006).
[11] M. Bender, K. Rutz, P.-G. Reinhard, J. A. Maruhn, and W. Greiner, Phys. Rev. C 60, 034304 (1999).
[12] S. Ćwiok, J. Dobaczewski, P.-H. Heenen, P. Magierski, and W. Nazarewicz, Nucl. Phys. A 611, 211 (1996).
[13] A. T. Kruppa, M. Bender, W. Nazarewicz, P.-G. Rein-hard, T. Vertse, and S. Ćwiok, Phys. Rev. C 61, 034313 (2000).
[14] K. Rutz, M. Bender, T. Bürvenich, T. Schilling, P.-G. Reinhard, J. A. Maruhn, and W. Greiner, Phys. Rev. C 56, 238 (1997).
[15] C. E. Düllmann and M. Block, Sci. Am. 318, 46 (2018).
[16] S. E. Agbemava, A. V. Afanasjev, T. Nakatsukasa, and P. Ring, Phys. Rev. C 92, 054310 (2015).
[17] M. Bender, W. Nazarewicz, and P.-G. Reinhard, Phys. Lett. B 515, 42 (2001).
[18] P. Jerabek, B. Schuetrumpf, P. Schwerdtfeger, and W. Nazarewicz, Phys. Rev. Lett. 120, 053001 (2018).
[19] Y. Oganessian, A. Demin, A. Iljinov, S. Tretyakova, A. Pleve, Y. Penionzhkevich, M. Ivanov, and Y. Tretyakov, Nucl. Phys. A 239, 157 (1975).
[20] Y. Oganessian, J. Phys. G 34, R165 (2007).
[21] R. Barber, P. J. Karol, H. Nakahara, E.Vardaci, and E. W. Vogt, Pure Appl. Chem. 83, 1485 (2011).
[22] P. J. Karol, R. C. Barber, B. M. Sherrill, V. Emanuele, and Y. Toshimitsu, Pure Appl. Chem. 88, 139 (2016).
[23] P. J. Karol, R. C. Barber, B. M. Sherrill, V. Emanuele, and Y. Toshimitsu, Pure Appl. Chem. 88, 155 (2016).
[24] C. E. Düllmann, GSI Scientific Report 2011 GSI Report 2012-1, 206 (2011).
[25] C. E. Düllmann, EPJ Web Conf. 163, 00015 (2017).
[26] Y. T. Oganessian et al., Phys. Rev. C 79, 024603 (2009).
[27] S. Hofmann, F. Heßberger, D. Ackermann, S. An-talic, P. Cagarda, B. Kindler, P. Kuusiniemi, M. Leino, B. Lommel, O. Malyshev, R. Mann, G. MuÂĺnzenberg, A. Popeko, S. S'aro, B. Streicher, and A. Yeremin, Nucl. Phys. A 734, 93 (2004).
[28] Y. T. Oganessian and V. K. Utyonkov, Rep. Prog. Phys. 78, 036301 (2015).
[29] Y. T. Oganessian, A. Sobiczewski, and G. M. Ter-Akopian, Phys. Scr. 92, 023003 (2017).
[30] Y. T. Oganessian et al., Phys. Rev. C 74, 044602 (2006).
[31] M. Bender, Phys. Rev. C 61, 031302 (2000).
[32] J. F. Berger, D. Hirata, M. Girod, and J. Decharge, Int. J. Mod. Phys. E 13, 79 (2004).
[33] S. Ćwiok, W. Nazarewicz, and P.-H. Heenen, Phys. Rev. Lett. 83, 1108 (1999).
[34] J. Erler, K. Langanke, H. P. Loens, G. Martínez-Pinedo, and P.-G. Reinhard, Phys. Rev. C 85, 025802 (2012).
[35] Y. K. Gambhir, A. Bhagwat, and M. Gupta, Phys. Rev. C 71, 037301 (2005).
[36] P.-H. Heenen, J. Skalski, A. Staszczak, and D. Vretenar, Nucl. Phys. A 944, 415 (2015).
[37] P. Jachimowicz, M. Kowal, and J. Skalski, Phys. Rev. C 89, 024304 (2014).
[38] I. Muntian, S. Hofmann, Z. Patyk, and A. Sobiczewski, Acta Phys. Pol. 34, 2073 (2003).
[39] A. Sobiczewski and K. Pomorski, Prog. Part. Nucl. Phys. 58, 292 (2007).
[40] S. V. Tolokonnikov, Y. S. Lutostansky, and E. E. Saper-stein, Phys. At. Nucl. 76, 708 (2013).
[41] S. V. Tolokonnikov, I. N. Borzov, M. Kortelainen, Y. S. Lutostansky, and E. E. Saperstein, Eur. Phys. J. A 53, 33 (2017).
[42] S. Typel and B. A. Brown, Phys. Rev. C 67, 034313 (2003).
[43] M. Warda and J. L. Egido, Phys. Rev. C 86, 014322 (2012).
[44] B. A. Brown, Phys. Rev. C 46, 811 (1992).
[45] A. Budaca, R. Budaca, and I. Silisteanu, Nucl. Phys. A 951, 60 (2016).
[46] P. R. Chowdhury, C. Samanta, and D. N. Basu, Phys. Rev. C 77, 044603 (2008).
[47] J. Dong, W. Zuo, and W. Scheid, Nucl. Phys. A 861, 1 (2011).
[48] H. Koura, J. Nucl. Sci. Technol. 49, 816 (2012).
[49] A. Parkhomenko and A. Sobiczewski, Acta Phys. Pol. B 36, 3095 (2005).
[50] G. Royer and H. F. Zhang, Phys. Rev. C 77, 037602 (2008).
[51] D. E. Ward, B. G. Carlsson, and S. Åberg, Phys. Rev. C 92, 014314 (2015).
[52] M. Bender, P.-H. Heenen, and P.-G. Reinhard, Rev. Mod. Phys. 75, 121 (2003).
[53] M. Stoitsov, J. Dobaczewski, W. Nazarewicz, and P. Bo-rycki, Int. J. Mass Spectrom. 251, 243 (2006).
[54] P. Klüpfel, P.-G. Reinhard, T. J. Bürvenich, and J. A. Maruhn, Phys. Rev. C 79, 034310 (2009).
[55] M. Kortelainen, T. Lesinski, J. Moré, W. Nazarewicz, J. Sarich, N. Schunck, M. V. Stoitsov, and S. Wild, Phys. Rev. C 82, 024313 (2010).
[56] T. Skyrme, Nucl. Phys. 9, 615 (1958).
[57] D. Vautherin and D. M. Brink, Phys. Rev. C 5, 626 (1972).
[58] J. Dobaczewski, W. Nazarewicz, and M. Stoitsov, Eur. Phys. J. A 15, 21 (2002).
[59] J. Bartel, P. Quentin, M. Brack, C. Guet, and H.-B. Håkansson, Nucl. Phys. A 386, 79 (1982).
[60] E. Chabanat, P. Bonche, P. Haensel, J. Meyer, and R. Schaeffer, Nucl. Phys. A 635, 231 (1998).
[61] M. Kortelainen, J. McDonnell, W. Nazarewicz, P.-G. Reinhard, J. Sarich, N. Schunck, M. V. Stoitsov, and S. M. Wild, Phys. Rev. C 85, 024304 (2012).
[62] M. Kortelainen, J. McDonnell, W. Nazarewicz, E. Olsen, P.-G. Reinhard, J. Sarich, N. Schunck, S. M. Wild, D. Davesne, J. Erler, and A. Pastore, Phys. Rev. C 89, 054314 (2014).
[63] Y. Shi, J. Dobaczewski, and P. T. Greenlees, Phys. Rev. C 89, 034309 (2014).
[64] J. Erler, N. Birge, M. Kortelainen, W. Nazarewicz, E. Olsen, A. M. Perhac, and M. Stoitsov, Nature 486, 509 (2012).
[65] P. Ring and P. Schuck, The Nuclear Many-Body Problem (Springer, New York, 1980).
[66] A. V. Afanasjev, J. Phys. G. 42, 034002 (2015).
[67] L. Bonneau, P. Quentin, and P. Möller, Phys. Rev. C 76, 024320 (2007).
[68] N. Schunck, J. Dobaczewski, J. McDonnell, J. Moré, W. Nazarewicz, J. Sarich, and M. V. Stoitsov, Phys. Rev. C 81, 024316 (2010).
[69] R. N. Perez, N. Schunck, R.-D. Lasseri, C. Zhang, and J. Sarich, Comp. Phys. Comm. 220, 363 (2017).
[70] S. Ćwiok, P.-H. Heenen, and W. Nazarewicz, Nature 433, 705 (2005), review Article.
[71] P. Möller, R. Bengtsson, B. Carlsson, P. Olivius, T. Ichikawa, H. Sagawa, and A. Iwamoto, At. Data Nucl. Data Tables 94, 758 (2008).
[72] M. V. Stoitsov, J. Dobaczewski, R. Kirchner, W. Nazarewicz, and J. Terasaki, Phys. Rev. C 76, 014308 (2007).
[73] http://massexplorer.frib.msu.edu.
[74] M. Wang, G. Audi, F. Kondev, W. Huang, S. Naimi, and X. Xu, Chin. Phys. C 41, 030003 (2017).
[75] N. T. Brewer et al., Phys. Rev. C 98, 024317 (2018).
[76] S. Ćwiok, V. Pashkevich, J. Dudek, and W. Nazarewicz, Nucl. Phys. A 410, 254 (1983).
[77] P. Möller and J. R. Nix, J. Phys. G 20, 1681 (1994).
[78] J. M. Bernardo and A. F. M. Smith, Bayesian Theory (Wiley, New Jersey, 1994).
[79] J. M. Gates et al., Phys. Rev. C 92, 021301 (2015).
[80] J. M. Gates, EPJ Web Conf. 131, 08003 (2016).
[81] J. M. Gates et al., Phys. Rev. Lett. (2018).
[82] C. E. Düllmann, EPJ Web Conf. 131, 08004 (2016).
[83] F. P. Heßberger and D. Ackermann, Eur. Phys. J. A 53, 123 (2017).
[84] S. Hofmann et al., Eur. Phys. J. A 52, 180 (2016).
[85] J. B. Roberto and K. P. Rykaczewski, Sep. Sci. Technol. 53, 1813 (2018).
[86] J. A. Hoeting, D. Madigan, A. E. Raftery, and C. T. Volinsky, Statist. Sci. 14, 382 (1999).
[87] L. Wasserman, J. Math. Psych. 44, 92 (2000).
[88] L. Neufcourt, Y. Cao, W. Nazarewicz, and F. Viens, Phys. Rev. C 98, 034318 (2018).

[[TRANSLATION_END]]