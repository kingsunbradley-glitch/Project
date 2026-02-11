# 翻译文稿: PhysRevC.110.054314.pdf

*翻译时间: 2026-02-04 17:29:36*
*模型: gemini-2.5-flash | temperature=0.2*

---

PHYSICAL REVIEW C 110, 054314 (2024)

**变形超重核的稳定性：一项壳修正分析**

郭雨鑫 (Yuxin Guo) 和 石跃 (Yue Shi)*
哈尔滨工业大学物理系，哈尔滨 150001，中国

(2024年6月19日收到；2024年9月25日修订；2024年10月25日接受；2024年11月14日发表)

**背景**：各种平均场计算表明，大多数质子数 $Z$ 介于 102 和 120 之间的超重核 (SHN) 是变形的。超重核的存在归因于壳稳定性。虽然球形超重核的壳稳定性效应已得到充分研究，但基于自洽平均场模型的变形超重核研究较少。
**目的**：本工作的目的是通过从核密度泛函理论 (DFT) 获得的单粒子 (s.p.) 能量中提取壳修正能，研究质子数 $Z$ 介于 98 和 120 之间的变形偶偶超重核的稳定性。
**方法**：我们使用基于混合基的无对称性限制非相对论平均场模型计算单粒子能量。斯特鲁金斯基壳修正能是利用这些单粒子能量获得的。对于对力道，我们使用利普金-野上 (LN) 对力模型。
**结果**：我们使用四种不同的能量密度泛函 (EDFs) 对质子数 $Z$ 介于 98 和 120 之间、中子数 $N$ 介于 140 和 190 之间的偶偶核进行了系统计算。提取的壳修正能被绘制为质子数和中子数的函数。所有四种 EDF 都预测核 $^{278}Hs_{162}$ 处存在最低极小值。对于 SLy4、SkM* 和 UNEDF1 EDF，下一个极小值位于 $^{284}Lv_{168}$。对于 UNEDF1$^{SO}$ EDF，极小值位于 $^{276}Cn_{164}$。计算得到的基态总能量使我们能够提取双粒子分离能 ($\delta_{2n,2p}$)、$\alpha$ 衰变能 ($Q_{\alpha}$) 和 $\alpha$ 衰变半衰期。SLy4、UNEDF1 和 UNEDF1$^{SO}$ EDF 的计算结果合理地再现了实验数据。对双粒子分离能和 $\alpha$ 衰变能的分析表明，同位素链和同中子素链中 $\delta_{2n,2p}$ 和 $Q_{\alpha}$ 值的变化可归因于相关核子数处变形壳隙的出现。
**结论**：我们的计算结果表明，高于 $^{100}Fm_{152}$ 的下一个变形壳闭合出现在 $^{108}Hs_{162}$。超重核区域中更重的变形壳闭合预计位于 $^{284}Lv_{168}$ 或 $^{112}Cn_{164}$。

DOI: 10.1103/PhysRevC.110.054314

**I. 引言**

超重核 (SHN) 一词通常指质子数 $Z > 104$ 的超锕元素。超重核观察到的主要衰变模式是 $\alpha$ 衰变和自发裂变 (SF) [1,2]。因此，超重核通常通过 $\alpha$ 衰变链分析来识别 [3]。迄今为止，最重的超重核已观测到 $Og$ ($Z = 118$)，完成了元素周期表的第七周期 [4]。这些最重核是通过涉及稀有 $^{48}Ca$ 束流和锕系靶的“热熔合”反应产生的。最近，新同位素 $^{286}Mc$ ($Z = 115$) 在 $^{48}Ca + ^{243}Am$ 反应中被合成 [5]。已知超重核主要位于 $\beta$ 稳定性谷的缺中子一侧 [6,7]。由于没有明显的方法来合成富中子超重体系，所有关于这些核的信息都必须通过理论建模来推断，这通常涉及大量的外推。

理论上，对于超重核，液滴模型预测库仑相互作用变得如此之强，以至于原子核不再对自发裂变稳定。然而，原子核是量子体系，其特征在于潜在的单粒子 (s.p.) 轨道。结果表明，费米面附近单粒子轨道的不均匀性可以为超重核贡献额外的结合能。

斯特鲁金斯基开发了一种从壳效应产生的涨落能中提取能量的方法 [8-11]。总能量的这部分涨落被称为壳修正能。通过将液滴能量与从核平均场模型中提取的壳修正能相结合，统一模型在预测基态和一些低能性质方面取得了定量成功 [12-14]。使用这种模型，即宏观-微观 (MM) 模型 [15-18]，计算得到的壳修正能作为质子数和中子数的函数，在球形超重核的 $Z = 114$ 和 $N = 184$ 处达到峰值，表明这些核子数处存在壳隙。当允许变形时 [19]，壳隙出现在 $Z = 108$ 和 $N = 162$ 处，对应于 $^{270}Hs$，正如 MM 模型所预测的。

与 MM 模型相比，自洽平均场 (SCMF) 方法能够从有效的核子-核子相互作用或拉格朗日量出发计算核基态性质 [20,21]。对于仅限于球形变形的计算，壳修正能已从使用以下方法计算的单粒子能量中提取：

*Contact author: yueshi@hit.edu.cn

2469-9985/2024/110(5)/054314(22) 054314-1 ©2024 American Physical Society

**II. 模型**

SCMF 模型被广泛用于描述超重核的近基态性质。因此，我们的壳修正分析基于从核密度泛函理论 (DFT) 获得的单粒子能级。由于所研究的超重核是束缚核，我们采用利普金-野上 (LN) 模型而不是哈特里-福克-玻戈留波夫 (HFB) 理论来处理对力相互作用。单极对力强度通过平滑能级密度确定。因此，我们将在斯特鲁金斯基形式主义之后讨论对力，其中已引入平滑能级密度。

鉴于上述考虑，我们首先在 II A 节中简要描述核 DFT 及其数值解。在 II B 节中，我们提供了斯特鲁金斯基壳修正方法的详细形式。II C 节描述了本工作中采用的 LN 对力形式。

**A. 核密度泛函理论**

核 DFT 的基本组成部分是一组单粒子波函数 $\varphi_{\mu,q}(\mathbf{r}, \sigma)$ 和占据振幅 $v_{\mu,q}$，

$$
\{\varphi_{\mu,q}(\mathbf{r}, \sigma), v_{\mu,q}, u_{\mu,q}, \mu = 1, \dots, N_{wf}\}, \quad (1)
$$

其中指标 $q$ 标记中子 ($q = n$) 或质子 ($q = p$) 量。由于 (1) 我们采用 LN 对力，并且 (2) 计算的原子核是束缚体系，因此包含的单粒子波函数数量 $N_{wf}$ 被取为核子数的两倍。具体来说，我们取中子的 $N_{wf} = 2 \times N$ 和质子的 $N_{wf} = 2 \times Z$。占据 ($v_{\mu,q}$) 和非占据 ($u_{\mu,q}$) 振幅被选择为正实数。它们通过 $u_{\mu,q}^2 + v_{\mu,q}^2 = 1$ 连接。$v_{\mu,q}$ 的值通过求解 LN 问题确定 (参见 II C 节)。

进入 Skyrme 能量密度泛函 (EDF) 的量是局域密度和流。这些包括局域粒子数密度 $\rho_q$、动能密度 $t_q$、对力能量密度 $\tilde{\rho}_q$ 和自旋流密度 $\mathbf{J}_q$，它们用波函数 $\varphi_{\mu,q}$ 和 (非) 占据概率 $v_{\mu,q}$ ($u_{\mu,q}$) 表示。关于这些量的详细信息可在参考文献 [42] 中找到。

在核 DFT 中，核总能量 $E_{tot}$ 由动能、Skyrme、对力 和 库仑项组成：

$$
E_{tot} = E_{Skyrme} + E_{kin} + E_{Coul} + E_{pair} + E_{cm} \\
= \int d^3r[K(\mathbf{r}) + \mathcal{E}_{Skyrme}(\mathbf{r}) + \mathcal{E}_{pair}(\mathbf{r}) + \mathcal{E}_{Coul}(\mathbf{r})]. \quad (2)
$$

能量密度 $K(\mathbf{r})$、$\mathcal{E}_{Skyrme}(\mathbf{r})$、$\mathcal{E}_{pair}(\mathbf{r})$、$\mathcal{E}_{Coul}(\mathbf{r})$ 是局域密度和流的泛函。各种能量密度的具体形式可在参考文献 [42] 中找到。

将总能量 (2) 对 $\rho$ 最小化，得到 SHF 方程：

$$
\begin{pmatrix} h_{\uparrow\uparrow} & h_{\uparrow\downarrow} \\ h_{\downarrow\uparrow} & h_{\downarrow\downarrow} \end{pmatrix} \begin{pmatrix} \Phi_{\uparrow,\mu} \\ \Phi_{\downarrow,\mu} \end{pmatrix} = \epsilon_{\mu} \begin{pmatrix} \Phi_{\uparrow,\mu} \\ \Phi_{\downarrow,\mu} \end{pmatrix}, \quad (3)
$$

其中箭头 $\uparrow$ 和 $\downarrow$ 分别表示核子自旋 $\sigma = 1/2$ 和 $\sigma = -1/2$。在方程 (3) 中，我们以旋量形式写出波函数；即 $\varphi_{\uparrow,\mu} = \varphi_{\mu}(\mathbf{r}, \sigma = 1/2)$。为简化起见，我们省略了下标 $q$。

方程 (3) 中的哈密顿量为

$$
h_q^{\sigma\sigma'}(\mathbf{r}) = -\frac{\hbar^2}{2m^*(\mathbf{r})}\nabla^2 + U_q(\mathbf{r}) + U_{Coul}(\mathbf{r})\delta_{\sigma,\sigma'} \\
- [i\mathbf{B}_q \cdot (\nabla \times \mathbf{\sigma})]_{\sigma\sigma'}. \quad (4)
$$

有效质量 [$m^*(\mathbf{r})$]、有效自旋密度 [$\mathbf{B}_q(\mathbf{r})$]、核势 [$U_q(\mathbf{r})$] 和质子的库仑势 [$U_{Coul}(\mathbf{r})$] 的定义可在参考文献 [42] 中找到。

054314-2

**变形超重核的稳定性：一项壳修正分析**

在本工作中，我们通过首先使用混合基展开 SHF 哈密顿量 [方程 (3)]，然后对其进行对角化 [43] 来求解 SHF 方程 (4)。该基组包含 $x$ 和 $y$ 方向上的两个谐振子 (HO) 基，以及 $z$ 方向上的笛卡尔基。我们使用两个数 $N_{max}$ 和 $N_z$ 来表征基组的维度 $n_x n_y i_z)$。HO 的基数通过要求 $n_x + n_y \le N_{max}$ 来确定。$N_z$ 值指定 $z$ 方向上的网格点数，网格点间距设置为 $d_z = 1.0 \text{ fm}$。

在我们的计算中，HF 计算的维度为 $(N_{max} + 1) \times (N_{max} + 2) \times N_z \times 2$。对于我们的超重核计算，我们设置 $N_{max} = 14$ 和 $N_z = 13$。计算得到的单粒子能量与使用 HFODD 代码 [44] 获得的能量非常相似 (差异小于 $100 \text{ keV}$)。

在实际计算中，我们首先使用一组单粒子波函数 $\varphi_{\mu,q}(\mathbf{r}, \sigma)$，这些波函数是使用 Woods-Saxon 势计算的。连同初始占据振幅 $v_{\mu}$，可以计算密度和势。对单粒子哈密顿量 [方程 (3)] 进行对角化，得到一组新的单粒子波函数及其对应的单粒子能量 $\epsilon_{\mu}$。随后，使用 $\epsilon_{\mu}$ 执行斯特鲁金斯基程序 (II B 节) 以获得平滑密度 (11)。然后，这些密度用于 LN 处理以获得占据振幅 $v_{\mu}$ 和非占据振幅 $u_{\mu}$。这完成了当前的 HF 迭代。该方法在本工作中被称为 SHF$_{mix}$ + LN 方法。

**B. 斯特鲁金斯基壳修正方法**

通过 SHF$_{mix}$ + LN 计算，可以确定超重核的基态性质。可以直接与实验数据进行比较，例如结合能、电荷半径、$\alpha$ 衰变能等。此外，还可以对尚未测量的超重核进行预测。然而，受量子效应影响的壳修正能，正如引言中提到的，在超重核的稳定性中起着关键作用。以前对变形超重核的研究依赖于检查 Nilsson 图来确定变形壳隙的位置。本工作旨在利用斯特鲁金斯基方法 [10] 提取这些变形核的壳修正能。据推测，这种定量研究为实验研究提供了更可靠的路线图。

本节描述了斯特鲁金斯基壳修正方法 [11] 的实现。使用 SHF 方法计算的原子核总能量可以分解为平滑部分和涨落部分，

$$
E_{sh} = \sum_{\mu,q=n,p} \epsilon_{\mu,q} \quad (5)
$$

$$
= \sum_{q=n,p} (\tilde{E}_{sh}^q + \delta E_{sh}^q), \quad (6)
$$

其中 $q$ 标记中子 ($q = n$) 和质子 ($q = p$) 量。由于两种核子的处理方式相同，本节的其余部分我们只讨论中子的情况，并省略 $q$ 指标。在 SHF 方法中，$E_{sh}$ (6) 与方程 (2) 中的总能量 $E_{tot}$ 相同。斯特鲁金斯基 [8,9] 开发了一种方法，用于使用单粒子谱 $\epsilon_{\mu}$ 提取能量的平滑部分 $\tilde{E}_{sh}$。首先，方程 (5) 中 $\epsilon_{\mu}$ 的求和可以用能级密度表示，

$$
E_{sh} = \int_{-\infty}^{\lambda} \epsilon g(\epsilon)d\epsilon, \quad (7)
$$

其中

$$
g(\epsilon) = \sum_{\mu} \delta(\epsilon - \epsilon_{\mu}) \quad (8)
$$

是给定单粒子谱 $\epsilon_{\mu}$ 的能级密度。方程 (7) 的上限 $\lambda$ 被选择为满足

$$
N = \int_{-\infty}^{\lambda} g(\epsilon)d\epsilon. \quad (9)
$$

由于束缚态的单粒子谱是离散的，费米能 $\lambda$ 可以选择在最后一个填充的单粒子能级和第一个未填充的单粒子能级之间。

能量的平滑部分定义为

$$
\tilde{E}_{sh} = \int_{-\infty}^{\tilde{\lambda}} \epsilon \tilde{g}(\epsilon)d\epsilon. \quad (10)
$$

平均单粒子能级密度通过折叠程序获得：

$$
\tilde{g}(\epsilon) = \frac{1}{\gamma_s} \int_{-\infty}^{+\infty} d\epsilon' g(\epsilon') f\left(\frac{\epsilon - \epsilon'}{\gamma_s}\right) \\
= \sum_{\mu} \frac{1}{\gamma_s} f\left(\frac{\epsilon - \epsilon_{\mu}}{\gamma_s}\right), \quad (11)
$$

其中 $\gamma_s$ 表示平滑宽度。折叠函数 $f(x)$ 表示为乘积：

$$
f(x) = w(x)P_{p_s}(x),
$$

其中权重函数 $w(x)$ 采用高斯形式：

$$
w(x) = \frac{1}{\sqrt{\pi}} e^{-x^2}.
$$

曲率修正多项式 $P_{p_s}(x)$ 由广义拉盖尔多项式定义

$$
P_{p_s}(x) = \sum_{n=0}^{p_s} C_n H_{2n}(x) = L_{p_s}^{(1/2)}(x^2),
$$

其中

$$
C_n = (-1)^n (2^{2n} n!)^{-1},
$$

其中 $H_{2n}(x)$ 是厄米多项式，$p_s$ 是平滑函数的阶数。对于接近滴线的原子核 [45-49]，需要从平滑能级密度 (11) 中减去连续谱的贡献。然而，对于本研究中考虑的超重核，连续谱的贡献可以忽略不计。在图 1 中，我们展示了超重核 $^{270}Hs$ 的 $\tilde{g}(\epsilon)$ 的典型行为。

054314-3

**物理评论 C 110, 054314 (2024)**

**[图 1]**：计算得到的 $^{270}Hs$ 的 $\tilde{g}(\epsilon)$ 作为 $\epsilon$ 的函数 ($\gamma_s = 1.44\hbar\omega_0, p_s = 8$)。使用了 UNEDF1 EDF。离散线是单粒子能量 $\epsilon_{\mu}$。
*该图展示了 $^{270}Hs$ 的平滑能级密度 $\tilde{g}(\epsilon)$ 随能量 $\epsilon$ 的变化。在超重核区域，由于中子数大和扁长变形，$\tilde{g}(\epsilon)$ 呈现连续增加的趋势。*

它表明，由于中子数大和扁长变形导致的高能级密度，$\tilde{g}$ 持续增加。

在方程 (10) 中，平滑费米能 $\tilde{\lambda}$ 由粒子数方程唯一确定

$$
N = \int_{-\infty}^{\tilde{\lambda}} \tilde{g}(\epsilon)d\epsilon \\
= \sum_{\mu} \left[ \frac{1}{2} + \frac{1}{2} \text{erf}(t_{\mu}) \right] e^{-t_{\mu}^2} \sum_{l=0}^{p_s} C_l H_{2l}(t_{\mu}), \quad (12)
$$

其中 $t_{\mu} = (\tilde{\lambda} - \epsilon_{\mu})/\gamma_s$，$\text{erf}(t_{\mu})$ 是误差函数：

$$
\text{erf}(x) = \frac{1}{\sqrt{\pi}} \int_{-x}^{x} e^{-t^2} dt.
$$

通过获得的 $\tilde{\lambda}$，总能量 (10) 的平滑部分计算如下：

$$
\tilde{E}_{sh} = \sum_{\mu} \frac{1}{2} \left[ 1 + \text{erf}(t_{\mu}) \right] - \frac{\gamma_s}{\sqrt{\pi}} e^{-t_{\mu}^2} \left[ \frac{1}{2} C_{p_s} H_{2p_s+1}(t_{\mu}) - \frac{\epsilon_{\mu} - \tilde{\lambda}}{\gamma_s} C_{p_s} H_{2p_s+1}(t_{\mu}) \right] \\
- \frac{2p_s \gamma_s}{\sqrt{\pi}} e^{-t_{\mu}^2} C_{p_s} H_{2p_s-1}(t_{\mu}).
$$

最后，壳修正能由下式给出

$$
\delta E_{sh}(\gamma_s, p_s) = E_{sh} - \tilde{E}_{sh}(\gamma_s, p_s). \quad (13)
$$

现在，对于给定的谱 $\epsilon_{\mu}$，壳修正能是平滑宽度 $\gamma_s$ 和曲率修正多项式阶数 $p_s$ 的函数。

为了确定两个参数 $\gamma_s$ 和 $p_s$ 的值，获得的壳修正能应满足平台条件

**物理评论 C 110, 054314 (2024)**

**[图 2]**：使用 SLy4 EDF 计算的 $^{270}Hs$ 的 $\delta E_{sh}$ 作为 $\gamma_s$ 和 $p_s$ 参数的函数。垂直线表示选择的 $\gamma_s = 1.44$ 值。
*该图展示了 $^{270}Hs$ 的壳修正能 $\delta E_{sh}$ 随平滑宽度 $\gamma_s$ 和曲率修正阶数 $p_s$ 的变化。当 $p_s > 4$ 时，$\delta E_{sh}$ 值在 $1.0 < \gamma_s < 2.0$ 范围内几乎保持不变，表明平台条件得到合理满足。*

$$
\frac{\partial (\delta E_{sh})}{\partial \gamma_s} = 0, \quad \frac{\partial (\delta E_{sh})}{\partial p_s} = 0. \quad (14)
$$

除了 Nilsson 势之外，这个条件通常不满足。我们选择平滑宽度为 $\gamma_s = 1.44$ (以振荡器频率 $\hbar\omega_0 = 41/A^{1/3} \text{ MeV}$ 为单位)，并在壳修正计算中选择 $p_s = 8$。

在图 2 中，我们绘制了 $^{270}Hs$ 的 $\delta E_{sh}$ 随 $\gamma_s$ 变化的曲线，对应于不同的 $p_s$ 值。可以看出，对于 $p_s > 4$，$\delta E_{sh}$ 值在 $1.0 < \gamma_s < 2.0$ 范围内几乎保持不变。具体来说，当 $\gamma_s = 1.44$ 时，对应于 $p_s = 6, 8, 10$ 的曲线重叠。这表明平台条件 (14) 对于 $\gamma_s = 1.44$ 和 $p_s = 8$ 得到了合理满足。选择的 $\gamma_s$ 和 $p_s$ 值与参考文献 [24] 中球形超重核计算所使用的值非常相似。

**C. 利普金-野上对力模型**

当前研究的主要目标是调查超重核中变形壳修正能。对力相互作用的纳入对于确定这些重核的变形至关重要。在本工作中，我们使用 LN 对力模型来描述对力相互作用 [50]。LN 模型相对于 BCS 模型的优势在于前者考虑了粒子数涨落效应。LN 模型增加了粒子数涨落常数 $\lambda_2$，这避免了 BCS 模型在计算中子数 $N = 142-152$ 的锕系核时出现的对力塌缩问题 [51]。

054314-4

**变形超重核的稳定性：一项壳修正分析**

对力强度 $G$ 使用 Möller-Nix 方案 [50] 计算：

$$
\frac{1}{G} = \frac{1}{2} \int_{y_1}^{y_2} \frac{dx}{\sqrt{x^2 + \tilde{\Delta}^2}} \\
= \frac{1}{2} \left[ \ln(y_2 + \sqrt{y_2^2 + \tilde{\Delta}^2}) - \ln(y_1 + \sqrt{y_1^2 + \tilde{\Delta}^2}) \right], \quad (15)
$$

其中 $\tilde{\rho}$ 是平均能级密度，定义为 $\tilde{\rho} = \tilde{g}(\tilde{\lambda})$。平滑能级密度 $\tilde{g}(\tilde{\lambda})$ 使用方程 (11) 计算，平滑费米能 $\tilde{\lambda}$ 使用方程 (12) 确定。

LN 模型中使用的有效相互作用对力隙 $\tilde{\Delta}$ 选择为 $\tilde{\Delta} = rB_s/N^{1/3} \text{ MeV}$，其中 $r = 3.3 \text{ MeV}$ 和 $B_s = 1$，如参考文献 [50] 中选择。方程 (15) 中积分限 $y_1$ 和 $y_2$ 定义为

$$
y_1 = \tilde{\lambda} - \frac{\tilde{N} + N_1 - 1}{\tilde{\rho}} \\
y_2 = \tilde{\lambda} - \frac{\tilde{N} + N_2}{\tilde{\rho}}
$$

参数 $N_1$ 和 $N_2$ 定义了相互作用区间。它始于费米能级以下 $N_1$ 处的能级，结束于费米能级以上 $N_2$ 处的能级。相互作用区间 $N_1$ 和 $N_2$ 定义为 [49]

$$
N_1 = \begin{cases} N - \tilde{g}(\tilde{\lambda})\Lambda & \text{for } N > \tilde{g}(\tilde{\lambda})\Lambda \\ \tilde{N} & \text{for } N \le \tilde{g}(\tilde{\lambda})\Lambda \end{cases},
$$

和

$$
N_2 = N + \tilde{g}(\tilde{\lambda})\Lambda,
$$

其中 $\Lambda = 1.2\hbar\omega_0$ 是截断参数。

下一步是迭代求解以下 LN 方程，以获得对力能隙 $\Delta$、费米能 $\lambda_{LN}$、粒子数涨落常数 $\lambda_2$、占据概率 $v_{\mu}$ 和改变的单粒子能量 $\tilde{\epsilon}_{\mu}$。LN 方程是 2(N$_2$ - N$_1$) + 5 个耦合非线性方程，由参考文献 [50] 给出

$$
N = 2 \sum_{\mu=N_1}^{N_2} v_{\mu}^2 + 2(N_1 - 1), \quad (16)
$$

$$
\frac{2}{G} = \sum_{\mu=N_1}^{N_2} \frac{1}{\sqrt{(\epsilon_{\mu} - \lambda_{LN})^2 + \Delta^2}}, \quad (17)
$$

$$
v_{\mu}^2 = \frac{1}{2} \left[ 1 - \frac{\epsilon_{\mu} - \lambda_{LN}}{\sqrt{(\epsilon_{\mu} - \lambda_{LN})^2 + \Delta^2}} \right], \\
\mu = N_1, N_1+1, \dots, N_2,
$$

$$
\tilde{\epsilon}_{\mu} = \epsilon_{\mu} + (4\lambda_2 - G)v_{\mu}^2, \\
\mu = N_1, N_1+1, \dots, N_2,
$$

$$
\Delta = G \sum_{\mu=N_1}^{N_2} u_{\mu}v_{\mu}, \quad (18)
$$

**物理评论 C 110, 054314 (2024)**

和

$$
\lambda_2 = \frac{G \left( \sum_{\mu=N_1}^{N_2} u_{\mu}^2 v_{\mu}^2 \right) \left( \sum_{\mu=N_1}^{N_2} u_{\mu}^2 v_{\mu}^2 \right) - \sum_{\mu=N_1}^{N_2} u_{\mu}^4 v_{\mu}^4}{\left( \sum_{\mu=N_1}^{N_2} u_{\mu}^2 v_{\mu}^2 \right)^2 - \sum_{\mu=N_1}^{N_2} u_{\mu}^4 v_{\mu}^4},
$$

其中

$$
u_{\mu}^2 = 1 - v_{\mu}^2, \quad \mu = N_1, N_1+1, \dots, N_2.
$$

在 LN 模型中，$\Delta + \lambda_2$ 被认为是奇偶质量差。

**III. 结果与讨论**

本节讨论我们使用 II 节中描述的 SHF$_{mix}$ + LN 程序计算的质子数 $Z \approx 98-120$ 和中子数 $N \approx 140-190$ 的偶偶核的结果。计算中使用了 SLy4 [52]、SkM* [53]、UNEDF1 [54] 和 UNEDF1$^{SO}$ [40] EDF。研究中包含锕系核 ($Z = 100-104$) 的动机是评估理论计算的性能。这些核的实验数据与计算结果的一致性为超锕核的预测提供了更多信心，因为超锕核的实验数据更为稀缺。

在 III A 节中，我们比较了计算得到的对力隙与实验数据。III B 节展示了平衡变形作为质子数和中子数函数的等高线图。在 III C 节中，我们展示了壳修正能的结果，并用它们预测超重核区域的变形壳闭合。此外，我们比较了不同 EDF 计算的壳修正能差异，并讨论了这些差异的根本原因。另一个用于分析壳结构的量是双粒子壳隙，在 III D 节中给出。我们还计算了几条同位素链的 $\alpha$ 衰变能 (III E 节) 和 $\alpha$ 衰变半衰期 (III F 节)，并与可用的实验数据进行比较。

**A. 对力隙**

由于费米面附近单粒子轨道密度高，超重核倾向于变形。因此，对力相互作用在描述这些重核的低能谱中起着关键作用。对力相互作用反映在奇偶质量交错中，这可以与计算得到的对力隙相关联 [56-58]。

本节中，我们比较了 SCMF 计算得到的对力隙与实验五点质量交错。中子五点质量交错由下式给出

$$
\Delta_n^{(5)}(N, Z) = \frac{(-1)^N}{8} [E(N+2, Z) \\
- 4E(N+1, Z) + 6E(N, Z) \\
- 4E(N-1, Z) + E(N-2, Z)], \quad (19)
$$

054314-5

**物理评论 C 110, 054314 (2024)**

**[图 3]**：计算得到的 (a) No ($Z = 102$) 同位素和 (b) $N = 150$ 同中子素的 (a) 中子和 (b) 质子对力隙 (18)。可用的 (a) 中子和 (b) 质子五点质量交错数据也包含在图中以供比较 [55]。
*该图展示了对力隙随中子数和质子数的变化。计算结果与实验五点质量交错数据在量级和趋势上基本一致，表明 LN 对力模型适用于预测超重核的对力性质。*

其中 $E(N, Z)$ 是质量数为 $A = N + Z$ 的原子核的质量。对于质子，$\Delta_p^{(5)}(N, Z)$ 由下式给出

$$
\Delta_p^{(5)}(N, Z) = \frac{(-1)^Z}{8} [E(N, Z+2) \\
- 4E(N, Z+1) + 6E(N, Z) \\
- 4E(N, Z-1) + E(N, Z-2)]. \quad (20)
$$

五点质量交错可以与量 $\Delta + \lambda_2$ (参见 II C 节) 进行比较。

图 3 显示了 No ($Z = 102$) 同位素的中子和 $N = 150$ 同中子素的质子计算对力隙 ($\Delta + \lambda_2$)。实验五点质量交错 ($\Delta^{(5)}$) 也包含在图中以供比较。可以看出，计算得到的对力隙 $\Delta + \lambda_2$ 的量级再现了 No ($Z = 102$) 同位素的 $\Delta^{(5)}$。总的来说，除了 $Z = 102$ 的情况外，对于 $N = 150$ 同中子素的 $\Delta^{(5)}$ 也得到了再现，其中 $\Delta^{(5)}$ 值被高估了约 $150 \text{ keV}$。特别是，在下图中，我们可以看到我们的计算再现了 $\Delta^{(5)}$ 值的下降趋势。这种一致性表明我们工作中使用的 LN 对力模型 (II C 节) 适用于预测超重核的对力性质。

**B. 平衡变形**

本节展示了本工作中考虑的锕系和超锕核的平衡变形。由于大的表面能项和库仑能项对变形的敏感性，结合能通常受核变形的影响。除了降低核变形能力的对力相互作用外，核变形还与中子-质子比、激发能和自旋值等各种因素有关 [59]。与许多早期对超重核壳结构进行球形限制的自洽平均场研究不同，本研究进行了无对称性限制的非相对论平均场计算。

图 4 显示了由下式定义的平衡四极变形

$$
\beta_2 = \frac{4\pi}{3r_0 A^{5/3}} Q_{20}, \quad (21)
$$

其中

$$
Q_{20} = \sum_{q=n,p} \int \rho_q(\mathbf{r}) (2z^2 - x^2 - y^2) d^3r.
$$

我们还使用 III A 节中使用的其他参数 (SLy4, SkM*, 和 UNEDF1) 计算了基态变形。结果与 UNEDF1$^{SO}$ 的结果相似，因此我们在此不展示它们。

如图 4 所示，质子数 $Z < 114$ 和中子数 $N \le 180$ 的大多数超重核基态具有扁长变形。$\beta_2$ 值随中子数的增加而逐渐减小，从 $N = 138$ 处的 $\beta_2 \approx 0.3$ 到 $N = 180$ 处的近球形变形 $\beta_2 \approx 0.1$，因为中子数接近球形幻数 $N \approx 184$。对于 No ($Z = 102$) 同位素，我们看到 $\beta_2 = 0.3$ 的等高线延伸到 $N = 154$ 更远。这可能是由于核 $^{254}No$ 在 $\beta_2 = 0.3$ 处打开了一个大的壳隙。值得注意的是，在预测的壳闭合附近，变形没有显著变化，壳闭合在 $Z = 108$ 和 $N = 162$ 处打开，如 III C 节所述。

对于缺中子同位素 ($N \le 180$)，质子数 $Z = 116, 118$ 和 $120$ 的超重核基态由于势能曲线中扁长极小值的降低而变为扁长形。扁长变形在 $Z = 120$ 和 $N \le 168$ 处达到最大值，$\beta_2 = -0.45$。对于 $N \ge 180$ 的同位素，基态变形接近近球形。对于 Lv ($Z = 116$) 同位素，$\beta_2 = 0.3$ 的等高线延伸到 $N = 168$。在 III C 节中，我们展示了壳闭合在 $Z = 116$ 和 $N = 168$ 处打开。

参考文献 [28] 中回顾了使用各种平均场模型计算的超重核平衡变形。使用 SkM* 和 SLy4 EDF 获得的结果与我们的结果相似。实验上，核 $^{254}No_{152}$ 附近的基态变形已被证实约为 $\beta_2 \approx 0.3$ [38,60]。

图 5(a) 和 5(b) 显示了总能量，图 5(c) 和 5(d) 显示了总壳修正能作为四极变形 $\beta_2$ 的函数，对应于 $^{256}Rf$ 和 $^{270}Hs$。

054314-6

**物理评论 C 110, 054314 (2024)**

**[图 4]**：使用 UNEDF1$^{SO}$ EDF 计算的偶偶超重核基态四极变形 $\beta_2$。两个黑点表示 $^{256}Rf$ 和 $^{280}Hs$ 核，其势能曲线绘制在图 5 中。
*该图展示了超重核的基态四极变形 $\beta_2$ 随质子数 $Z$ 和中子数 $N$ 的变化。对于 $Z < 114$ 和 $N \le 180$ 的核，主要表现为扁长变形，而对于 $Z \ge 116$ 的核，则表现为扁球变形。变形值随核子数的变化而平滑过渡。*

**[图 5]**：使用 UNEDF1 和 UNEDF1$^{SO}$ EDF 计算的 $^{256}Rf$ 和 $^{270}Hs$ (在图 4 中标记) 的总能量和总壳修正能 (以 MeV 为单位) 作为四极变形 $\beta_2$ 的函数。
*该图展示了 $^{256}Rf$ 和 $^{270}Hs$ 的总能量和壳修正能随变形 $\beta_2$ 的变化。对于这些核，扁长极小值在能量上低于扁球极小值。总能量和壳修正能的极小值出现在相似的变形处，表明壳修正对稳定性的重要性。*

这两个核位于或接近四种 EDF 计算的总壳修正能等高线图中的极小值。

图 5(a) 和 5(b) 显示了扁长和扁球两侧两个极小值的共存。对于这些核，扁长极小值在能量上低于扁球极小值。从图 4 可以推断，随着质子数的增加，在 $Z \approx 116$ (对于 $N \le 180$) 处从扁长变形到扁球变形的转变。在图 5(c) 和 5(d) 中，总壳修正能达到接近总能量曲线极小值的变形处的极小值。从 MM 模型的角度来看，这是预期的，因为变形解的形成归因于给定变形处 Nilsson 能级的稀疏性。较低的单粒子能级密度导致较低的壳修正能。

特别是，反射不对称变形可能在评估超重核的稳定性中发挥作用 [61]。例如，八极变形似乎存在于某些 Th 同位素的基态中 [62]。然而，在本工作中考虑的原子核中，我们的无对称性限制计算中没有出现反射不对称变形。

**C. 壳修正能**

我们将壳修正能视为稳定性的指标。具体来说，具有较大负壳修正能的原子核应该比具有较小壳修正能的原子核更稳定。这尤其适用于超锕核，其中壳修正能提供了额外的结合能。本节试图通过研究壳修正能的等高线图来阐明这些变形超锕核的稳定性。实验学家似乎也采用了类似的壳修正能态度 [41]。

图 6-9 显示了偶偶核 $Z = 98-120$ 和 $N = 140-190$ 的总壳修正能 $\delta E_{sh}^{n+p} = \delta E_{sh}^n + \delta E_{sh}^p$ 的结果。计算使用 SLy4、SkM*、UNEDF1 和 UNEDF1$^{SO}$ EDF。我们首先讨论对 $^{254}No$ 附近核的预测，这些核有更多的实验数据。之后，我们讨论实验上知之甚少的更重核。

在图 6 中，SLy4 EDF 预测在 $^{256}Rf$ 处存在一个极小值 (约 $-10.5 \text{ MeV}$)，对应于最强的稳定性。如图 4 所示，$^{256}Rf$ 附近的核具有强烈的扁长变形，$\beta_2 \approx 0.30$。SkM* 和 UNEDF1 EDF 在锘区域给出相同的壳，如图 7 和 8 所示。使用 SkM* 和 UNEDF1 EDF 计算的极小值比 SLy4 EDF 的极小值浅。在图 9 中，

054314-7

**物理评论 C 110, 054314 (2024)**

**[图 6]**：使用 SLy4 EDF 计算的偶偶超重核的总壳修正能 (以 MeV 为单位) 作为中子数和质子数的函数。图中标记了实验上可用的偶偶 $\alpha$ 衰变链。$\alpha$ 衰变半衰期取自参考文献 [63]。
*该图展示了 SLy4 EDF 预测的壳修正能等高线图。图中显示在 $^{270}Hs_{162}$ 处存在一个深极小值，表明该核具有较高的稳定性。实验数据（$\alpha$ 衰变半衰期）也在此图中标记，与理论预测的稳定性区域一致。*

我们看到 UNEDF1$^{SO}$ EDF 在 $^{252}Fm$ 处给出了一个壳闭合，其 ($Z, N$) = (100, 152)。

实验上，$^{254}No_{152}$ 附近的核的能谱信息相对丰富 [38]。奇 A 体系的实验准粒子能量表明，变形壳隙在 ($Z, N$) = (100, 152) 处打开。这得到了 Fm 同位素和 $^{254}No$ 具有比其邻近核更长半衰期的事实的支持。此外，No 和 Fm 区域集体转动和振动态的系统调查似乎支持了 Woods-Saxon 势预测的单粒子结构，其中变形子壳在 ($Z, N$) = (100, 152) 处打开 [64,65]。

UNEDF1$^{SO}$ EDF 成功地再现了 $^{252}Fm$ 附近的壳闭合。这是预期的，因为该 EDF 是通过稍微增强原始 UNEDF1 EDF [54] 的自旋-轨道项以拟合附近奇 A 核 $^{251}Cf$ 和 $^{249}Bk$ 的实验准粒子能量而获得的。这些壳闭合的出现反映在 $\delta E_{sh}^{n+p}$ 等高线图中。有趣的是，看看这样一套微调的参数集是否能对更重超重核的基态性质做出合理预测。

**[图 7]**：与图 6 相同，但计算使用 SkM* EDF。
*该图展示了 SkM* EDF 预测的壳修正能等高线图。与 SLy4 EDF 类似，SkM* EDF 也在 $^{270}Hs_{162}$ 处显示出深极小值，表明该核的稳定性。*

054314-8

**物理评论 C 110, 054314 (2024)**

**[图 8]**：与图 6 相同，但计算使用 UNEDF1 EDF。
*该图展示了 UNEDF1 EDF 预测的壳修正能等高线图。图中显示在 $^{270}Hs_{162}$ 处存在一个深极小值，表明该核具有较高的稳定性。*

在超重核区域，壳修正能的等高线图显示，对于所有四种 EDF (参见图 6-9)，在 $^{108}Hs_{162}$ 处存在更深的极小值。这种情况与球形计算不同，球形计算中壳闭合在不同核子数处打开，如引言中所述。与 $^{270}Hs$ 处的极小值相比，$^{254}No$ (或 $^{250}Fm$) 附近的极小值浅约 $2-3 \text{ MeV}$。使用 Woods-Saxon 势的 MM 计算预测了相似的等高线结构，即在 $^{270}Hs$ 和 $^{252}Fm$ 核处存在两个极小值 [66]。然而，$^{252}Fm$ 处极小值的深度与 $^{270}Hs$ 附近的深度相似 [41]，这与我们的 SHF$_{mix}$ + LN 计算不同。

实验上，$^{270}Hs$ 的衰变性质已在 $^{26}Mg + ^{248}Cm$ [67] 和 $^{48}Ca + ^{226}Ra$ [68] 反应中进行了研究。从这些研究中获得的 $^{270}Hs$ 的产生截面、$\alpha$ 衰变能、$\alpha$ 衰变半衰期和自发裂变半衰期证实了该核的更高稳定性，这与 MM 和当前计算结果一致。图 6 中标记了实验上已知的偶偶 $\alpha$ 衰变核及其 $\alpha$ 衰变半衰期。Fm ($Z = 100$) 和 Hs ($Z = 108$) 同位素链的更长半衰期以及更多实验数据的可用性表明，这些核比其他同位素更稳定。

**[图 9]**：与图 6 相同，但计算使用 UNEDF1$^{SO}$ EDF。图中标记了包含 Cn 同位素的实验上可用的奇偶 $\alpha$ 衰变链。$\alpha$ 衰变半衰期取自参考文献 [63]。
*该图展示了 UNEDF1$^{SO}$ EDF 预测的壳修正能等高线图。图中显示在 $^{270}Hs_{162}$ 处存在一个深极小值，并且在 $^{276}Cn_{164}$ 处也存在一个极小值。实验数据（$\alpha$ 衰变半衰期）也在此图中标记，与理论预测的稳定性区域一致。*

054314-9

**物理评论 C 110, 054314 (2024)**

**[图 10]**：计算得到的壳修正能，对应于 (a) Rf ($Z = 104$) 同位素 [红色方块表示使用 UNEDF1$^{SO}$ EDF 计算的 Fm ($Z = 100$) 同位素]，(b) Hs ($Z = 108$) 同位素，(c) Lv ($Z = 116$) 同位素，(d) $N = 152$ 同中子素，(e) $N = 162$ 同中子素，以及 (f) $N = 168$ 同中子素。
*该图展示了不同 EDFs 预测的壳修正能随中子数和质子数的变化。对于 Rf 和 Hs 同位素，在 $N = 152$ 和 $N = 162$ 处存在极小值。对于 $N = 162$ 同中子素，在 $Z = 108$ 处存在显著极小值。这些极小值表明了这些核子数处的壳闭合效应。*

高于 $^{270}Hs$，我们的计算表明，下一个变形壳闭合出现在 SLy4 (图 6)、SkM* (图 7) 和 UNEDF1 (图 8) EDF 的 $^{284}Lv$ 处。对于 UNEDF1$^{SO}$ (图 9)，极小值出现在 ($Z, N$) = (112, 164) 处。对于 SLy4，近球形极小值出现在 ($Z, N$) = (118, 174) 处。如图 6 所示，这个极小值与已知最重 $\alpha$ 衰变链 ($Z > 112$) 的极小值非常接近。等高线图中的这些极小值深度与 $^{270}Hs$ 处的极小值深度相当。

实验上，Cn ($Z = 112$) 同位素链没有偶偶 $\alpha$ 衰变数据。在图 9 中，我们展示了实验上可用的奇 A 核的 $\alpha$ 衰变链和半衰期。有趣的是，起始于 Cn 同位素的 $\alpha$ 衰变链，正如 UNEDF1$^{SO}$ EDF 预测的，位于壳修正等高线图的谷附近，如图 9 所示。这条链中奇 A 核的长半衰期似乎证实了 UNEDF1$^{SO}$ EDF 预测的增强稳定性。必须指出的是，对这些核稳定性的定量理解需要对竞争衰变模式（如 $\alpha$ 衰变 [69] 和自发裂变 [70]）的半衰期进行微观计算。然而，这些计算超出了当前研究的范围。

为了更深入地理解极小值的形成，我们绘制了 Rf [图 10(a)]、Hs [图 10(b)] 和 Lv [图 10(c)] 同位素的中子壳修正能。我们看到，对于 Rf 同位素，各种 EDF 在 $N = 152$ 和 $162$ 处显示极小值，其中 $N = 162$ 处的极小值深约 $1 \text{ MeV}$。UNEDF1 和 UNEDF1$^{SO}$ EDF 的绝对壳修正能小于 SLy4 和 SkM* EDF 计算的能量。这可能是由于 UNEDF1 和 UNEDF1$^{SO}$ EDF 在费米面附近具有更高的能级密度。Hs 同位素在 $N = 152, 162$ 处也出现类似的极小值。

图 10 的下半部分显示了 $N = 152$ [图 10(d)]、$N = 162$ [图 10(e)] 和 $N = 168$ [图 10(f)] 同中子素的质子壳修正能。对于 $N = 152$ 同中子素，SLy4 和 UNEDF1 EDF 在 $Z = 104$ 处显示浅极小值。UNEDF1$^{SO}$ EDF 在 Fm ($Z = 100$) 同位素处显示一个凹陷。$\delta E_{sh}$ 值在 $Z = 100-108$ 范围内几乎保持不变，并在 $Z > 108$ 时迅速上升。对于 $N = 162$ 同中子素，四种 EDF 显示在 $Z = 108$ 处具有显著极小值的结果。对于 $N = 168$ 同中子素，四种 EDF 的结果显示在 $Z = 108, 116$ 处存在极小值。此外，UNEDF1$^{SO}$ EDF 也在 $Z = 112$ 处具有极小值。

图 10 中极小值的出现可能与费米面附近较低的单粒子能级密度有关。例如，在图 11 中，我们绘制了使用 SLy4 EDF 计算的 $Z = 108$ 同位素链和 $N = 162$ 同中子素链的单粒子能级。

与 $^{256}Rf$ 和 $^{270}Hs$ 附近的壳闭合相比，$^{284}Lv$ 附近极小值形成的原因不太清楚 [参见图 10(c)]。对于 SLy4 EDF，$^{284}Lv$ 处的极小值是质子和中子壳闭合在 $Z = 116$ 和 $N = 168$ 处共同贡献的结果。对于 SkM* 和 UNEDF1 EDF，在 ($Z, N$) = (116, 168) 处出现次极小值主要是由于质子壳闭合，而中子壳修正能 $\delta E_{sh}^n$ 在 $N = 162-184$ 范围内几乎保持不变。

图 11(a) 显示了使用 SLy4 EDF 计算的 $Z = 108$ 同位素的中子单粒子轨道。可以看出，随着中子数的增加，费米面达到 $N = 150, 152$ 和 $162$ 处的壳闭合。这与图 10(b) 一致，其中中子壳修正能显示在 $N = 150, 152$ 处存在极小值，最低点在 $N = 162$ 处。

图 11(b) 显示了使用 SLy4 EDF 计算的 $N = 162$ 同中子素的质子单粒子轨道。在 $Z = 104, 108$ 和 $116$ 处出现三个壳闭合。对于图 10(d)，壳闭合仅在 $N = 162$ 的 $Z = 108$ 处可见。$Z = 104$ 壳闭合的影响反映在 $Z = 104$ 处壳修正能斜率的变化中。

054314-10

**物理评论 C 110, 054314 (2024)**

**[图 11]**：使用 SLy4 EDF 计算的 (a) $Z = 108$ 同位素的中子单粒子能量，以及 (b) $N = 162$ 同中子素的质子单粒子能量。
*该图展示了单粒子能级随中子数和质子数的变化。图中清晰可见在特定核子数处（如 $N=150, 152, 162$ 和 $Z=104, 108, 116$）存在能隙，这些能隙对应于壳闭合效应。*

在讨论了 $Z < 120$ 的超重核变形壳闭合后，我们尝试定位 $Z \ge 120$ 的球形壳闭合。由于参考文献 [25] 中已对球形计算进行了更广泛的搜索，我们仅限于有限的计算。具体来说，在图 12 中，我们展示了 $Z = 124$ 同位素的中子壳修正能 $\delta E_{sh}^n$ 和 $N = 184$ 同中子素的质子壳修正能 $\delta E_{sh}^p$。

对于 SLy4 EDF，在 $N = 184$ 和 $Z = 124, 126$ 处观察到深极小值，表明这些核子数处存在壳闭合。SLy4 的这些计算结果与参考文献 [23-25] 中报道的结果一致。然而，对于

**[图 12]**：计算得到的 (a) $Z = 124$ 同位素的 $\delta E_{sh}^n$ (灰色线表示扁长变形下的结果) 和 (b) $N = 184$ 同中子素的 $\delta E_{sh}^p$。
*该图展示了壳修正能随中子数和质子数的变化。对于 SLy4 EDF，在 $N=184$ 和 $Z=124, 126$ 处存在深极小值，表明这些核子数处存在壳闭合。而 UNEDF1$^{SO}$ EDF 的极小值较浅，且在 $N$ 介于 114 到 124 之间以及 $Z$ 介于 178 到 184 之间呈现平坦底部。*

UNEDF1$^{SO}$ EDF，极小值变得更浅。如图 12 所示，壳修正能显示一个平坦的底部，其中 $N$ 范围从 114 到 124， $Z$ 范围从 178 到 184。以等高线图的形式，可以预期一个核区域将共享相似的低壳修正能。参考文献 [25] 中已注意到 $Z = 82$ 和 $N = 126$ 之外的幻数模式缺乏。

**D. 经验壳隙**

另一个可用于分析单粒子结构的量是经验壳隙参数 $\delta_{2n,2p}(N, Z)$，其定义为

$$
\delta_{2n}(N, Z) = S_{2n}(N, Z) - S_{2n}(N+2, Z) \\
= -B(N-2, Z) + 2B(N, Z) - B(N+2, Z),
$$

$$
\delta_{2p}(N, Z) = S_{2p}(N, Z) - S_{2p}(N, Z+2) \\
= -B(N, Z-2) + 2B(N, Z) - B(N, Z+2),
$$

其中 $S_{2n,2p}(N, Z)$ 是双中子和双质子分离能，$B(N, Z)$ 是结合能。量 $\delta_{2n/2p}(N, Z)$ 与结合能作为中子或质子数的二阶导数相关，是单粒子能级密度的敏感指标。这个量放大了弱壳效应，被称为双中子或双质子壳隙 [32,39]。

将双粒子壳隙仅与单粒子轨道间距关联起来有些误导，因为单粒子能量随不同核的变形而变化。因此，$\delta_{2n/2p}$ 值包含了变形效应。然而，对于计算中使用的 EDF，核的变形随核子数平滑变化 (参见图 4)。因此，双粒子壳隙在某个核子数处的突然增加可以被视为该核子数处单粒子壳闭合出现的指标。

图 13 的上半部分显示了使用几种 EDF 计算的 No ($Z = 102$) 和 Hs ($Z = 108$) 同位素链的 $\delta_{2n}$ 值。图 13 的下半部分显示了 $N = 152$ 和 $N = 162$ 同中子素链的 $\delta_{2p}$ 值。可用数据是根据参考文献 [55] 中的实验结合能计算的。实验数据显示在 $N = 152, 162$ 处存在大的双中子壳隙，在 $Z = 100, 108$ 处存在双质子壳隙。

尽管计算结果与实验值在细节上有所不同，但它们合理地再现了峰值。对于 No ($Z = 102$) 同位素 [图 13(a)]，所有四种参数在 $N = 152$ 处都出现峰值。对于 SLy4 EDF，在 $N = 150$ 处出现一个额外的峰值。这是由于使用 SLy4 EDF 计算的单粒子能级在 $N = 150$ 处存在大的能隙。对于 $N = 152$ 同中子素 [图 13(c)]，除了 UNEDF1 EDF 外，所有参数在 $Z = 100$ 处都出现峰值。UNEDF1$^{SO}$ 如预期般再现了数据的总体峰值结构。

对于 Hs ($Z = 108$) 同位素 [图 13(b)]，所有四种参数在 $N = 152$ 和 $162$ 处都出现强峰值。这可以归因于这些中子数处的壳闭合，如图 11(a) 中的单粒子能级所示。对于 $N = 162$

054314-11

**物理评论 C 110, 054314 (2024)**

**[图 13]**：计算得到的 $\delta_{2n}$ 值作为中子数函数，以及 $\delta_{2p}$ 值作为质子数函数，对应于 No ($Z = 102$) 和 Hs ($Z = 108$) 同位素链，以及 $N = 152$ 和 $N = 162$ 同中子素链。计算使用了四种 EDF。五角星表示实验值，叉号表示根据相邻核的系统趋势估计的值。实验和估计值均取自参考文献 [55]。
*该图展示了双中子分离能 $\delta_{2n}$ 和双质子分离能 $\delta_{2p}$ 随核子数的变化。图中显示在 $N=152, 162$ 和 $Z=100, 108$ 处存在明显的峰值，这些峰值与实验数据基本一致，表明这些核子数处存在壳闭合。*

同中子素 [图 13(d)]，所有 EDF 都预测在 $Z = 108$ 处存在闭合。UNEDF1$^{SO}$ 预测在 $Z = 100$ 处存在一个壳闭合，这与实验数据一致。UNEDF1$^{SO}$ 的 $\delta_{2n}$ 值在 $N = 152$ 处比 UNEDF1 EDF 的值更深。这导致 UNEDF1$^{SO}$ 在 $^{270}Hs$ 核东北方向的壳修正能始终低于 UNEDF1 EDF，参见图 9。

**E. $\alpha$ 衰变能**

$\alpha$ 衰变能 $Q_{\alpha}$ 是超重核单粒子结构的敏感指标。它定义为

$$
Q_{\alpha}(N, Z) = E(N, Z) - E(N-2, Z-2) - E(2, 2), \quad (22)
$$

其中 $E(N, Z)$ 是质量数为 $A = N + Z$ 的原子核的结合能，$E(2, 2) = -28.5 \text{ MeV}$ [71] 是 $\alpha$ 粒子的结合能。一个重要的观测结果是 $Q_{\alpha}$ 值可以用于识别核壳结构。核壳结构反映在 $Q_{\alpha}$ 值随中子数或质子数的变化中，这在参考文献 [72] 中有讨论。

图 14 显示了 No ($Z = 102$) 到 Og ($Z = 120$) 同位素链中偶偶核的 $Q_{\alpha}$ 值，计算使用了四种 EDF。为了比较，实验数据和估计值也包含在图中。实验数据取自参考文献 [55]。可以看出，我们的计算结果与实验数据在量级和峰值位置上基本一致。对于 SLy4 EDF，在 $N = 162$ 处存在一个深谷，这可能是由于 SLy4 EDF 计算的单粒子能级在 $N = 162$ 处存在大的能隙。UNEDF1 和 UNEDF1$^{SO}$ EDF 的计算结果低估了 $N = 162$ 处谷的深度，这可能是由于它们在费米面附近具有更大的单粒子能级密度。

图 15 显示了包含 $^{270}Hs$ 的 $\alpha$ 衰变链的单粒子能级。图 15(a) 和 15(b) 显示了使用 UNEDF1 EDF 计算的 $^{270}Hs$ 的中子和质子单粒子能级。图 15(c) 和 15(d) 显示了使用 UNEDF1$^{SO}$ EDF 计算的 $^{270}Hs$ 的中子和质子单粒子能级。

一般来说，随着中子数的增加，超重核的 $Q_{\alpha}$ 值会下降。然而，在 $N = 162$ 处，所有 EDF 都预测 $Q_{\alpha}$ 值会突然增加，这表明在 $N = 162$ 处存在一个壳闭合。对于 $Z = 108$ 同位素，SLy4 EDF 预测在 $N = 162$ 处存在一个深谷，这与图 10(b) 中壳修正能的极小值一致。对于 $Z = 116$ 同位素，SLy4 EDF 预测在 $N = 168$ 处存在一个深谷，这与图 10(c) 中壳修正能的极小值一致。

**物理评论 C 110, 054314 (2024)**

**[图 14]**：计算得到的 No ($Z = 102$) 到 Og ($Z = 120$) 同位素链中偶偶核的 $Q_{\alpha}$ 值。计算使用了四种 EDF。星号表示实验 $\alpha$ 衰变能值，叉号表示根据相邻核的系统趋势估计的值。实验和估计值均取自参考文献 [55]。
*该图展示了 $\alpha$ 衰变能 $Q_{\alpha}$ 随中子数 $N$ 的变化。图中显示 $Q_{\alpha}$ 值通常随 $N$ 的增加而下降，但在 $N=162$ 处出现一个明显的上升，表明在该中子数处存在壳闭合。不同 EDFs 的计算结果与实验数据在趋势上基本一致。*

$\alpha$ 衰变能的超重核单粒子结构效应可以系统地再现。一般来说，随着中子数的增加，$\alpha$ 衰变能 $Q_{\alpha}$ 值会下降 [73]。

更相似的是，实验单粒子能级在 Skyrme EDF 的 $Z = 108$ 和 $N = 162$ 处显示出超重核的壳效应。我们可以系统地再现超重核的 $\alpha$ 衰变能。

一般来说，随着中子数的增加，预测的 $Q_{\alpha}$ 值会下降。然而，在某些核子数处，例如 $N = 150, 152, 162, 168, 174, 184$，我们观察到 $Q_{\alpha}$ 值发生突然变化。这些变化是由于这些核子数处壳闭合引起的。图 14 中 $N = 162$ 处的壳闭合与 SLy4 EDF 的大能隙有关。对于 $N = 150$ 处的壳闭合，SLy4 EDF 计算的单粒子能级中存在大的能隙，如图 11(a) 所示。对于 UNEDF1 和 UNEDF1$^{SO}$ EDF，这些能隙较浅，且在 $N = 152$ 处比 SLy4 EDF 的能隙浅。这可能是由于 UNEDF1 和 UNEDF1$^{SO}$ EDF 在费米面附近具有更大的单粒子能级密度。

图 14 中显示了 UNEDF1 和 UNEDF1$^{SO}$ EDF 的 $Q_{\alpha}$ 值。UNEDF1$^{SO}$ EDF 的 $Q_{\alpha}$ 值在 $N = 162$ 处比 UNEDF1 EDF 的值更深。这导致 UNEDF1$^{SO}$ 在 $^{270}Hs$ 核东北方向的壳修正能始终低于 UNEDF1 EDF，参见图 9。

对于实验 $\alpha$ 衰变链，起始于 $^{270}Ds$ 和 $^{264}Hs$ [图 16(b) 和 16(c)]， $Q_{\alpha}$ 值随中子数的减少而平滑下降。除了 SkM* EDF，计算结果很好地再现了数据。尽管壳修正强烈影响 $Q_{\alpha}$ 值，但壳效应的大小与 $\alpha$ 衰变能之间没有简单的关联。这些 $\alpha$ 衰变链被观测到锘核的事实可能归因于相关核位于壳修正等高线图的谷中，如图 6 所示。

为了量化壳修正能对 $\alpha$ 衰变能的贡献，我们定义了两个新量，$\delta Q_{\alpha}$ 和 $\delta Q_{\alpha}^{Expt}$。$\delta Q_{\alpha}$ 定义如下：

$$
\delta Q_{\alpha} = Q_{\alpha}(N, Z) - E_{Mac}(N, Z) \\
+ E_{Mac}(N-2, Z-2) + E(2, 2). \quad (23)
$$

宏观能量 $E_{Mac}(N, Z)$ 由下式给出

$$
E_{Mac}(N, Z) = E(N, Z) - \delta E_{sh}^{n+p}(N, Z), \quad (24)
$$

其中 $E(N, Z)$ 是使用方程 (2) 计算的总能量。壳修正能 $\delta E_{sh}^{n+p}$ 是使用方程 (13) 计算的。“$n+p$”表示质子和中子贡献的总和。

将方程 (22) 和 (24) 代入方程 (23)，我们得到

$$
\delta Q_{\alpha} = \delta E_{sh}^{n+p}(N, Z) - \delta E_{sh}^{n+p}(N-2, Z-2). \quad (25)
$$

从方程 (25) 可以清楚地看出，量 $\delta Q_{\alpha}$ 是 $\alpha$ 衰变链中两个相邻核之间 $\delta E_{sh}^{n+p}$ 的差值。

054314-14

**物理评论 C 110, 054314 (2024)**

**[图 16]**：计算得到的 $^{282}Fl$、$^{278}Fl$ 和 $^{276}Fl$ 核的 $\alpha$ 衰变链的 $Q_{\alpha}$ 值。计算使用了几种 EDF。红色星号表示实验 $\alpha$ 衰变能值，叉号表示根据相邻核的系统趋势估计的值。实验和估计值均取自参考文献 [55]。
*该图展示了 $\alpha$ 衰变能 $Q_{\alpha}$ 随质量数 $A$ 的变化。图中显示 $Q_{\alpha}$ 值通常随 $A$ 的增加而下降，但在特定质量数处出现局部极小值或极大值，反映了壳效应。计算结果与实验数据基本一致。*

其中

$$
\delta Q_{\alpha}^{Expt} = Q_{\alpha}^{Expt}(N, Z) - E_{Mac}(N, Z) \\
+ E_{Mac}(N-2, Z-2) + E(2, 2). \quad (26)
$$

类似地，通过从实验质量中减去理论 $\tilde{E}_{sh}$，可以获得实验壳修正能等类似量 [11]。

两个量 $\delta Q_{\alpha}$ 和 $\delta Q_{\alpha}^{Expt}$ 对应于单粒子能级不均匀性引起的涨落贡献。通过减去宏观能量，不同同位素链的 $\delta Q_{\alpha}$ 和 $\delta Q_{\alpha}^{Expt}$ 值使我们能够更清楚地将 $Q_{\alpha}$ 图 (图 15) 中的扭结归因于变形壳闭合，与 $Q_{\alpha}$ 对 $N$ 的图 (图 14) 相比。

图 17 显示了使用 SLy4 (左图) 和 UNEDF1$^{SO}$ (右图) EDF 计算的 No 到 $Z = 120$ 同位素偶偶核的 $\delta Q_{\alpha}$ 值。用于计算 $\delta Q_{\alpha}^{Expt}$ 的实验 $Q_{\alpha}$ 值取自参考文献 [55]。用于计算这些图中 $\delta Q_{\alpha}$ 和 $\delta Q_{\alpha}^{Expt}$ 值的 $E_{Mac}$ 和 $\delta E_{sh}^{n+p}$ 值在附录 A 中给出。新的 $\delta Q_{\alpha}^{Expt}$ 值可以通过从新测量的 $Q_{\alpha}$ 值中减去表 I 中提供的 $E_{Mac}$ 来获得。

$\delta Q_{\alpha}$ 值作为中子数函数的曲线极小值比 $Q_{\alpha}$ 曲线更清晰。这表明 $\delta Q_{\alpha}$ 确实是比 $Q_{\alpha}$ 值更敏感的壳隙指标。在图 17 左图中使用 SLy4 EDF，我们看到 No、Rf 和 Sg 同位素在 $N = 150, 152$ 处存在浅极小值。所有同位素链在 $N = 162$ 处显示深极小值。除了 No 和 Rf 同位素，其他同位素在 $N = 168$ 处显示极小值。对于 Cn、Fl、Lv 和 Og 同位素，在 $N = 174$ 处可以看到极小值。实验 $\delta Q_{\alpha}^{Expt}$ 数据得到了相当好的再现。

图 17 的右图显示了使用 UNEDF1$^{SO}$ EDF 计算的 $\delta Q_{\alpha}$ 值作为中子数的函数。将这些图与图 14(d) 比较，我们看到了绘制 $\delta Q_{\alpha}$ 而不是 $Q_{\alpha}$ 的优势，对于不同的同位素链。事实上，图 14(d) 中斜率的细微变化在图 17(g)-17(l) 中变成了极小值。

对于较轻的同位素 (No、Rf 和 Sg)，可以看到 $N \approx 152$ 处的极小值。极小值的出现显示出同位素依赖性。实验 $\delta Q_{\alpha}^{Expt}$ 值显示，对于较重的同位素 (Hs、Ds、Cn、Fl、Lv 和 Og)，在 $N = 162$ 处存在强极小值。对于 Og 和 $Z = 120$ 同位素，在 $N = 184$ 处可以看到极小值。这些极小值被 UNEDF1$^{SO}$ 的计算低估了。

**F. $\alpha$ 衰变半衰期**

$\alpha$ 衰变半衰期可以使用以下 Viola-Seaborg 型公式 [74] 计算：

$$
\log_{10} T_{\alpha} = a_Z Q_{\alpha}^{-1/2} + b_Z, \quad (27)
$$

054314-15

**物理评论 C 110, 054314 (2024)**

**[图 17]**：$\delta Q_{\alpha}$ 值作为中子数函数，对应于 No ($Z = 102$) 到 Cn ($Z = 120$) 偶偶核，使用 SLy4 (左图) 和 UNEDF1$^{SO}$ (右图) EDF。实验数据 $\delta Q_{\alpha}^{Expt}$ 用五角星表示。叉号表示根据相邻核的系统趋势估计的值。两种实验数据均根据参考文献 [55] 中采用的 $Q_{\alpha}$ 计算。
*该图展示了 $\delta Q_{\alpha}$ 值随中子数的变化。图中显示在 $N=152, 162, 168, 174, 184$ 处存在明显的极小值，表明这些中子数处存在壳闭合。计算结果与实验数据基本一致，且 $\delta Q_{\alpha}$ 比 $Q_{\alpha}$ 更敏感地揭示了壳隙。*

其中

$$
a_Z = 1.66175Z - 8.5166, \\
b_Z = -0.20228Z - 33.9069,
$$

其中 $Q_{\alpha}$ 是 $\alpha$ 衰变能 (以 MeV 为单位)，$Z$ 是母核的质子数。这是一个基于 Wenzel-Kramers-Brillouin (WKB) 近似现象学公式，其参数在参考文献 [75] 中给出。使用方程 (27) 计算的半衰期 ($T_{\alpha}$) 以秒为单位。

图 18 绘制了图 16 中所示 $\alpha$ 衰变核的 $\log_{10} T_{\alpha}$ 值。星号表示从参考文献 [76] 中提取的实验 $\log_{10} T_{\alpha}$ 值。除了 SkM* EDF，我们的计算很好地再现了起始于 $^{276,278}Fl$ 的 $\alpha$ 衰变链的实验 $\log_{10} T_{\alpha}$ 值，如图 18(b) 和 18(c) 所示。

**物理评论 C 110, 054314 (2024)**

**[图 18]**：起始于 $^{282}Fl$ 和 $^{276,278}Fl$ 核的 $\alpha$ 衰变链的 $\log_{10} T_{\alpha}$ 值。星号表示参考文献 [76] 中的实验数据。
*该图展示了 $\alpha$ 衰变半衰期 $\log_{10} T_{\alpha}$ 随质量数 $A$ 的变化。图中显示 $\log_{10} T_{\alpha}$ 值通常随 $A$ 的增加而下降，但在特定质量数处出现局部极小值或极大值，反映了壳效应。计算结果与实验数据基本一致。*

关注包含 $^{270}Hs$ 的衰变链，我们从图 18(a) 中看到，$\alpha$ 衰变链中较轻的核 ($^{266}Sg, ^{262}Rf, ^{258}No$) 被计算为具有更长的半衰期。据推测，实验应该能够观察到这些核的 $\alpha$ 衰变链的延伸。我们尝试解释它们的缺失。检查图 6-9 表明，$\delta E_{sh}$ 值在 $^{270}Hs$ 处达到极小值。$\delta E_{sh}$ 值沿 $\alpha$ 衰变链 (从 $^{270}Hs$ 到 $^{266}Sg$) 急剧增加。这表明较轻核 ($^{266}Sg, ^{262}Rf, ^{258}No$) 的稳定性下降。相比之下，对于图 18(b) 和 18(c) 中较重的 $\alpha$ 衰变链，相关核位于 $\delta E_{sh}$ 图的谷中，表明 $\alpha$ 衰变链中核的稳定性增强。

对于 SLy4 EDF，$^{274}Ds$ 的 $\log_{10} T_{\alpha}$ 值被计算为显著小于 $^{270}Hs$ 的值，如图 18(a) 所示。然而，对于 UNEDF1 和 UNEDF1$^{SO}$，我们注意到比 SLy4 EDF 更长的半衰期。结合 UNEDF1$^{SO}$ 预测的 $^{270}Hs$ 东北方向核的更大绝对 $\delta E_{sh}$ 值 (参见图 9)，母核 $^{274}Ds$ 到 $^{270}Hs$ 的 $\alpha$ 衰变可能由于更长的半衰期而在实验中观测到。

**IV. 总结**

大多数质子数 $Z$ 介于 102 和 120 之间的超重核 (SHN) 被计算为变形的。由于破坏性的库仑相互作用，超重核的存在依赖于壳稳定性效应。虽然球形超重核的壳稳定性效应已得到充分研究，但变形超重核在自洽平均场模型方面的研究较少。

在本工作中，我们旨在通过从核密度泛函理论获得的单粒子 (s.p.) 能量中提取壳修正能，研究质子数 $Z \approx 98-120$ 的变形偶偶超重核的稳定性。Fm ($Z = 100$) 和 No ($Z = 102$) 同位素被纳入研究，以评估理论对这些实验上更知名的核的性能。

斯特鲁金斯基形式主义中使用的单粒子能量是使用基于混合基的无对称性限制核密度泛函理论计算的。能量密度泛函使用了四组参数 (SLy4、SkM*、UNEDF1 和 UNEDF1$^{SO}$)。特别是，我们包含了 UNEDF1 和 UNEDF1$^{SO}$ 的结果，后者基于前者，但增强了自旋-轨道耦合常数以提高锘区域奇 A 核低激发态的性能。

对于质子数 $Z = 98-124$、中子数 $N = 140-190$ 的偶偶核，计算预测 $Z < 114$ 的核具有扁长基态形状，而 $Z > 114$ 的核具有扁球基态形状，变形随 $Z$ 接近 114 和 $N$ 接近 180 而平滑减小。对于对力性质，我们的利普金-野上处理给出的对力隙合理地再现了 $Z = 102$ 同位素 (质子) 和 $N = 150$ 同中子素 (中子) 的实验奇偶质量差。

计算得到的总壳修正能被绘制为质子数和中子数的函数。使用四种 EDF 计算的等高线图显示在 $^{152}No_{152}$ 附近存在次极小值。对于 UNEDF1$^{SO}$ EDF，极小值位于 $N = 152, Z = 100$ 处。而其他 EDF 的等高线图在 $N = 152, Z = 104$ 处存在次极小值。这可以用单粒子能级中壳隙或闭合在相关核子数处打开来解释。等高线图显示，对于所有四种 EDF，在 $^{108}Hs_{162}$ 处存在绝对极小值，表明该变形核具有双幻数。高于 $^{270}Hs$，SLy4、SkM* 和 UNEDF1 EDF 预测变形次极小值出现在 $^{184}Lv_{168}$ 处。对于 UNEDF1$^{SO}$，次极小值转移到 $^{272}Cn_{164}$。UNEDF1$^{SO}$ EDF 与 UNEDF1 相比，自旋-轨道耦合常数略有增强，使其能够正确再现 $N = 152$ 和 $Z = 100$ 处的实验壳闭合。UNEDF1$^{SO}$ 预测 $^{270}Hs$ 处存在与其他三种 EDF 相同的变形壳；然而，它预测下一个变形壳位于 $^{276}Cn$，而其他 EDF 预测它位于 $^{284}Lv$。未来的实验有望验证这些预测。

使用计算得到的总能量，可以与实验数据进行比较。首先，我们从总能量中提取经验壳隙参数。实验数据显示，对于 $Z = 102$ 同位素，在 $N = 152$ 处存在峰值；对于 $N = 152$ 同中子素，在 $Z = 100$ 处存在峰值。数据还显示，对于 $Z = 108$ 同位素，在 $N = 162$ 处存在显著峰值；对于 $N = 162$ 同中子素，在 $Z = 108$ 处存在显著峰值。我们的计算合理地再现了这些特征。

054314-17

**物理评论 C 110, 054314 (2024)**

**[表 I]**：图 17 中 $\alpha$ 衰变核的宏观能量 $E_{Mac}$ (25)、总壳修正能 $\delta E_{sh}$ (13)、基态四极变形 $\beta_2$ (22) 和电荷半径 $r_{ch}$ (fm)。壳修正能 $\delta E_{sh}$ 使用 SLy4 和 UNEDF1$^{SO}$ EDF 计算。
*该表列出了 SLy4 和 UNEDF1$^{SO}$ EDF 计算的核性质，包括宏观能量、壳修正能、四极变形和电荷半径。数据显示不同 EDFs 对这些性质的预测存在差异，尤其是在壳修正能方面。*

**物理评论 C 110, 054314 (2024)**

**[表 I]**：(续)

**物理评论 C 110, 054314 (2024)**

**[表 I]**：(续)

其次，我们提取 $Q_{\alpha}$ 值，并将其绘制为中子数与相同质子数的关系图。我们使用 SLy4 的结果很好地再现了实验数据，因为其量级和峰值位置与数据一致。使用 UNEDF1 和 UNEDF1$^{SO}$ 的结果低估了 $N = 162$ 处凹陷的深度，因为它们的单粒子能级密度较大。我们检查了起始于 $^{282,278,276}Fl$ 的 $\alpha$ 衰变链的 $Q_{\alpha}$ 值，这些链分别穿过 $^{270,266,264}Hs$。计算得到的 $Q_{\alpha}$ 值合理地再现了数据。包含 $^{270}Hs$ 的 $\alpha$ 衰变链的 $Q_{\alpha}$ 值显示从 $^{270}Hs$ 到 $^{274}Cn$ 突然增加，这是由于壳修正能贡献的增加。为了量化壳效应在 $\alpha$ 衰变能中的贡献，我们引入了两个新量 $\delta Q_{\alpha}$ (23) 和 $\delta Q_{\alpha}^{Expt}$ (26)。$\delta Q_{\alpha}$ 值作为中子数的函数是中子壳隙的敏感指标。使用 $\delta Q_{\alpha}$， $Q_{\alpha}$ 值作为中子数函数的斜率的微小变化变得更加明显。

054314-20

**物理评论 C 110, 054314 (2024)**

最后，我们使用 Viola-Seaborg 公式计算了选定 $\alpha$ 衰变链的 $\alpha$ 衰变半衰期。计算得到的 $\log_{10} T_{\alpha}$ 值与可用的实验数据一致。$^{270}Hs$ $\alpha$ 衰变链中 $Q_{\alpha}$ 值的突然变化也发生在 $\log_{10} T_{\alpha}$ 中，随着 $\alpha$ 衰变链中中子数的增加，$\log_{10} T_{\alpha}$ 急剧增加。穿过 $^{270}Hs$ 的 $\alpha$ 衰变链中的核具有更长的半衰期。结合更深的壳修正能，其他核在衰变链中的 $\alpha$ 衰变可能在实验中观测到。因此，当前工作建议未来的实验可以集中研究 $^{270}Hs$、$^{276}Cn$ 和 $^{284}Lv$ 附近的 $\alpha$ 衰变链中的核。

**致谢**

作者感谢 A. T. Kruppa 的有益讨论。当前工作得到了国家自然科学基金 (项目号 12075068 和 11705038) 和中央高校基本科研业务费 (HIT.BRET.2021003) 的支持。

**附录**

表 I 列出了图 17 中 $\alpha$ 衰变核的 $E_{Mac}$、$\delta E_{sh}$、$\beta_2$ 和电荷半径。通过 $E_{Mac}$ (24)，可以提取新测量的 $Q_{\alpha}$ 值对应的 $Q_{\alpha}^{Expt}$。

[1] Y. T. Oganessian and V. K. Utyonkov, Rep. Prog. Phys. 78, 036301 (2015).
[2] Y. T. Oganessian, A. Sobiczewski, and G. M. Ter-Akopian, Phys. Scr. 92, 023003 (2017).
[3] S. Hofmann and G. Münzenberg, Rev. Mod. Phys. 72, 733 (2000).
[4] Y. T. Oganessian, V. K. Utyonkov, Y. V. Lobanov, F. S. Abdullin, A. N. Polyakov, R. N. Sagaidak, I. V. Shirokovsky, Y. S. Tsyganov, A. A. Voinov, G. G. Gulbekian, S. L. Bogomolov, B. N. Gikal, A. N. Mezentsev, S. Iliev, V. G. Subbotin, A. M. Sukhov, K. Subotic, V. I. Zagrebaev, G. K. Vostokin, M. G. Itkis et al., Phys. Rev. C 74, 044602 (2006).
[5] Y. T. Oganessian, V. K. Utyonkov, N. D. Kovrizhnykh, F. S. Abdullin, S. N. Dmitriev, A. A. Dzhioev, D. Ibadullayev, M. G. Itkis, A. V. Karpov, D. A. Kuznetsov, O. V. Petrushkin, A. V. Podshibiakin, A. N. Polyakov, A. G. Popeko, I. S. Rogov, R. N. Sagaidak, L. Schlattauer, V. D. Shubin, M. V. Shumeiko, D. I. Solovyev et al., Phys. Rev. C 106, 064306 (2022).
[6] S. A. Giuliani, Z. Matheson, W. Nazarewicz, E. Olsen, P.-G. Reinhard, J. Sadhukhan, B. Schuetrumpf, N. Schunck, and P. Schwerdtfeger, Rev. Mod. Phys. 91, 011001 (2019).
[7] G. G. Adamian, N. V. Antonenko, H. Lenske, L. A. Malov, and S.-G. Zhou, Eur. Phys. J. A 57, 89 (2021).
[8] V. Strutinsky, Nucl. Phys. A 95, 420 (1967).
[9] V. Strutinsky, Nucl. Phys. A 122, 1 (1968).
[10] M. Brack, J. Damgaard, A. Jensen, H. Pauli, V. Strutinsky, and C. Wong, Rev. Mod. Phys. 44, 320 (1972).
[11] P. Ring and P. Schuck, The Nuclear Many-Body Problem (Springer-Verlag, New York, 1980).
[12] A. Sobiczewski and K. Pomorski, Prog. Part. Nucl. Phys. 58, 292 (2007).
[13] A. Sobiczewski, Radiochim. Acta 99, 395 (2011).
[14] B. Nerlo-Pomorska, K. Pomorski, and J. Bartel, Phys. Rev. C 84, 044310 (2011).
[15] S. G. Nilsson, C. F. Tsang, A. Sobiczewski, Z. Szymański, S. Wycech, C. Gustafson, I.-L. Lamm, P. Möller, and B. Nilsson, Nucl. Phys. A 131, 1 (1969).
[16] P. Möller and J. R. Nix, J. Phys. G 20, 1681 (1994).
[17] S. Ćwiok, W. Dudek, P. Kaszyński, and W. Nazarewicz, Eur. Phys. J. A 23, 387 (2005).
[18] J. Zhang and H. F. Zhang, Phys. Rev. C 102, 044308 (2020).
[19] S. G. Nilsson and I. Ragnarsson, Shapes and Shells in Nuclear Structure (Cambridge University Press, Cambridge, 1995).
[20] M. Bender, P.-H. Heenen, and P.-G. Reinhard, Rev. Mod. Phys. 75, 121 (2003).
[21] J. D. Walecka, Theoretical Nuclear and Subnuclear Physics, 2nd ed. (Imperical College Press and World Scientific Publishing Co. Pte. Ltd. 2004).
[22] S. Ćwiok, J. Dobaczewski, P.-H. Heenen, P. Magierski, and W. Nazarewicz, Nucl. Phys. A 611, 211 (1996).
[23] M. Bender, K. Rutz, P.-G. Reinhard, J. A. Maruhn, and W. Greiner, Phys. Rev. C 60, 034304 (1999).
[24] A. T. Kruppa, M. Bender, W. Nazarewicz, P.-G. Reinhard, T. Vertse, and S. Ćwiok, Phys. Rev. C 61, 034313 (2000).
[25] M. Bender, in Fusion Dynamics at the Extremes (World Scientific Publishing Company, Singapore, 2001), pp. 51–64.
[26] W. Nazarewicz, M. Bender, S. Ćwiok, P.-H. Heenen, A. Kruppa, P.-G. Reinhard, and T. Vertse, Nucl. Phys. A 701, 165 (2002).
[27] N. Nikolov, N. Schunck, W. Nazarewicz, M. Bender, and J. Pei, Phys. Rev. C 83, 034305 (2011).
[28] P.-H. Heenen, J. Skalski, A. Staszczak, and D. Vretenar, Nucl. Phys. A 944, 415 (2015).
[29] T. Bürvenich, K. Rutz, M. Bender, P.-G. Reinhard, J. A. Maruhn, and W. Greiner, Eur. Phys. J. A 3, 139 (1998).
[30] W. Zhang, J. Meng, S. Zhang, L. Geng, and H. Toki, Nucl. Phys. A 753, 106 (2005).
[31] A. V. Afanasjev and O. Abdurazakov, Phys. Rev. C 88, 014320 (2013).
[32] S. E. Agbemava, A. V. Afanasjev, T. Nakatsukasa, and P. Ring, Phys. Rev. C 92, 054310 (2015).
[33] S. Raeder, D. Ackermann, H. Backe, R. Beerwerth, J. C. Berengut, M. Block, A. Borschevsky, B. Cheal, P. Chhetri, C. E. Düllmann, V. A. Dzuba, E. Eliav, J. Even, R. Ferrer, V. V. Flambaum, S. Fritzsche, F. Giacoppo, S. Götz, F. P. Heßberger, M. Huyse et al., Phys. Rev. Lett. 120, 232503 (2018).
[34] P. Jiang, Z. M. Niu, Y. F. Niu, and W. H. Long, Phys. Rev. C 98, 064323 (2018).
[35] M. Bender, W. Nazarewicz, and P.-G. Reinhard, Phys. Lett. B 515, 42 (2001).
[36] M. Bender and P.-H. Heenen, J. Phys.: Conf. Ser., 420, 012002 (2013).
[37] J. Dobaczewski, A. Afanasjev, M. Bender, L. Robledo, and Y. Shi, Nucl. Phys. A 944, 388 (2015).
[38] R.-D. Herzberg and P. Greenlees, Prog. Part. Nucl. Phys. 61, 674 (2008).
[39] M. Block, Nucl. Phys. A 944, 471 (2015).

054314-21

**物理评论 C 110, 054314 (2024)**

[40] Y. Shi, J. Dobaczewski, and P. T. Greenlees, Phys. Rev. C 89, 034309 (2014).
[41] S. Hofmann, S. Heinz, R. Mann, J. Maurer, G. Münzenberg, S. Antalic, W. Barth, H. G. Burkhard, L. Dahl, K. Eberhardt, R. Grzywacz, J. H. Hamilton, R. A. Henderson, J. M. Kenneally, B. Kindler, I. Kojouharov, R. Lang, B. Lommel, K. Miernik, D. Miller et al., Eur. Phys. J. A 52, 180 (2016).
[42] Y. Shi, Phys. Rev. C 98, 014329 (2018).
[43] Y. Shi and N. Hinohara, arXiv:2111.00144.
[44] N. Schunck, J. Dobaczewski, J. McDonnell, W. Satuła, J. Sheikh, A. Staszczak, M. Stoitsov, and P. Toivanen, Comput. Phys. Commun. 183, 166 (2012).
[45] W. Nazarewicz, T. R. Werner, and J. Dobaczewski, Phys. Rev. C 50, 2860 (1994).
[46] A. Kruppa, Phys. Lett. B 431, 237 (1998).
[47] T. Vertse, A. T. Kruppa, R. J. Liotta, W. Nazarewicz, N. Sandulescu, and T. R. Werner, Phys. Rev. C 57, 3089 (1998).
[48] P. Salamon, A. T. Kruppa, and T. Vertse, Phys. Rev. C 81, 064322 (2010).
[49] N. Tajima, Y. R. Shimizu, and S. Takahara, Phys. Rev. C 82, 034316 (2010).
[50] P. Möller and J. Nix, Nucl. Phys. A 536, 20 (1992).
[51] P. Möller, A. J. Sierk, T. Ichikawa, and H. Sagawa, At. Data Nucl. Data Tables 109-110, 1 (2016).
[52] E. Chabanat, P. Bonche, P. Haensel, J. Meyer, and R. Schaeffer, Nucl. Phys. A 635, 231 (1998).
[53] J. Bartel, P. Quentin, M. Brack, C. Guet, and H.-B. Håkansson, Nucl. Phys. A 386, 79 (1982).
[54] M. Kortelainen, J. McDonnell, W. Nazarewicz, P.-G. Reinhard, J. Sarich, N. Schunck, M. V. Stoitsov, and S. M. Wild, Phys. Rev. C 85, 024304 (2012).
[55] M. Wang, W. Huang, F. Kondev, G. Audi, and S. Naimi, Chin. Phys. C 45, 030003 (2021).
[56] F. Xu, P. Walker, J. Sheikh, and R. Wyss, Phys. Lett. B 435, 257 (1998).
[57] M. Bender, K. Rutz, P. G. Reinhard, and J. A. Maruhn, Eur. Phys. J. A 8, 59 (2000).
[58] G. F. Bertsch, C. A. Bertulani, W. Nazarewicz, N. Schunck, and M. V. Stoitsov, Phys. Rev. C 79, 034306 (2009).
[59] S. E. Agbemava, A. V. Afanasjev, A. Taninah, and A. Gyawali, Phys. Rev. C 99, 034316 (2019).
[60] P. T. Greenlees, J. Rubert, J. Piot, B. J. P. Gall, L. L. Andersson, M. Asai, Z. Asfari, D. M. Cox, F. Dechery, O. Dorvaux, T.
Grahn, K. Hauschild, G. Henning, A. Herzan, R.-D. Herzberg, F. P. Heßberger, U. Jakobsson, P. Jones, R. Julin, S. Juutinen et al., Phys. Rev. Lett. 109, 012501 (2012).
[61] J. Chi, Y. Qiang, C. Gao, and J. Pei, Nucl. Phys. A 1032, 122626 (2023).
[62] S. E. Agbemava, A. V. Afanasjev, and P. Ring, Phys. Rev. C 93, 044304 (2016).
[63] G. Audi, O. Bersillon, J. Blachot, and A. Wapstra, Nucl. Phys. A 729, 3 (2003).
[64] R. R. Chasman, I. Ahmad, A. M. Friedman, and J. R. Erskine, Rev. Mod. Phys. 49, 833 (1977).
[65] P. T. Greenlees, R.-D. Herzberg, S. Ketelhut, P. A. Butler, P. Chowdhury, T. Grahn, C. Gray-Jones, G. D. Jones, P. Jones, R. Julin, S. Juutinen, T.-L. Khoo, M. Leino, S. Moon, M. Nyman, J. Pakarinen, P. Rahkila, D. Rostron, J. Sarén, C. Scholey et al., Phys. Rev. C 78, 021303(R) (2008).
[66] A. Sobiczewski, I. Muntian, and Z. Patyk, Phys. Rev. C 63, 034306 (2001).
[67] J. Dvorak, W. Brüchle, M. Chelnokov, R. Dressler, C. Düllmann, K. Eberhardt, V. Gorshkov, E. Jäger, R. Krücken, A. Kuznetsov, Y. Nagame, F. Nebel, Z. Novackova, Z. Qin, M. Schädel, B. Schausten, E. Schimpf, A. Semchenkov, P. Thörle, A. Türler et al., Phys. Rev. Lett. 97, 242501 (2006).
[68] Y. T. Oganessian, V. K. Utyonkov, F. S. Abdullin, S. N. Dmitriev, R. Graeger, R. A. Henderson, M. G. Itkis, Y. V. Lobanov, A. N. Mezentsev, K. J. Moody, S. L. Nelson, A. N. Polyakov, M. A. Ryabinin, R. N. Sagaidak, D. A. Shaughnessy, I. V. Shirokovsky, M. A. Stoyer, N. J. Stoyer, V. G. Subbotin, K. Subotic et al., Phys. Rev. C 87, 034605 (2013).
[69] J. C. Pei, F. R. Xu, Z. J. Lin, and E. G. Zhao, Phys. Rev. C 76, 044326 (2007).
[70] A. Staszczak, A. Baran, and W. Nazarewicz, Phys. Rev. C 87, 024320 (2013).
[71] M. Wang, G. Audi, A. Wapstra, F. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, Chin. Phys. C 36, 1603 (2012).
[72] E. Olsen and W. Nazarewicz, Phys. Rev. C 99, 014317 (2019).
[73] D. Guan and J. Pei, Phys. Lett. B 851, 138578 (2024).
[74] V. E. Viola and G. Seaborg, J. Inorg. Nucl. Chem. 28, 741 (1966).
[75] A. Sobiczewski, Z. Patyk, and S. Ćwiok, Phys. Lett. B 224, 1 (1989).
[76] G. Audi, F. Kondev, M. Wang, B. Pfeiffer, X. Sun, J. Blachot, and M. MacCormick, Chin. Phys. C 36, 1157 (2012).

054314-22
[[TRANSLATION_END]]