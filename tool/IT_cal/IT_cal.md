# IT_cal：超重核低能电磁跃迁寿命的 Weisskopf + BrIcc 估算方法

## 1. 目的与基本思路

对于低能同质异能态（isomeric state, i.s.）的电磁退激，若只做数量级估算，可以把计算拆成两部分：

1. **核跃迁部分**：用 Weisskopf 单粒子估计得到纯 $\gamma$ 辐射跃迁率 $\lambda_\gamma$；
2. **原子电子部分**：用 BrIcc 得到内部转换系数 $\alpha_{\rm IC}$，从而把 $\gamma$ 衰变率修正为总电磁退激率。

最终：

$$
\boxed{
\lambda_{\rm IT}
=
\lambda_\gamma(1+\alpha_{\rm tot})
}
$$

以及

$$
\boxed{
T_{1/2}^{\rm IT}
=
\frac{\ln 2}{\lambda_\gamma(1+\alpha_{\rm tot})}
=
\frac{T_{1/2}^{\gamma}}{1+\alpha_{\rm tot}}
}
$$

其中：

- $T_{1/2}^{\gamma}$：只考虑光子发射的部分半衰期；
- $\alpha_{\rm tot}$：总内部转换系数；
- $T_{1/2}^{\rm IT}$：考虑 $\gamma$ 与内部转换后的总电磁退激半衰期。

---

## 2. 第一步：确定跃迁多极性

由初末态自旋和宇称决定允许的多极性。

选择定则：

$$
|J_i-J_f|\le L\le J_i+J_f
$$

电跃迁 $EL$ 的宇称关系：

$$
\pi_f=\pi_i(-1)^L
$$

磁跃迁 $ML$ 的宇称关系：

$$
\pi_f=\pi_i(-1)^{L+1}
$$

因此，应先选择满足 $\Delta J$ 和宇称变化的**最低允许多极性**。最低阶通常给出最大的跃迁率，因此最适合作 Weisskopf 数量级估算。

---

## 3. 第二步：Weisskopf 纯 $\gamma$ 跃迁率

一般电磁跃迁率可写成

$$
\lambda_\gamma(\pi L)
=
\frac{8\pi(L+1)}
{L[(2L+1)!!]^2}
\frac{1}{\hbar}
\left(
\frac{E_\gamma}{\hbar c}
\right)^{2L+1}
B(\pi L).
$$

Weisskopf 估计令

$$
B(\pi L)=1~{\rm W.u.}
$$

从而得到适合直接计算的数值公式。

### 本文例题需要的两个公式

对于 $E5$：

$$
\boxed{
\lambda_W(E5)
=
2.395\times10^{-12}
A^{10/3}E_\gamma^{11}
\quad{\rm s^{-1}}
}
$$

对于 $M3$：

$$
\boxed{
\lambda_W(M3)
=
10.415
A^{4/3}E_\gamma^7
\quad{\rm s^{-1}}
}
$$

这里必须注意：

$$
\boxed{E_\gamma\ {\rm 的单位必须是~MeV}}
$$

纯 $\gamma$ 半衰期为

$$
\boxed{
T_{1/2}^{\gamma}
=
\frac{\ln2}{\lambda_W}
}
$$

平均寿命为

$$
\tau_\gamma=\frac{1}{\lambda_W}.
$$

---

# 4. 例题一：$^{273}$Ds 的 $E5$ 跃迁

考虑

$$
^{273}{\rm Ds}:
\qquad
3/2^+[611]
\rightarrow
13/2^-[716].
$$

假设能级差为

$$
E_\gamma\simeq300~{\rm keV}
=0.300~{\rm MeV}.
$$

## 4.1 多极性

有

$$
\Delta J=
\left|
\frac{13}{2}-\frac{3}{2}
\right|
=5,
$$

且宇称发生变化：

$$
+\rightarrow-.
$$

最低允许多极性为

$$
\boxed{E5}.
$$

## 4.2 Weisskopf 跃迁率

取

$$
A=273,
$$

则

$$
\lambda_W(E5)
=
2.395\times10^{-12}
\times273^{10/3}
\times0.300^{11}.
$$

计算得到

$$
\boxed{
\lambda_W(E5)
\simeq5.60\times10^{-10}~{\rm s^{-1}}
}
$$

因此

$$
T_{1/2}^{\gamma}
=
\frac{0.693}{5.60\times10^{-10}}
\simeq1.24\times10^9~{\rm s}.
$$

换算为年：

$$
\boxed{
T_{1/2}^{\gamma}
\simeq39.2~{\rm yr}
}
$$

所以在 $1$ W.u. Weisskopf 基准下，

$$
\boxed{
^{273}{\rm Ds}:
\quad
3/2^+[611]
\xrightarrow[300~{\rm keV}]{E5}
13/2^-[716]
\quad
T_{1/2}^{\gamma}\approx39~{\rm yr}
}
$$

这是极慢的直接高多极 $\gamma$ 跃迁。

> 本例按当前讨论只给出纯 Weisskopf $\gamma$ 寿命；若还需要考虑内部转换，可继续使用第 5 节的 BrIcc 方法。

---

# 5. 内部转换：$\alpha_{\rm tot}$ 如何定义

BrIcc 中某一个电子亚壳层 $i$ 的内部转换系数定义为

$$
\boxed{
\alpha_i
=
\frac{\lambda_{{\rm IC},i}}
{\lambda_\gamma}
}
$$

因此

$$
\lambda_{{\rm IC},i}
=
\alpha_i\lambda_\gamma.
$$

对所有允许发生内转换的电子亚壳层求和：

$$
\lambda_{\rm IC}
=
\sum_i\lambda_{{\rm IC},i}
=
\lambda_\gamma\sum_i\alpha_i.
$$

所以

$$
\boxed{
\alpha_{\rm tot}
=
\sum_{\text{open subshells}}\alpha_i
}
$$

并有

$$
\boxed{
\lambda_{\rm IT}
=
\lambda_\gamma+\lambda_{\rm IC}
=
\lambda_\gamma(1+\alpha_{\rm tot})
}
$$

---

## 5.1 哪些电子壳层可以发生内部转换？

对某一个亚壳层 $i$，必须满足

$$
\boxed{
E_\gamma>E_{b,i}
}
$$

其中 $E_{b,i}$ 是该电子亚壳层结合能。

转换电子动能约为

$$
E_e=E_\gamma-E_{b,i}.
$$

如果

$$
E_\gamma<E_{b,i},
$$

则该亚壳层不能放出转换电子，因此

$$
\alpha_i=0.
$$

---

## 5.2 总 ICC 如何求和

例如若 $L_1,L_2,L_3$ 均打开：

$$
\alpha_L
=
\alpha_{L1}
+\alpha_{L2}
+\alpha_{L3}.
$$

如果 M 壳也全部打开：

$$
\alpha_M
=
\alpha_{M1}
+\alpha_{M2}
+\alpha_{M3}
+\alpha_{M4}
+\alpha_{M5}.
$$

总系数则为

$$
\boxed{
\alpha_{\rm tot}
=
\alpha_K+\alpha_L+\alpha_M+\alpha_N+\alpha_O+\cdots
}
$$

只对**能量上开放的亚壳层**求和。

### 注意

BrIcc 给出的 $\alpha_{L1}$、$\alpha_{L2}$ 等已经是完整亚壳层的 ICC，电子占据数已包含在原子计算中，因此：

$$
\boxed{\text{不要再额外乘 }2j+1\text{ 或亚壳层电子数}}
$$

---

# 6. BrIcc 中的参数从哪里找？

推荐使用：
https://bricc.anu.edu.au/
**T. Kibédi et al., Nucl. Instrum. Methods Phys. Res. A 589 (2008) 202–229**

BrIccFO 使用 relativistic Dirac–Fock 原子计算和 **Frozen Orbital** approximation，默认数据覆盖：

- $Z=5$–110；
- 所有主要原子壳层；
- $E1$–$E5$ 与 $M1$–$M5$；
- 低能到高能的广泛跃迁能区。

对于一个确定跃迁，需要输入：

| 参数 | 来源 |
|---|---|
| $Z$ | 元素原子序数 |
| $E_\gamma$ | 核能级差 |
| multipolarity | 由 $J^\pi_i\rightarrow J^\pi_f$ 选择定则得到 |
| data set | 推荐 BrIccFO |
| $\alpha_i,\alpha_{\rm tot}$ | BrIcc 输出 |
| 电子结合能 | BrIcc 原子参数表 / Appendix B.1 |

---

# 7. 例题二：$^{269}$Hs 的低能 $M3$ 跃迁

考虑

$$
^{269}{\rm Hs}:
\qquad
9/2^+[604]
\rightarrow
3/2^+[622].
$$

这里

$$
\Delta J=3,
\qquad
\Delta\pi=0.
$$

所以最低允许多极性为

$$
\boxed{M3}
$$

下一阶为 $E4$，因此可以写成

$$
\boxed{M3(E4)}.
$$

在没有实验 mixing ratio 的情况下，本例按纯 $M3$ 的 $1$ W.u. Weisskopf 基准估算，因为最低阶 $M3$ 的跃迁率远大于 $E4$。

取

$$
A=269,\qquad Z=108.
$$

---

## 7.1 Hs 的电子结合能

BrIcc 的 Hs 原子参数给出近似：

$$
E_K=173.200~{\rm keV},
$$

$$
E_{L1}=34.750~{\rm keV},
$$

$$
E_{L2}=33.740~{\rm keV},
$$

$$
E_{L3}=24.920~{\rm keV}.
$$

M 壳结合能约为

$$
9.270,\ 8.790,\ 6.760,\ 5.995,\ 5.620~{\rm keV},
$$

更外层结合能更低。

---

# 8. $^{269}$Hs：$E_\gamma=30$ keV

## 8.1 纯 $\gamma$ Weisskopf 寿命

$$
E_\gamma=30~{\rm keV}=0.030~{\rm MeV}.
$$

代入

$$
\lambda_W(M3)
=
10.415
A^{4/3}E_\gamma^7,
$$

得到

$$
\lambda_W(M3)
=
10.415
\times269^{4/3}
\times0.030^7
$$

$$
\boxed{
\lambda_W(M3)
\simeq3.96\times10^{-7}~{\rm s^{-1}}
}
$$

所以

$$
T_{1/2}^{\gamma}
=
\frac{0.693}{3.96\times10^{-7}}
\simeq1.75\times10^6~{\rm s},
$$

即

$$
\boxed{
T_{1/2}^{\gamma}(30~{\rm keV})
\simeq20.3~{\rm d}
}
$$

---

## 8.2 哪些电子壳层打开？

因为

$$
30<34.750=E_{L1},
$$

$$
30<33.740=E_{L2},
$$

而

$$
30>24.920=E_{L3},
$$

所以

$$
\boxed{
K,\ L_1,\ L_2\ {\rm closed};
\qquad
L_3,\ M,\ N,\ O,\ldots\ {\rm open}.
}
$$

因此

$$
\boxed{
\alpha_{\rm tot}(30)
=
\alpha_{L3}
+\sum_j\alpha_{Mj}
+\sum_j\alpha_{Nj}
+\sum_j\alpha_{Oj}
+\cdots
}
$$

不能把 $L_1$ 和 $L_2$ 加进去。

---

## 8.3 加入内部转换

正式计算时应从 BrIccFO 直接取得

$$
\alpha_{\rm tot}(Z=108,E_\gamma=30~{\rm keV},M3).
$$

作为当前数量级例题，可采用约

$$
\alpha_{\rm tot}\sim2.2\times10^2.
$$

于是

$$
T_{1/2}^{\rm IT}
=
\frac{1.75\times10^6}
{1+2.2\times10^2}
$$

得到约

$$
\boxed{
T_{1/2}^{\rm IT}(30~{\rm keV})
\sim8\times10^3~{\rm s}
\sim2.2~{\rm h}
}
$$

> 这里的 $\alpha_{\rm tot}\sim2.2\times10^2$ 是当前用于数量级说明的估算值。论文最终数值应以 BrIccFO 对 $Z=108$ 的直接输出为准。

---

# 9. $^{269}$Hs：$E_\gamma=50$ keV

## 9.1 纯 $\gamma$ Weisskopf寿命

$$
E_\gamma=50~{\rm keV}=0.050~{\rm MeV}.
$$

$$
\lambda_W(M3)
=
10.415
\times269^{4/3}
\times0.050^7.
$$

得到

$$
\boxed{
\lambda_W(M3)
\simeq1.41\times10^{-5}~{\rm s^{-1}}
}
$$

从而

$$
T_{1/2}^{\gamma}
=
\frac{0.693}{1.41\times10^{-5}}
\simeq4.91\times10^4~{\rm s},
$$

即

$$
\boxed{
T_{1/2}^{\gamma}(50~{\rm keV})
\simeq13.6~{\rm h}
}
$$

---

## 9.2 哪些电子壳层打开？

此时

$$
50>E_{L1},E_{L2},E_{L3},
$$

但

$$
50<E_K.
$$

因此

$$
\boxed{
K\ {\rm closed};
\qquad
L_1,L_2,L_3,M,N,O,\ldots\ {\rm open}.
}
$$

所以

$$
\boxed{
\begin{aligned}
\alpha_{\rm tot}(50)
={}&
\alpha_{L1}
+\alpha_{L2}
+\alpha_{L3}\\
&+\sum_j\alpha_{Mj}
+\sum_j\alpha_{Nj}
+\sum_j\alpha_{Oj}
+\cdots .
\end{aligned}
}
$$

---

## 9.3 加入内部转换

正式计算应使用

$$
\alpha_{\rm tot}(Z=108,E_\gamma=50~{\rm keV},M3)
$$

的 BrIccFO 直接输出。

当前数量级例题取

$$
\alpha_{\rm tot}\sim1.9\times10^2,
$$

则

$$
T_{1/2}^{\rm IT}
=
\frac{4.91\times10^4}
{1+1.9\times10^2}
$$

得到约

$$
\boxed{
T_{1/2}^{\rm IT}(50~{\rm keV})
\sim2.6\times10^2~{\rm s}
\sim4.3~{\rm min}
}
$$

---

# 10. 为什么 30 keV 与 50 keV 不能简单用幂律缩放 ICC？

虽然纯 $\gamma$ 的 Weisskopf 跃迁率严格满足

$$
\lambda_\gamma(M3)\propto E_\gamma^7,
$$

但内部转换系数 $\alpha_{\rm IC}$ 并不在整个低能区服从单一简单幂律。

原因是随着 $E_\gamma$ 增加，不同电子亚壳层依次越过结合能阈值：

$$
E_\gamma=E_{b,i}.
$$

例如 Hs：

- 30 keV：仅 $L_3$ 打开，$L_1,L_2$ 尚未打开；
- 50 keV：$L_1,L_2,L_3$ 全部打开。

因此 $\alpha_{\rm tot}$ 会受到电子壳层阈值的明显影响，必须使用 BrIcc 对各亚壳层进行计算并求和，不能仅从一个高能点按 $E^{-n}$ 外推。

---

# 11. 如果存在混合多极性

对于

$$
M3+E4
$$

若已知混合比

$$
\delta
=
\frac{\gamma(E4)}{\gamma(M3)},
$$

则 BrIcc 使用的有效 ICC 为

$$
\boxed{
\alpha_{\rm mix}
=
\frac{
\alpha(M3)+\delta^2\alpha(E4)
}{
1+\delta^2
}
}
$$

如果 $\delta$ 未知，则不能严格计算混合跃迁的 ICC。

对于当前 $^{269}$Hs 的数量级估算，由于 Weisskopf $M3$ 远快于 $E4$，采用

$$
\boxed{M3(E4)\approx M3}
$$

作为第一步估算。

---

# 12. 两个例题的结果汇总

| Nucleus | Transition | $E_\gamma$ | Multipolarity | $T_{1/2}^{\gamma}$ | IC treatment | $T_{1/2}^{\rm IT}$ |
|---|---|---:|---|---:|---|---:|
| $^{273}$Ds | $3/2^+[611]\rightarrow13/2^-[716]$ | 300 keV | $E5$ | $\approx39.2$ yr | 本例未加入 | — |
| $^{269}$Hs | $9/2^+[604]\rightarrow3/2^+[622]$ | 30 keV | $M3(E4)\approx M3$ | $\approx20.3$ d | BrIccFO，$\alpha_{\rm tot}\sim2.2\times10^2$ | $\sim2.2$ h |
| $^{269}$Hs | $9/2^+[604]\rightarrow3/2^+[622]$ | 50 keV | $M3(E4)\approx M3$ | $\approx13.6$ h | BrIccFO，$\alpha_{\rm tot}\sim1.9\times10^2$ | $\sim4.3$ min |

---

# 13. 实际计算时的标准流程

$$
\boxed{
\begin{array}{c}
J_i^\pi,\ J_f^\pi
\\
\Downarrow\\
\text{确定最低允许多极性 }\pi L
\\
\Downarrow\\
A,\ E_\gamma
\\
\Downarrow\\
\lambda_W(\pi L)
\\
\Downarrow\\
T_{1/2}^{\gamma}
=
\ln2/\lambda_W
\\
\Downarrow\\
Z,\ E_\gamma,\ \pi L
\ \text{输入 BrIccFO}
\\
\Downarrow\\
\text{检查电子结合能，确定 open subshells}
\\
\Downarrow\\
\alpha_{\rm tot}
=
\sum_{\rm open}\alpha_i
\\
\Downarrow\\
T_{1/2}^{\rm IT}
=
T_{1/2}^{\gamma}/(1+\alpha_{\rm tot})
\end{array}
}
$$

---

# 14. 参数清单

每算一个新跃迁，只需准备：

| 参数 | 如何获得 |
|---|---|
| $A$ | 核素质量数 |
| $Z$ | 元素原子序数 |
| $J_i^\pi$ | 初态能级指认 |
| $J_f^\pi$ | 末态能级指认 |
| $E_\gamma$ | 两能级能量差 |
| $\pi L$ | 由角动量与宇称选择定则判断 |
| $B(\pi L)$ | Weisskopf 估算时取 $1$ W.u. |
| $E_{b,i}$ | BrIcc 原子结合能表 |
| $\alpha_i$ | BrIccFO 各亚壳层输出 |
| $\alpha_{\rm tot}$ | 所有开放亚壳层 $\alpha_i$ 求和 |
| $\delta$ | 若为混合跃迁，需要实验或理论 mixing ratio |

---

# 15. 参考文献

1. T. Kibédi, T. W. Burrows, M. B. Trzhaskovskaya, P. M. Davidson, and C. W. Nestor Jr.,  
   *Evaluation of theoretical conversion coefficients using BrIcc*,  
   **Nucl. Instrum. Methods Phys. Res. A 589, 202–229 (2008)**.  
   DOI: 10.1016/j.nima.2008.02.051.

2. T. Kibédi, M. B. Trzhaskovskaya, M. Gupta, and A. E. Stuchbery,  
   *Conversion coefficients for superheavy elements*,  
   **At. Data Nucl. Data Tables 98, 313–355 (2012)**.

---

## 一句话总结

$$
\boxed{
\text{Weisskopf 给 }\lambda_\gamma;
\qquad
\text{BrIcc 给 }\alpha_{\rm tot};
\qquad
T_{1/2}^{\rm IT}
=
\frac{\ln2}
{\lambda_\gamma(1+\alpha_{\rm tot})}.
}
$$
