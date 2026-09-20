# 205Fr 传输效率误差计算

## 1. 误差传递公式、符号含义及计算方法

在本计算中，传输效率与计数、ACCT 束流强度、反应截面、靶厚和 $\alpha$ 分支比的关系写为

$$
\epsilon_{\mathrm{tran}}\propto
\frac{N}{\bar I\,\sigma\,t\,B_\alpha}.
$$

假定各误差来源相互独立，并采用一阶高斯误差传播，则

$$
\left(
\frac{u_{\epsilon_{\mathrm{tran}}}}
{\epsilon_{\mathrm{tran}}}
\right)^2
=
\left(\frac{u_N}{N}\right)^2
+\left(\frac{u_I}{\bar I}\right)^2
+\left(\frac{u_\sigma}{\sigma}\right)^2
+\left(\frac{u_t}{t}\right)^2
+\left(\frac{u_{B_\alpha}}{B_\alpha}\right)^2.
$$

各符号的含义和计算方法如下。

| 符号 | 含义 | 本计算中的取值或计算方法 |
|---|---|---|
| $\epsilon_{\mathrm{tran}}$ | 传输效率中心值 | 取 `205Fr.csv` 中的 `TranEff` |
| $u_{\epsilon_{\mathrm{tran}}}$ | 传输效率的合成标准不确定度 | 由五项相对误差的平方和计算 |
| $N$ | 目标核计数 | 取 `205Fr.csv` 中的 `Counts` |
| $u_N$ | 计数统计误差 | 假定计数服从泊松分布，$u_N=\sqrt{N}$ |
| $I_j$ | ACCT 文件中的第 $j$ 个原始测量值 | 从对应 run 的时间窗内读取 |
| $j$ | ACCT 采样点序号 | $j=1,2,\ldots,m$ |
| $m$ | 时间窗内 ACCT 采样点数 | 统计满足 $\tau_0\leq\tau<\tau_0+T_{\mathrm{sustain}}$ 的数据点数 |
| $\bar I$ | 对应时间窗内的 ACCT 平均值 | $\bar I=\frac{1}{m}\sum_{j=1}^{m}I_j$ |
| $s_I^2$ | ACCT 样本方差 | $s_I^2=\frac{1}{m-1}\sum_{j=1}^{m}(I_j-\bar I)^2$ |
| $u_I$ | ACCT 误差 | 采用样本标准差，$u_I=s_I=\sqrt{s_I^2}$ |
| $\sigma$ | 反应截面 | $\sigma=545\ \mathrm{\mu b}$ |
| $u_\sigma$ | 反应截面误差 | 相对误差为 1.2%，即 $u_\sigma=545\times0.012=6.54\ \mathrm{\mu b}$ |
| $t$ | 靶厚 | $t=0.5\ \mathrm{mg/cm^2}$ |
| $u_t$ | 靶厚误差 | 暂取 10%，即 $u_t=0.5\times0.10=0.05\ \mathrm{mg/cm^2}$ |
| $B_\alpha$ | $^{205}\mathrm{Fr}$ 的 $\alpha$ 分支比 | $B_\alpha=98.5\%=0.985$ |
| $u_{B_\alpha}$ | $\alpha$ 分支比误差 | $B_\alpha(^{205}\mathrm{Fr})=98.5(4)\%$ 表示绝对误差为 0.4 个百分点，即 $u_{B_\alpha}=0.004$ |
| $\tau$ | ACCT 数据的采样时刻 | 取 ACCT CSV 中的时间戳 |
| $\tau_0$ | run 的起始时刻 | 由 `date` 和 `Tsatrt` 共同确定 |
| $T_{\mathrm{sustain}}$ | run 的持续时间 | 取 `Tsustain/s`，单位为秒 |

因此，实际采用的无量纲相对误差公式为

$$
\frac{u_{\epsilon_{\mathrm{tran}}}}
{\epsilon_{\mathrm{tran}}}
=
\sqrt{
\frac{1}{N}
+\frac{s_I^2}{\bar I^2}
+0.012^2
+0.10^2
+\left(\frac{0.004}{0.985}\right)^2
}.
$$

计算过程中所有相对误差均使用无量纲小数，最后再乘以 100% 转换为百分数。若 `TranEff` 的百分数数值记为 $\epsilon_{\mathrm{tran},\%}$，则传输效率的绝对误差为

$$
u_{\epsilon_{\mathrm{tran}},\mathrm{pp}}
=
\epsilon_{\mathrm{tran},\%}
\frac{u_{\epsilon_{\mathrm{tran}}}}
{\epsilon_{\mathrm{tran}}},
$$

其中 `pp` 表示百分点。`DetectEff`、束流占空比和靶占空比在本次计算中视为无误差。

## 2. run00068 的计算过程

run00068 的日期为 2026-05-11，起始时刻为 17:28:00，持续时间为 600 s，因此 ACCT 数据选取范围为

$$
\mathrm{2026{-}05{-}11\ 17{:}28{:}00}
\leq \tau<
\mathrm{2026{-}05{-}11\ 17{:}38{:}00}.
$$

该时间窗内共有 $m=3437$ 个 ACCT 采样点。由原始数据计算得到

$$
\bar I=5.5873591864,
$$

$$
s_I^2=0.1642543380.
$$

首先将 ACCT 方差开平方，得到 ACCT 误差：

$$
u_I=\sqrt{s_I^2}
=\sqrt{0.1642543380}
=0.4052830344.
$$

ACCT 的无量纲相对误差为

$$
\frac{u_I}{\bar I}
=\frac{0.4052830344}{5.5873591864}
=0.07253570.
$$

run00068 的目标核计数为 $N=109027$，泊松计数误差为

$$
u_N=\sqrt{109027}=330.1923681735,
$$

对应的无量纲相对误差为

$$
\frac{u_N}{N}
=\frac{1}{\sqrt{109027}}
=0.00302854.
$$

反应截面、靶厚和 $\alpha$ 分支比的无量纲相对误差分别为

$$
\frac{u_\sigma}{\sigma}=0.012,
\qquad
\frac{u_t}{t}=0.10.
$$

$$
\frac{u_{B_\alpha}}{B_\alpha}
=\frac{0.4}{98.5}
=\frac{0.004}{0.985}
=0.00406091.
$$

将五项相对误差代入误差传播公式：

$$
\frac{u_{\epsilon_{\mathrm{tran}}}}
{\epsilon_{\mathrm{tran}}}
=
\sqrt{
0.00302854^2
+0.07253570^2
+0.012^2
+0.10^2
+0.00406091^2
}
=0.12422194.
$$

最后将总相对误差转换为百分数：

$$
0.12422194\times100\%=12.422194\%.
$$

run00068 的传输效率中心值为 $\epsilon_{\mathrm{tran}}=23.17\%$，因此绝对误差为

$$
u_{\epsilon_{\mathrm{tran}},\mathrm{pp}}
=23.17\times0.12422194
=2.8782\%.
$$

最终结果为

$$
\boxed{
\epsilon_{\mathrm{tran}}
=(23.17\pm2.88)\%
}.
$$

对应的对称误差区间为

$$
20.29\%\leq\epsilon_{\mathrm{tran}}\leq26.05\%.
$$
