# calcBass — LISE++ 风格的核融合入口道势垒计算器

本项目是**独立实现**，不是 LISE++ 官方源代码，也不能保证数值逐位等同。它基于公开的 Bass77、Bass80 核势公式，以及 LISE++ 2015 版公开教程展示的 Woods–Saxon 形式，实现 **球形核、固定取向近似、逐部分波的入口道势垒** 计算。

## 安装

Python >= 3.9，安装依赖：

```bash
python -m pip install -r requirements.txt
```

## 最简单的运行方式

```bash
python calc_bass.py
```

默认使用 `40Ar (A=40,Z=18) + 238U (A=238,Z=92)`，并分析三个总实验室系束流能量 `217 220 225 MeV`。默认比较 Bass77、Bass80、Woods–Saxon 三种模型，输出三个单独的 PNG 图及 CSV 数据。

显式指定参数：

```bash
python calc_bass.py --projectile 40 18 --target 238 92 \
    --energies 217 220 225 --model all --lmax 140 \
    --show-l 0 20 40 60 --out-dir results_ArU
```

选择一套模型：

```bash
python calc_bass.py --model bass80
python calc_bass.py --model bass77
python calc_bass.py --model woods-saxon --ws-depth 105 --ws-r0 1.12 --ws-a 0.75
```

额外计算**近似**的经典 / Hill–Wheeler 部分波俘获截面（仅作说明，**不是** LISE++ 的完整截面）：

```bash
python calc_bass.py --capture --hw 3.0 --lmax 140
```

如无 matplotlib，可运行 `--no-plots`，仍需 numpy、scipy。

## 输出目录及物理意义

| 文件 | 内容 |
|---|---|
| `summary.csv` | 核半径、`R_B`、`V_B`（c.m. / lab）、近似 `l_crit` |
| `beam_energies.csv` | 每个输入能量的 `E_cm` 和 `E_cm - V_B` |
| `potential_vs_r.csv` | 每个模型逐个半径的 `V_C`、`V_N`、`V_total (l=0)` |
| `barrier_vs_l.csv` | 每个 `l` 的外侧势垒位置、势垒高度和内侧极小值 |
| `potentials_l0.png` | 不同模型的 `V(r,l=0)` 总势能比较 |
| `potential_components.png` | 各项核力、库仑势、总势能的贡献 |
| `potentials_vs_l.png` | 单一模型中不同角动量的势能曲线 |
| `barrier_vs_l.png` | 势垒 `B_l` 的角动量依赖 |
| `capture_estimate.csv` | 仅当选择 `--capture` 时输出；简化的俘获截面 |

PNG 均为 **每幅图一个 Figure**，没有多子图。

## 计算公式

### 总相互作用势

`V(r,l) = V_C(r) + V_N(r) + hbar^2*l*(l+1)/(2*mu*r^2)`

- 取两个核的中心距离 `r`（fm），采用点电荷库仑势 `V_C = (1.4399764 MeV·fm) * Z1*Z2/r`。
- 约化质量近似 `mu*c^2 = (A1*A2/(A1+A2)) * 931.49410242 MeV`。
- 在势能曲线上寻找**外侧局部极大点** `dV/dr=0`，要求在其内侧存在局部极小点。对于 Bass 模型，内侧极小点的位置可能落在势模型未经验证的重叠区，只作为数值判据，不宜直接视为真实核动力学势阱。

### Bass77

`R_i = 1.16*A_i^(1/3) - 1.39*A_i^(-1/3)` fm

`Phi(s) = 1 / [0.030*exp(s/3.30) + 0.0061*exp(s/0.65)]`

### Bass80

`R_si = 1.28*A_i^(1/3) - 0.76 + 0.8*A_i^(-1/3)` fm

`R_i = R_si - 0.98/R_si` fm

`Phi(s) = 1 / [0.033*exp(s/3.50) + 0.007*exp(s/0.65)]`

对两个 Bass 模型：

`V_N(r) = -R1*R2/(R1+R2) * Phi(r-R1-R2)`

### Woods–Saxon

`V_N(r) = -V0/[1+exp((r-r0*(A1^(1/3)+A2^(1/3)))/a)]`

默认 `V0 = 105 MeV`, `r0 = 1.12 fm`, `a = 0.75 fm`，为 LISE++ 2015 教程中示例参数，不是 Bass77/80 的更新系数。

### 动能换算

对于靶静止时（非相对论、质量数近似）：

`E_cm = E_lab * A_target/(A_projectile+A_target)`

输入 `--energies` 的单位是束流**总动能 MeV**，不是 MeV/u。

### 可选俘获截面

采用 Hill–Wheeler 透射率：

`T_l(E) = [1+exp(2*pi*(B_l-E)/hbar_omega)]^-1`

`Sigma_capture = pi/k^2 * sum_l[(2*l+1)*T_l]`，`k²=2*mu*E/hbar²`，计算中将 fm² 转成 mb。

只对程序检测到数值势阱的部分波求和；这属于简化的入射道俘获估算，忽略了靶核形变、耦合道、准裂变、`P_CN`、蒸发竞争和 LISE++ 的私有实现细节，**不要把计算的截面直接用作 273Ds 蒸发余核截面或论文中的精确 LISE++ 结果**。

## 参考结果

对 `40Ar + 238U`（`l=0`）应得到：

| 模型 | `R_B` (fm) | `V_B` (MeV, c.m.) | `V_B` (MeV, lab) |
|---|---:|---:|---:|
| Bass77 | 12.488 | 177.815 | 207.700 |
| Bass80 | 12.661 | 175.485 | 204.978 |
| WS 示例参数 | 12.182 | 181.861 | 212.425 |

可用官方 LISE++ 2015 例子 `12C + 238U` 作对照：

- LISE++ 官方截图 **64.22 MeV**（`Barrier position 11.55 fm`）。
- 该脚本按公开 Bass77 核势得 **64.512 MeV**（`R_B 11.504 fm`）。
- 该脚本按公开 Bass80 核势得 **63.218 MeV**（`R_B 11.761 fm`）。

差别说明不能以标准 Bass77/80 直接冒充 LISE++ 的内部 Bass 融合算法。

## 单元测试

```bash
python -m unittest -v test_calc_bass.py
```

## 使用与局限

1. **`V_B` 是入口道、球形近似、`l=0`（或指定 `l`）的外侧势垒**。并非复合核裂变势垒，也不是蒸发余核截面。
2. **`l_crit` 与内侧势阱对内核重叠区模型的外推十分敏感**；特别是 Bass77、Bass80 核势最初适用于表面相互作用，不要将脚本显示的 `l_crit` 当成 LISE++ 所采用的严格值。
3. 对于强形变 `238U`，尖端/侧面碰撞取向的势垒不同。该脚本没有加入形变与耦合道，也没有复现 LISE++ 的 `P_CN`、裂变势垒或 `LisFus` 蒸发模型。
4. 使用核质量数近似代替精确原子核质量，能量换算使用非相对论公式；高能束流应换成精确质量/相对论运动学。

## 文献

- O. B. Tarasov, *Update of Fusion reaction mechanism in LISE++* (2015), [官方 PDF](https://lise.frib.msu.edu/9_10/9_10_Fusion.pdf)（第 6、8、14、18 页）。
- O. B. Tarasov and D. Bazin, *Development of the program LISE: application to fusion–evaporation*, Nucl. Instrum. Methods Phys. Res. B **204**, 174–178 (2003), [PDF](https://lise.frib.msu.edu/doc/lise_fusion.pdf)。
- X. Liu et al., *Systematic study of cluster radioactivity in trans-lead nuclei with various versions of proximity potential formalisms*, Chinese Physics C **48**, 054101 (2024), [公式(25)–(28)](https://cpc.ihep.ac.cn/article/doi/10.1088/1674-1137/ad260b)。

**注意：**不提供任何 LISE++ 的软件复制、反编译代码或版权资源；仅对公开核模型进行从零实现。
