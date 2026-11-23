# Codex 任务：使用 LaTeX/TikZ 绘制 261Rf、265Sg、269Hs 理论能级示意图

## 一、任务目标

根据我提供的理论能级表，参考对话中最后上传的**横向短线表示能级、按核素分列**的论文示意图，绘制一张清晰、可用于论文或学术报告的矢量能级图。

**只绘制以下三个核素，顺序从左到右：**

1. $^{261}\mathrm{Rf}$（$Z=104,\ N=157$）
2. $^{265}\mathrm{Sg}$（$Z=106,\ N=159$）
3. $^{269}\mathrm{Hs}$（$Z=108,\ N=161$）

### 画面内容（严格要求）

- 每个核素占一列，在**下方**标出相应的核素名称（如 $^{261}\mathrm{Rf}$）。
- 每个理论能级画成**一条短水平实线**；线的垂直位置由其相对基态的激发能决定。
- 三个核素使用**相同的垂直能量比例尺**；每个核素自身的基态都在 $E_x=0$，三条基态水平对齐。
- **主图只显示核素名称与能级短横线。** 不显示激发能的数值、不显示能量轴或刻度、不显示 $J^\pi$ / Nilsson 轨道标签、不显示衰变箭头、不显示网格或图例。
- **不要人为把相邻能级拉开**，要保留真实的相对高度；即使两个能级只相差 10 keV，也要在同一垂直比例尺下表示。
- 最终采用黑色线条、白色背景、尽量接近给定示意图的简洁学术风格。

> 注意：最后上传的六核素示例图**只作为绘图样式参考**。这张新图使用下述三个核素的**理论能级数值**，不要从示例图中量取能级高度，也不要把其他核素加进图内。

## 二、必须使用的数据

下表的数据是**绘图坐标输入**，而不是需要在图片上标注的文字。激发能统一用 keV；相对于各自的基态，故每一列的最低能级都是 0 keV。

| 核素 | 激发能 $E_x$ (keV) | 主导 Nilsson 组态（只用于辨认数据，不绘制到主图） |
|---|---:|---|
| $^{261}\mathrm{Rf}$ | 0 | $\nu 11/2^{-}[725]$ |
| $^{261}\mathrm{Rf}$ | 50 | $\nu 3/2^{+}[622]$ |
| $^{261}\mathrm{Rf}$ | 60 | $\nu 7/2^{+}[613]$ |
| $^{261}\mathrm{Rf}$ | 80 | $\nu 1/2^{+}[620]$ |
| $^{261}\mathrm{Rf}$ | 280 | $\nu 9/2^{+}[615]$ |
| $^{265}\mathrm{Sg}$ | 0 | $\nu 7/2^{+}[613]$ |
| $^{265}\mathrm{Sg}$ | 50 | $\nu 9/2^{+}[615]$ |
| $^{265}\mathrm{Sg}$ | 210 | $\nu 11/2^{-}[725]$ |
| $^{265}\mathrm{Sg}$ | 280 | $\nu 3/2^{+}[622]$ |
| $^{265}\mathrm{Sg}$ | 310 | $\nu 1/2^{+}[620]$ |
| $^{269}\mathrm{Hs}$ | 0 | $\nu 9/2^{+}[615]$ |
| $^{269}\mathrm{Hs}$ | 50 | $\nu 7/2^{+}[613]$ |
| $^{269}\mathrm{Hs}$ | 340 | $\nu 3/2^{+}[622]$ |
| $^{269}\mathrm{Hs}$ | 350 | $\nu 11/2^{-}[725]$ |
| $^{269}\mathrm{Hs}$ | 360 | $\nu 1/2^{+}[620]$ |

请核对每列**恰好五条能级短线**（共 15 条）。特别注意 $^{269}\mathrm{Hs}$ 顶部 340、350、360 keV 是三条**不同的能级**，不要合并；$^{261}\mathrm{Rf}$ 的 50、60 keV 也是两条不同的线。

## 三、布局与美术规范

- 横向三列间距应一致。建议每条水平线长约 1.6–1.8 cm、列距约 3.5 cm。
- 纵坐标从 0 到至少 360 keV，共用同一转换比例。建议用 **20 cm/MeV**，即 **0.02 cm/keV**，因此相隔 10 keV 的两线在物理画幅中相距 0.2 cm，能够被清晰辨认。
- 三列 0 keV 的水平线严格齐平。
- 核素名在对应线段中心下方，数学模式排版；不显示核素的 $Z,N$ 数字（核素自身的质量数除外）。
- 能级线宽建议 0.7–0.9 pt，无圆点、无箭头、无装饰；避免把 10 keV 间隔的线画成一条粗线。
- 整幅图以 PDF 矢量为主；需要预览时另外输出 PNG。不要使用截图作为最终成图。

## 四、推荐可直接编译的 TikZ 基础代码

请在工作目录创建 `nuclear_levels.tex`，以下代码可以作为最终实现的起点。只要没有新的明确要求，**不要增加主图文字或轴线**。

```tex
\documentclass[tikz,border=9pt]{standalone}
\usepackage{amsmath}
\begin{document}
\begin{tikzpicture}[x=1cm,y=0.02cm]

  % All y coordinates are excitation energies in keV.
  % x coordinate represents columns, not energy.
  \tikzset{level/.style={line width=0.8pt, black}}

  % Horizontal layout (in cm)
  \def\levelwidth{1.70}
  \def\xRf{0}
  \def\xSg{3.5}
  \def\xHs{7.0}

  % 261Rf (N=157): 0, 50, 60, 80, 280 keV
  \foreach \energy in {0,50,60,80,280} {
    \draw[level] (\xRf,\energy) -- ++(\levelwidth,0);
  }

  % 265Sg (N=159): 0, 50, 210, 280, 310 keV
  \foreach \energy in {0,50,210,280,310} {
    \draw[level] (\xSg,\energy) -- ++(\levelwidth,0);
  }

  % 269Hs (N=161): 0, 50, 340, 350, 360 keV
  \foreach \energy in {0,50,340,350,360} {
    \draw[level] (\xHs,\energy) -- ++(\levelwidth,0);
  }

  % Only nuclide names are printed. No numerical energies / orbital labels.
  \node[anchor=north, font=\large] at (0.85,-25) {$^{261}\mathrm{Rf}$};
  \node[anchor=north, font=\large] at (4.35,-25) {$^{265}\mathrm{Sg}$};
  \node[anchor=north, font=\large] at (7.85,-25) {$^{269}\mathrm{Hs}$};

\end{tikzpicture}
\end{document}
```

### 编译示例

```bash
pdflatex -interaction=nonstopmode -halt-on-error nuclear_levels.tex
# Optional preview:
pdftoppm -f 1 -singlefile -r 300 -png nuclear_levels.pdf nuclear_levels
```

若使用其他 TeX 引擎也可，但应保证 PDF 为真正的矢量绘图，核素上标位置正确。

## 五、交付与验收

最终请生成：

1. `nuclear_levels.tex`：可编辑、可独立编译的源代码。
2. `nuclear_levels.pdf`：矢量最终图片。
3. `nuclear_levels.png`：高分辨率预览（如已具备 PDF→PNG 工具）。

请检查：共 3 列、每列 5 条线，共 15 条；最低水平线高度一致；**绝不显示激发能数值或轨道标签**；每个核素内部能级的高度与给定数据相符；顶部紧邻的 10 keV 级差仍能区分。

如需未来在能级旁加上 Nilsson 标记，请把它做成可选的第二版本，不要改变当前这份极简主图。
