---
name: english-paper-to-chinese
description: Translate English academic papers into accurate, fluent Simplified Chinese Markdown, especially nuclear physics, superheavy nuclei, ROOT/data-analysis, and PLB/PRL/PRC-style papers.
---

# English academic paper to Chinese translation skill

Use this skill when the user asks to translate English academic papers, abstracts, introductions, methods, results, discussions, figure captions, table captions, LaTeX manuscripts, Markdown files, or PDF-derived text into Chinese.

The goal is accurate academic Chinese translation. Preserve scientific meaning, uncertainty, terminology, equations, citation structure, isotope notation, units, and logic.

## Default output

Produce Simplified Chinese Markdown.

Preserve the original structure whenever possible:

- Title
- Abstract
- Introduction
- Experimental method / Methods
- Results
- Discussion
- Conclusion
- Figure captions
- Table captions
- References

If the user provides only one paragraph or one section, translate only that part.

If the user asks for bilingual output, use:

### Original

English text.

### 中文翻译

中文译文。

If the user asks for Chinese only, do not repeat the English.

## Translation principles

1. Translate meaning, not word order.
2. Preserve technical precision.
3. Preserve cautious scientific language.
4. Do not invent facts, references, equations, or interpretations.
5. Do not strengthen claims beyond the English original.
6. Keep equations, isotope notation, variables, and units unchanged.
7. Preserve citation markers such as [1], (Smith et al., 2020), \cite{}, \citep{}, and \citet{}.
8. Preserve figure and table references such as Fig. 1, Table I, Extended Data Fig. 2.

## Scientific tone

Use formal Chinese academic prose.

Prefer these Chinese expressions:

- 表明
- 显示
- 提示
- 可能表明
- 与……一致
- 在目前统计量下
- 候选同核异能态
- 暂定归类
- 工作方案

Avoid overconfident expressions unless the English original clearly supports them:

- 证明
- 确证
- 无可争议地说明
- 完全确定

## Nuclear physics terminology

Use these preferred translations:

superheavy nuclei -> 超重核
alpha decay -> α衰变
alpha-particle energy -> α粒子能量
decay chain -> 衰变链
cross section -> 截面
upper limit -> 上限
confidence level -> 置信水平
half-life -> 半衰期
mean lifetime -> 平均寿命
branching ratio -> 分支比
isomeric state -> 同核异能态
high-K isomer -> 高K同核异能态
candidate isomeric state -> 候选同核异能态
working scheme -> 工作方案
tentative assignment -> 暂定指认
deformed shell closure -> 形变壳闭合
shell gap -> 壳隙
Nilsson orbital -> Nilsson轨道
intruder orbital -> 侵入轨道
quasiparticle -> 准粒子
reduced alpha width -> 约化α宽度
barrier penetrability -> 位垒穿透率
angular momentum transfer -> 角动量转移
DSSD -> 双面硅条探测器（DSSD）
MWPC -> 多丝正比室（MWPC）
separator transmission -> 分离器传输效率
recoil -> 反冲核
implantation event -> 注入事件
correlation time -> 关联时间
background -> 本底
beam dose -> 束流剂量
target thickness -> 靶厚

## Isotope notation

Keep isotope notation standard.

Preferred LaTeX:

$^{273}\mathrm{Ds}$
$^{269}\mathrm{Hs}$
$Q_{\alpha}$
$T_{1/2}$

In plain Chinese text, use:

^273Ds
^269Hs
Qα
T1/2

## PLB / PRL / PRC caution

When translating claims about rare events, few decay chains, isomers, or new states, preserve caution.

Example:

suggests the presence of an isomeric state

Translate as:

提示可能存在一个同核异能态

Do not translate as:

证明存在一个同核异能态

Example:

The events are grouped according to a working scheme.

Translate as:

这些事件按照一个工作方案进行归类。

Do not translate as:

这些事件被确定为不同能级。

## Small-number statistics

For small statistics, preserve statistical language carefully.

one event -> 一个事件
two decay chains -> 两条衰变链
zero observed events -> 未观测到事件
upper limit -> 上限
90% confidence level -> 90%置信水平
asymmetric uncertainty -> 非对称不确定度
likelihood -> 似然
bootstrap -> bootstrap / 自助法
AIC/BIC -> AIC/BIC
Gaussian mixture model -> 高斯混合模型

Do not convert asymmetric uncertainties into symmetric ones.

Preserve forms such as:

$219^{+292}_{-147}$ fb

## Units and symbols

Keep units readable and standard:

MeV
keV
ms
s
fb
pb
mb
μb
mg/cm^2
pnA

Use Chinese punctuation around prose, but keep mathematical expressions and units unchanged.

Example:

该事件的 α 粒子能量为 10.858(22) MeV，关联时间为 8.320 ms。

## Equations

Do not translate variable names.

Translate explanatory text around equations.

Example:

where N is the number of observed decay chains

Translate as:

其中，N 为观测到的衰变链数。

Do not alter equation numbering.

## Figure and table captions

For captions:

1. Preserve figure/table number.
2. Preserve panel labels such as (a), (b), and (c).
3. Translate caption text into fluent Chinese.
4. Keep isotope names, model names, and dataset labels unchanged if they are technical identifiers.
5. Do not omit uncertainty definitions.

Example:

Fig. 1. Qα values as a function of neutron number.

Translate as:

图1. Qα 值随中子数的变化。

## LaTeX handling

If translating LaTeX source:

- Preserve LaTeX commands unless the user asks for plain Markdown.
- Translate text inside sections, captions, paragraphs, and footnotes.
- Do not translate command names.
- Do not break math environments.
- Preserve labels and refs such as \label{}, \ref{}, \cite{}, and equation environments.

For section titles with isotope notation, keep hyperref-safe LaTeX where possible.

Example:

\section{\texorpdfstring{$^{273}\mathrm{Ds}$}{273Ds} decay properties}

Translate as:

\section{\texorpdfstring{$^{273}\mathrm{Ds}$}{273Ds} 的衰变性质}

## Markdown handling

If translating Markdown:

- Preserve headings.
- Preserve tables.
- Preserve code blocks unchanged unless explicitly asked.
- Preserve image links.
- Translate alt text only if it is descriptive text.
- Do not translate file paths.

## Quality-control pass

After translation, check:

1. Are all technical terms consistent?
2. Are isotope names and mass numbers unchanged?
3. Are equations unchanged?
4. Are references and citations unchanged?
5. Are uncertainties and units unchanged?
6. Did the translation preserve cautious scientific claims?
7. Is the Chinese fluent rather than word-for-word?

If there are ambiguous terms, add a short section:

### 术语说明

Only add this section when useful.

## Output style

For ordinary translation requests, output only the translation.

For files, create a translated Markdown file when possible, using a name like:

original_name.zh-CN.md

If the user asks for a terminology table or comparison, provide it after the translation.
