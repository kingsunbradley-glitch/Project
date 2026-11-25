# PLB 双栏通栏图

当前文件：figures/excitation_functions_origin.pdf（矢量主文件），另有 SVG、600 dpi PNG。
最终尺寸 190 × 120 mm，依据 Elsevier 通用双栏图尺寸基准；PLB 在线作者指南访问返回 403，未核实其专门例外。
字体为已嵌入的 Times New Roman，按用户最新要求加大加粗：刻度 13 pt、轴标题 16 pt，曲线标签 17 pt、势垒标签 16 pt。
曲线 1.2 pt，图框 0.8 pt。PDF 字号检查和碰撞检查通过。
保持无图例、4n 黑虚线、5n 红实线、测量点红色、上限黑色、纵轴下限 1.1e-3 pb。
按用户要求省略 50 MeV 的 4n 上限；42、45 MeV 两个道重合的上限共用符号。原始数据不删除。

## Suggested caption

Calculated excitation functions for the 4n (black dashed line) and 5n (red solid line)
evaporation channels of the 278Ds compound system, compared with experimental results.
Black downward arrows denote the coincident 4n and 5n cross-section upper limits at
E*=42 and 45 MeV; the red point denotes the measured 5n cross section at E*=50 MeV.
Horizontal error bars extend by ±1 MeV for the upper limits and by −1/+1.2 MeV for
the measured point. The vertical error bar reproduces the reported cross-section
uncertainty of −0.08/+0.18 pb. The confidence levels are not specified in the supplied
table. Each dataset uses its supplied excitation energies. The 4n upper limit at
50 MeV is omitted from this display.

## Placement

Use a figure* environment to span both columns and include the PDF at the manuscript's
text width. If this width is smaller than 190 mm, all text scales by width/190;
check final superscripts again (6.3 pt reaches 6 pt at approximately 181 mm).
If the actual manuscript width is below 181 mm, regenerate at that width instead
of simply shrinking this figure. Keep the caption definitions of lines and symbols
because the artwork intentionally has no legend.

Official reference: https://www.elsevier.com/en-au/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing

最新标注：在 E*=43.5 MeV 处加入向下小箭头及 V_BASS，位置由用户指定，未另行计算 BASS 势垒；4n、5n 直接标在曲线旁。以上字体按图像可读性放大，不沿用此前正文式的小字号。
