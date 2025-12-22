// =======================================================

// 文件: draw_tdiff.C

// 修正后的代码结构

// =======================================================

#include <TROOT.h>

#include <TCanvas.h>

#include <TTree.h>

#include <TFile.h>

#include <TString.h> // 注意：您代码中的 TStr 可能是 TString 的笔误，TString 是 ROOT 中常用的字符串类

#include <cstdio> // 为了使用 printf

  

// -------------------------------------------------------

// 1. 这是一个“辅助函数”，专门负责绘制一张图

// -------------------------------------------------------

void draw_single_plot_SSD(TTree* tree, int ch_min, int ch_max) {

if (!tree) {

printf("错误: 传入的TTree指针无效！\n");

return;

}

  

// 动态生成唯一的直方图名称，避免冲突

TString hist_name = Form("SSD_tdiff_%d_%d", ch_min, ch_max);

  

// 动态生成筛选条件

TString cut_condition = Form("DSSDX_mul==1 && DSSDY_mul==1 && SSD_mul==1 && (SSD_Ch[0]>=%d && SSD_Ch[0]<=%d)", ch_min, ch_max);

  

// 构建完整的Draw命令

TString draw_command = Form("DSSDX_Ts[0]-SSD_Ts[0]>>%s", hist_name.Data());

  

// 创建一个新的画布来显示新图

TString canvas_name = Form("c_%d_%d", ch_min, ch_max);

TString canvas_title = Form("Time Difference for SSD Ch %d-%d", ch_min, ch_max);

new TCanvas(canvas_name, canvas_title, 800, 600);

  

// 执行绘图

printf("正在为通道 %d-%d 绘图...\n", ch_min, ch_max);

tree->Draw(draw_command, cut_condition, "h");

}

  
  

// -------------------------------------------------------

// 2. 这是这个宏的“主函数”

// 当你在ROOT中执行 .x draw_tdiff.C 时，这个函数会被调用

// -------------------------------------------------------

void draw_tdiffSSD() {

// 假设 TTree 对象 tr_map 已经存在于内存中

// 如果不存在，你需要先打开文件

// TFile *f = TFile::Open("/path/to/your/file.root");

// TTree *tr_map = (TTree*)f->Get("tr_map");

  

TTree *tr_map = (TTree*)gROOT->FindObject("tr_map");

if (!tr_map) {

printf("错误：在内存中找不到名为 'tr_map' 的 TTree！请先用 TFile::Open 打开文件。\n");

return;

}

  

// 在主函数内部，按顺序调用辅助函数来完成你的任务

// 这个循环将自动绘制6张图：(0-7), (8-15), (16-23), (24-31), (32-39)

for (int i = 0; i < 5; ++i) // 0-4 对应通道范围 (0-7), (8-15), (16-23), (24-31), (32-39),(39-47);分别对应六块八条的SSD通道

{ draw_single_plot_SSD(tr_map, i * 8, i * 8 + 7); }

// 您可以继续添加更多调用

// draw_single_plot(tr_map, 24, 31);

}