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

void draw_single_plot_Veto(TTree* tree, int ch) {

if (!tree) {

printf("错误: 传入的TTree指针无效！\n");

return;

}

  

// 动态生成唯一的直方图名称，避免冲突

TString hist_name = Form("Veto_tdiff_%d", ch);

  

// 动态生成筛选条件

TString cut_condition = Form("DSSDX_mul==1 && DSSDY_mul==1 && Veto_mul==1 && Veto_Ch[0] ==%d", ch);

  

// 构建完整的Draw命令

TString draw_command = Form("DSSDX_Ts[0]-Veto_Ts[0]>>%s(100,-200,200)", hist_name.Data());

  

// 创建一个新的画布来显示新图

TString canvas_name = Form("c_%d", ch);

TString canvas_title = Form("Time Difference for Veto Ch %d", ch);

new TCanvas(canvas_name, canvas_title, 800, 600);

  

// 执行绘图

printf("正在为通道 %d 绘图...\n", ch);

tree->Draw(draw_command, cut_condition, "h");

}

  
  

// -------------------------------------------------------

// 2. 这是这个宏的“主函数”

// 当你在ROOT中执行 .x draw_tdiff.C 时，这个函数会被调用

// -------------------------------------------------------

void draw_tdiffVeto() {

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

draw_single_plot_Veto(tr_map,0);

draw_single_plot_Veto(tr_map,1);

draw_single_plot_Veto(tr_map,2);

// 您可以继续添加更多调用

// draw_single_plot(tr_map, 24, 31);

}