#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TString.h"
#include <iostream>

void draw_lego() {
    // 1. 设置更美观的绘图风格
    gStyle->SetOptStat(0);       // 去掉统计框
    gStyle->SetPalette(kBird);   // 现代颜色 (kBird)
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetNumberContours(255); // 让颜色渐变更平滑

    // 2. 创建 TTree 并读取数据
    TTree *t = new TTree("t", "Data Tree");
    
    // 关键修正：显式指定列格式 "D1/D:Q1/D:Q2/D:N/D"
    // /D 表示 Double (浮点数)。这样 ROOT 就不会在读完第一列后停止了。
    // processed_result.txt 是文件名
    Long64_t nlines = t->ReadFile("processed_result.txt", "D1/D:Q1/D:Q2/D:N/D");

    if (nlines <= 0) {
        std::cout << "Error: File is empty or not found!" << std::endl;
        return;
    }
    std::cout << "Successfully read " << nlines << " lines." << std::endl;

    // 3. 提取 D1 = 0.347 的数据
    // 使用 TTree::Draw 将数据提取到内存数组中 (GetV1, GetV2, GetV3)
    // 过滤条件: abs(D1 - 0.347) < 0.001 (浮点数比较)
    t->Draw("Q1:Q2:N", "abs(D1 - 0.347) < 0.001", "goff");
    
    Int_t n = t->GetSelectedRows();
    if (n == 0) {
        std::cout << "No data found for D1 = 0.347" << std::endl;
        return;
    }
    std::cout << "Found " << n << " points for D1 = 0.347" << std::endl;

    // 4. 使用 TGraph2D 绘制 LEGO 图
    TCanvas *c1 = new TCanvas("c1", "Lego Plot", 900, 700);
    
    // 调整边距以便放下 Z 轴色标
    gPad->SetRightMargin(0.15); 
    gPad->SetLeftMargin(0.12);

    // 创建 Graph2D (利用 ROOT 的自动插值功能生成网格)
    // V1 是 Q1 (X), V2 是 Q2 (Y), V3 是 N (Z)
    TGraph2D *g = new TGraph2D(n, t->GetV1(), t->GetV2(), t->GetV3());
    
    g->SetTitle("Lego Plot (D1 = 0.347);Q_{1};Q_{2};N (Average)");
    
    // 设置插值网格密度 (决定了柱子的多少)
    // 数值越大，柱子越细密。40x40 通常效果不错。
    g->SetNpx(40);
    g->SetNpy(40);
    
    // 绘制命令
    // LEGO2: 彩色柱状图
    // Z: 显示色标
    g->Draw("LEGO2Z"); 

    // 调整视角 (Theta: 倾斜角, Phi: 旋转角)
    gPad->SetTheta(30);
    gPad->SetPhi(30);

    // 5. 保存图片
    c1->SaveAs("lego_0347.png");
    // c1->SaveAs("lego_0347.pdf"); // 如果需要矢量图
}