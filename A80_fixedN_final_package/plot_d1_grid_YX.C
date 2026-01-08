#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>
#include <map>

void plot_d1_grid_YX() {
    // 1. 设置绘图风格
    gStyle->SetOptStat(0);      // 隐藏统计框
    gStyle->SetPalette(kBird);  // 设置色阶
    gStyle->SetTitleSize(0.06, "t"); // 稍微调大标题字体，看得清楚点

    // 2. 打开文件
    TFile *file = TFile::Open("qa_plots.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: 无法打开 qa_plots.root" << std::endl;
        return;
    }

    // 3. 定义 Run Number 和 D1 值的映射 (根据你提供的图片)
    std::map<int, double> d1_map;
    /*
    d1_map[186] = 1.13;
    d1_map[157] = 1.11;
    d1_map[99] = 1.09;
    d1_map[98] = 1.07;
    //(0.347,Q1,1.13)*/
    /*
    d1_map[188] = 1.13;
    d1_map[159] = 1.11;
    d1_map[90] = 1.09;
    d1_map[123] = 1.07;
    //(0.347,Q1,1.13)*/
    // 4. 定义画布布局顺序 (2列 x 3行)
    // Pad顺序: 1(左上), 2(右上), 3(左中), 4(右中), 5(左下), 6(右下)

    int run_order[] = {71,82,84,86,99,90,92,94,96,101,115,117};
    d1_map[run_order[0]] = 1.19;
    d1_map[run_order[1]] = 1.18;
    d1_map[run_order[2]] = 1.17;
    d1_map[run_order[3]] = 1.15;
    d1_map[run_order[4]] = 1.13;
    d1_map[run_order[5]] = 1.11;
    d1_map[run_order[6]] = 1.09;
    d1_map[run_order[7]] = 1.07;
    d1_map[run_order[8]] = 1.05;
    d1_map[run_order[9]] = 1.03;
    d1_map[run_order[10]] = 1.00;
    d1_map[run_order[11]] = 0.96;

    // 5. 创建画布
    int x = 4;
    int y = 3;
    TCanvas *c1 = new TCanvas("c1", "D1 Field Scan", x*400, y*300);
    c1->Divide(x, y); // 2列 x 3行

    // 6. 循环绘制
    for (int i = 0; i < x * y; i++) {
        int run_num = run_order[i];
        
        // 切换 Pad
        c1->cd(i + 1);
        gPad->SetRightMargin(0.12); // 给色标留位置

        // 获取直方图
        TString histPath = TString::Format("run%05d/hXY_p3", run_num);
        TH2F *h2 = (TH2F*)file->Get(histPath);

        if (h2) {
            // 获取对应的 D1 值
            double d1_val = d1_map[run_num];

            // 设置新标题
            // %.3f 强制保留3位小数 (0.35 -> 0.350)
            // #cdot 在 ROOT 里会被渲染成中间的点 (·)
            TString newTitle = TString::Format("Q2=%.2f T / m", d1_val);
            h2->SetTitle(newTitle);
            
            // 绘制
            h2->Draw("COLZ");
        } else {
            std::cout << "Warning: 没找到 " << histPath << std::endl;
        }
    }

    c1->Update();
}