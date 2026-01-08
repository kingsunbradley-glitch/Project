#include "TFile.h"

#include "TH2.h"

#include "TCanvas.h"

#include "TStyle.h"

#include "TString.h"

#include <iostream>

#include <map>



void plot_d1_grid_XY() {

    // 1. 设置绘图风格

    gStyle->SetOptStat(0);      

    gStyle->SetPalette(kBird);  

    gStyle->SetTitleSize(0.06, "t");

    gStyle->SetTitleFont(132, "t"); // 使用 Times New Roman 字体，数学符号支持更好



    // 2. 打开文件

    TFile *file = TFile::Open("qa_plots.root");

    if (!file || file->IsZombie()) {

        std::cout << "Error: 无法打开 qa_plots.root" << std::endl;

        return;

    }



    // 3. 定义新数据的映射关系 (Run -> D1)

    std::map<int, double> d1_map;

    //（0.343,0.58,Q2）

    d1_map[335] = 1.06;

    d1_map[347] = 1.05;

    d1_map[332] = 1.04;

    d1_map[345] = 1.03;

    d1_map[326] = 1.02;

    d1_map[328] = 1.01;

    d1_map[330] = 1.00;
    int run_order[] = {335,347,332,345,326,328,330};
    

    // 4. 定义绘图顺序 (按照你图片中的列表顺序)

    // 只有 5 个数据，Pad 6 将会留空

    /*
    //（0.343,Q1,1.02）=
    d1_map[318] = 0.60;

    d1_map[326] = 0.58;

    d1_map[361] = 0.57;

    d1_map[310] = 0.56;

    d1_map[296] = 0.54;

    d1_map[276] = 0.52;

    d1_map[260] = 0.50;

    int run_order[] = {260,276,296,310,361,326,318};
    */

    int total_runs = 7;



    // 5. 创建画布
    TCanvas *c1 = new TCanvas("c1", "D1 Field Scan New", 3*1200, 2*900);
    c1->Divide(4, 2);
    // 6. 循环绘制
    for (int i = 0; i < total_runs; i++) {

        int run_num = run_order[i];

       

        // 切换 Pad (从 1 开始)

        c1->cd(i + 1);

        gPad->SetRightMargin(0.12);



        TString histPath = TString::Format("run%05d/hXY_p3", run_num);

        TH2F *h2 = (TH2F*)file->Get(histPath);



        if (h2) {

            double d1_val = d1_map[run_num];



            // 使用 #bullet (实心圆点) 来代替 #cdot，彻底解决乱码问题

            TString newTitle = TString::Format("Q2=%.2f T / m", d1_val);

            h2->SetTitle(newTitle);

           

            h2->Draw("COLZ");

        } else {

            std::cout << "Warning: 没找到 " << histPath << std::endl;

        }

    }



    c1->Update();

}
/*
#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>
#include <map>

void plot_d1_grid_XY() {
    // 1. 设置绘图风格
    gStyle->SetOptStat(0);      
    gStyle->SetPalette(kBird);  
    gStyle->SetTitleSize(0.06, "t"); 
    gStyle->SetTitleFont(132, "t"); 

    // 2. 打开文件
    TFile *file = TFile::Open("qa_plots.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: 无法打开 qa_plots.root" << std::endl;
        return;
    }

    // 3. 定义映射关系 (Run -> Q1/D1 value)
    std::map<int, double> d1_map;
    d1_map[318] = 0.60;
    d1_map[326] = 0.58;
    d1_map[361] = 0.57;
    d1_map[310] = 0.56;
    d1_map[296] = 0.54;
    d1_map[276] = 0.52;
    d1_map[260] = 0.50;

    // 4. 定义绘图顺序
    int run_order[] = {260, 276, 296, 310, 361, 326, 318};
    int total_runs = 7; 

    // 5. 创建画布
    // 宽x3, 高x2 的尺寸
    TCanvas *c1 = new TCanvas("c1", "Q1 Field Scan", 3600, 1800);
    c1->Divide(4, 2); // 切分成 4列 x 2行 = 8个格子

    // 6. 循环绘制
    for (int i = 0; i < total_runs; i++) {
        int run_num = run_order[i];
        
        // ============ 关键修改逻辑 ============
        // 原始逻辑是 pad_index = i + 1
        // 我们需要跳过 4。所以当 i+1 >= 4 时，我们让它再去 +1
        
        int pad_index = i + 1;
        if (pad_index >= 4) {
            pad_index = pad_index + 1; 
        }
        // ======================================

        // 切换到计算好的 Pad
        c1->cd(pad_index);
        gPad->SetRightMargin(0.12); 

        TString histPath = TString::Format("run%05d/hXY_p3", run_num);
        TH2F *h2 = (TH2F*)file->Get(histPath);

        if (h2) {
            double d1_val = d1_map[run_num];

            // 【修正】TString::Format 参数匹配
            // 原代码传了两个变量但只有一个%f。
            // 这里改成显示 "Run XXX (Q1=0.xxx T/m)" 比较清晰
            TString newTitle = TString::Format("Q1=%.3f T/m", d1_val);
            
            // 如果你只想要 Q1=... 不想要 Run Number，就用下面这行：
            // TString newTitle = TString::Format("Q1=%.3f T/m", d1_val);

            h2->SetTitle(newTitle);
            h2->Draw("COLZ");
        } else {
            std::cout << "Warning: 没找到 " << histPath << std::endl;
        }
    }

    c1->Update();
}
    */
