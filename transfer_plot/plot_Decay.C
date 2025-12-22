// draw_plot.C
// 该脚本用于绘制二维图，并在图上添加辅助线,用以验证转移产物衰变链。

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "TStyle.h"
#include <iostream>

// 定义主函数，并为所有参数提供默认值
void draw_plot(
    
    float x0 = 8090, float y0 = 7430,
    const char* filePath = "/home/evalie2/Project/document/273Ds/inter2/SS032_inter2_chain.root"
) { 
    float xmin = x0-50; float xmax = x0+50;
    float ymin = y0-50; float ymax = y0+50;
    // --- 1. 打开文件并获取TTree ---
    TFile *file = TFile::Open(filePath);
    if (!file || file->IsZombie()) {
        std::cerr << "错误：无法打开文件 " << filePath << std::endl;
        return;
    }

    TTree *tr_chain = (TTree*)file->Get("tr_chain");
    if (!tr_chain) {
        std::cerr << "错误：在文件中找不到名为 'tr_chain' 的TTree" << std::endl;
        file->Close();
        return;
    }
    std::cout << "文件和TTree成功加载." << std::endl;

    // --- 2. 创建画布和二维直方图 ---
    // 使用 gStyle 美化图像
    gStyle->SetOptStat(1111); // 显示统计框
    gStyle->SetPadGridX(true);  // 显示X网格
    gStyle->SetPadGridY(true);  // 显示Y网格

    TCanvas *c1 = new TCanvas("c1", "DSSD E1 vs E2", 800, 600);

    // ============================ 修改点 1 ============================
    // 交换这里的坐标轴标题，使其与您的要求一致 (X轴是E[1], Y轴是E[2])
    // 原来是: "DSSD E2 vs E1;DSSD_E[2] (keV);DSSD_E[1] (keV)"
    TH2F *h2 = new TH2F("h2", 
                        "DSSD E1 vs E2;DSSD_E[1] (keV);DSSD_E[2] (keV)", 
                        100, xmin, xmax, 
                        100, ymin, ymax);
    // =================================================================
    
    // 设置标记点样式（如果使用 "SCAT" 绘图选项）
    h2->SetMarkerStyle(kFullCircle);
    h2->SetMarkerSize(1.0);

    /*// --- 3. 填充并绘制直方图 ---
    // 定义筛选条件
    const char* cut = "Delta_Ts[1]>0 && SSD_E==0";

    // ============================ 修改点 2 ============================
    // 交换这里的绘图变量，使其符合 "Y:X" 的格式
    // 原来是: "DSSD_E[2]:DSSD_E[1] >> h2"
    // 现在Y轴是DSSD_E[2], X轴是DSSD_E[1]
    TString drawCmd = TString::Format("DSSD_E[2]:DSSD_E[1] >> h2");
    // =================================================================

    tr_chain->Draw(drawCmd, cut, "COLZ"); // "COLZ" 绘制彩色图并带图例
    std::cout << "二维图绘制完成." << std::endl;*/
    // --- 3. 绘制散点图 ---
    // 首先，绘制空的直方图 h2 来设置坐标轴范围和标题
    h2->Draw(); // 这一步会画出坐标轴框架

    // 定义筛选条件
    const char* cut = "Delta_Ts[1]>0 && SSD_E==0";

    // 将TTree的标记点样式设置为实心圆点
    // 注意：这里我们直接设置TTree的属性，而不是h2的
    tr_chain->SetMarkerStyle(kFullCircle);
    tr_chain->SetMarkerSize(1.0);

    // 接下来，将 TTree 中的数据以散点图的形式绘制在同一个画布上
    // 选项 "SAME" 告诉 ROOT 在已有的坐标轴上进行绘制，而不是新建一个
    tr_chain->Draw("DSSD_E[2]:DSSD_E[1]", cut, "SAME");

    std::cout << "散点图绘制完成." << std::endl;


    // --- 4. 创建并绘制辅助线 ---
    // 使用传入的参数 x0 和 y0
    TLine *vert_line = new TLine(x0, ymin, x0, ymax);
    TLine *horz_line = new TLine(xmin, y0, xmax, y0);

    // 设置线条样式
    vert_line->SetLineColor(kBlack);
    vert_line->SetLineStyle(kDashed);
    vert_line->SetLineWidth(2);

    horz_line->SetLineColor(kBlack);
    horz_line->SetLineStyle(kDashed);
    horz_line->SetLineWidth(2);

    // 在已有的图上叠加绘制线条
    vert_line->Draw("SAME");
    horz_line->Draw("SAME");
    std::cout << "辅助线绘制完成." << std::endl;
    
    // 文件指针不再需要，可以关闭 (注意：TTree会随文件关闭而失效)
    // 如果后续还需要操作TTree，可以暂时不关闭文件
    // file->Close(); 
    // delete file;
}