/******************************************************************************
 * *
 * 一个用于分析多个ROOT文件的ROOT宏。                                           *
 * 功能：                                                                    *
 * 1. ... (同前) ...                                                         *
 * 5. 将详细的峰分析参数（能量，FWHM等）输出到一个制表符分隔的 .txt 文件中， *
 * 便于后续导入Excel等软件。                                              *
 * *
 ******************************************************************************/
//$ g++ c.cpp $(root-config --cflags --libs) -lSpectrum -o c

#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <fstream> // **新增**：包含文件流操作的头文件
#include <TH1.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLine.h>

// 定义主分析函数
void analyze_and_find_peaks( int startRun = 1, int endRun = 182,const char* expName = "SS032") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ... (前面链接文件、画图、扣本底的代码保持不变) ...
    // 1. 链接文件
    TChain* tr_chain = new TChain("tr_chain");
    const char* dataPath = "/home/evalie2/Project/document/273Ds/inter2/";
    for (int i = startRun; i <= endRun; ++i) {
        tr_chain->Add(TString::Format("%s%s%05d_chain.root", dataPath, expName, i));
    }
    if (tr_chain->GetNtrees() == 0) {
        std::cerr << "错误：没有文件被添加到TChain中。" << std::endl; return;
    }
    std::cout << "文件添加完成。总事件数: " << tr_chain->GetEntries() << std::endl;

    // 2. "在内存中"绘制原始一维谱
    tr_chain->Draw("DSSD_E>>h1(1000, 5000, 10000)", "Delta_Ts>0", "goff");
    TH1F* h1 = (TH1F*)gROOT->FindObject("h1");
    if (!h1) { std::cerr << "错误: 找不到直方图 h1！" << std::endl; return; }
    h1->SetTitle("Original Spectrum & Background");
    h1->GetXaxis()->SetTitle("Energy (channel)");
    h1->GetYaxis()->SetTitle("Counts");

    // 3. 估算并扣除本底
    TSpectrum *s = new TSpectrum();
    TH1 *h_bkg = s->Background(h1, 20, "");
    TH1F *h1_sub = (TH1F*)h1->Clone("h1_sub");
    h1_sub->Add(h_bkg, -1);
    h1_sub->SetTitle("Peaks on Background-Subtracted Spectrum");
    h1_sub->GetXaxis()->SetTitle("Energy (channel)");
    h1_sub->GetYaxis()->SetTitle("Net Counts");

    // 4. 生成图一：对比图
    TCanvas *c_comparison = new TCanvas("c_comparison", "Background Subtraction Comparison", 900, 600);
    h1->SetLineColor(kBlack);
    h_bkg->SetLineColor(kRed);
    h1_sub->SetLineColor(kBlue);
    h1->Draw("hist");
    h_bkg->Draw("same hist");
    h1_sub->Draw("same hist");
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "Original Spectrum", "l");
    leg->AddEntry(h_bkg, "Estimated Background", "l");
    leg->AddEntry(h1_sub, "Spectrum after Subtraction", "l");
    leg->Draw();
    
    // 5. 生成图二：寻峰、拟合并标注
    TCanvas *c_peaks = new TCanvas("c_peaks", "Peak Search and Fit Results", 1200, 800);
    h1_sub->Draw("hist");

    // 6. 寻峰
    Int_t nfound = s->Search(h1_sub, 2, "goff", 0.05);
    Double_t *xpeaks = s->GetPositionX();
    
    // **新增**：创建并打开一个用于输出的文本文件
    std::ofstream outfile("peak_analysis_results.txt");
    // **新增**：写入表头，使用制表符(\t)分隔，方便Excel导入
    outfile << "Peak_ID\tEnergy_Mean\tEnergy_Error\tSigma\tSigma_Error\tFWHM\tFWHM_Error\tN_Events" << std::endl;

    std::cout << "\n**************************************************" << std::endl;
    std::cout << "           Peak Analysis Results" << std::endl;
    std::cout << "**************************************************" << std::endl;

    // 7. 对找到的每个峰进行高斯拟合和标注
    // 7. 对找到的每个峰进行高斯拟合和标注
    for (int p = 0; p < nfound; p++) {
    Double_t xp = xpeaks[p];
    Int_t bin = h1_sub->GetXaxis()->FindBin(xp);
    Float_t yp = h1_sub->GetBinContent(bin);
    Double_t fit_range = 30.0;
    Double_t fit_min = xp - fit_range;
    Double_t fit_max = xp + fit_range;
    TF1 *fitFunc = new TF1("fitFunc", "gaus", fit_min, fit_max);
    fitFunc->SetParameters(yp, xp, 5);
    fitFunc->SetParNames("Constant", "Mean", "Sigma");
    h1_sub->Fit(fitFunc, "QR+");

    Double_t fit_const = fitFunc->GetParameter(0);
    Double_t fit_mean = fitFunc->GetParameter(1);
    Double_t fit_mean_err = fitFunc->GetParError(1);
    Double_t fit_sigma = fitFunc->GetParameter(2);
    Double_t fit_sigma_err = fitFunc->GetParError(2);
    Double_t fwhm = 2.355 * fit_sigma;
    Double_t fwhm_err = 2.355 * fit_sigma_err;

    // vvvvvvvvvv 【修改点 2：计算峰内事件数 n】 vvvvvvvvvv
    // 我们在峰中心左右 3*sigma 的范围内进行积分，这覆盖了约99.7%的事件
    Double_t integration_min = fit_mean - 3.0 * fit_sigma;
    Double_t integration_max = fit_mean + 3.0 * fit_sigma;
    Int_t bin_min = h1_sub->GetXaxis()->FindBin(integration_min);
    Int_t bin_max = h1_sub->GetXaxis()->FindBin(integration_max);
    // 使用 Integral() 函数计算指定 bin 范围内的计数总和
    Double_t n_events = h1_sub->Integral(bin_min, bin_max);
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // --- 在终端打印详细信息 ---
    std::cout << "--- Peak " << std::setw(2) << p + 1 << " ---" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  - Fit Mean (Energy) : " << std::setw(8) << fit_mean << " +/- " << fit_mean_err << std::endl;
    std::cout << "  - Fit Sigma (Width) : " << std::setw(8) << fit_sigma << " +/- " << fit_sigma_err << std::endl;
    std::cout << "  - FWHM              : " << std::setw(8) << fwhm << " +/- " << fwhm_err << std::endl;
    // vvvvvvvvvv 【修改点 3：在终端打印事件数】 vvvvvvvvvv
    std::cout << "  - N_Events (+-3sig) : " << std::setw(8) << n_events << std::endl;
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    std::cout << "------------------------------------------" << std::endl;

    // **新增**：将结果写入文件，同样使用制表符分隔
    outfile << std::fixed << std::setprecision(4); // 文件中可以保存更高精度
    // vvvvvvvvvv 【修改点 4：将事件数写入文件】 vvvvvvvvvv
    outfile << p + 1 << "\t"
            << fit_mean << "\t"
            << fit_mean_err << "\t"
            << fit_sigma << "\t"
            << fit_sigma_err << "\t"
            << fwhm << "\t"
            << fwhm_err << "\t"
            << n_events << std::endl;

        // ... (绘制引出线和标签的代码保持不变) ...
        /* 
        Double_t callout_length = fit_const * 0.2 + 200;
        const Double_t max_callout_length = 1500.0;
        if (callout_length > max_callout_length) {
            callout_length = max_callout_length;
        }
        Double_t line_x0 = fit_mean;
        Double_t line_y0 = fit_const;
        Double_t line_y1 = fit_const + callout_length;
        TLine *line = new TLine(line_x0, line_y0, line_x0, line_y1);
        line->SetLineColor(kRed);
        line->SetLineStyle(2);
        line->Draw();
        Double_t label_x = line_x0 + 20;
        Double_t label_y = line_y0 + callout_length * 0.5;
        TLatex *label = new TLatex(label_x, label_y, TString::Format("%.0f", fit_mean));
        label->SetTextAngle(90);
        label->SetTextAlign(22);
        label->SetTextColor(kRed);
        label->SetTextSize(0.025);
        label->Draw();*/
    }

    // **新增**：关闭文件并给出提示
    outfile.close();
    std::cout << "\n分析结果已成功保存到 peak_analysis_results.txt 文件中。" << std::endl;
    std::cout << "**************************************************" << std::endl;
}