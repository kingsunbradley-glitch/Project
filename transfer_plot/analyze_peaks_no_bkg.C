/******************************************************************************
 * *
 * 一个用于分析多个ROOT文件的ROOT宏。
 *
 * 功能：
 * 1. 链接 TChain 中的多个 .root 文件。
 * 2. 绘制原始一维能谱 (h1)。
 * 3. 【已修改】不在原始谱上进行本底扣除。
 * 4. 在原始谱 (h1) 上自动寻峰。
 * 5. 使用 "高斯 + 恒定本底" (gaus+pol0) 函数拟合找到的峰。
 * 6. 将详细的峰分析参数（能量，FWHM等）输出到一个制表符分隔的 .txt 文件中。
 * 7. 专门用于分析本底谱
 *
 ******************************************************************************/

#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <fstream> // 用于文件流操作
#include <TH1.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>
#include <TLine.h>

// 定义主分析函数
//void analyze_peaks_no_bkg( int startRun = 1, int endRun = 182,const char* expName = "SS032") {
void analyze_peaks_no_bkg( ) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0); // 默认不显示拟合参数框

    // 1. 链接文件
    // !!! 关键：请确保 "tr_PV" 是你 _map.root 文件中 TTree 的真实名字 !!!
    TChain* tr_map = new TChain("tr_map");
    
    // --- 方案A：手动添加文件 ---
    // (已修复：从 chain.Add 改为 tr_chain->Add)
    tr_map->Add("/home/evalie2/Project/document/273Ds/inter_map/SS03200075_map.root");
    tr_map->Add("/home/evalie2/Project/document/273Ds/inter_map/SS03200080_map.root");
    tr_map->Add("/home/evalie2/Project/document/273Ds/inter_map/SS03200115_map.root");

    /* // --- 方案B：循环添加文件 (如果需要) ---
    const char* dataPath = "/home/evalie2/Project/document/273Ds/inter_map/";
    for (int i = startRun; i <= endRun; ++i) {
        tr_chain->Add(TString::Format("%s%s%05d_map.root", dataPath, expName, i));
    }
    */

    // 检查文件是否成功添加
    if (tr_map->GetNtrees() == 0) {
        std::cerr << "错误：没有文件被添加到TChain中。" << std::endl;
        std::cerr << "请检查 TChain 的名字 (当前是 'tr_map') 和文件路径是否正确。" << std::endl;
        return;
    }
    std::cout << "文件添加完成。总事件数: " << tr_map->GetEntries() << std::endl;
    
    // 2. "在内存中"绘制原始一维谱
    tr_map->Draw("DSSDX_E>>h1(1000, 5000, 10000)", "DSSDX_mul == 1 && DSSDY_mul == 1 && MWPC_mul == 0 && Veto_mul == 0", "goff");
    TH1F* h1 = (TH1F*)gROOT->FindObject("h1");
    if (!h1) { 
        std::cerr << "错误: 找不到直方图 h1！" << std::endl;
        std::cerr << "TChain::Draw 命令可能失败 (TTree是空的或分支名'DSSD_E'/'Delta_Ts'错误)。" << std::endl; 
        return; 
    }
    h1->SetTitle("Original Spectrum with Peak Fits");
    h1->GetXaxis()->SetTitle("Energy (channel)");
    h1->GetYaxis()->SetTitle("Counts");

    // 3. 【跳过本底扣除】
    // (本底扣除、克隆 h1_sub、绘制对比图 c_comparison 的相关代码已全部移除)
    TSpectrum *s = new TSpectrum();

    // 4. 生成图：寻峰、拟合并标注 (直接在 h1 上操作)
    TCanvas *c_peaks = new TCanvas("c_peaks", "Peak Search and Fit Results", 1200, 800);
    h1->Draw("hist"); // 直接绘制原始谱

    // 5. 寻峰 (在 h1 上)
    // (参数 2 = sigma, 0.05 = threshold)
    Int_t nfound = s->Search(h1, 2, "goff", 0.05); 
    Double_t *xpeaks = s->GetPositionX();
    
    // 6. 创建并打开一个用于输出的文本文件
    std::ofstream outfile("peak_analysis_results_noBKG.txt");
    // 写入表头，使用制表符(\t)分隔
    outfile << "Peak_ID\tEnergy_Mean\tEnergy_Error\tSigma\tSigma_Error\tFWHM\tFWHM_Error\tBkg_Level" << std::endl;

    std::cout << "\n**************************************************" << std::endl;
    std::cout << "           Peak Analysis Results (No Bkg Sub)" << std::endl;
    std::cout << "**************************************************" << std::endl;

    // 7. 对找到的每个峰进行高斯拟合和标注 (在 h1 上)
    for (int p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = h1->GetXaxis()->FindBin(xp);
        Float_t yp = h1->GetBinContent(bin);
        Double_t fit_range = 30.0;
        Double_t fit_min = xp - fit_range;
        Double_t fit_max = xp + fit_range;

        // 【修改】使用 "gaus + pol0" 拟合
        // (使用 TString::Format 确保每个函数名唯一)
        TF1 *fitFunc = new TF1(TString::Format("fitFunc_%d", p), "gaus(0) + pol0(3)", fit_min, fit_max);
        
        // 【修改】设置初始参数 (4个参数)
        // gaus(0) -> 参数 0, 1, 2 (Constant, Mean, Sigma)
        // pol0(3) -> 参数 3 (Background)
        
        // 估计本底值 (取拟合范围最左边的bin)
        Double_t bg_estimate = h1->GetBinContent(h1->GetXaxis()->FindBin(fit_min));
        // 估计峰高 (总高度 - 本底)
        Double_t peak_height_estimate = yp - bg_estimate;
        if (peak_height_estimate < 0) peak_height_estimate = yp; // 安全检查

        fitFunc->SetParameters(peak_height_estimate, xp, 5, bg_estimate);
        fitFunc->SetParNames("Constant", "Mean", "Sigma", "Background");
        
        // 执行拟合 (在 h1 上)
        h1->Fit(fitFunc, "QR+"); // "R" = 在范围内拟合, "Q" = 安静模式, "+" = 添加到直方图

        // 提取高斯部分的参数 (本底参数是 [3])
        Double_t fit_const = fitFunc->GetParameter(0);
        Double_t fit_mean = fitFunc->GetParameter(1);
        Double_t fit_mean_err = fitFunc->GetParError(1);
        Double_t fit_sigma = fitFunc->GetParameter(2);
        Double_t fit_sigma_err = fitFunc->GetParError(2);
        Double_t fit_bkg = fitFunc->GetParameter(3); // 提取本底参数
        Double_t fwhm = 2.355 * fit_sigma;
        Double_t fwhm_err = 2.355 * fit_sigma_err;
        
        // --- 在终端打印详细信息 ---
        std::cout << "--- Peak " << std::setw(2) << p + 1 << " ---" << std::endl;
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "   - Fit Mean (Energy) : " << std::setw(8) << fit_mean << " +/- " << fit_mean_err << std::endl;
        std::cout << "  - Fit Sigma (Width) : " << std::setw(8) << fit_sigma << " +/- " << fit_sigma_err << std::endl;
        std::cout << "  - FWHM              : " << std::setw(8) << fwhm << " +/- " << fwhm_err << std::endl;
        std::cout << "  - Background Level  : " << std::setw(8) << fit_bkg << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        // --- 将结果写入文件 ---
        outfile << std::fixed << std::setprecision(4); // 文件中保存更高精度
        outfile << p + 1 << "\t"
                << fit_mean << "\t"
                << fit_mean_err << "\t"
                << fit_sigma << "\t"
                << fit_sigma_err << "\t"
                << fwhm << "\t"
                << fwhm_err << "\t"
                << fit_bkg << std::endl;

        // (绘制引出线和标签的代码保持不变，如果你需要它们，取消注释即可)
        /* Double_t callout_length = fit_const * 0.2 + 200;
        const Double_t max_callout_length = 1500.0;
        if (callout_length > max_callout_length) {
            callout_length = max_callout_length;
        }
        Double_t line_x0 = fit_mean;
        Double_t line_y0 = fit_const + fit_bkg; // 峰的顶点 = 峰高 + 本底
        Double_t line_y1 = line_y0 + callout_length;
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

    // 8. 关闭文件并给出提示
    outfile.close();
    std::cout << "\n分析结果已成功保存到 peak_analysis_results_noBKG.txt 文件中。" << std::endl;
    std::cout << "**************************************************" << std::endl;
}