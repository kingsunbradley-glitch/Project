/******************************************************************************
 * *
 * 一个用于分析多个ROOT文件的ROOT宏。                                           *
 * 功能：                                                                    *
 * 1. ... (同前) ...                                                         *
 * 5. 将详细的峰分析参数（能量，FWHM, N_Events等）输出到 .txt 文件。          *
 * 6. 计算每个峰与给定参考能量的偏差 (Delta_E_Cal)，并计算总平均偏差。        *
 * 7. 【新】计算 Delta_E_Cal 值的统计方差 (Variance) 和 标准差 (StdDev)。     *
 * *
 ******************************************************************************/

#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>     // **【新增点 1】**：包含 vector
#include <cmath>      // **【新增点 2】**：包含 cmath
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

    // 【新增点 3：定义参考能量和匹配阈值】
    std::vector<Double_t> ref_energies = {
        5304,6118,6341,6558,6622,7066,7129,7312,7386,7686,8026,8089,8784
    }; // 示例值，请替换为您实际的参考值
    const Double_t match_threshold =15 ; // 示例值，请按需调整

    // ... (链接文件、画图、扣本底的代码保持不变) ...
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

    // 2. 绘制原始谱 (代码不变)
    tr_chain->Draw("DSSD_E>>h1(1000, 5000, 10000)", "Delta_Ts>0", "goff");
    TH1F* h1 = (TH1F*)gROOT->FindObject("h1");
    if (!h1) { std::cerr << "错误: 找不到直方图 h1！" << std::endl; return; }
    h1->SetTitle("Original Spectrum & Background"); h1->GetXaxis()->SetTitle("Energy (channel)"); h1->GetYaxis()->SetTitle("Counts");

    // 3. 扣本底 (代码不变)
    TSpectrum *s = new TSpectrum();
    TH1 *h_bkg = s->Background(h1, 20, "");
    TH1F *h1_sub = (TH1F*)h1->Clone("h1_sub");
    h1_sub->Add(h_bkg, -1);
    h1_sub->SetTitle("Peaks on Background-Subtracted Spectrum"); h1->GetXaxis()->SetTitle("Energy (channel)"); h1->GetYaxis()->SetTitle("Net Counts");

    // 4. 生成图一 (代码不变)
    TCanvas *c_comparison = new TCanvas("c_comparison", "Background Subtraction Comparison", 900, 600);
    h1->SetLineColor(kBlack); h_bkg->SetLineColor(kRed); h1_sub->SetLineColor(kBlue);
    h1->Draw("hist"); h_bkg->Draw("same hist"); h1_sub->Draw("same hist");
    TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "Original Spectrum", "l"); leg->AddEntry(h_bkg, "Estimated Background", "l"); leg->AddEntry(h1_sub, "Spectrum after Subtraction", "l");
    leg->Draw();
    
    // 5. 生成图二 (代码不变)
    TCanvas *c_peaks = new TCanvas("c_peaks", "Peak Search and Fit Results", 1200, 800);
    h1_sub->Draw("hist");

    // 6. 寻峰 (代码不变)
    Int_t nfound = s->Search(h1_sub, 2, "goff", 0.05);
    Double_t *xpeaks = s->GetPositionX();
    
    // 打开输出文件 (修改点 1：表头不变)
    std::ofstream outfile("peak_analysis_results.txt");
    outfile << "Peak_ID\tEnergy_Mean\tEnergy_Error\tSigma\tSigma_Error\tFWHM\tFWHM_Error\tN_Events\tDelta_E_Cal\tReference_Value" << std::endl;

    std::cout << "\n**************************************************" << std::endl;
    std::cout << "           Peak Analysis Results" << std::endl;
    std::cout << "**************************************************" << std::endl;

    // vvvvvvvvvv 【新增点 4：初始化总偏差计算器 和 存储器】 vvvvvvvvvv
    Double_t total_Delta_E_Cal = 0.0;
    Int_t matched_peaks_count = 0;
    std::vector<Double_t> delta_e_cal_values; // 用于存储所有的 Delta_E_Cal
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // 7. 对找到的每个峰进行高斯拟合和标注
    for (int p = 0; p < nfound; p++) {
        // ... (拟合代码不变) ...
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

        // ... (计算峰内事件数 n (代码不变)) ...
        Double_t integration_min = fit_mean - 3.0 * fit_sigma;
        Double_t integration_max = fit_mean + 3.0 * fit_sigma;
        Int_t bin_min = h1_sub->GetXaxis()->FindBin(integration_min);
        Int_t bin_max = h1_sub->GetXaxis()->FindBin(integration_max);
        Double_t n_events = h1_sub->Integral(bin_min, bin_max);

        // ... (【新增点 5】查找最接近的参考值 (代码不变)) ...
        Double_t min_diff = 1e18;
        Double_t best_A_ref = -1.0;
        for (size_t i = 0; i < ref_energies.size(); ++i) {
            Double_t diff = std::fabs(fit_mean - ref_energies[i]);
            if (diff < min_diff) {
                min_diff = diff;
                best_A_ref = ref_energies[i];
            }
        }
        
        Double_t current_Delta_E_Cal = 0.0;
        TString ref_val_str = "None";

        if (min_diff < match_threshold) {
            current_Delta_E_Cal = fit_mean - best_A_ref;
            total_Delta_E_Cal += current_Delta_E_Cal;
            matched_peaks_count++;
            ref_val_str.Form("%.2f", best_A_ref);
            
            // vvvvvvvvvv 【修改点 2：存储 Delta_E_Cal 值】 vvvvvvvvvv
            delta_e_cal_values.push_back(current_Delta_E_Cal);
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        }

        // ... (在终端打印 (代码不变)) ...
        std::cout << "--- Peak " << std::setw(2) << p + 1 << " ---" << std::endl;
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  - Fit Mean (Energy) : " << std::setw(8) << fit_mean << " +/- " << fit_mean_err << std::endl;
        //... (省略)
        std::cout << "  - Reference Value   : " << std::setw(8) << ref_val_str.Data() << std::endl;
        std::cout << "  - Delta_E_Cal       : " << std::setw(8) << current_Delta_E_Cal << std::endl;
        std::cout << "------------------------------------------" << std::endl;

        // ... (写入文件 (代码不变)) ...
        outfile << std::fixed << std::setprecision(4);
        outfile << p + 1 << "\t"
                << fit_mean << "\t" << fit_mean_err << "\t"
                << fit_sigma << "\t" << fit_sigma_err << "\t"
                << fwhm << "\t" << fwhm_err << "\t"
                << n_events << "\t"
                << current_Delta_E_Cal << "\t"
                << ref_val_str.Data() << std::endl;

        // ... (绘制引出线和标签的代码保持不变) ...
    }

    // vvvvvvvvvv 【修改点 3：计算并写入总结（平均值、方差、标准差）】 vvvvvvvvvv
    outfile << "\n# --- Summary ---" << std::endl;
    std::cout << "\n------------------------------------------" << std::endl;
    std::cout << "            Summary Statistics" << std::endl;
    std::cout << "------------------------------------------" << std::endl;

    if (matched_peaks_count > 0) {
        // 1. 计算平均值
        Double_t avg_Delta_E_Cal = total_Delta_E_Cal / matched_peaks_count;
        
        // 2. 计算方差 (Variance) 和 标准差 (Standard Deviation)
        Double_t sum_sq_diff = 0.0;
        for (size_t i = 0; i < delta_e_cal_values.size(); ++i) {
            sum_sq_diff += std::pow(delta_e_cal_values[i] - avg_Delta_E_Cal, 2);
        }
        
        Double_t variance_delta_e = 0.0;
        Double_t stddev_delta_e = 0.0;

        // (样本方差 N-1，在N>1时才有意义)
        if (matched_peaks_count > 1) { 
            variance_delta_e = sum_sq_diff / (matched_peaks_count - 1);
            stddev_delta_e = std::sqrt(variance_delta_e);
        }

        // 3. 写入文件
        outfile << std::fixed << std::setprecision(6);
        outfile << "Matched_Peaks_Count\t" << matched_peaks_count << std::endl;
        outfile << "Average_Delta_E_Cal\t" << avg_Delta_E_Cal << std::endl;
        outfile << "Variance_Delta_E_Cal\t" << variance_delta_e << std::endl;
        outfile << "StdDev_Delta_E_Cal\t" << stddev_delta_e << std::endl;

        // 4. 打印到终端
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "  Matched Peaks Count : " << matched_peaks_count << std::endl;
        std::cout << "  Average_Delta_E_Cal : " << std::setw(10) << avg_Delta_E_Cal << " (平均偏差)" << std::endl;
        std::cout << "  Variance_Delta_E_Cal: " << std::setw(10) << variance_delta_e << " (方差)" << std::endl;
        std::cout << "  StdDev_Delta_E_Cal  : " << std::setw(10) << stddev_delta_e << " (标准差)" << std::endl;

    } else {
        outfile << "Matched_Peaks_Count\t0" << std::endl;
        outfile << "Average_Delta_E_Cal\tNaN" << std::endl;
        outfile << "Variance_Delta_E_Cal\tNaN" << std::endl;
        outfile << "StdDev_Delta_E_Cal\tNaN" << std::endl;
        
        std::cout << "  没有找到与参考值匹配的峰。" << std::endl;
    }
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    // 关闭文件
    outfile.close();
    std::cout << "\n分析结果已成功保存到 peak_analysis_results.txt 文件中。" << std::endl;
    std::cout << "**************************************************" << std::endl;
}