/******************************************************************************
 * *
 * 高效版：多个ROOT文件的能谱分析 + ΔE_Cal计算（基于直方图积分法）            *
 * 功能：                                                                    *
 * 1. 自动扣背景并寻峰                                                       *
 * 2. 高斯拟合提取参数（Mean, Sigma, FWHM, N_Events）                        *
 * 3. 与参考能量匹配，计算每个峰的 ΔE_Cal                                    *
 * 4. 输出每峰结果与总平均 ΔE_Cal                                            *
 * *
 ******************************************************************************/

//g++ plot_Calerro.cpp $(root-config --cflags --libs ) -lSpectrum -o plot_Calerro


#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <TH1.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TF1.h>

void analyze_and_find_peaks(int startRun = 1, int endRun = 182, const char* expName = "SS032") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ------------------- 用户可修改部分 -------------------
    std::vector<Double_t> ref_energies = {
        5304,6118,6341,6558,6622,7066,7129,7312,7386,7686,8026,8089,8784
    }; // 参考能量 (keV)
    const Double_t match_threshold = 10.0; // 匹配阈值 (keV)
    const char* dataPath = "/home/evalie2/Project/document/273Ds/inter2/";
    // ------------------------------------------------------

    // 1. 构建 TChain
    TChain* tr_chain = new TChain("tr_chain");
    for (int i = startRun; i <= endRun; ++i)
        tr_chain->Add(TString::Format("%s%s%05d_chain.root", dataPath, expName, i));

    if (tr_chain->GetNtrees() == 0) {
        std::cerr << "错误：未找到可用的 ROOT 文件。" << std::endl;
        return;
    }
    std::cout << "文件加载成功，总事件数: " << tr_chain->GetEntries() << std::endl;

    // 2. 绘制并获取能量谱
    tr_chain->Draw("DSSD_E>>h1(1000,5000,10000)", "Delta_Ts>0", "goff");
    TH1F* h1 = (TH1F*)gROOT->FindObject("h1");
    if (!h1) { std::cerr << "错误: 找不到直方图 h1！" << std::endl; return; }

    // 3. 扣背景
    TSpectrum* s = new TSpectrum();
    TH1* h_bkg = s->Background(h1, 20, "");
    TH1F* h1_sub = (TH1F*)h1->Clone("h1_sub");
    h1_sub->Add(h_bkg, -1);

    // 4. 绘制比较图
    TCanvas* c1 = new TCanvas("c_comparison", "Background Subtraction", 900, 600);
    h1->SetLineColor(kBlack);
    h_bkg->SetLineColor(kRed);
    h1_sub->SetLineColor(kBlue);
    h1->Draw("hist");
    h_bkg->Draw("same hist");
    h1_sub->Draw("same hist");
    auto leg = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "Original", "l");
    leg->AddEntry(h_bkg, "Background", "l");
    leg->AddEntry(h1_sub, "Subtracted", "l");
    leg->Draw();

    // 5. 寻峰
    Int_t nfound = s->Search(h1_sub, 2, "goff", 0.05);
    Double_t* xpeaks = s->GetPositionX();

    // 6. 输出文件
    std::ofstream outfile("peak_analysis_results.txt");
    outfile << "Peak_ID\tEnergy_Mean\tEnergy_Error\tSigma\tSigma_Error\tFWHM\tFWHM_Error\tN_Events\tDelta_E_Cal\tReference_Value\n";

    Double_t total_Delta_E_Cal = 0.0;
    Int_t matched_peaks_count = 0;

    // 7. 对每个峰进行拟合与 ΔE_Cal 计算
    for (int p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        Int_t bin = h1_sub->GetXaxis()->FindBin(xp);
        Float_t yp = h1_sub->GetBinContent(bin);
        Double_t fit_range = 30.0;
        Double_t fit_min = xp - fit_range;
        Double_t fit_max = xp + fit_range;

        TF1* fitFunc = new TF1("fitFunc", "gaus", fit_min, fit_max);
        fitFunc->SetParameters(yp, xp, 5);
        h1_sub->Fit(fitFunc, "QR+");

        Double_t mean = fitFunc->GetParameter(1);
        Double_t mean_err = fitFunc->GetParError(1);
        Double_t sigma = fitFunc->GetParameter(2);
        Double_t sigma_err = fitFunc->GetParError(2);
        Double_t fwhm = 2.355 * sigma;
        Double_t fwhm_err = 2.355 * sigma_err;

        Double_t integration_min = mean - 3*sigma;
        Double_t integration_max = mean + 3*sigma;
        Int_t bin_min = h1_sub->GetXaxis()->FindBin(integration_min);
        Int_t bin_max = h1_sub->GetXaxis()->FindBin(integration_max);

        Double_t n_events = h1_sub->Integral(bin_min, bin_max);

        // ---------- 计算 ΔE_Cal（基于直方图积分） ----------
        Double_t min_diff = 1e18;
        Double_t best_ref = -1.0;
        for (auto ref : ref_energies) {
            Double_t diff = std::fabs(mean - ref);
            if (diff < min_diff) { min_diff = diff; best_ref = ref; }
        }

        Double_t current_Delta_E_Cal = 0.0;
        TString ref_val_str = "None";

        if (min_diff < match_threshold) {
            Double_t sum_energy = 0.0;
            Double_t sum_counts = 0.0;
            for (int b = bin_min; b <= bin_max; ++b) {
                Double_t E = h1_sub->GetXaxis()->GetBinCenter(b);
                Double_t counts = h1_sub->GetBinContent(b);
                sum_energy += E * counts;
                sum_counts += counts;
            }
            if (sum_counts > 0) {
                Double_t mean_energy = sum_energy / sum_counts;
                current_Delta_E_Cal = mean_energy - best_ref;
                total_Delta_E_Cal += current_Delta_E_Cal;
                matched_peaks_count++;
                ref_val_str.Form("%.2f", best_ref);
            }
        }
        // ----------------------------------------------------

        // 打印输出
        
        /**/
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Peak " << p + 1 << ":\n"
                  << "  Mean = " << mean << " keV\n"
                  << "  Sigma = " << sigma << "\n"
                  << "  N_Events = " << n_events << "\n"
                  << "  Ref = " << ref_val_str << "\n"
                  << "  ΔE_Cal = " << current_Delta_E_Cal << "\n"
                  << "-----------------------------------\n";

        outfile << std::fixed << std::setprecision(4)
                << p + 1 << "\t" << mean << "\t" << mean_err << "\t"
                << sigma << "\t" << sigma_err << "\t"
                << fwhm << "\t" << fwhm_err << "\t"
                << n_events << "\t" << current_Delta_E_Cal << "\t"
                << ref_val_str.Data() << "\n";
    }


    

    // 8. 总平均 ΔE_Cal
    outfile << "\n# --- Summary ---\n";
    if (matched_peaks_count > 0) {
        Double_t avg_Delta_E_Cal = total_Delta_E_Cal / matched_peaks_count;
        outfile << "Average_Delta_E_Cal\t" << avg_Delta_E_Cal << "\n";
        outfile << "Matched_Peaks_Count\t" << matched_peaks_count << "\n";
        std::cout << "\n平均 ΔE_Cal = " << avg_Delta_E_Cal
                  << " (" << matched_peaks_count << " 个匹配峰)\n";
    } else {
        outfile << "Average_Delta_E_Cal\tNaN\nMatched_Peaks_Count\t0\n";
        std::cout << "\n没有匹配到参考能量的峰。\n";
    }

    outfile.close();
    std::cout << "\n结果已保存到 peak_analysis_results.txt\n";
}
int main() {
    analyze_and_find_peaks(1, 182, "SS032");
    return 0;
}