#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <TString.h>

// ---------------------------------------------------------------
// 改进版 Cal error 分析程序 (v2.0 - 健壮版)
// 修正：
//   1. 忽略背景扣除后负 bin
//   2. 输出拟合质心与积分质心
//   3. 同时计算 signed / abs / RMS ΔE
//   4. 输出更详细的调试信息
//   5.用于计算刻度误差
// ---------------------------------------------------------------

void analyze_and_find_peaks(int startRun = 1, int endRun = 182, const char* expName = "SS032") {
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // ------------------ 输入配置 ------------------
    const char* dataPath = "/home/evalie2/Project/document/273Ds/inter2/"; // 假设的路径，请根据实际情况修改
    TString tree_name = "tr_chain"; // 与原版文件保持一致
    TString branch_name = "DSSD_E";

    // ------------------ 构建 TChain ------------------
    TChain *tr_chain = new TChain(tree_name);
    for (int i = startRun; i <= endRun; ++i) {
        tr_chain->Add(TString::Format("%s%s%05d_chain.root", dataPath, expName, i));
    }

    if (tr_chain->GetNtrees() == 0) {
        std::cerr << "错误：未找到可用的 ROOT 文件。" << std::endl;
        delete tr_chain;
        return;
    }
    std::cout << "文件加载成功，总事件数: " << tr_chain->GetEntries() << std::endl;

    // ------------------ 参考能量表 (keV) ------------------
    // 注意：这里的参考能量列表沿用了原 fixed.cpp 中的列表，与 plot_Calerro.cpp 略有不同。
    std::vector<Double_t> ref_energies = {
        5304,6118,6341,6558,6622,7066,7129,7312,7386,7686,8026,8089,8784
    };
    
    // ------------------ 创建直方图 ------------------
    TH1D *h1 = new TH1D("h1", "Energy Spectrum", 4000, 5000, 10000);
    // 使用原版中的 Delta_Ts>0 条件进行绘图
    tr_chain->Draw(Form("%s >> h1", branch_name.Data()), "Delta_Ts>0", "goff");
    
    if (h1->GetEntries() == 0) {
        std::cerr << "错误：直方图 h1 中没有数据点，请检查输入文件路径和分支名。" << std::endl;
        delete h1; delete tr_chain;
        return;
    }

    // 背景拟合与扣除
    TSpectrum* s = new TSpectrum();
    TH1* h_bkg = s->Background(h1, 20, ""); // 使用 TSpectrum 扣背景
    TH1D *h1_sub = (TH1D *)h1->Clone("h1_sub");
    h1_sub->Add(h_bkg, -1.0);

    // ------------------ 寻峰 ------------------
    TSpectrum *sp = new TSpectrum(50);
    int nfound = sp->Search(h1_sub, 2, "goff", 0.05);
    std::cout << "共找到 " << nfound << " 个峰.\n";

    Double_t *xpeaks = sp->GetPositionX();

    // ------------------ ΔE 统计变量 ------------------
    Double_t total_Delta_E_Cal = 0.0;
    Double_t total_abs_Delta = 0.0;
    Double_t sum_sq_Delta = 0.0;
    Int_t matched_peaks_count = 0;
    
    // ------------------ 主循环 ------------------
    for (int p = 0; p < nfound; ++p) {
        Double_t peak_pos = xpeaks[p];
        Int_t bin = h1_sub->GetXaxis()->FindBin(peak_pos);
        Double_t peak_height = h1_sub->GetBinContent(bin);

        // 高斯拟合窗口
        Double_t fit_range = 30.0;
        Double_t fit_min = peak_pos - fit_range;
        Double_t fit_max = peak_pos + fit_range;
        
        // 确保拟合函数在循环内部创建，避免重名
        TF1 *fitf = new TF1(Form("fitf_%d", p), "gaus", fit_min, fit_max);
        fitf->SetParameters(peak_height, peak_pos, 10);
        h1_sub->Fit(fitf, "QNR");

        Double_t mean_fit = fitf->GetParameter(1);
        Double_t sigma = fitf->GetParameter(2);

        // 积分计算质心（忽略负 bin），窗口 ±2σ
        Int_t bin_min = h1_sub->GetXaxis()->FindBin(mean_fit - 2 * sigma);
        Int_t bin_max = h1_sub->GetXaxis()->FindBin(mean_fit + 2 * sigma);
        Double_t sum_energy = 0.0, sum_counts = 0.0;
        for (int b = bin_min; b <= bin_max; ++b) {
            Double_t E = h1_sub->GetXaxis()->GetBinCenter(b);
            Double_t counts = h1_sub->GetBinContent(b);
            if (counts < 0) counts = 0; // 忽略负 bin
            sum_energy += E * counts;
            sum_counts += counts;
        }
        Double_t centroid_int = (sum_counts > 0) ? (sum_energy / sum_counts) : mean_fit;

        // 匹配参考能量
        Double_t best_ref = -1;
        Double_t min_diff = 1e9;
        for (auto ref : ref_energies) {
            Double_t diff = fabs(mean_fit - ref);
            if (diff < min_diff) {
                min_diff = diff;
                best_ref = ref;
            }
        }

        Double_t match_threshold = 10.0;
        if (min_diff < match_threshold) {
            Double_t dE_fit = mean_fit - best_ref;
            Double_t dE_int = centroid_int - best_ref;
            
            // 统计基于拟合质心的 ΔE
            total_Delta_E_Cal += dE_fit;
            total_abs_Delta += fabs(dE_fit);
            sum_sq_Delta += dE_fit * dE_fit;
            matched_peaks_count++;

            std::cout << std::fixed << std::setprecision(2);
            std::cout << "\n峰 #" << std::setw(2) << p + 1
                      << "   Ref = " << best_ref << " keV"
                      << "   FitMean = " << mean_fit
                      << "   IntCentroid = " << centroid_int
                      << "   ΔE_fit = " << dE_fit
                      << "   ΔE_int = " << dE_int
                      << std::endl;
        }
        delete fitf; // 清理拟合函数
    }

    // ------------------ 总结 ------------------
    if (matched_peaks_count > 0) {
        Double_t mean_signed = total_Delta_E_Cal / matched_peaks_count;
        Double_t mean_abs = total_abs_Delta / matched_peaks_count;
        Double_t rms = sqrt(sum_sq_Delta / matched_peaks_count);

        std::cout << "\n------------------ 结果统计 ------------------\n";
        std::cout << "匹配峰数量: " << matched_peaks_count << std::endl;
        std::cout << "平均 signed ΔE_Cal (Fit) = " << mean_signed << " keV\n";
        std::cout << "平均 |ΔE_Cal| (Fit) = " << mean_abs << " keV\n";
        std::cout << "RMS(ΔE_Cal) (Fit) = " << rms << " keV\n";
        std::cout << "---------------------------------------------\n";
    } else {
        std::cout << "⚠️ 未匹配到任何参考峰！\n";
    }

    // 输出图像
    TCanvas *c1 = new TCanvas("c1", "Cal Error", 1000, 700);
    h1_sub->Draw();
    // 增加一个标记，表明这是由哪个范围生成的结果
    TString output_name = TString::Format("Cal_error_%s_%d_to_%d.png", expName, startRun, endRun);
    c1->SaveAs(output_name);

    // 清理 ROOT 对象
    delete c1;
    delete h1_sub;
    delete h_bkg;
    delete s;
    delete h1;
    delete tr_chain;
    delete sp;
}

int main() {
    // 默认调用方式，与原版 main 函数一致
    analyze_and_find_peaks(1, 182, "SS032");
    
    // 示例：可以调用其他范围
    // analyze_and_find_peaks(183, 200, "SS032");
    
    return 0;
}