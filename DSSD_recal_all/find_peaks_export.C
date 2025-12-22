#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPaveText.h"

void find_peaks_export(){
    // ================== 参数配置 ==================
    const char* fileName  = "Diagnose_Calibration.root";
    const char* histName  = "h_calib_X";
    const char* outTxt    = "peaks_sigma_output.txt"; // 输出文件名

    Double_t sigma_est = 4.0;    // 寻峰时的预估 Sigma
    Double_t threshold = 0.05;   // 阈值 (1%)
    Double_t fitWindow = 3.0;    // 拟合窗口: +/- 3 * sigma
    // =============================================

    // 1. 读取数据
    TFile *f = TFile::Open(fileName);
    if (!f || f->IsZombie()) { printf("Error: File not found\n"); return; }
    TH1D *h1 = (TH1D*)f->Get(histName);
    if (!h1) { printf("Error: Hist not found\n"); return; }

    // 2. 准备画布
    TCanvas *c1 = new TCanvas("c1", "Fit Sigma Output", 1400, 800);
    h1->Draw();
    gStyle->SetOptFit(0); 

    // 3. 寻峰
    TSpectrum *s = new TSpectrum(50);
    Int_t nfound = s->Search(h1, sigma_est, "nobackground", threshold);

    Double_t *xpeaks = s->GetPositionX();
    std::vector<Double_t> peaks_sorted(xpeaks, xpeaks + nfound);
    std::sort(peaks_sorted.begin(), peaks_sorted.end());

    std::cout << "\n>>> Found " << nfound << " peaks. Fitting..." << std::endl;

    // 4. 导出文件准备 (添加了 Sigma 和 Sigma_Err 列)
    std::ofstream outFile(outTxt);
    outFile << "ID\tMean\tSigma\tSigma_Err\tFWHM\tChi2_NDF" << std::endl;
    
    // 终端表头
    printf("%-3s | %-8s | %-8s | %-8s | %-8s | %-6s\n", 
           "ID", "Mean", "Sigma", "SigErr", "FWHM", "Chi2");
    printf("------------------------------------------------------------\n");

    // 5. 循环拟合
    for (int i = 0; i < nfound; i++) {
        Double_t xp = peaks_sorted[i];
        Double_t xmin = xp - fitWindow * sigma_est;
        Double_t xmax = xp + fitWindow * sigma_est;

        // 定义函数: gaus(0) + pol1(3) -> [0]=H, [1]=Mean, [2]=Sigma, [3]=b, [4]=k
        TF1 *ffit = new TF1(Form("fit_%d", i), "gaus(0) + pol1(3)", xmin, xmax);
        
        // --- 初始化参数 ---
        Double_t y_peak = h1->GetBinContent(h1->GetXaxis()->FindBin(xp));
        Double_t y_min  = h1->GetBinContent(h1->GetXaxis()->FindBin(xmin));
        Double_t y_max  = h1->GetBinContent(h1->GetXaxis()->FindBin(xmax));
        Double_t slope = (y_max - y_min) / (xmax - xmin);
        Double_t intercept = y_min - slope * xmin;
        Double_t height = y_peak - (intercept + slope * xp);
        if (height < 0) height = y_peak * 0.2;

        ffit->SetParameters(height, xp, sigma_est, intercept, slope);
        ffit->SetParLimits(2, 0.1, sigma_est * 4.0); // 限制 Sigma 范围
        ffit->SetParLimits(1, xmin, xmax);           // 限制 Mean 范围

        // 绘图设置
        ffit->SetLineColor(kRed);
        ffit->SetLineWidth(2);

        // 执行拟合
        h1->Fit(ffit, "RQ+"); 

        // --- 获取结果 (重点在这里) ---
        Double_t mean_val  = ffit->GetParameter(1);
        
        Double_t sigma_val = TMath::Abs(ffit->GetParameter(2)); // 获取 Sigma
        Double_t sigma_err = ffit->GetParError(2);              // 获取 Sigma 的误差
        
        Double_t fwhm_val  = 2.35482 * sigma_val;
        Double_t chi2      = ffit->GetChisquare() / (ffit->GetNDF()>0 ? ffit->GetNDF() : 1);

        // --- 输出到终端 ---
        printf("%3d | %8.2f | %8.2f | %8.2f | %8.2f | %6.2f\n", 
               i+1, mean_val, sigma_val, sigma_err, fwhm_val, chi2);

        // --- 输出到TXT ---
        outFile << i+1 << "\t" 
                << mean_val << "\t" 
                << sigma_val << "\t" 
                << sigma_err << "\t" 
                << fwhm_val << "\t" 
                << chi2 << std::endl;
    }
    
    outFile.close();
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "Results saved to: " << outTxt << std::endl;
    
    c1->Update();
}