#include "TCanvas.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TMath.h"
#include "TCut.h"
#include "TBox.h"

//autoFitDecay(&tree, 6500, 1000, 1200, 500, 1);
//1为输出Latex标签，0为不输出
//绘制二维能谱的局部放大图，并对母核和子核的半衰期进行Schmidt函数拟合
// ==========================================================
// 辅助工具: 智能单位转换与格式化
// ==========================================================
TString FormatResult(double val, double err) {
    if (val <= 0) return "0";
    
    // 定义单位阈值和除数
    const char* units[] = {"ns", "#mus", "ms", "s"};
    double divisors[] = {1.0, 1.0e3, 1.0e6, 1.0e9};
    int idx = 0;

    // 自动判断量级
    if (val >= 1.0e9)      idx = 3; // > 1s
    else if (val >= 1.0e6) idx = 2; // > 1ms
    else if (val >= 1.0e3) idx = 1; // > 1us
    else                   idx = 0; // ns

    // 执行单位转换
    double val_scaled = val / divisors[idx];
    double err_scaled = err / divisors[idx];
    
    // 返回格式化字符串: Value +/- Error Unit
    return Form("%.2f #pm %.2f %s", val_scaled, err_scaled, units[idx]);
}

// ==========================================================
// Schmidt 函数: 用于拟合 ln(t) 谱
// ==========================================================
Double_t SchmidtFunc(Double_t *x, Double_t *par) {
    Double_t t_log = x[0]; 
    Double_t norm = par[0];
    Double_t mu = par[1];  // ln(tau)
    Double_t bg = par[2];  // 常数本底
    
    return norm * TMath::Exp( (t_log - mu) - TMath::Exp(t_log - mu) ) + bg;
}

// ==========================================================
// 主函数
// 参数修改说明:
// int showLabel: 1 代表显示 (True)，0 代表不显示 (False)。默认值为 1。
// ==========================================================
void autoFitDecay(TTree *tr, double e1, double t1_guess, 
                    double e2, double t2_guess, 
                    int showLabel = 1) 
{
    
    if (!tr) {
        printf("Error: Tree/Chain is null!\n");
        return;
    }
    
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0); // 隐藏默认统计框
    
    // ================= 配置参数 =================
    double view_range = 150.0; // 二维图绘图范围
    double cut_range  = 30.0;  // 实际开窗范围
    
    TCut cut_physics = ""; 
    
    // 能量截断
    TCut cut_energy = Form("abs(DSSD_E[1]-%f)<%f && abs(DSSD_E[2]-%f)<%f", 
                           e1, cut_range, e2, cut_range);
    
    TCut final_cut = cut_energy && cut_physics;

    // ================= 绘图部分 =================
    TCanvas *c1 = new TCanvas("c_decay", "Decay Analysis", 1200, 450);
    c1->Divide(3, 1);

    // --- Pad 1: 能量关联图 ---
    c1->cd(1);
    gPad->SetRightMargin(0.12); gPad->SetLeftMargin(0.12);
    
    TH2F *h2 = new TH2F("h_check", "Energy Selection;Parent E (keV);Daughter E (keV)", 
                        100, e1 - view_range, e1 + view_range, 
                        100, e2 - view_range, e2 + view_range);
    
    tr->Draw("DSSD_E[2]:DSSD_E[1]>>h_check", cut_physics, "colz");
    
    // 绘制红框
    TBox *box = new TBox(e1 - cut_range, e2 - cut_range, e1 + cut_range, e2 + cut_range);
    box->SetFillStyle(0); box->SetLineColor(kRed); box->SetLineWidth(3); 
    box->Draw("same");
    
    // [控制点 1] 是否显示红框说明文字
    // C++ 逻辑: if(1) 为真，if(0) 为假
    if (showLabel) {
        TLatex lat_info; lat_info.SetNDC(); lat_info.SetTextSize(0.04);
        lat_info.DrawLatex(0.15, 0.85, Form("#color[2]{Box: #pm %.0f keV}", cut_range));
    }

    // --- Pad 2: 母核拟合 ---
    c1->cd(2);
    double ln_min1 = TMath::Log(t1_guess * 0.05);
    double ln_max1 = TMath::Log(t1_guess * 20.0);
    
    TH1F *h1 = new TH1F("h1", "Parent Decay;ln(t / ns);Counts", 50, ln_min1, ln_max1);
    tr->Draw("log(Delta_Ts[1])>>h1", final_cut && "Delta_Ts[1]>0");
    
    if(h1->GetEntries() > 0) {
        TF1 *f1 = new TF1("f1", SchmidtFunc, ln_min1, ln_max1, 3);
        f1->SetParameters(h1->GetMaximum(), TMath::Log(t1_guess/0.693), 0.1);
        f1->SetParLimits(2, 0, 100); 
        f1->SetLineColor(kRed);
        
        h1->Fit("f1", "L Q"); 
        
        // [控制点 2] 是否显示母核半衰期
        if (showLabel) {
            double tau = TMath::Exp(f1->GetParameter(1));
            double t_half = 0.693147 * tau;
            double t_err = t_half * f1->GetParError(1); 
            
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.05); lat.SetTextColor(kRed);
            lat.DrawLatex(0.13, 0.8, Form("T_{1/2} = %s", FormatResult(t_half, t_err).Data()));
        }
    }

    // --- Pad 3: 子核拟合 ---
    c1->cd(3);
    double ln_min2 = TMath::Log(t2_guess * 0.05);
    double ln_max2 = TMath::Log(t2_guess * 20.0);
    
    TH1F *h2t = new TH1F("h2t", "Daughter Decay;ln(t / ns);Counts", 50, ln_min2, ln_max2);
    tr->Draw("log(Delta_Ts[2])>>h2t", final_cut && "Delta_Ts[2]>0");
    
    if(h2t->GetEntries() > 0) {
        TF1 *f2 = new TF1("f2", SchmidtFunc, ln_min2, ln_max2, 3);
        f2->SetParameters(h2t->GetMaximum(), TMath::Log(t2_guess/0.693), 0.1);
        f2->SetParLimits(2, 0, 100);
        f2->SetLineColor(kRed);
        
        h2t->Fit("f2", "L Q");
        
        // [控制点 3] 是否显示子核半衰期
        if (showLabel) {
            double tau = TMath::Exp(f2->GetParameter(1));
            double t_half = 0.693147 * tau;
            double t_err = t_half * f2->GetParError(1);
            
            TLatex lat; lat.SetNDC(); lat.SetTextSize(0.05); lat.SetTextColor(kRed);
            lat.DrawLatex(0.13, 0.8, Form("T_{1/2} = %s", FormatResult(t_half, t_err).Data()));
        }
    }
    
    c1->Update();
}