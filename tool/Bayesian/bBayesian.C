#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLegendEntry.h"

// 主分析函数
void bBayesian() {
    // --- 1. 定义实验参数 (单位: ms) ---
    const double t1 = 8.319;       // 第一个事件的衰变时间
    const double t2 = 41.703;      // 第二个事件的衰变时间
    const int n_events = 2;        // 事件总数
    const double T0 = 0.0005;      // 观测窗口开始时间 (0.5 us)
    const double T1 = 100000.0;    // 观测窗口结束时间 (100 s)
    const double LN2 = TMath::Log(2);

    // =================================================================
    // Section A: Bayesian Analysis
    // =================================================================

    // --- 2. 定义后验概率分布函数 P(t_1/2 | data) ---
    auto posterior_t_half_unnormalized = [&](double* x, double* /*p*/) {
        double t_half = x[0];
        if (t_half <= 1e-3) return 0.0;
        double lambda = LN2 / t_half;
        double term_lambda_part = lambda * TMath::Exp(-lambda * (t1 + t2));
        double denominator = TMath::Exp(-lambda * T0) - TMath::Exp(-lambda * T1);
        double jacobian = LN2 / (t_half * t_half);
        return term_lambda_part / (denominator * denominator) * jacobian;
    };

    // --- 3. 创建 TF1 对象并进行数值分析 ---
    double t_min = 1.0;
    double t_max = 150.0;
    TF1 *f_posterior = new TF1("f_posterior", posterior_t_half_unnormalized, t_min, t_max, 0);
    f_posterior->SetNpx(2000);

    // --- 4. 计算贝叶斯分析的关键统计量 ---
    double integral = f_posterior->Integral(t_min, t_max);
    if (integral < 1e-9) {
        std::cout << "Error: PDF integral is zero. Check function or range." << std::endl;
        return;
    }
    auto posterior_t_half_normalized = [&](double* x, double* /*p*/) {
        return posterior_t_half_unnormalized(x, nullptr) / integral;
    };
    TF1 *f_normalized = new TF1("f_normalized", posterior_t_half_normalized, t_min, t_max, 0);
    f_normalized->SetNpx(2000);

    double map_estimate = f_posterior->GetMaximumX(t_min, t_max);
    double quantiles_needed[] = {0.1585, 0.5, 0.8415};
    double values[3];
    f_normalized->GetQuantiles(3, values, quantiles_needed);
    
    double bayes_lower_bound_68 = values[0];
    double bayes_median_estimate = values[1];
    double bayes_upper_bound_68 = values[2];
    
    double bayes_error_high = bayes_upper_bound_68 - bayes_median_estimate;
    double bayes_error_low = bayes_median_estimate - bayes_lower_bound_68;

    // =================================================================
    // Section B: MLE Analysis (Schmidt et al., 1984)
    // =================================================================

    // --- 4b. 计算 MLE 分析的关键统计量 ---
    double t_mean = (t1 + t2) / n_events; // 平均寿命 MLE
    
    // 从论文 Table 1, n=2, Exponential distr. 中获取系数
    const double factor_lower = 0.606;
    const double factor_upper = 2.82;

    // 计算平均寿命的置信区间
    double tau_lower = t_mean * factor_lower;
    double tau_upper = t_mean * factor_upper;

    // 转换为半衰期
    double mle_t_half_best = t_mean * LN2;
    double mle_lower_bound_68 = tau_lower * LN2;
    double mle_upper_bound_68 = tau_upper * LN2;

    // 计算不对称误差
    double mle_error_high = mle_upper_bound_68 - mle_t_half_best;
    double mle_error_low = mle_t_half_best - mle_lower_bound_68;


    // --- 5. 打印所有结果 ---
    printf("\n\n");
    std::cout << "######################################################" << std::endl;
    std::cout << "#             Combined Analysis Results              #" << std::endl;
    std::cout << "######################################################" << std::endl;
    
    std::cout << "\n=====================================================" << std::endl;
    std::cout << "           Bayesian Analysis (Jeffreys Prior)        " << std::endl;
    std::cout << "=====================================================" << std::endl;
    printf("Best Estimates:\n");
    printf("  - Mode (MAP)                 : %.2f ms (Probability peak)\n", map_estimate);
    printf("  - Median                     : %.2f ms (Probability center)\n", bayes_median_estimate);
    std::cout << "-----------------------------------------------------" << std::endl;
    printf("68.3%% Credible Interval (1-sigma):\n");
    printf("  - Range                      : [%.2f ms, %.2f ms]\n", bayes_lower_bound_68, bayes_upper_bound_68);
    std::cout << "-----------------------------------------------------" << std::endl;
    printf("Final Result in Standard Format (Median-based):\n");
    printf("  t_1/2 = %.2f (+%.2f, -%.2f) ms\n", bayes_median_estimate, bayes_error_high, bayes_error_low);
    std::cout << "=====================================================" << std::endl;

    std::cout << "\n=====================================================" << std::endl;
    std::cout << "      MLE Analysis (Schmidt et al., Z. Phys. A 316)    " << std::endl;
    std::cout << "=====================================================" << std::endl;
    printf("Best Estimate (from mean lifetime):\n");
    printf("  - t_1/2                      : %.2f ms\n", mle_t_half_best);
    std::cout << "-----------------------------------------------------" << std::endl;
    printf("68.3%% Confidence Interval (from Table 1, n=2):\n");
    printf("  - Range                      : [%.2f ms, %.2f ms]\n", mle_lower_bound_68, mle_upper_bound_68);
    std::cout << "-----------------------------------------------------" << std::endl;
    printf("Final Result in Standard Format (MLE-based):\n");
    printf("  t_1/2 = %.2f (+%.2f, -%.2f) ms\n", mle_t_half_best, mle_error_high, mle_error_low);
    std::cout << "=====================================================" << std::endl;


    // --- 6. 可视化 (仅展示贝叶斯后验分布) ---
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1", "Bayesian Posterior PDF for Half-life", 900, 700);
    c1->SetGrid();
    c1->SetLeftMargin(0.12);

    f_posterior->SetTitle("Posterior PDF for Half-life;Half-life t_{1/2} (ms);Probability Density");
    f_posterior->GetXaxis()->SetTitleSize(0.04);
    f_posterior->GetYaxis()->SetTitleSize(0.04);
    f_posterior->GetXaxis()->SetTitleOffset(1.1);
    f_posterior->GetYaxis()->SetTitleOffset(1.3);
    f_posterior->SetLineColor(kAzure+2);
    f_posterior->SetLineWidth(3);
    f_posterior->Draw("L");

    // ... (绘图部分的代码与之前相同) ...
    TLine line;
    line.SetLineWidth(2);
    line.SetLineColor(kRed);
    line.SetLineStyle(2);
    line.DrawLine(map_estimate, 0, map_estimate, f_posterior->Eval(map_estimate));
    line.SetLineColor(kGreen+2);
    line.SetLineStyle(7);
    line.DrawLine(bayes_median_estimate, 0, bayes_median_estimate, f_posterior->Eval(bayes_median_estimate));
    line.SetLineColor(kGray+2);
    line.SetLineStyle(3);
    line.DrawLine(bayes_lower_bound_68, 0, bayes_lower_bound_68, f_posterior->Eval(bayes_lower_bound_68));
    line.DrawLine(bayes_upper_bound_68, 0, bayes_upper_bound_68, f_posterior->Eval(bayes_upper_bound_68));
    
    TLegend *leg = new TLegend(0.55, 0.60, 0.88, 0.88);
    leg->SetHeader("Bayesian Analysis Summary", "C");
    leg->SetBorderSize(1);
    leg->AddEntry(f_posterior, "Posterior PDF", "l");
    auto entry_map = leg->AddEntry((TObject*)0, Form("Mode (MAP) = %.2f ms", map_estimate), "l");
    entry_map->SetLineColor(kRed);
    entry_map->SetLineStyle(2);
    entry_map->SetLineWidth(2);
    auto entry_median = leg->AddEntry((TObject*)0, Form("Median = %.2f ms", bayes_median_estimate), "l");
    entry_median->SetLineColor(kGreen+2);
    entry_median->SetLineStyle(7);
    entry_median->SetLineWidth(2);
    auto entry_ci = leg->AddEntry((TObject*)0, Form("68.3%% C.I. [%.2f, %.2f] ms", bayes_lower_bound_68, bayes_upper_bound_68), "l");
    entry_ci->SetLineColor(kGray+2);
    entry_ci->SetLineStyle(3);
    entry_ci->SetLineWidth(2);
    leg->Draw();
    
    c1->Update();
}