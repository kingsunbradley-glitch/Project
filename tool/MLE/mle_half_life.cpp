// mle_half_life.cpp
// 编译: 
//g++ -O2 `root-config --cflags --libs` -o mle_half_life mle_half_life.cpp
// 用法: ./mle_half_life [input_file]
// 如果没有 input_file，会从 stdin 读取，或可以在命令行中直接输入数字（参考下面的读取逻辑）
//./mle_half_life 0.107 0.184 0.52 0.0339 0.373 0.31 0.11 类似这种

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "TApplication.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"

// Solve for x in: g(x) = ln x - (x-1) + delta_over_n = 0
// We find one root in (eps,1) and one root in (1, high).
double g_of_x(double x, double delta_over_n) {
    return std::log(x) - (x - 1.0) + delta_over_n;
}

double find_root_bisect(double a, double b, double delta_over_n, int max_iter=100, double tol=1e-12) {
    double fa = g_of_x(a, delta_over_n);
    double fb = g_of_x(b, delta_over_n);
    if (fa * fb > 0) {
        return NAN;
    }
    for (int iter=0; iter<max_iter; ++iter) {
        double m = 0.5*(a+b);
        double fm = g_of_x(m, delta_over_n);
        if (std::abs(fm) < tol) return m;
        if (fa * fm <= 0) {
            b = m;
            fb = fm;
        } else {
            a = m;
            fa = fm;
        }
        if (std::abs(b-a) < tol) return 0.5*(a+b);
    }
    return 0.5*(a+b);
}

// Wrapper that finds the two solutions x_left in (0,1) and x_right in (1, +inf)
bool find_x_brackets(double delta_over_n, double &x_left, double &x_right) {
    // left root search in (eps, 1)
    double a = 1e-12; // small positive
    double b = 1.0;
    double fa = g_of_x(a, delta_over_n);
    double fb = g_of_x(b, delta_over_n);
    if (fa*fb > 0) {
        return false;
    }
    x_left = find_root_bisect(a, b, delta_over_n);
    if (!std::isfinite(x_left)) return false;

    // right root: find b such that g(b) < 0; start with b = 1 + step and increase
    double lo = 1.0;
    double hi = 2.0;
    double glo = g_of_x(lo, delta_over_n);
    double ghi = g_of_x(hi, delta_over_n);
    int tries = 0;
    while (glo * ghi > 0 && tries < 200) {
        hi *= 2.0;
        ghi = g_of_x(hi, delta_over_n);
        ++tries;
    }
    if (glo * ghi > 0) return false;
    x_right = find_root_bisect(lo, hi, delta_over_n);
    if (!std::isfinite(x_right)) return false;
    return true;
}

// log-likelihood (up to additive const) as function of lambda: ell = n ln lambda - lambda S
double logL_lambda(double lambda, int n, double S) {
    if (lambda <= 0) return -1e300;
    return n * std::log(lambda) - lambda * S;
}

// main MLE + plotting routine
int main(int argc, char **argv) {
    // 关键修改 1: 创建 TApplication 实例以处理图形窗口
    TApplication theApp("App", &argc, argv);

    std::vector<double> times;
    // read input:
    if (argc >= 2) {
        std::ifstream ifs(argv[1]);
        if (ifs) {
            std::string line;
            while (std::getline(ifs, line)) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream iss(line);
                double t;
                if (iss >> t) times.push_back(t);
            }
        } else {
            for (int i=1;i<argc;++i) {
                double t;
                if (std::istringstream(argv[i]) >> t) times.push_back(t);
            }
        }
    } else {
        std::cout << "Reading times from stdin (one time per line). Ctrl-D (Unix) or Ctrl-Z (Windows) to end.\n";
        std::string line;
        while (std::getline(std::cin, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            double t;
            if (iss >> t) times.push_back(t);
        }
    }

    if (times.empty()) {
        std::cerr << "No times provided. Exiting.\n";
        return 1;
    }

    int n = times.size();
    double S = 0.0;
    for (double t: times) S += t;

    double lambda_hat = double(n) / S;
    double T_half_hat = std::log(2.0) / lambda_hat;

    std::cout.setf(std::ios::fixed);
    std::cout.precision(8);
    std::cout << "n = " << n << "\n";
    std::cout << "Sum t_i = " << S << "\n";
    std::cout << "MLE lambda = " << lambda_hat << "\n";
    std::cout << "MLE T1/2 = " << T_half_hat << "\n";

    const double deltas[3] = {0.5, 2.0, 4.5}; // 对应 1σ, 2σ, 3σ
    double T_lo[3], T_hi[3];
    double err_minus[3], err_plus[3];

    for (int k=0;k<3;++k) {
        double Delta = deltas[k];
        double delta_over_n = Delta / double(n);
        double xl, xr;
        bool ok = find_x_brackets(delta_over_n, xl, xr);
        if (!ok) {
            std::cerr << "Failed to find x roots for sigma level k="<<k<<"\n";
            T_lo[k] = T_hi[k] = NAN;
            continue;
        }
        double lam_left = xl * double(n) / S;
        double lam_right = xr * double(n) / S;
        T_hi[k] = std::log(2.0) / lam_left;
        T_lo[k] = std::log(2.0) / lam_right;

        // *** 新增功能：计算非对称误差 ***
        err_minus[k] = T_half_hat - T_lo[k];
        err_plus[k]  = T_hi[k] - T_half_hat;

        std::cout << k+1 << " sigma interval: T_half = " << T_half_hat << " ( -" << err_minus[k] << " / +" << err_plus[k] << " ) s\n";
        std::cout << "                -> [" << T_lo[k] << ", " << T_hi[k] << "] s\n";
    }
    

    // =================================================================
    // 第 1 部分：绘制主要的对数似然函数图 (在 c1 上)
    // =================================================================
    TCanvas *c1 = new TCanvas("c1","Profile log-likelihood for T1/2 (normalized)",1200,800);
    c1->cd(); // 激活第一个画布
    gStyle->SetOptStat(0);

    double Tcenter = T_half_hat;
    double xmin = Tcenter * 0.02;
    double xmax = Tcenter * 5.0;
    for (int k=0;k<3;++k) {
        if (std::isfinite(T_lo[k]) && std::isfinite(T_hi[k])) {
            xmin = std::min(xmin, T_lo[k]*0.5);
            xmax = std::max(xmax, T_hi[k]*2.0);
        }
    }
    if (xmin <= 0) xmin = Tcenter*1e-4;

    const int NPOINTS = 800;
    std::vector<double> xs(NPOINTS), ys(NPOINTS);
    double logLmax = logL_lambda(lambda_hat, n, S);

    for (int i=0;i<NPOINTS;++i) {
        double T = xmin + (xmax - xmin) * double(i) / double(NPOINTS-1);
        double lambda = std::log(2.0) / T;
        xs[i] = T;
        ys[i] = logL_lambda(lambda, n, S) - logLmax;
    }
    
    TGraph *gr = new TGraph(NPOINTS, &xs[0], &ys[0]);
    gr->SetTitle("Profile log-likelihood (normalized);T_{1/2};#Delta ln L");
    gr->GetXaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetTitleSize(0.045);
    gr->GetYaxis()->SetRangeUser( - (deltas[2] + 1.0), 0.1 );
    gr->SetLineWidth(2);
    gr->Draw("AL");

    TLine *h1 = new TLine(xmin, -deltas[0], xmax, -deltas[0]);
    TLine *h2 = new TLine(xmin, -deltas[1], xmax, -deltas[1]);
    TLine *h3 = new TLine(xmin, -deltas[2], xmax, -deltas[2]);
    h1->SetLineStyle(2); h2->SetLineStyle(3); h3->SetLineStyle(5);
    h1->SetLineWidth(1); h2->SetLineWidth(1); h3->SetLineWidth(1);
    h1->Draw("same"); h2->Draw("same"); h3->Draw("same");

    TLine *vMLE = new TLine(T_half_hat, gr->GetYaxis()->GetXmin(), T_half_hat, 0.0);
    vMLE->SetLineColor(kRed); vMLE->SetLineWidth(2); vMLE->SetLineStyle(1);
    vMLE->Draw("same");

    int colors[3] = {kBlue+1, kGreen+2, kMagenta+2};
    for (int k=0;k<3;++k) {
        if (std::isfinite(T_lo[k]) && std::isfinite(T_hi[k])) {
            TLine *vlo = new TLine(T_lo[k], gr->GetYaxis()->GetXmin(), T_lo[k], 0.0);
            TLine *vhi = new TLine(T_hi[k], gr->GetYaxis()->GetXmin(), T_hi[k], 0.0);
            vlo->SetLineColor(colors[k]); vlo->SetLineStyle(7 - k); vlo->SetLineWidth(2);
            vhi->SetLineColor(colors[k]); vhi->SetLineStyle(7 - k); vhi->SetLineWidth(2);
            vlo->Draw("same"); vhi->Draw("same");
        }
    }

    TLegend *leg = new TLegend(0.65,0.6,0.89,0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(gr, "Profile ln L (normalized)", "l");
    leg->AddEntry(vMLE, "MLE T_{1/2}", "l");
    TLine d1, d2, d3;
    d1.SetLineColor(colors[0]); d2.SetLineColor(colors[1]); d3.SetLineColor(colors[2]);
    d1.SetLineWidth(2); d2.SetLineWidth(2); d3.SetLineWidth(2);
    leg->AddEntry(&d1, "1#sigma interval", "l");
    leg->AddEntry(&d2, "2#sigma interval", "l");
    leg->AddEntry(&d3, "3#sigma interval", "l");
    leg->Draw();
    c1->Update();

    // =================================================================
    // 第 2 部分：创建新画布，并单独绘制文本摘要 (在 c2 上)
    // =================================================================
    TCanvas *c2 = new TCanvas("c2", "Summary Data", 600, 400);
    c2->cd();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.05);
    double y_start = 0.9;
    double x_pos = 0.1;
    double y_step = 0.08;

    tex.DrawLatex(x_pos, y_start, Form("n = %d", n));
    tex.DrawLatex(x_pos, y_start - y_step, Form("Sum t_i = %.6g s", S));
    tex.DrawLatex(x_pos, y_start - 2*y_step, Form("MLE T_{1/2} = %.6g s", T_half_hat));
    
    TLine line(0.05, y_start - 2.5*y_step, 0.95, y_start - 2.5*y_step);
    line.Draw();

    for (int k=0;k<3;++k) {
        if (std::isfinite(T_lo[k]) && std::isfinite(T_hi[k])) {
            tex.DrawLatex(x_pos, y_start - (3.5+k)*y_step, Form("%d#sigma: [%.6g, %.6g] s", k+1, T_lo[k], T_hi[k]));
        } else {
            tex.DrawLatex(x_pos, y_start - (3.5+k)*y_step, Form("%d#sigma: failed", k+1));
        }
    }
    c2->Update();

    /* 自动保存画布为图片文件
    c1->SaveAs("mle_plot.png");
    c2->SaveAs("mle_summary.png");
    std::cout << "\nSaved plot to mle_plot.png and summary to mle_summary.png\n";
    */
    // 关键修改 2: 启动事件循环，使窗口保持显示
    std::cout << "Canvas windows created. Close any window to exit.\n";
    theApp.Run();

    return 0;
}
