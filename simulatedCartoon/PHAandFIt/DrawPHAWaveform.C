// DrawPHAWaveform.C
// root -l DrawPHAWaveform.C
//
// 输出：
//   PHA_waveform_algorithm_demo.pdf
//   PHA_waveform_algorithm_demo.png

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TF1.h>
#include <TStyle.h>
#include <TAxis.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

using std::vector;

// ============================================================
// 全局参数
// ============================================================
const double T_START = 0.0;       // us
const double T_END   = 30.0;      // us
const double DT      = 0.01;      // us, 10 ns/sample
const double PRETRIG = 5.0;       // us

const double TAU_R = 0.08;        // us, pulse rise time
const double TAU_D = 9.0;         // us, preamp decay time

// ============================================================
// 前放脉冲模型：快上升 + 慢衰减
// ============================================================
double PulseShape(double t, double t0, double A)
{
    if (t < t0) return 0.0;

    double x = t - t0;

    return A * (1.0 - std::exp(-x / TAU_R)) * std::exp(-x / TAU_D);
}

// ============================================================
// 双脉冲波形拟合函数
// f(t) = B + A1 P(t-t1) + A2 P(t-t2)
// ============================================================
double TwoPulseFitFunc(double *xx, double *par)
{
    double t = xx[0];

    double B  = par[0];
    double A1 = par[1];
    double t1 = par[2];
    double A2 = par[3];
    double t2 = par[4];

    return B + PulseShape(t, t1, A1) + PulseShape(t, t2, A2);
}

// ============================================================
// 生成原始波形，带 baseline offset 和轻微噪声
// ============================================================
TGraph* MakeRawWaveform(const vector<std::pair<double, double>>& pulses,
                        const char* name)
{
    int n = int((T_END - T_START) / DT) + 1;

    vector<double> x(n), y(n);

    const double baseline = 0.12;

    for (int i = 0; i < n; ++i) {
        double t = T_START + i * DT;
        double v = baseline;

        for (auto &p : pulses) {
            v += PulseShape(t, p.first, p.second);
        }

        // small noise / ripple
        v += 0.003 * std::sin(2.0 * M_PI * t / 0.37);

        x[i] = t;
        y[i] = v;
    }

    TGraph *g = new TGraph(n, x.data(), y.data());
    g->SetName(name);
    g->SetLineWidth(2);
    return g;
}

// ============================================================
// 基零修正：
// 用触发前 0--4.5 us 的平均值作为 baseline
// ============================================================
TGraph* BaselineCorrection(TGraph *graw,
                           double bmin,
                           double bmax,
                           const char* name)
{
    int n = graw->GetN();

    double *x = graw->GetX();
    double *y = graw->GetY();

    double sum = 0.0;
    int cnt = 0;

    for (int i = 0; i < n; ++i) {
        if (x[i] >= bmin && x[i] <= bmax) {
            sum += y[i];
            cnt++;
        }
    }

    double baseline = (cnt > 0) ? sum / cnt : 0.0;

    vector<double> xx(n), yy(n);

    for (int i = 0; i < n; ++i) {
        xx[i] = x[i];
        yy[i] = y[i] - baseline;
    }

    TGraph *g = new TGraph(n, xx.data(), yy.data());
    g->SetName(name);
    g->SetLineWidth(2);

    return g;
}

// ============================================================
// 三角算法示意：
// 用对称差分近似 triangle filter output。
// 实际系统中三角滤波通常用于前沿时间提取。
// ============================================================
TGraph* TriangleFilter(TGraph *g,
                       double width_us,
                       const char* name)
{
    int n = g->GetN();

    double *x = g->GetX();
    double *y = g->GetY();

    int W = std::max(1, int(width_us / DT));

    vector<double> xx;
    vector<double> yy;

    for (int i = W; i < n - W; ++i) {
        double tri = y[i + W] - y[i - W];

        xx.push_back(x[i]);
        yy.push_back(tri);
    }

    // 归一化到方便和波形画在一起
    double maxabs = 0.0;
    for (auto v : yy) {
        maxabs = std::max(maxabs, std::abs(v));
    }

    if (maxabs > 0) {
        for (auto &v : yy) {
            v = 0.8 * v / maxabs;
        }
    }

    TGraph *gt = new TGraph((int)xx.size(), xx.data(), yy.data());
    gt->SetName(name);
    gt->SetLineWidth(2);
    gt->SetLineStyle(2);

    return gt;
}

// ============================================================
// 梯形算法示意：
// late moving average - early moving average
// 平台高度对应 PHA energy
// ============================================================
TGraph* TrapezoidFilter(TGraph *g,
                        double rise_us,
                        double gap_us,
                        const char* name)
{
    int n = g->GetN();

    double *x = g->GetX();
    double *y = g->GetY();

    int L = std::max(1, int(rise_us / DT));
    int G = std::max(0, int(gap_us  / DT));

    vector<double> xx;
    vector<double> yy;

    for (int i = 2 * L + G; i < n; ++i) {
        double sum_early = 0.0;
        double sum_late  = 0.0;

        for (int k = i - 2 * L - G; k < i - L - G; ++k) {
            sum_early += y[k];
        }

        for (int k = i - L; k < i; ++k) {
            sum_late += y[k];
        }

        double trap = sum_late / L - sum_early / L;

        xx.push_back(x[i]);
        yy.push_back(trap);
    }

    // 归一化到方便和波形画在一起
    double maxabs = 0.0;
    for (auto v : yy) {
        maxabs = std::max(maxabs, std::abs(v));
    }

    if (maxabs > 0) {
        for (auto &v : yy) {
            v = 0.8 * v / maxabs;
        }
    }

    TGraph *gt = new TGraph((int)xx.size(), xx.data(), yy.data());
    gt->SetName(name);
    gt->SetLineWidth(2);
    gt->SetLineStyle(3);

    return gt;
}

// ============================================================
// 找三角滤波峰位置作为 timestamp 示意
// ============================================================
double FindPeakTime(TGraph *g, double xmin, double xmax)
{
    int n = g->GetN();

    double *x = g->GetX();
    double *y = g->GetY();

    double best_t = xmin;
    double best_y = -1e99;

    for (int i = 0; i < n; ++i) {
        if (x[i] < xmin || x[i] > xmax) continue;

        if (y[i] > best_y) {
            best_y = y[i];
            best_t = x[i];
        }
    }

    return best_t;
}

// ============================================================
// 画竖线
// ============================================================
void DrawVLine(double x, double ymin, double ymax, int color, int style = 2)
{
    TLine *l = new TLine(x, ymin, x, ymax);
    l->SetLineColor(color);
    l->SetLineStyle(style);
    l->SetLineWidth(2);
    l->Draw();
}

// ============================================================
// 每个 panel 的通用格式
// ============================================================
void FormatFrame(TGraph *g,
                 const char *title,
                 double ymin,
                 double ymax)
{
    g->SetTitle(title);
    g->GetXaxis()->SetRangeUser(0, 30);
    g->GetYaxis()->SetRangeUser(ymin, ymax);
    g->GetXaxis()->SetTitle("Time [#mus]");
    g->GetYaxis()->SetTitle("Amplitude [a.u.]");
}

// ============================================================
// 主程序
// ============================================================
void DrawPHAWaveform()
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetLabelSize(0.040, "XYZ");

    TCanvas *c = new TCanvas("c", "PHA and waveform algorithm demo", 1500, 1100);
    c->Divide(2, 2);

    // ========================================================
    // panel 1: DSSD single pulse
    // ========================================================
    c->cd(1);

    TGraph *raw1 = MakeRawWaveform({{5.0, 1.0}}, "raw1");
    TGraph *wf1  = BaselineCorrection(raw1, 0.0, 4.5, "wf1");

    TGraph *tri1  = TriangleFilter(wf1, 0.25, "tri1");
    TGraph *trap1 = TrapezoidFilter(wf1, 0.8, 0.5, "trap1");

    double ts1 = FindPeakTime(tri1, 4.5, 6.0);

    wf1->SetLineColor(kBlack);
    tri1->SetLineColor(kBlue + 1);
    trap1->SetLineColor(kRed + 1);

    FormatFrame(wf1,
                "(1) DSSD single pulse",
                -0.20,
                1.20);

    wf1->Draw("AL");
    tri1->Draw("L SAME");
    trap1->Draw("L SAME");

    DrawVLine(5.0, -0.20, 1.20, kGray + 2, 3);
    DrawVLine(ts1, -0.20, 1.20, kBlue + 1, 2);

    TLegend *leg1 = new TLegend(0.55, 0.68, 0.88, 0.88);
    leg1->SetBorderSize(0);
    leg1->AddEntry(wf1,  "baseline-corrected waveform", "l");
    leg1->AddEntry(tri1, "triangle filter", "l");
    leg1->AddEntry(trap1, "trapezoid filter", "l");
    leg1->AddEntry((TObject*)0, "vertical line: timestamp", "");
    leg1->Draw();

    // ========================================================
    // panel 2: 5 us + 8 us pile-up, separable
    // ========================================================
    c->cd(2);

    TGraph *raw2 = MakeRawWaveform({{5.0, 1.0}, {8.0, 0.75}}, "raw2");
    TGraph *wf2  = BaselineCorrection(raw2, 0.0, 4.5, "wf2");

    TGraph *tri2  = TriangleFilter(wf2, 0.25, "tri2");
    TGraph *trap2 = TrapezoidFilter(wf2, 0.45, 0.20, "trap2");

    double ts2a = FindPeakTime(tri2, 4.5, 6.0);
    double ts2b = FindPeakTime(tri2, 7.5, 9.0);

    wf2->SetLineColor(kBlack);
    tri2->SetLineColor(kBlue + 1);
    trap2->SetLineColor(kRed + 1);

    FormatFrame(wf2,
                "(2) Pile-up: 5 #mus and 8 #mus",
                -0.25,
                1.50);

    wf2->Draw("AL");
    tri2->Draw("L SAME");
    trap2->Draw("L SAME");

    DrawVLine(5.0, -0.25, 1.50, kGray + 2, 3);
    DrawVLine(8.0, -0.25, 1.50, kGray + 2, 3);
    DrawVLine(ts2a, -0.25, 1.50, kBlue + 1, 2);
    DrawVLine(ts2b, -0.25, 1.50, kBlue + 1, 2);

    TLegend *leg2 = new TLegend(0.55, 0.68, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->AddEntry(wf2,  "baseline-corrected waveform", "l");
    leg2->AddEntry(tri2, "triangle filter", "l");
    leg2->AddEntry(trap2, "fine trapezoid filter", "l");
    leg2->AddEntry((TObject*)0, "vertical line: timestamp", "");
    leg2->Draw();

    // ========================================================
    // panel 3: 5 us + 5.5 us pile-up, trapezoid fails
    // ========================================================
    c->cd(3);

    TGraph *raw3 = MakeRawWaveform({{5.0, 1.0}, {5.5, 0.75}}, "raw3");
    TGraph *wf3  = BaselineCorrection(raw3, 0.0, 4.5, "wf3");

    TGraph *tri3  = TriangleFilter(wf3, 0.25, "tri3");
    TGraph *trap3 = TrapezoidFilter(wf3, 0.8, 0.5, "trap3");

    double ts3a = FindPeakTime(tri3, 4.5, 5.3);
    double ts3b = FindPeakTime(tri3, 5.3, 6.2);

    wf3->SetLineColor(kBlack);
    tri3->SetLineColor(kBlue + 1);
    trap3->SetLineColor(kRed + 1);

    FormatFrame(wf3,
                "(3) Close pile-up: 5 #mus and 5.5 #mus",
                -0.25,
                1.70);

    wf3->Draw("AL");
    tri3->Draw("L SAME");
    trap3->Draw("L SAME");

    DrawVLine(5.0, -0.25, 1.70, kGray + 2, 3);
    DrawVLine(5.5, -0.25, 1.70, kGray + 2, 3);
    DrawVLine(ts3a, -0.25, 1.70, kBlue + 1, 2);
    DrawVLine(ts3b, -0.25, 1.70, kBlue + 1, 2);

    TLegend *leg3 = new TLegend(0.55, 0.68, 0.88, 0.88);
    leg3->SetBorderSize(0);
    leg3->AddEntry(wf3,  "baseline-corrected waveform", "l");
    leg3->AddEntry(tri3, "triangle filter", "l");
    leg3->AddEntry(trap3, "overlapped trapezoid filter", "l");
    leg3->AddEntry((TObject*)0, "vertical line: timestamp", "");
    leg3->Draw();

    // ========================================================
    // panel 4: close pile-up + waveform fitting
    // ========================================================
    c->cd(4);

    TGraph *raw4 = MakeRawWaveform({{5.0, 1.0}, {5.5, 0.75}}, "raw4");
    TGraph *wf4  = BaselineCorrection(raw4, 0.0, 4.5, "wf4");

    TGraph *tri4 = TriangleFilter(wf4, 0.25, "tri4");

    TF1 *fit4 = new TF1("fit4", TwoPulseFitFunc, 4.2, 12.0, 5);

    // 因为 wf4 已经做过 baseline correction，所以 B 初值为 0
    fit4->SetParameters(0.0, 1.0, 5.0, 0.75, 5.5);
    fit4->SetParNames("B", "A1", "t1", "A2", "t2");

    fit4->SetParLimits(0, -0.05, 0.05);
    fit4->SetParLimits(1, 0.1, 2.0);
    fit4->SetParLimits(2, 4.8, 5.2);
    fit4->SetParLimits(3, 0.1, 2.0);
    fit4->SetParLimits(4, 5.3, 5.8);

    wf4->Fit(fit4, "RQ0");

    double ts4a = FindPeakTime(tri4, 4.5, 5.3);
    double ts4b = FindPeakTime(tri4, 5.3, 6.2);

    wf4->SetLineColor(kBlack);
    tri4->SetLineColor(kBlue + 1);

    fit4->SetLineColor(kGreen + 2);
    fit4->SetLineWidth(3);

    FormatFrame(wf4,
                "(4) f(t)=B+A_{1}P(t-t_{1})+A_{2}P(t-t_{2})",
                -0.25,
                1.70);

    wf4->Draw("AL");
    tri4->Draw("L SAME");
    fit4->Draw("SAME");

    DrawVLine(5.0, -0.25, 1.70, kGray + 2, 3);
    DrawVLine(5.5, -0.25, 1.70, kGray + 2, 3);
    DrawVLine(ts4a, -0.25, 1.70, kBlue + 1, 2);
    DrawVLine(ts4b, -0.25, 1.70, kBlue + 1, 2);

    TLegend *leg4 = new TLegend(0.53, 0.66, 0.88, 0.88);
    leg4->SetBorderSize(0);
    leg4->AddEntry(wf4,  "baseline-corrected waveform", "l");
    leg4->AddEntry(tri4, "triangle filter", "l");
    leg4->AddEntry(fit4, "waveform fit", "l");
    leg4->AddEntry((TObject*)0, "vertical line: timestamp", "");
    leg4->Draw();

    c->SaveAs("PHA_waveform_algorithm_demo.pdf");
    c->SaveAs("PHA_waveform_algorithm_demo.png");

    std::cout << "\nSaved:\n"
              << "  PHA_waveform_algorithm_demo.pdf\n"
              << "  PHA_waveform_algorithm_demo.png\n\n";
}