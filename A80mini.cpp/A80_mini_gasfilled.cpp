// Build:
//   g++ -O2 -std=c++17 A39_mini.cpp -o A39_mini $(root-config --cflags --glibs)
//
// Run:
//   ./A39_mini <runnum> [E_loss_keV]
//   ./A39_mini <filename.root> [E_loss_keV]
//
// Examples:
//   ./A39_mini 362
//   ./A39_mini 362 20
//   ./A39_mini run00362_map.root 20
// Why 39?
// because the probability within $1\sigma$ of a 2D Gaussian distribution is 39%, and the probability within $n\sigma$ is given by $1 - \exp(-n^2/2)$

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <string>
#include <cctype>
#include <cstdio>

// ROOT Includes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TRootCanvas.h>    // GUI
#include <TCanvasImp.h>     // GUI

// ==========================================
//  1.Config.h
// ==========================================
const int NX = 128;
const int NY = 48;
const double XLO = -0.5, XHI = 127.5;
const double YLO = -0.5, YHI = 47.5;

// 固定抽样数 (N0)
const int RESERVOIR_SIZE = 5000;
// 覆盖率
const double TARGET_FRAC = 0.39;

// ==========================================
// 2.*.cpp
// ==========================================

struct AeqResult {
    long long n_pix = 0;
    double    n_eq  = 0.0;
    double    cov_at_npix = 0.0;
};

// 旋转：把 TH1 横着画成 “Counts-x / Channel-y” 的阶梯线
TGraph* CreateRotatedGraph(TH1* h) {
    if (!h) return nullptr;

    TGraph* gr = new TGraph();
    int nbins = h->GetNbinsX();

    double y_min = h->GetXaxis()->GetBinLowEdge(1);
    gr->SetPoint(0, 0, y_min);

    int pointIndex = 1;
    for (int i = 1; i <= nbins; ++i) {
        double content  = h->GetBinContent(i);
        double edgeLow  = h->GetBinLowEdge(i);
        double edgeHigh = h->GetBinLowEdge(i + 1);

        gr->SetPoint(pointIndex++, content, edgeLow);
        gr->SetPoint(pointIndex++, content, edgeHigh);
    }

    double y_max = h->GetXaxis()->GetBinLowEdge(nbins + 1);
    gr->SetPoint(pointIndex, 0, y_max);

    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    return gr;
}

// Metrics.cpp
AeqResult CalcAeqFromCounts(std::vector<double> counts, double targetFrac) {
    AeqResult r;
    if (counts.empty()) return r;

    double tot = 0.0;
    for (double c : counts) tot += c;
    if (tot <= 0.0) return r;

    std::sort(counts.begin(), counts.end(), std::greater<double>());

    const double target = targetFrac * tot;
    double acc = 0.0, accPrev = 0.0;
    long long nFull = 0;

    for (double c : counts) {
        accPrev = acc;
        acc += c;
        ++nFull;
        if (acc >= target) {
            r.n_pix = nFull;
            r.cov_at_npix = acc / tot;

            const double need = target - accPrev;
            double f = (c > 0.0) ? (need / c) : 1.0;
            if (f < 0.0) f = 0.0;
            if (f > 1.0) f = 1.0;

            r.n_eq = (nFull - 1) + f;
            return r;
        }
    }

    r.n_pix = (long long)counts.size();
    r.n_eq = (double)counts.size();
    r.cov_at_npix = 1.0;
    return r;
}

// 从蓄水池样本中重建计数并计算 A39
AeqResult CalculateReservoirA39(const std::vector<int>& binIdx) {
    int nbin = NX * NY;
    std::vector<int> cnt(nbin, 0);

    for (int idx : binIdx) {
        if (idx >= 0 && idx < nbin) cnt[idx]++;
    }

    std::vector<double> counts;
    counts.reserve(nbin);
    for (int c : cnt) {
        if (c > 0) counts.push_back((double)c);
    }
    return CalcAeqFromCounts(counts, TARGET_FRAC);
}

// ==========================================
// 3. main.cpp
// ==========================================
int main(int argc, char** argv) {
    auto IsAllDigits = [](const std::string& s) -> bool {
        if (s.empty()) return false;
        for (unsigned char ch : s) {
            if (!std::isdigit(ch)) return false;
        }
        return true;
    };

    if (argc < 2) {
        std::cerr << "Usage:\n"
                  << "  1) ./A39_mini <runnum> [E_loss_keV]\n"
                  << "     -> opens default file: run%05d_map.root\n"
                  << "  2) ./A39_mini <filename.root> [E_loss_keV]\n"
                  << "Examples:\n"
                  << "  ./A39_mini 362\n"
                  << "  ./A39_mini 362 20\n"
                  << "  ./A39_mini run00362_map.root 20\n";
        return 1;
    }

    // 解析 E_loss（keV），默认 0
    double E_loss = 0.0;
    if (argc >= 3) {
        try {
            E_loss = std::stod(argv[2]);
        } catch (...) {
            std::cerr << "Error: E_loss must be a number (keV). Got: " << argv[2] << std::endl;
            return 1;
        }
    }

    // argv[1]：run 号或文件名
    std::string arg1 = argv[1];
    std::string fileName;
    if (IsAllDigits(arg1)) {
        int run = std::stoi(arg1);
        char buf[512];
        std::snprintf(buf, sizeof(buf), "run%05d_map.root", run);
        fileName = buf;
    } else {
        fileName = arg1;
    }

    // ROOT app（必须在任何 GUI / canvas 前）
    TApplication app("app", &argc, argv);

    TFile fin(fileName.c_str(), "READ");
    if (fin.IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return 1;
    }

    TTree* tr = (TTree*)fin.Get("tr_map");
    if (!tr) {
        std::cerr << "Error: tr_map not found in file." << std::endl;
        return 1;
    }

    UShort_t mx = 0, my = 0;
    Double_t Xch[256]{0}, Ych[256]{0};
    Double_t XE[256]{0};

    tr->SetBranchAddress("DSSDX_mul", &mx);
    tr->SetBranchAddress("DSSDY_mul", &my);
    tr->SetBranchAddress("DSSDX_Ch", Xch);
    tr->SetBranchAddress("DSSDY_Ch", Ych);
    tr->SetBranchAddress("DSSDX_E", XE);

    // ====== 能量参数（名义值） ======
    const double SPEC_NOM_LO = 5000.0;
    const double SPEC_NOM_HI = 6000.0;
    const double ROI_NOM_LO  = 5700.0;
    const double ROI_NOM_HI  = 5900.0;

    // ✅ 你的要求：所有能量参数都减去 E_loss
    const double SPEC_LO = SPEC_NOM_LO - E_loss;
    const double SPEC_HI = SPEC_NOM_HI - E_loss;
    const double ROI_LO  = ROI_NOM_LO  - E_loss;
    const double ROI_HI  = ROI_NOM_HI  - E_loss;

    // 定义直方图：范围使用“减去 E_loss 后”的窗口
    TH1D hE("hE", Form("Energy Spectrum (%.0f-%.0f keV);Energy [keV];Counts", SPEC_LO, SPEC_HI),
            500, SPEC_LO, SPEC_HI);

    TH2D hXY("hXY", Form("Spatial Dist (%.0f-%.0f keV);X [ch];Y [ch]", ROI_LO, ROI_HI),
             NX, XLO, XHI, NY, YLO, YHI);

    // 蓄水池抽样变量
    std::vector<int> reservoir;
    reservoir.reserve(RESERVOIR_SIZE);
    long long seenCount = 0;
    std::mt19937_64 rng(12345);

    Long64_t nEntries = tr->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tr->GetEntry(i);

        // 基本 Cut: 多重性为 1
        if (mx != 1 || my != 1) continue;

        // 事件能量：保持原始 XE[0]，不再做 E = E_raw - E_loss
        double E = XE[0];

        // 全能谱窗口（参数已减 E_loss）
        if (E >= SPEC_LO && E <= SPEC_HI) {
            hE.Fill(E);
        }

        // ROI（参数已减 E_loss）
        if (E >= ROI_LO && E <= ROI_HI) {
            double x = Xch[0];
            double y = Ych[0];

            hXY.Fill(x, y);

            // 蓄水池抽样
            int ix = (int)std::round(x);
            int iy = (int)std::round(y);

            if (ix >= 0 && ix < NX && iy >= 0 && iy < NY) {
                int binIdx = iy * NX + ix;
                seenCount++;

                if ((int)reservoir.size() < RESERVOIR_SIZE) {
                    reservoir.push_back(binIdx);
                } else {
                    std::uniform_int_distribution<long long> dist(0, seenCount - 1);
                    long long r = dist(rng);
                    if (r < RESERVOIR_SIZE) {
                        reservoir[(size_t)r] = binIdx;
                    }
                }
            }
        }
    }

    // ==========================================
    // 4. 计算结果
    // ==========================================
    AeqResult res;
    if (seenCount >= RESERVOIR_SIZE) {
        res = CalculateReservoirA39(reservoir);
    } else {
        std::cout << "[Warning] Total events in ROI (" << seenCount
                  << ") < Reservoir Size (" << RESERVOIR_SIZE
                  << "). Using all events.\n";
        res = CalculateReservoirA39(reservoir);
    }

    // 终端输出
    std::cout << "------------------------------------------------\n";
    std::cout << "File:                 " << fileName << "\n";
    std::cout << "E_loss (keV):          " << E_loss << "\n";
    std::cout << "Spectrum window (keV): " << SPEC_LO << " - " << SPEC_HI << "  (nominal "
              << SPEC_NOM_LO << "-" << SPEC_NOM_HI << ")\n";
    std::cout << "ROI window (keV):      " << ROI_LO << " - " << ROI_HI << "  (nominal "
              << ROI_NOM_LO  << "-" << ROI_NOM_HI  << ")\n";
    std::cout << "Events in ROI:         " << seenCount << "\n";
    std::cout << "Reservoir N0:          " << RESERVOIR_SIZE << "\n";
    std::cout << "------------------------------------------------\n";
    std::cout << "Fixed-N A39 (EqPix):   " << res.n_eq << "\n";
    std::cout << "Raw Pixels used:       " << res.n_pix << "\n";
    std::cout << "Coverage achieved:     " << res.cov_at_npix * 100.0 << " %\n";
    std::cout << "------------------------------------------------\n";

    // ==========================================
    // 5. 绘图与保存
    // ==========================================
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "A39 Summary", 1200, 900);
    c1->Divide(2, 2);

    // Pad 1: XY 2D Map
    c1->cd(1);
    gPad->SetRightMargin(0.12);
    hXY.SetTitle(Form("XY Distribution (%.0f-%.0f keV), N=%lld,Eloss: %.1f keV", ROI_LO, ROI_HI, seenCount, E_loss));
    hXY.Draw("COLZ");

    // Pad 2: Y Projection 
    c1->cd(2);
    TH1D* hY_temp = hXY.ProjectionY("hY_temp");
    TGraph* gY = CreateRotatedGraph(hY_temp);
    gY->SetTitle("Y Projection;Counts;Channel");
    gY->SetLineColor(kBlue);
    gY->Draw("AL");
    gPad->RedrawAxis();

    // Pad 3: X Projection
    c1->cd(3);
    TH1D* hX = hXY.ProjectionX("hX");
    hX->SetTitle("X Projection;Channel;Counts");
    hX->Draw();

    // Pad 4: Energy Spectrum
    c1->cd(4);
    hE.SetTitle(Form("Global Energy Spectrum (%.0f-%.0f keV);Energy [keV];Counts", SPEC_LO, SPEC_HI));
    hE.Draw();

    // ROI 线
    c1->cd(4);
    TLine l; l.SetLineColor(kRed); l.SetLineStyle(2); l.SetLineWidth(2);
    double ymax = hE.GetMaximum();
    l.DrawLine(ROI_LO, 0, ROI_LO, ymax);
    l.DrawLine(ROI_HI, 0, ROI_HI, ymax);

    // 标注结果
    c1->cd(1);
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.04);
    lat.SetTextColor(kRed);
    lat.DrawLatex(0.15, 0.85, Form("A39(FixN): %.2f", res.n_eq));
    lat.SetTextSize(0.03);
    lat.SetTextColor(kBlack);
    //lat.DrawLatex(0.15, 0.79, Form("E_{loss}: %.1f keV", E_loss));

    // 保存图片
    std::string outImg = fileName + "_A39.pdf";
    c1->SaveAs(outImg.c_str());
    std::cout << "Plot saved to: " << outImg << std::endl;

    // 关闭窗口自动退出
    TRootCanvas* rc = dynamic_cast<TRootCanvas*>(c1->GetCanvasImp());
    if (rc) {
        rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    }

    c1->Update();
    app.Run(kTRUE);

    fin.Close();
    return 0;
}
