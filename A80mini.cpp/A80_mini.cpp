//g++ A80_mini.cpp -o A80_mini $(root-config --cflags --glibs)
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <string>

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
#include <TRootCanvas.h>    //GUI
#include <TCanvasImp.h>     //GUI


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
const double TARGET_FRAC = 0.80;

// ==========================================
// 2.*.cpp
// ==========================================

struct AeqResult {
    long long n_pix = 0;
    double    n_eq  = 0.0;
    double    cov_at_npix = 0.0;
};
// 旋转
TGraph* CreateRotatedGraph(TH1* h) {
    if (!h) return nullptr;
    
    TGraph* gr = new TGraph();
    int nbins = h->GetNbinsX();
    
    
    double y_min = h->GetXaxis()->GetBinLowEdge(1);
    gr->SetPoint(0, 0, y_min);
    
    int pointIndex = 1; // 从 1 开始，因为 0 被占用了
    
    for (int i = 1; i <= nbins; ++i) {
        double content = h->GetBinContent(i);
        double edgeLow = h->GetBinLowEdge(i);
        double edgeHigh = h->GetBinLowEdge(i+1); 
        
        // 交换 X 和 Y
        gr->SetPoint(pointIndex++, content, edgeLow);
        gr->SetPoint(pointIndex++, content, edgeHigh);
    }
    
   
    double y_max = h->GetXaxis()->GetBinLowEdge(nbins + 1);
    gr->SetPoint(pointIndex, 0, y_max);
    
    
    gr->SetLineWidth(h->GetLineWidth());
    gr->SetLineColor(h->GetLineColor());
    
    return gr;
}
//  Metrics.cpp 
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
            if (f < 0.0) f = 0.0; if (f > 1.0) f = 1.0;
            r.n_eq = (nFull - 1) + f;
            return r;
        }
    }
    r.n_pix = (long long)counts.size();
    r.n_eq = (double)counts.size();
    r.cov_at_npix = 1.0;
    return r;
}

// 从蓄水池样本中重建计数并计算 A80
AeqResult CalculateReservoirA80(const std::vector<int>& binIdx) {
    int nbin = NX * NY;
    std::vector<int> cnt(nbin, 0);
    
    // 统计样本落在哪些 bin
    for (int idx : binIdx) {
        if (idx >= 0 && idx < nbin) cnt[idx]++;
    }

    // 转为 double vector 传入通用算法
    std::vector<double> counts;
    for (int c : cnt) {
        if (c > 0) counts.push_back((double)c);
    }
    return CalcAeqFromCounts(counts, TARGET_FRAC);
}

// ==========================================
// 3. main.cpp
// ==========================================
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: ./A80 [filename.root]" << std::endl;
        return 1;
    }

    std::string fileName = argv[1];
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
    Double_t Xch[256]{0}, Ych[256]{0}; // 数组大小根据实际情况，256通常够用
    Double_t XE[256]{0};

    tr->SetBranchAddress("DSSDX_mul", &mx);
    tr->SetBranchAddress("DSSDY_mul", &my);
    tr->SetBranchAddress("DSSDX_Ch", Xch);
    tr->SetBranchAddress("DSSDY_Ch", Ych);
    tr->SetBranchAddress("DSSDX_E", XE);

    // 定义直方图
    // 1. 全能谱 (5000 - 6000)
    TH1D hE("hE", "Energy Spectrum (5000-6000 keV);Energy [keV];Counts", 500, 5000, 6000);
    
    // 2. 目标 ROI 的空间分布 (5700 - 5900)
    TH2D hXY("hXY", "Spatial Dist (5700-5900 keV);X [ch];Y [ch]", NX, XLO, XHI, NY, YLO, YHI);

    // 蓄水池抽样变量
    std::vector<int> reservoir;
    reservoir.reserve(RESERVOIR_SIZE);
    long long seenCount = 0;
    std::mt19937_64 rng(12345); // 固定种子以保证结果可复现

    Long64_t nEntries = tr->GetEntries();
    // std::cout << "Processing " << nEntries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tr->GetEntry(i);

        // 基本 Cut: 多重性为 1
        if (mx != 1 || my != 1) continue;
        
        double E = XE[0];
        
        // 填充全能谱
        if (E >= 5000 && E <= 6000) {
            hE.Fill(E);
        }

        // 仅针对目标区间 5700-5900 进行空间统计和抽样
        if (E >= 5700 && E <= 5900) {
            double x = Xch[0];
            double y = Ych[0];

            // 填充直方图 (用于绘图)
            hXY.Fill(x, y);

            // --- 蓄水池抽样逻辑 ---
            int ix = (int)std::round(x);
            int iy = (int)std::round(y);

            // 确保在范围内
            if (ix >= 0 && ix < NX && iy >= 0 && iy < NY) {
                int binIdx = iy * NX + ix;
                seenCount++;

                if (reservoir.size() < RESERVOIR_SIZE) {
                    reservoir.push_back(binIdx);
                } else {
                    // 随机替换
                    std::uniform_int_distribution<long long> dist(0, seenCount - 1);
                    long long r = dist(rng);
                    if (r < RESERVOIR_SIZE) {
                        reservoir[r] = binIdx;
                    }
                }
            }
        }
    }

    // ==========================================
    // 4. 计算结果
    // ==========================================
    AeqResult res;
    // 只有当收集到的样本数达到 N0 时，计算才有意义（或者你也可以允许小于N0直接算）
    if (seenCount >= RESERVOIR_SIZE) {
        res = CalculateReservoirA80(reservoir);
    } else {
        std::cout << "[Warning] Total events in ROI (" << seenCount << ") < Reservoir Size (" << RESERVOIR_SIZE << "). Using all events." << std::endl;
        res = CalculateReservoirA80(reservoir); // 样本不足时，reservoir里就是所有事件
    }

    // 终端输出
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "File: " << fileName << std::endl;
    std::cout << "ROI:  5700 - 5900 keV" << std::endl;
    std::cout << "Events in ROI: " << seenCount << std::endl;
    std::cout << "Reservoir N0:  " << RESERVOIR_SIZE << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "Fixed-N A80 (Eq. Pixels): " << res.n_eq << std::endl;
    std::cout << "Raw Pixels used:          " << res.n_pix << std::endl;
    std::cout << "Coverage achieved:        " << res.cov_at_npix * 100.0 << " %" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

    // ==========================================
    // 5. 绘图与保存
    // ==========================================
    gStyle->SetOptStat(0);
    
    // 【修改点 1】必须使用指针 (new TCanvas)，否则 Close 信号无法正确终止程序
    TCanvas *c1 = new TCanvas("c1", "A80 Summary", 1200, 900);
    c1->Divide(2, 2); // 注意：所有的 . 都要变成 ->

    // Pad 1: XY 2D Map
    c1->cd(1);
    gPad->SetRightMargin(0.12);
    hXY.SetTitle(Form("XY Distribution (5700-5900 keV), N=%lld", seenCount));
    hXY.Draw("COLZ");

    // Pad 2: Y Projection (右上角)
    c1->cd(2);
    TH1D* hY_temp = hXY.ProjectionY("hY_temp");
    TGraph* gY = CreateRotatedGraph(hY_temp);
    gY->SetTitle("Y Projection;Counts;Channel");
    // gY->SetFillColor(kBlue-7); 
    gY->SetLineColor(kBlue);   
    gY->Draw("ALF");
    gPad->RedrawAxis();
    
    // Pad 3: X Projection
    c1->cd(3);
    TH1D* hX = hXY.ProjectionX("hX");
    hX->SetTitle("X Projection;Channel;Counts");
    hX->Draw();

    // Pad 4: Energy Spectrum
    c1->cd(4);
    hE.SetTitle("Global Energy Spectrum (5000-6000 keV);Energy [keV];Counts");
    hE.Draw();
    
    c1->cd(4);
    TLine l; l.SetLineColor(kRed); l.SetLineStyle(2); l.SetLineWidth(2);
    double ymax = hE.GetMaximum();
    l.DrawLine(5700, 0, 5700, ymax);
    l.DrawLine(5900, 0, 5900, ymax);

    // 标注结果
    c1->cd(1);
    TLatex lat; 
    lat.SetNDC(); lat.SetTextSize(0.04); lat.SetTextColor(kRed);
    lat.DrawLatex(0.15, 0.85, Form("A80(FixN): %.2f", res.n_eq));

    // 保存图片
    std::string outImg = fileName + "_A80.pdf";
    c1->SaveAs(outImg.c_str());
    std::cout << "Plot saved to: " << outImg << std::endl;

    TRootCanvas *rc = dynamic_cast<TRootCanvas*>(c1->GetCanvasImp());
    if (rc) {
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
    }

    
    c1->Update(); 
    app.Run(kTRUE);   // 程序在这里暂停，直到窗口关闭触发 Terminate()
    
    // 清理资源（可选，因为程序都要退出了）
    fin.Close();
    // delete c1; // 指针可以手动 delete，但通常操作系统会回收
    return 0;
}