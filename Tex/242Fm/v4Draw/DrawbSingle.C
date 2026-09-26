#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <vector>

void DrawbSingle()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const TString macroDirectory = gSystem->DirName(__FILE__);
    const TString inputPath = macroDirectory + "/plot_data.root";

    TFile *input = TFile::Open(inputPath, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "[Error] Cannot open " << inputPath << std::endl;
        delete input;
        return;
    }

    // Use the GREEN low-energy data, corresponding to the original 7--9.5 MeV plot.
    TTree *tree = nullptr;
    input->GetObject("tree_low_green", tree);

    if (!tree) {
        std::cerr << "[Error] Cannot find tree_low_green in plot_data.root"
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    // Same low-energy X range as the original single-green plot.
    const double xMin = 7.0;      // MeV
    const double xMax = 9.5;      // MeV
    const int nBinsX = 125;       // 20 keV/bin

    // Requested logarithmic time range: 1 us -- 10 s.
    const double yMin = 1.0e-6;   // s
    const double yMax = 10.0;     // s
    const int nBinsY = 100;

    // Logarithmic binning in Y so the COLZ cells are uniform on a log-Y canvas.
    std::vector<double> yBins(nBinsY + 1);
    const double logYMin = std::log10(yMin);
    const double logYMax = std::log10(yMax);

    for (int i = 0; i <= nBinsY; ++i) {
        yBins[i] = std::pow(
            10.0,
            logYMin + (logYMax - logYMin) * i / nBinsY);
    }

    TCanvas *c = new TCanvas("c_b_single", "", 1100, 600);
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.14);   // room for COLZ palette
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.05);
    c->SetTickx(1);
    c->SetTicky(1);
    c->SetLogy(1);

    TH2D *h2 = new TH2D(
        "h2_green",
        "",
        nBinsX, xMin, xMax,
        nBinsY, yBins.data());

    h2->SetDirectory(nullptr);
    h2->GetXaxis()->SetTitle("Energy (MeV)");
    h2->GetYaxis()->SetTitle("DeltaT_s (s)");
    h2->GetZaxis()->SetTitle("Counts");

    // Explicitly keep only 7.0--9.5 MeV and 1 us--10 s.
    TString cut =
        "SumE>=7000 && SumE<=9500"
        " && DeltaT_s>=1e-6 && DeltaT_s<=10";

    tree->Draw("DeltaT_s:SumE/1000.0>>h2_green", cut, "COLZ");

    // Force the displayed range as well, so ROOT cannot auto-expand it.
    h2->GetXaxis()->SetRangeUser(xMin, xMax);
    h2->GetYaxis()->SetRangeUser(yMin, yMax);

    c->Modified();
    c->Update();

    // No SaveAs(): keep everything interactive in ROOT Editor.
}
