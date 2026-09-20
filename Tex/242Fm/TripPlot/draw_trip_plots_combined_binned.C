#include <TCanvas.h>
#include <TColor.h>
#include <TCut.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2C.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>

namespace {

constexpr int kFont = 133;
constexpr double kCommonTextSize = 44.0;
constexpr double kPanelLabelSize = 48.0;

void StyleAxis(TAxis *axis)
{
    axis->SetTitleFont(kFont);
    axis->SetLabelFont(kFont);
    axis->SetTitleSize(kCommonTextSize);
    axis->SetLabelSize(kCommonTextSize);
    axis->SetAxisColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
    axis->CenterTitle(true);
}

void DrawPanelLabel(TPad *pad, const char *label)
{
    pad->cd();

    const double frameLeft = pad->GetLeftMargin();
    const double frameRight = 1.0 - pad->GetRightMargin();
    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();

    const double x = frameLeft + 0.97 * (frameRight - frameLeft);
    const double y = frameBottom + 0.94 * (frameTop - frameBottom);

    TLatex *panelLabel = new TLatex();
    panelLabel->SetNDC();
    panelLabel->SetTextAlign(33);
    panelLabel->SetTextFont(kFont);
    panelLabel->SetTextSize(kPanelLabelSize);
    panelLabel->DrawLatex(x, y, label);
}

void DrawYTitle(TPad *pad, const char *title, double xNDC)
{
    pad->cd();

    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();
    const double frameCenterY = 0.5 * (frameBottom + frameTop);

    TLatex *yTitle = new TLatex();
    yTitle->SetNDC();
    yTitle->SetTextAlign(21);
    yTitle->SetTextAngle(90.0);
    yTitle->SetTextFont(kFont);
    yTitle->SetTextSize(kCommonTextSize);
    yTitle->DrawLatex(xNDC, frameCenterY, title);
}

TGraph *OccupiedBinsToGraph(const TH2C *occupancy, const char *name)
{
    int occupiedBins = 0;
    for (int xBin = 1; xBin <= occupancy->GetNbinsX(); ++xBin) {
        for (int yBin = 1; yBin <= occupancy->GetNbinsY(); ++yBin) {
            if (occupancy->GetBinContent(xBin, yBin) > 0.0)
                ++occupiedBins;
        }
    }

    TGraph *graph = new TGraph(occupiedBins);
    graph->SetName(name);

    int point = 0;
    for (int xBin = 1; xBin <= occupancy->GetNbinsX(); ++xBin) {
        const double x = occupancy->GetXaxis()->GetBinCenter(xBin);
        for (int yBin = 1; yBin <= occupancy->GetNbinsY(); ++yBin) {
            if (occupancy->GetBinContent(xBin, yBin) <= 0.0)
                continue;

            const double y = occupancy->GetYaxis()->GetBinCenter(yBin);
            graph->SetPoint(point++, x, y);
        }
    }

    return graph;
}

}  // namespace

void draw_trip_plots_combined_binned(bool save_pdf = false)
{
    TCut cut_base =
        "chain_mul > 1"
        " && Delta_Ts[1] > 0"
        " && DSSD_E[0] > 3000"
        " && DSSD_E[0] < 20000";

    TCut cut_10s =
        cut_base &&
        "Delta_Ts[1] < 10e9";

    const TString macroDir = gSystem->DirName(__FILE__);
    const TString dataPath = macroDir + "/plot_data.root";

    TFile *input = TFile::Open(dataPath, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "[Error] Cannot open " << dataPath << std::endl;
        delete input;
        return;
    }

    auto *tree4_11 = dynamic_cast<TTree *>(input->Get("tree_run4_11"));
    auto *tree12_40 = dynamic_cast<TTree *>(input->Get("tree_run12_40"));
    if (!tree4_11 || !tree12_40) {
        std::cerr << "[Error] plot_data.root is missing a required tree."
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    std::cout << "cut_base: " << cut_base.GetTitle() << '\n'
              << "cut_10s : " << cut_10s.GetTitle() << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(kFont);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    TGaxis::SetMaxDigits(4);

    constexpr double energyMinMeV = 0.0;
    constexpr double energyMaxMeV = 250.0;
    constexpr double logTMin = 1.0;
    constexpr double logTMax = 11.2;

    // Panel (a) retains its exact 20 keV energy binning.
    constexpr double spectrumBinWidthMeV = 0.020;
    constexpr int nSpectrumBins = static_cast<int>(
        (energyMaxMeV - energyMinMeV) / spectrumBinWidthMeV);

    // Panel (b) uses a fine occupancy grid.  This is already finer than the
    // final plotting frame in both directions, but caps the number of plotted
    // objects at 1.1 million instead of one object per original event.
    constexpr int nCorrelationXBins = 1250;
    constexpr int nCorrelationYBins = 880;

    gROOT->cd();
    auto *h4_11 = new TH1D(
        "hBinnedSumE_4_11_20keV",
        "",
        nSpectrumBins,
        energyMinMeV,
        energyMaxMeV);
    auto *h12_40 = new TH1D(
        "hBinnedSumE_12_40_20keV",
        "",
        nSpectrumBins,
        energyMinMeV,
        energyMaxMeV);

    auto *occupancy4_11 = new TH2C(
        "hOccupancy_4_11",
        "",
        nCorrelationXBins,
        energyMinMeV,
        energyMaxMeV,
        nCorrelationYBins,
        logTMin,
        logTMax);
    auto *occupancy12_40 = new TH2C(
        "hOccupancy_12_40",
        "",
        nCorrelationXBins,
        energyMinMeV,
        energyMaxMeV,
        nCorrelationYBins,
        logTMin,
        logTMax);

    tree4_11->Draw(
        "SumE/1000.0>>hBinnedSumE_4_11_20keV", "pass_10s", "goff");
    tree12_40->Draw(
        "SumE/1000.0>>hBinnedSumE_12_40_20keV", "pass_10s", "goff");

    tree4_11->Draw(
        "log10_DeltaT_ns:SumE/1000.0>>hOccupancy_4_11", "", "goff");
    tree12_40->Draw(
        "log10_DeltaT_ns:SumE/1000.0>>hOccupancy_12_40", "", "goff");

    const Long64_t events4_11 = tree4_11->GetEntries();
    const Long64_t events12_40 = tree12_40->GetEntries();

    h4_11->SetDirectory(nullptr);
    h12_40->SetDirectory(nullptr);
    occupancy4_11->SetDirectory(nullptr);
    occupancy12_40->SetDirectory(nullptr);

    input->Close();
    delete input;

    h4_11->SetLineColor(kGreen + 2);
    h4_11->SetLineWidth(2);
    h12_40->SetLineColor(kRed + 1);
    h12_40->SetLineWidth(2);

    const double spectrumMaximum =
        std::max(h4_11->GetMaximum(), h12_40->GetMaximum());
    const double spectrumYMax =
        TMath::Power(10.0, TMath::Ceil(TMath::Log10(2.0 * spectrumMaximum)));

    auto *graph4_11 =
        OccupiedBinsToGraph(occupancy4_11, "gBinnedDeltaT_4_11");
    auto *graph12_40 =
        OccupiedBinsToGraph(occupancy12_40, "gBinnedDeltaT_12_40");

    std::cout << "Panel (b): " << events4_11 << " + "
              << events12_40 << " events -> "
              << graph4_11->GetN() << " + " << graph12_40->GetN()
              << " occupied bins" << std::endl;

    delete occupancy4_11;
    delete occupancy12_40;

    // Use explicit opaque square markers.  Marker style 7 can be rendered as
    // a round, antialiased dot by PNG/PDF backends even when it looks square
    // in an interactive ROOT window.
    graph4_11->SetMarkerStyle(21);
    graph4_11->SetMarkerSize(0.35);
    graph4_11->SetMarkerColor(kGreen + 2);
    graph12_40->SetMarkerStyle(21);
    graph12_40->SetMarkerSize(0.35);
    graph12_40->SetMarkerColor(kRed + 1);

    constexpr double panelBoundary = 0.5367;
    auto *canvas = new TCanvas(
        "cTripPlotsCombinedBinned", "Combined TripPlot (binned)", 1100, 1250);
    auto *upperPad =
        new TPad("pTripUpperBinned", "", 0.0, panelBoundary, 1.0, 1.0);
    auto *lowerPad =
        new TPad("pTripLowerBinned", "", 0.0, 0.0, 1.0, panelBoundary);

    upperPad->SetLeftMargin(0.12);
    upperPad->SetRightMargin(0.04);
    upperPad->SetTopMargin(0.05);
    upperPad->SetBottomMargin(0.00);
    upperPad->SetLogy();
    upperPad->SetTickx(1);
    upperPad->SetTicky(1);
    upperPad->SetFrameLineWidth(2);

    lowerPad->SetLeftMargin(0.12);
    lowerPad->SetRightMargin(0.04);
    lowerPad->SetTopMargin(0.00);
    lowerPad->SetBottomMargin(0.18);
    lowerPad->SetTickx(1);
    lowerPad->SetTicky(1);
    lowerPad->SetFrameLineWidth(2);

    upperPad->Draw();
    lowerPad->Draw();

    upperPad->cd();
    h12_40->SetMinimum(0.5);
    h12_40->SetMaximum(spectrumYMax);
    h12_40->GetXaxis()->SetTitle("");
    h12_40->GetXaxis()->SetLabelSize(0.0);
    h12_40->GetXaxis()->SetLabelOffset(0.0);
    h12_40->GetXaxis()->SetNdivisions(510);
    h12_40->GetYaxis()->SetTitle("");
    h12_40->GetYaxis()->SetNdivisions(510);
    h12_40->GetYaxis()->SetLabelOffset(0.015);
    StyleAxis(h12_40->GetXaxis());
    StyleAxis(h12_40->GetYaxis());

    h12_40->Draw("HIST");
    h4_11->Draw("HIST SAME");

    TLatex latexCut;
    latexCut.SetNDC();
    latexCut.SetTextAlign(13);
    latexCut.SetTextFont(kFont);
    latexCut.SetTextSize(kCommonTextSize);
    latexCut.DrawLatex(
        0.15,
        0.893,
        "#Delta#it{t}(ER-#alpha/SF)<10 s");

    DrawYTitle(upperPad, "Counts / 20 keV", 0.040);
    DrawPanelLabel(upperPad, "(a)");
    upperPad->RedrawAxis();

    lowerPad->cd();
    auto *correlationFrame = new TH2D(
        "hBinnedSumELogTFrame",
        "",
        10,
        energyMinMeV,
        energyMaxMeV,
        10,
        logTMin,
        logTMax);
    correlationFrame->SetDirectory(nullptr);
    correlationFrame->GetXaxis()->SetTitle("Energy (MeV)");
    correlationFrame->GetYaxis()->SetTitle("");
    correlationFrame->GetXaxis()->SetNdivisions(510);
    correlationFrame->GetYaxis()->SetNdivisions(510);
    correlationFrame->GetXaxis()->SetTitleOffset(1.35);
    correlationFrame->GetXaxis()->SetLabelOffset(0.020);
    correlationFrame->GetYaxis()->SetLabelOffset(0.015);
    StyleAxis(correlationFrame->GetXaxis());
    StyleAxis(correlationFrame->GetYaxis());
    correlationFrame->Draw();

    graph12_40->Draw("P SAME");
    graph4_11->Draw("P SAME");

    DrawYTitle(lowerPad, "log_{10} #Delta#it{t}(ER-#alpha/SF)", 0.046);
    DrawPanelLabel(lowerPad, "(b)");
    lowerPad->Modified();
    lowerPad->Update();
    lowerPad->RedrawAxis();

    canvas->cd();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(macroDir + "/SumE_combined_binned.png");
    if (save_pdf)
        canvas->SaveAs(macroDir + "/SumE_combined_binned.pdf");

    std::cout << "[Done] Wrote " << macroDir
              << "/SumE_combined_binned.png" << std::endl;
}
