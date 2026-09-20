#include <TCanvas.h>
#include <TColor.h>
#include <TCut.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TMath.h>
#include <TNamed.h>
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

void DrawYTitle(TPad *pad, const char *title)
{
    pad->cd();

    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();
    const double frameCenterY = 0.5 * (frameBottom + frameTop);

    TLatex *yTitle = new TLatex();
    yTitle->SetNDC();
    // Align on the font baseline (not the full bounding box), so the "10"
    // subscript in log_{10} does not shift the main title horizontally.
    yTitle->SetTextAlign(21);
    yTitle->SetTextAngle(90.0);
    yTitle->SetTextFont(kFont);
    yTitle->SetTextSize(kCommonTextSize);
    yTitle->DrawLatex(0.046, frameCenterY, title);
}

}  // namespace

void draw_trip_plots_combined(bool save_pdf = false)
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
    auto *graph4_11In = dynamic_cast<TGraph *>(input->Get("g_DeltaT_run4_11"));
    auto *graph12_40In = dynamic_cast<TGraph *>(input->Get("g_DeltaT_run12_40"));

    if (!tree4_11 || !tree12_40 || !graph4_11In || !graph12_40In) {
        std::cerr << "[Error] plot_data.root is missing one or more required "
                     "trees/graphs."
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
    constexpr double binWidthMeV = 0.020;
    constexpr int nEnergyBins =
        static_cast<int>((energyMaxMeV - energyMinMeV) / binWidthMeV);

    gROOT->cd();
    auto *h4_11 = new TH1D(
        "hCombinedSumE_4_11_20keV",
        "",
        nEnergyBins,
        energyMinMeV,
        energyMaxMeV);
    auto *h12_40 = new TH1D(
        "hCombinedSumE_12_40_20keV",
        "",
        nEnergyBins,
        energyMinMeV,
        energyMaxMeV);

    tree4_11->Draw(
        "SumE/1000.0>>hCombinedSumE_4_11_20keV", "pass_10s", "goff");
    tree12_40->Draw(
        "SumE/1000.0>>hCombinedSumE_12_40_20keV", "pass_10s", "goff");

    h4_11->SetDirectory(nullptr);
    h12_40->SetDirectory(nullptr);
    h4_11->SetLineColor(kGreen + 2);
    h4_11->SetLineWidth(2);
    h12_40->SetLineColor(kRed + 1);
    h12_40->SetLineWidth(2);

    const double spectrumMaximum =
        std::max(h4_11->GetMaximum(), h12_40->GetMaximum());
    const double spectrumYMax =
        TMath::Power(10.0, TMath::Ceil(TMath::Log10(2.0 * spectrumMaximum)));

    auto *graph4_11 =
        dynamic_cast<TGraph *>(graph4_11In->Clone("gCombinedDeltaT_4_11_MeV"));
    auto *graph12_40 =
        dynamic_cast<TGraph *>(graph12_40In->Clone("gCombinedDeltaT_12_40_MeV"));

    for (int i = 0; i < graph4_11->GetN(); ++i)
        graph4_11->GetX()[i] /= 1000.0;
    for (int i = 0; i < graph12_40->GetN(); ++i)
        graph12_40->GetX()[i] /= 1000.0;

    graph4_11->SetMarkerStyle(7);
    graph4_11->SetMarkerSize(1.0);
    graph4_11->SetMarkerColorAlpha(kGreen + 2, 0.50);
    graph12_40->SetMarkerStyle(7);
    graph12_40->SetMarkerSize(1.0);
    graph12_40->SetMarkerColorAlpha(kRed + 1, 0.28);

    input->Close();
    delete input;

    // This boundary gives the upper and lower plotting frames the same
    // physical height with the margins below, as in Combined_Spectrum_modified.C.
    constexpr double panelBoundary = 0.5367;

    auto *canvas =
        new TCanvas("cTripPlotsCombined", "Combined TripPlot", 1100, 1250);

    auto *upperPad =
        new TPad("pTripUpper", "", 0.0, panelBoundary, 1.0, 1.0);
    auto *lowerPad =
        new TPad("pTripLower", "", 0.0, 0.0, 1.0, panelBoundary);

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

    // ------------------------------------------------------------------
    // Panel (a): spectrum.  Its X labels/title are hidden because the two
    // panels share the single X axis displayed below panel (b).
    // ------------------------------------------------------------------
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

    DrawYTitle(upperPad, "Counts / 20 keV");
    DrawPanelLabel(upperPad, "(a)");
    upperPad->RedrawAxis();

    // ------------------------------------------------------------------
    // Panel (b): event-by-event energy-time correlation.
    // ------------------------------------------------------------------
    lowerPad->cd();
    auto *correlationFrame = new TH2D(
        "hCombinedSumELogTFrame",
        "",
        10,
        energyMinMeV,
        energyMaxMeV,
        10,
        1.0,
        11.2);
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

    DrawYTitle(lowerPad, "log_{10} #Delta#it{t}(ER-#alpha/SF)");
    DrawPanelLabel(lowerPad, "(b)");
    lowerPad->Modified();
    lowerPad->Update();
    lowerPad->RedrawAxis();

    canvas->cd();
    canvas->Modified();
    canvas->Update();
    canvas->SaveAs(macroDir + "/SumE_combined.png");
    if (save_pdf)
        canvas->SaveAs(macroDir + "/SumE_combined.pdf");

    std::cout << "[Done] Wrote " << macroDir << "/SumE_combined.png"
              << std::endl;
}
