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
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

constexpr int kFont = 133;  // Times New Roman, pixel-sized text
constexpr double kCommonTextSize = 44.0;

void StyleAxis(TAxis *axis, double titleSize, double labelSize)
{
    axis->SetTitleFont(kFont);
    axis->SetLabelFont(kFont);
    axis->SetTitleSize(titleSize);
    axis->SetLabelSize(labelSize);
    axis->CenterTitle(true);
    axis->SetAxisColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
}

void StyleCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.12);
    canvas->SetRightMargin(0.04);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.18);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetFrameLineWidth(2);
}

}  // namespace

void draw_trip_plots(bool save_pdf = false)
{
    // These are the cuts used when plot_data.root was produced.  The compact
    // trees in that file already satisfy cut_base; pass_10s stores the result
    // of applying the additional condition in cut_10s.
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
    auto *selectionInfo = dynamic_cast<TNamed *>(input->Get("selection_info"));

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
    if (selectionInfo)
        std::cout << "file selection_info: " << selectionInfo->GetTitle()
                  << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(kFont);
    gStyle->SetLineWidth(2);
    TGaxis::SetMaxDigits(4);

    // ------------------------------------------------------------------
    // Figure 1: SumE spectrum after cut_10s.
    // The stored trees are already selected with cut_base, so pass_10s is
    // precisely the remaining Delta_Ts[1] < 10e9 condition.
    // ------------------------------------------------------------------
    constexpr double energyMinMeV = 0.0;
    constexpr double energyMaxMeV = 250.0;
    constexpr double binWidthMeV = 0.020;  // 20 keV
    constexpr int nEnergyBins =
        static_cast<int>((energyMaxMeV - energyMinMeV) / binWidthMeV);

    gROOT->cd();
    auto *h4_11 = new TH1D(
        "hSumE_4_11_20keV",
        "",
        nEnergyBins,
        energyMinMeV,
        energyMaxMeV);
    auto *h12_40 = new TH1D(
        "hSumE_12_40_20keV",
        "",
        nEnergyBins,
        energyMinMeV,
        energyMaxMeV);

    tree4_11->Draw("SumE/1000.0>>hSumE_4_11_20keV", "pass_10s", "goff");
    tree12_40->Draw("SumE/1000.0>>hSumE_12_40_20keV", "pass_10s", "goff");

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

    auto *spectrumCanvas =
        new TCanvas("cSumESpectrum", "SumE spectrum", 1100, 670);
    StyleCanvas(spectrumCanvas);
    spectrumCanvas->SetLogy();

    h12_40->SetMinimum(0.5);
    h12_40->SetMaximum(spectrumYMax);
    h12_40->GetXaxis()->SetTitle("Energy (MeV)");
    h12_40->GetYaxis()->SetTitle("Counts / 20 keV");
    h12_40->GetXaxis()->SetNdivisions(510);
    h12_40->GetYaxis()->SetNdivisions(510);
    h12_40->GetXaxis()->SetTitleOffset(1.35);
    h12_40->GetYaxis()->SetTitleOffset(1.35);
    h12_40->GetXaxis()->SetLabelOffset(0.020);
    h12_40->GetYaxis()->SetLabelOffset(0.015);
    StyleAxis(
        h12_40->GetXaxis(), kCommonTextSize, kCommonTextSize);
    StyleAxis(
        h12_40->GetYaxis(), kCommonTextSize, kCommonTextSize);

    // Draw the higher-statistics spectrum first so that the green spectrum
    // remains visible where the two distributions overlap.
    h12_40->Draw("HIST");
    h4_11->Draw("HIST SAME");

    TLatex latexCut;
    latexCut.SetNDC();
    latexCut.SetTextFont(kFont);
    latexCut.SetTextSize(kCommonTextSize);
    latexCut.DrawLatex(
        0.15,
        0.76,
        "#Delta#it{t}(ER-#alpha/SF)<10 s");

    spectrumCanvas->RedrawAxis();
    spectrumCanvas->SaveAs(macroDir + "/SumE_spectrum_10s.png");
    if (save_pdf)
        spectrumCanvas->SaveAs(macroDir + "/SumE_spectrum_10s.pdf");

    // ------------------------------------------------------------------
    // Figure 2: event-by-event SumE-log10(DeltaT/ns) correlation after
    // cut_base.  Clone before converting keV to MeV so the input file is
    // treated strictly as a read-only data container.
    // ------------------------------------------------------------------
    auto *graph4_11 =
        dynamic_cast<TGraph *>(graph4_11In->Clone("gDeltaT_4_11_MeV"));
    auto *graph12_40 =
        dynamic_cast<TGraph *>(graph12_40In->Clone("gDeltaT_12_40_MeV"));

    for (int i = 0; i < graph4_11->GetN(); ++i)
        graph4_11->GetX()[i] /= 1000.0;
    for (int i = 0; i < graph12_40->GetN(); ++i)
        graph12_40->GetX()[i] /= 1000.0;

    graph4_11->SetMarkerStyle(7);  // medium pixel dot; thicker and fast for 10^6 points
    graph4_11->SetMarkerSize(1.0);
    graph4_11->SetMarkerColorAlpha(kGreen + 2, 0.50);
    graph12_40->SetMarkerStyle(7);
    graph12_40->SetMarkerSize(1.0);
    graph12_40->SetMarkerColorAlpha(kRed + 1, 0.28);

    input->Close();
    delete input;

    auto *correlationCanvas =
        new TCanvas("cSumELogT", "SumE-logT correlation", 1100, 670);
    StyleCanvas(correlationCanvas);

    auto *correlationFrame = new TH2D(
        "hSumELogTFrame",
        "",
        10,
        energyMinMeV,
        energyMaxMeV,
        10,
        1.0,
        11.2);
    correlationFrame->SetDirectory(nullptr);
    correlationFrame->GetXaxis()->SetTitle("Energy (MeV)");
    correlationFrame->GetYaxis()->SetTitle(
        "log_{10} #Delta t (ER-#alpha/SF)");
    correlationFrame->GetXaxis()->SetNdivisions(510);
    correlationFrame->GetYaxis()->SetNdivisions(510);
    correlationFrame->GetXaxis()->SetTitleOffset(1.35);
    correlationFrame->GetYaxis()->SetTitleOffset(1.05);
    correlationFrame->GetXaxis()->SetLabelOffset(0.020);
    correlationFrame->GetYaxis()->SetLabelOffset(0.015);
    StyleAxis(
        correlationFrame->GetXaxis(), kCommonTextSize, kCommonTextSize);
    StyleAxis(
        correlationFrame->GetYaxis(), kCommonTextSize, kCommonTextSize);
    // Draw an empty 2-D histogram rather than AXIS-only.  This makes ROOT
    // retain both axis titles reliably when millions of markers are painted
    // afterwards in batch mode.
    correlationFrame->Draw();

    // Again draw the high-statistics group first.  Both objects remain true
    // scatter plots; no density conversion or point sampling is performed.
    graph12_40->Draw("P SAME");
    graph4_11->Draw("P SAME");

    correlationCanvas->Modified();
    correlationCanvas->Update();
    correlationCanvas->RedrawAxis();
    correlationCanvas->SaveAs(macroDir + "/SumE_logT_correlation.png");
    if (save_pdf)
        correlationCanvas->SaveAs(macroDir + "/SumE_logT_correlation.pdf");

    std::cout << "[Done] Wrote:\n  "
              << macroDir << "/SumE_spectrum_10s.png\n  "
              << macroDir << "/SumE_logT_correlation.png" << std::endl;
}
