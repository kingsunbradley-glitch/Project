#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

namespace {

// Match ../242Fm_plot/Combined_Spectrum_modified.C.
constexpr int kFontCode = 133;
constexpr double kTextSize = 44.0;
constexpr double kPanelLabelSize = 48.0;
constexpr int kFrameLineWidth = 2;

constexpr double kXMinMeV = 7.0;
constexpr double kXMaxMeV = 9.5;
constexpr double kYMinSeconds = 50.0e-9;
constexpr double kYMaxSeconds = 18.0;

void StyleAxis(TAxis *axis)
{
    axis->SetTitleFont(kFontCode);
    axis->SetLabelFont(kFontCode);
    axis->SetTitleSize(kTextSize);
    axis->SetLabelSize(kTextSize);
    axis->SetAxisColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitleColor(kBlack);
    axis->CenterTitle(true);
}

void StyleFrame(TH2D *frame, bool showXAxis)
{
    frame->SetStats(false);
    frame->SetLineWidth(kFrameLineWidth);

    frame->GetYaxis()->SetTitle("");
    frame->GetYaxis()->SetLabelOffset(0.015);
    frame->GetYaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetMoreLogLabels(false);
    frame->GetYaxis()->SetNoExponent(false);

    frame->GetXaxis()->SetNdivisions(505);
    if (showXAxis) {
        frame->GetXaxis()->SetTitle("Energy (MeV)");
        frame->GetXaxis()->SetTitleOffset(1.35);
        frame->GetXaxis()->SetLabelOffset(0.020);
        // Keep one decimal place on the shared X axis: 7.0, 7.5, 8.0, ...
        frame->GetXaxis()->SetDecimals(true);
    } else {
        frame->GetXaxis()->SetTitle("");
        frame->GetXaxis()->SetLabelSize(0.0);
        frame->GetXaxis()->SetLabelOffset(0.0);
    }

    StyleAxis(frame->GetXaxis());
    StyleAxis(frame->GetYaxis());
}

TGraph *MakeScatterGraph(TTree *tree, const char *name)
{
    if (!tree)
        return nullptr;

    double energyKeV = 0.0;
    double deltaTSeconds = 0.0;
    tree->SetBranchStatus("*", false);
    tree->SetBranchStatus("SumE", true);
    tree->SetBranchStatus("DeltaT_s", true);
    tree->SetBranchAddress("SumE", &energyKeV);
    tree->SetBranchAddress("DeltaT_s", &deltaTSeconds);

    const Long64_t entries = tree->GetEntries();
    TGraph *graph = new TGraph(static_cast<int>(entries));
    graph->SetName(name);

    int kept = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        if (tree->GetEntry(entry) <= 0)
            continue;

        const double energyMeV = energyKeV / 1000.0;
        if (energyMeV < kXMinMeV || energyMeV > kXMaxMeV)
            continue;
        if (deltaTSeconds < kYMinSeconds || deltaTSeconds > kYMaxSeconds)
            continue;

        graph->SetPoint(kept++, energyMeV, deltaTSeconds);
    }
    graph->Set(kept);

    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", true);

    std::cout << "[Info] " << name << ": kept " << kept << " / "
              << entries << " events (continuous Y values)" << std::endl;
    return graph;
}

void StyleGraph(TGraph *graph, Color_t color)
{
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(7);
    graph->SetMarkerSize(0.50);
}

TH2D *MakeFrame(const char *name)
{
    TH2D *frame = new TH2D(
        name,
        "",
        125,
        kXMinMeV,
        kXMaxMeV,
        100,
        kYMinSeconds,
        kYMaxSeconds);
    frame->SetDirectory(nullptr);
    return frame;
}

void DrawPanelLabel(TPad *pad, const char *label)
{
    pad->cd();
    const double frameLeft = pad->GetLeftMargin();
    const double frameRight = 1.0 - pad->GetRightMargin();
    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();

    TLatex *panelLabel = new TLatex();
    panelLabel->SetNDC();
    panelLabel->SetTextFont(kFontCode);
    panelLabel->SetTextSize(kPanelLabelSize);
    panelLabel->SetTextAlign(31);
    panelLabel->DrawLatex(
        frameLeft + 0.97 * (frameRight - frameLeft),
        frameBottom + 0.06 * (frameTop - frameBottom),
        label);
}

}  // namespace

void Draw_v4_combined()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetTextFont(kFontCode);
    gStyle->SetLabelFont(kFontCode, "XYZ");
    gStyle->SetTitleFont(kFontCode, "XYZ");
    gStyle->SetFrameLineWidth(kFrameLineWidth);
    gStyle->SetLineScalePS(1.0);

    const TString macroDirectory = gSystem->DirName(__FILE__);
    const TString inputPath = macroDirectory + "/plot_data.root";
    TFile *input = TFile::Open(inputPath, "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "[Error] Cannot open " << inputPath << std::endl;
        delete input;
        return;
    }

    TTree *lowGreen = nullptr;
    TTree *lowRed = nullptr;
    input->GetObject("tree_low_green", lowGreen);
    input->GetObject("tree_low_red", lowRed);
    if (!lowGreen || !lowRed) {
        std::cerr << "[Error] Cannot read tree_low_green/tree_low_red"
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    TGraph *greenGraph =
        MakeScatterGraph(lowGreen, "g_combined_green_7_9p5MeV");
    TGraph *redGraph =
        MakeScatterGraph(lowRed, "g_combined_red_7_9p5MeV");

    input->Close();
    delete input;
    if (!greenGraph || !redGraph)
        return;

    StyleGraph(greenGraph, kGreen + 2);
    StyleGraph(redGraph, kRed + 1);

    // This boundary makes the two physical plotting frames equal in height
    // for the margins below, as in Combined_Spectrum_modified.C.
    constexpr double panelBoundary = 0.5367;
    TCanvas *canvas = new TCanvas(
        "c_green_red_combined", "", 1100, 1250);
    TPad *upperPad = new TPad(
        "p_green_upper", "", 0.0, panelBoundary, 1.0, 1.0);
    TPad *lowerPad = new TPad(
        "p_red_lower", "", 0.0, 0.0, 1.0, panelBoundary);

    upperPad->SetLeftMargin(0.15);
    upperPad->SetRightMargin(0.04);
    upperPad->SetTopMargin(0.05);
    upperPad->SetBottomMargin(0.00);
    upperPad->SetFrameLineWidth(kFrameLineWidth);
    upperPad->SetTickx(1);
    upperPad->SetTicky(1);
    upperPad->SetLogy(true);

    lowerPad->SetLeftMargin(0.15);
    lowerPad->SetRightMargin(0.04);
    lowerPad->SetTopMargin(0.00);
    lowerPad->SetBottomMargin(0.18);
    lowerPad->SetFrameLineWidth(kFrameLineWidth);
    lowerPad->SetTickx(1);
    lowerPad->SetTicky(1);
    lowerPad->SetLogy(true);

    upperPad->Draw();
    lowerPad->Draw();

    upperPad->cd();
    TH2D *greenFrame = MakeFrame("h_combined_green_frame");
    StyleFrame(greenFrame, false);
    greenFrame->Draw();
    greenGraph->Draw("P SAME");
    DrawPanelLabel(upperPad, "(a)");
    upperPad->Modified();
    upperPad->Update();
    upperPad->RedrawAxis();

    lowerPad->cd();
    TH2D *redFrame = MakeFrame("h_combined_red_frame");
    StyleFrame(redFrame, true);
    redFrame->Draw();
    redGraph->Draw("P SAME");
    DrawPanelLabel(lowerPad, "(b)");
    lowerPad->Modified();
    lowerPad->Update();
    lowerPad->RedrawAxis();

    // Transparent full-canvas overlay for one vertically centered Y title.
    canvas->cd();
    TPad *overlay = new TPad("p_combined_overlay", "", 0, 0, 1, 1);
    overlay->SetFillStyle(4000);
    overlay->SetFrameFillStyle(0);
    overlay->SetFrameLineColor(0);
    overlay->SetBorderMode(0);
    overlay->SetBorderSize(0);
    overlay->SetMargin(0, 0, 0, 0);
    overlay->Draw();
    overlay->cd();

    TLatex *commonYTitle = new TLatex();
    commonYTitle->SetNDC();
    commonYTitle->SetTextFont(kFontCode);
    commonYTitle->SetTextSize(kTextSize);
    commonYTitle->SetTextAlign(22);
    commonYTitle->SetTextAngle(90.0);
    commonYTitle->DrawLatex(
        0.040, 0.50, "#Delta#it{t}(ER-#alpha/SF) (s)");

    canvas->cd();
    canvas->Modified();
    canvas->Update();

    const TString outputBase =
        macroDirectory + "/Figure12_combined_7_9p5MeV";
    for (const char *extension : {"pdf", "eps", "png"})
        canvas->SaveAs(outputBase + "." + extension);

    std::cout << "[Done] Saved PDF/EPS/PNG to " << outputBase
              << ".*" << std::endl;
}
