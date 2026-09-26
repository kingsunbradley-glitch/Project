#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
// c_all_combined_mac->SetCanvasSize(1100, 1600);
namespace {

constexpr int kFontCode = 133;
constexpr double kTextSize = 44.0;
constexpr double kPanelLabelSize = 48.0;
constexpr int kFrameLineWidth = 2;

constexpr double kLowXMinMeV = 7.0;
constexpr double kLowXMaxMeV = 9.5;
constexpr double kFullXMinMeV = 0.0;
constexpr double kFullXMaxMeV = 250.0;
constexpr double kYMinSeconds = 1.0e-8;
constexpr double kYMaxSeconds = 30.0;
// Equivalent to the original-data cut Delta_Ts[1] < 12e9 (Delta_Ts in ns).
constexpr double kDeltaTSelectionMaxSeconds = 12.0;
constexpr double kYMajorTickLength = 0.010;
constexpr double kYMinorTickLength = 0.006;
constexpr double kYLabelGap = 0.010;
constexpr double kLeftMargin = 0.13;

// macOS ROOT may clamp a tall TCanvas window to the screen height.  The
// dedicated macOS entry point below restores the requested drawable size with
// SetCanvasSize().  Its dense full-energy scatter plots are also reduced to
// one point per display cell: points in the same cell are visually
// indistinguishable, but every occupied part of the plot remains represented.
constexpr int kMacCanvasWidth = 1100;
constexpr int kMacCanvasHeight = 1600;
constexpr Long64_t kMacReductionThreshold = 100000;
constexpr int kMacDisplayXBins = 440;
constexpr int kMacDisplayYBins = 260;

// Canvas layout.  The full-energy panel is on top; the two low-energy panels
// share their X axis below it.  Pad heights are derived so all three frames
// have exactly the same physical height despite their different margins.
constexpr double kMiddleGap = 0.004;
constexpr double kTopPanelTopMargin = 0.05;
constexpr double kTopPanelBottomMargin = 0.25;
constexpr double kMiddlePanelTopMargin = 0.05;
constexpr double kMiddlePanelBottomMargin = 0.00;
constexpr double kBottomPanelTopMargin = 0.00;
constexpr double kBottomPanelBottomMargin = 0.27;
constexpr double kTopPanelFrameFraction =
    1.0 - kTopPanelTopMargin - kTopPanelBottomMargin;
constexpr double kMiddlePanelFrameFraction =
    1.0 - kMiddlePanelTopMargin - kMiddlePanelBottomMargin;
constexpr double kBottomPanelFrameFraction =
    1.0 - kBottomPanelTopMargin - kBottomPanelBottomMargin;
constexpr double kCommonFrameHeight =
    (1.0 - kMiddleGap) /
    (1.0 / kTopPanelFrameFraction +
     1.0 / kMiddlePanelFrameFraction +
     1.0 / kBottomPanelFrameFraction);
constexpr double kTopPanelHeight =
    kCommonFrameHeight / kTopPanelFrameFraction;
constexpr double kMiddlePanelHeight =
    kCommonFrameHeight / kMiddlePanelFrameFraction;
constexpr double kBottomPanelHeight =
    kCommonFrameHeight / kBottomPanelFrameFraction;
constexpr double kBottomPadTopY = kBottomPanelHeight;
constexpr double kMiddlePadTopY =
    kBottomPadTopY + kMiddlePanelHeight;
constexpr double kTopPadBottomY =
    kMiddlePadTopY + kMiddleGap;

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

void StylePad(
    TPad *pad,
    double topMargin,
    double bottomMargin)
{
    pad->SetLeftMargin(kLeftMargin);
    pad->SetRightMargin(0.04);
    pad->SetTopMargin(topMargin);
    pad->SetBottomMargin(bottomMargin);
    pad->SetFrameLineWidth(kFrameLineWidth);
    pad->SetTickx(1);
    pad->SetTicky(0);
    pad->SetLogy(true);
}

void StyleFrame(
    TH2D *frame,
    bool showXAxis,
    bool fixedOneDecimal)
{
    frame->SetStats(false);
    frame->SetLineWidth(kFrameLineWidth);

    StyleAxis(frame->GetXaxis());
    StyleAxis(frame->GetYaxis());

    frame->GetYaxis()->SetTitle("");
    frame->GetYaxis()->SetLabelSize(0.0);
    frame->GetYaxis()->SetLabelOffset(0.0);
    frame->GetYaxis()->SetTickLength(0.0);
    frame->GetYaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetMoreLogLabels(false);
    frame->GetYaxis()->SetNoExponent(false);

    frame->GetXaxis()->SetNdivisions(505);
    if (showXAxis) {
        frame->GetXaxis()->SetTitle("Energy (MeV)");
        frame->GetXaxis()->SetTitleOffset(1.35);
        frame->GetXaxis()->SetLabelOffset(0.020);
        frame->GetXaxis()->SetDecimals(fixedOneDecimal);
    } else {
        frame->GetXaxis()->SetTitle("");
        frame->GetXaxis()->SetLabelSize(0.0);
        frame->GetXaxis()->SetLabelOffset(0.0);
    }

}

TH2D *MakeFrame(
    const char *name,
    double xMinMeV,
    double xMaxMeV,
    int nXBins)
{
    TH2D *frame = new TH2D(
        name,
        "",
        nXBins,
        xMinMeV,
        xMaxMeV,
        100,
        kYMinSeconds,
        kYMaxSeconds);
    frame->SetDirectory(nullptr);
    return frame;
}

TGraph *MakeScatterGraph(
    TTree *tree,
    const char *name,
    double xMinMeV,
    double xMaxMeV,
    bool applyFigure3Cut,
    bool reduceDensePoints)
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
    const bool useDisplayReduction =
        reduceDensePoints && entries >= kMacReductionThreshold;
    TGraph *graph = useDisplayReduction
        ? new TGraph()
        : new TGraph(static_cast<int>(entries));
    graph->SetName(name);

    std::vector<unsigned char> occupiedDisplayCells;
    if (useDisplayReduction) {
        occupiedDisplayCells.assign(
            kMacDisplayXBins * kMacDisplayYBins, 0);
    }

    int kept = 0;
    int selected = 0;
    int removedByDisplayReduction = 0;
    int removedByDeltaTCut = 0;
    int removedByFigure3Cut = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        if (tree->GetEntry(entry) <= 0)
            continue;

        const double energyMeV = energyKeV / 1000.0;
        if (energyMeV < xMinMeV || energyMeV > xMaxMeV)
            continue;
        if (deltaTSeconds < kYMinSeconds)
            continue;
        if (deltaTSeconds >= kDeltaTSelectionMaxSeconds) {
            ++removedByDeltaTCut;
            continue;
        }

        if (applyFigure3Cut &&
            energyMeV >= 10.0 && energyMeV <= 20.0 &&
            deltaTSeconds < 1.0e-6) {
            ++removedByFigure3Cut;
            continue;
        }

        ++selected;
        if (useDisplayReduction) {
            int xBin = static_cast<int>(
                (energyMeV - xMinMeV) / (xMaxMeV - xMinMeV) *
                kMacDisplayXBins);
            int yBin = static_cast<int>(
                (std::log10(deltaTSeconds) - std::log10(kYMinSeconds)) /
                (std::log10(kYMaxSeconds) - std::log10(kYMinSeconds)) *
                kMacDisplayYBins);
            xBin = std::max(0, std::min(kMacDisplayXBins - 1, xBin));
            yBin = std::max(0, std::min(kMacDisplayYBins - 1, yBin));
            const int cell = yBin * kMacDisplayXBins + xBin;
            if (occupiedDisplayCells[cell]) {
                ++removedByDisplayReduction;
                continue;
            }
            occupiedDisplayCells[cell] = 1;
        }

        graph->SetPoint(kept++, energyMeV, deltaTSeconds);
    }
    graph->Set(kept);

    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", true);

    std::cout << "[Info] " << name << ": selected " << selected << " / "
              << entries << ", drawing " << kept
              << "; Delta_Ts[1] < 12e9 cut removed "
              << removedByDeltaTCut;
    if (applyFigure3Cut)
        std::cout << "; Figure 3 cut removed " << removedByFigure3Cut;
    if (useDisplayReduction)
        std::cout << "; display-overlap reduction removed "
                  << removedByDisplayReduction;
    std::cout << std::endl;
    return graph;
}

void StyleGraph(TGraph *graph, Color_t color)
{
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    // Marker style 7 has a fixed pixel size; use a scalable filled circle so
    // SetMarkerSize actually makes the scatter points larger.
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
}

double UserXToNDC(TPad *pad, double x)
{
    const double frameLeft = pad->GetLeftMargin();
    const double frameWidth =
        1.0 - pad->GetLeftMargin() - pad->GetRightMargin();
    const double padX = pad->XtoPad(x);
    return frameLeft + frameWidth *
        (padX - pad->GetUxmin()) /
        (pad->GetUxmax() - pad->GetUxmin());
}

double UserYToNDC(TPad *pad, double y)
{
    const double frameBottom = pad->GetBottomMargin();
    const double frameHeight =
        1.0 - pad->GetBottomMargin() - pad->GetTopMargin();
    // YtoPad applies log10 when the pad uses a logarithmic Y axis.
    const double padY = pad->YtoPad(y);
    return frameBottom + frameHeight *
        (padY - pad->GetUymin()) /
        (pad->GetUymax() - pad->GetUymin());
}

void DrawNuclideLabel(
    double xMeV,
    double ySeconds,
    const char *labelText)
{
    TPad *pad = static_cast<TPad *>(gPad);
    const double xNDC = UserXToNDC(pad, xMeV);
    const double yNDC = UserYToNDC(pad, ySeconds);

    // TLatex hit testing is unreliable in nested NDC pads in ROOT's macOS
    // GUI.  A transparent TPaveText preserves the same appearance while
    // providing a generous, stable selection box for the GUI Editor.
    constexpr double halfWidthNDC = 0.065;
    constexpr double halfHeightNDC = 0.055;
    TPaveText *label = new TPaveText(
        xNDC - halfWidthNDC,
        yNDC - halfHeightNDC,
        xNDC + halfWidthNDC,
        yNDC + halfHeightNDC,
        "NDC");
    label->SetName(TString::Format("nuclide_%s", labelText));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetLineWidth(0);
    label->SetMargin(0.0);
    label->SetTextFont(kFontCode);
    label->SetTextSize(kTextSize);
    label->SetTextAlign(22);
    label->AddText(labelText);
    label->Draw();
}

void DrawNuclideLabelWithArrow(
    double labelXMeV,
    double labelYSeconds,
    double targetXMeV,
    double targetYSeconds,
    const char *labelText)
{
    TPad *pad = static_cast<TPad *>(gPad);
    TArrow *arrow = new TArrow(
        UserXToNDC(pad, labelXMeV),
        UserYToNDC(pad, labelYSeconds),
        UserXToNDC(pad, targetXMeV),
        UserYToNDC(pad, targetYSeconds),
        0.012,
        "|>");
    arrow->SetNDC();
    arrow->SetLineColor(kBlack);
    arrow->SetFillColor(kBlack);
    arrow->SetLineWidth(kFrameLineWidth);
    arrow->Draw();
    DrawNuclideLabel(labelXMeV, labelYSeconds, labelText);
}

void DrawDashedEllipse(
    TPad *pad,
    double xMinMeV,
    double xMaxMeV,
    double yMinSeconds,
    double yMaxSeconds)
{
    const double xMinNDC = UserXToNDC(pad, xMinMeV);
    const double xMaxNDC = UserXToNDC(pad, xMaxMeV);
    const double yMinNDC = UserYToNDC(pad, yMinSeconds);
    const double yMaxNDC = UserYToNDC(pad, yMaxSeconds);
    const double centerX = 0.5 * (xMinNDC + xMaxNDC);
    const double centerY = 0.5 * (yMinNDC + yMaxNDC);
    const double radiusX = 0.5 * (xMaxNDC - xMinNDC);
    const double radiusY = 0.5 * (yMaxNDC - yMinNDC);

    constexpr int nEllipsePoints = 161;
    constexpr double twoPi = 6.2831853071795864769;
    TPolyLine *ellipse = new TPolyLine(nEllipsePoints);
    ellipse->SetNDC();
    for (int point = 0; point < nEllipsePoints; ++point) {
        const double angle =
            twoPi * point / static_cast<double>(nEllipsePoints - 1);
        ellipse->SetPoint(
            point,
            centerX + radiusX * std::cos(angle),
            centerY + radiusY * std::sin(angle));
    }
    ellipse->SetFillStyle(0);
    ellipse->SetLineColor(kBlack);
    ellipse->SetLineStyle(2);
    ellipse->SetLineWidth(kFrameLineWidth);
    ellipse->Draw();
}

void DrawRunLabel(TPad *pad, const char *labelText)
{
    pad->cd();
    const double frameLeft = pad->GetLeftMargin();
    const double frameRight = 1.0 - pad->GetRightMargin();
    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();

    TLatex *label = new TLatex();
    label->SetNDC();
    label->SetTextFont(kFontCode);
    label->SetTextSize(kTextSize);
    label->SetTextAlign(11);
    label->DrawLatex(
        frameLeft + 0.05 * (frameRight - frameLeft),
        frameBottom + 0.055 * (frameTop - frameBottom),
        labelText);
}

void DrawPanelLabel(TPad *pad, const char *label, bool upperRight)
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
    panelLabel->SetTextAlign(upperRight ? 33 : 31);
    panelLabel->DrawLatex(
        frameLeft + 0.97 * (frameRight - frameLeft),
        frameBottom + (upperRight ? 0.94 : 0.06) *
        (frameTop - frameBottom),
        label);
}

void DrawManualLogYTicks(TPad *pad)
{
    pad->cd();
    pad->Update();

    const double frameLeft = pad->GetLeftMargin();
    const double frameRight = 1.0 - pad->GetRightMargin();
    const double frameBottom = pad->GetBottomMargin();
    const double frameTop = 1.0 - pad->GetTopMargin();
    const double logMin = std::log10(kYMinSeconds);
    const double logMax = std::log10(kYMaxSeconds);
    const int decadeMin = static_cast<int>(std::floor(logMin));
    const int decadeMax = static_cast<int>(std::ceil(logMax));

    for (int decade = decadeMin; decade <= decadeMax; ++decade) {
        const double base = std::pow(10.0, decade);
        for (int multiplier = 1; multiplier <= 9; ++multiplier) {
            const double value = multiplier * base;
            if (value < kYMinSeconds || value > kYMaxSeconds)
                continue;

            const double fraction =
                (std::log10(value) - logMin) / (logMax - logMin);
            const double yNDC =
                frameBottom + fraction * (frameTop - frameBottom);
            const double tickLength =
                multiplier == 1 ? kYMajorTickLength : kYMinorTickLength;

            // Match Combined_Spectrum_modified.C: left ticks point outward,
            // while the mirrored right-side ticks point into the frame.
            TLine *leftTick = new TLine(
                frameLeft - tickLength, yNDC, frameLeft, yNDC);
            leftTick->SetNDC();
            leftTick->SetLineColor(kBlack);
            leftTick->SetLineWidth(kFrameLineWidth);
            leftTick->Draw();

            TLine *rightTick = new TLine(
                frameRight, yNDC, frameRight - tickLength, yNDC);
            rightTick->SetNDC();
            rightTick->SetLineColor(kBlack);
            rightTick->SetLineWidth(kFrameLineWidth);
            rightTick->Draw();
        }

    }
}

void DrawManualLogYLabels(
    double padBottomY,
    double padTopY,
    double topMargin,
    double bottomMargin)
{
    const double logMin = std::log10(kYMinSeconds);
    const double logMax = std::log10(kYMaxSeconds);
    const int decadeMin = static_cast<int>(std::floor(logMin));
    const int decadeMax = static_cast<int>(std::ceil(logMax));

    for (int decade = decadeMin; decade <= decadeMax; ++decade) {
        // Label every second decade: 10, 10^-1, 10^-3, ...; omit 1.
        const double labelValue = std::pow(10.0, decade);
        if (decade % 2 == 0 ||
            labelValue < kYMinSeconds || labelValue > kYMaxSeconds)
            continue;

        const double fraction =
            (std::log10(labelValue) - logMin) / (logMax - logMin);
        const double localY =
            bottomMargin +
            fraction * (1.0 - topMargin - bottomMargin);
        const double canvasY =
            padBottomY + localY * (padTopY - padBottomY);

        TLatex *tickLabel = new TLatex();
        tickLabel->SetNDC();
        tickLabel->SetTextFont(kFontCode);
        tickLabel->SetTextSize(kTextSize);
        tickLabel->SetTextAlign(32);
        const TString labelText =
            decade == 1 ? "10" : TString::Format("10^{%d}", decade);
        tickLabel->DrawLatex(
            kLeftMargin - kYMajorTickLength - kYLabelGap,
            canvasY,
            labelText);
    }
}

void DrawCommonYTitle(double centerY)
{
    TLatex *title = new TLatex();
    title->SetNDC();
    title->SetTextFont(kFontCode);
    title->SetTextSize(kTextSize);
    title->SetTextAlign(22);
    title->SetTextAngle(90.0);
    title->DrawLatex(
        0.022, centerY, "#Delta#it{t}(ER-#alpha/SF) (s)");
}

}  // namespace

void Draw_v4_all_combined_impl(bool macOptimized)
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
    TTree *fullGreen = nullptr;
    TTree *fullRed = nullptr;
    input->GetObject("tree_low_green", lowGreen);
    input->GetObject("tree_low_red", lowRed);
    input->GetObject("tree_full_green", fullGreen);
    input->GetObject("tree_full_red", fullRed);
    if (!lowGreen || !lowRed || !fullGreen || !fullRed) {
        std::cerr << "[Error] plot_data.root is missing a required tree."
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    TGraph *greenLow = MakeScatterGraph(
        lowGreen,
        "g_all_green_low",
        kLowXMinMeV,
        kLowXMaxMeV,
        false,
        macOptimized);
    TGraph *redLow = MakeScatterGraph(
        lowRed,
        "g_all_red_low",
        kLowXMinMeV,
        kLowXMaxMeV,
        false,
        macOptimized);
    TGraph *greenFull = MakeScatterGraph(
        fullGreen,
        "g_all_green_full",
        kFullXMinMeV,
        kFullXMaxMeV,
        true,
        macOptimized);
    TGraph *redFull = MakeScatterGraph(
        fullRed,
        "g_all_red_full",
        kFullXMinMeV,
        kFullXMaxMeV,
        true,
        macOptimized);

    input->Close();
    delete input;
    if (!greenLow || !redLow || !greenFull || !redFull)
        return;

    StyleGraph(greenLow, kGreen + 2);
    StyleGraph(redLow, kRed + 1);
    StyleGraph(greenFull, kGreen + 2);
    StyleGraph(redFull, kRed + 1);

    TCanvas *canvas = new TCanvas(
        macOptimized ? "c_all_combined_mac" : "c_all_combined",
        "",
        kMacCanvasWidth,
        kMacCanvasHeight);
    if (macOptimized) {
        // Cocoa limits the initial window to the visible screen.  This call
        // fixes the drawable area itself, so exported files retain the aspect
        // ratio for which the pixel-font sizes and pad margins were designed.
        canvas->SetCanvasSize(kMacCanvasWidth, kMacCanvasHeight);
    }

    TPad *panelA = new TPad(
        "p_all_a", "", 0.0, kTopPadBottomY, 1.0, 1.0);
    TPad *panelB = new TPad(
        "p_all_b", "", 0.0, kBottomPadTopY, 1.0, kMiddlePadTopY);
    TPad *panelC = new TPad(
        "p_all_c", "", 0.0, 0.0, 1.0, kBottomPadTopY);

    StylePad(panelA, kTopPanelTopMargin, kTopPanelBottomMargin);
    StylePad(panelB, kMiddlePanelTopMargin, kMiddlePanelBottomMargin);
    StylePad(panelC, kBottomPanelTopMargin, kBottomPanelBottomMargin);

    panelA->Draw();
    panelB->Draw();
    panelC->Draw();

    panelA->cd();
    TH2D *frameA = MakeFrame(
        "h_all_frame_a", kFullXMinMeV, kFullXMaxMeV, 500);
    StyleFrame(frameA, true, false);
    frameA->Draw();
    redFull->Draw("P SAME");
    greenFull->Draw("P SAME");
    panelA->Modified();
    panelA->Update();
    panelA->RedrawAxis();
    DrawManualLogYTicks(panelA);
    // Keep annotations at the end of the pad primitive list.  Besides drawing
    // them above the data and axes, this makes them straightforward to select
    // and move with ROOT's GUI Editor.
    // Highlight the delayed events between 100 and 230 MeV below 10^-3 s.
    // The lower edge is placed below the lowest selected point, so all events
    // in this cluster are enclosed without touching the plot frame.
    DrawDashedEllipse(panelA, 100.0, 210.0, 3.0e-7, 1.0e-4);
    DrawNuclideLabel(150.0, 0.1, "^{246}Fm");
    DrawNuclideLabel(230.0, 1.0e-5, "^{242}Fm");
    DrawRunLabel(panelA, "Run 1 + Run 2");
    DrawPanelLabel(panelA, "(a)", true);
    panelA->Modified();
    panelA->Update();

    panelB->cd();
    TH2D *frameB = MakeFrame(
        "h_all_frame_b", kLowXMinMeV, kLowXMaxMeV, 125);
    StyleFrame(frameB, false, false);
    frameB->Draw();
    greenLow->Draw("P SAME");
    panelB->Modified();
    panelB->Update();
    panelB->RedrawAxis();
    DrawManualLogYTicks(panelB);
    DrawNuclideLabel(7.3, 5.0e-4, "^{211}Po");
    DrawNuclideLabel(7.6, 1.0e-3, "^{212}At");
    DrawNuclideLabelWithArrow(
        8.2, 1.0e-3, 8.05, 5.0e-2, "^{213}Rn");
    DrawNuclideLabel(8.3, 1.0e-7, "^{213}Po");
    DrawNuclideLabel(8.05, 1.0e-7, "^{215}At");
    DrawNuclideLabel(8.55, 1.0e-6, "^{215}Rn");
    DrawNuclideLabelWithArrow(
        8.5, 1.0e-1, 8.2, 1.0e-1, "^{246}Fm");
    DrawNuclideLabel(9.0, 1.0e-6, "^{212}Po");
    DrawNuclideLabel(9.05, 1.0e-3, "^{216}Ac");
    DrawRunLabel(panelB, "Run 1");
    DrawPanelLabel(panelB, "(b)", true);
    panelB->Modified();
    panelB->Update();

    panelC->cd();
    TH2D *frameC = MakeFrame(
        "h_all_frame_c", kLowXMinMeV, kLowXMaxMeV, 125);
    StyleFrame(frameC, true, true);
    frameC->Draw();
    redLow->Draw("P SAME");
    panelC->Modified();
    panelC->Update();
    panelC->RedrawAxis();
    DrawManualLogYTicks(panelC);
    DrawNuclideLabelWithArrow(
        8.1, 1.0e-5, 8.05, 1.0e-4, "^{213}Rn");
    DrawRunLabel(panelC, "Run 2");
    DrawPanelLabel(panelC, "(c)", false);
    panelC->Modified();
    panelC->Update();

    // Full-canvas transparent layer keeps both Y titles aligned and fully
    // inside the output while leaving the middle whitespace untouched.
    canvas->cd();
    TPad *overlay = new TPad("p_all_overlay", "", 0, 0, 1, 1);
    overlay->SetFillStyle(4000);
    overlay->SetFrameFillStyle(0);
    overlay->SetFrameLineColor(0);
    overlay->SetBorderMode(0);
    overlay->SetBorderSize(0);
    overlay->SetMargin(0, 0, 0, 0);
    overlay->SetBit(TObject::kCannotPick);
    overlay->Draw();
    overlay->cd();

    // Draw all Y labels on the final overlay so every label stays centered on
    // its tick and no neighboring pad can cover it.
    DrawManualLogYLabels(
        kTopPadBottomY,
        1.0,
        kTopPanelTopMargin,
        kTopPanelBottomMargin);
    DrawManualLogYLabels(
        kBottomPadTopY,
        kMiddlePadTopY,
        kMiddlePanelTopMargin,
        kMiddlePanelBottomMargin);
    DrawManualLogYLabels(
        0.0,
        kBottomPadTopY,
        kBottomPanelTopMargin,
        kBottomPanelBottomMargin);

    const double topFrameCenter =
        kTopPadBottomY + kTopPanelHeight * 0.5 *
        (kTopPanelBottomMargin + 1.0 - kTopPanelTopMargin);
    DrawCommonYTitle(topFrameCenter);

    const double lowerFramesBottom =
        kBottomPanelHeight * kBottomPanelBottomMargin;
    const double lowerFramesTop =
        kBottomPadTopY +
        kMiddlePanelHeight * (1.0 - kMiddlePanelTopMargin);
    DrawCommonYTitle(0.5 * (lowerFramesBottom + lowerFramesTop));

    canvas->cd();
    canvas->Modified();
    canvas->Update();

    const TString outputBase =
        macroDirectory +
        (macOptimized
             ? "/Figure123_combined_mac"
             : "/Figure123_combined");
    for (const char *extension : {"pdf", "eps", "png"})
        canvas->SaveAs(outputBase + "." + extension);

    std::cout << "[Done] Saved PDF/EPS/PNG to " << outputBase
              << ".*" << std::endl;
}

void Draw_v4_all_combined()
{
    Draw_v4_all_combined_impl(false);
}
