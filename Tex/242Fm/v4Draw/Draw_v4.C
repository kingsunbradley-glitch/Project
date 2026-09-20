#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

namespace {

// Typography, proportions, and line weights follow
// ../242Fm_plot/Combined_Spectrum_modified.C.
constexpr int kFontCode = 133;
constexpr double kTextSize = 44.0;
constexpr int kFrameLineWidth = 2;

constexpr double kLowXMinMeV = 7.0;
constexpr double kLowXMaxMeV = 9.5;

constexpr double kFullXMinMeV = 0.0;
constexpr double kFullXMaxMeV = 250.0;

constexpr double kYMinSeconds = 50.0e-9;
constexpr double kYMaxSeconds = 18.0;

void SetCanvasStyle(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.04);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.18);
    canvas->SetFrameLineWidth(kFrameLineWidth);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->SetLogy(true);
}

void SetAxisStyle(TAxis *axis)
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

void SetFrameStyle(TH2D *frame)
{
    frame->SetStats(false);
    frame->SetLineWidth(kFrameLineWidth);

    frame->GetXaxis()->SetTitle("Energy (MeV)");
    frame->GetYaxis()->SetTitle("#Delta#it{t}(ER-#alpha/SF) (s)");
    frame->GetXaxis()->SetTitleOffset(1.35);
    frame->GetYaxis()->SetTitleOffset(1.85);
    frame->GetXaxis()->SetLabelOffset(0.020);
    frame->GetYaxis()->SetLabelOffset(0.015);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetMoreLogLabels(false);
    frame->GetYaxis()->SetNoExponent(false);

    SetAxisStyle(frame->GetXaxis());
    SetAxisStyle(frame->GetYaxis());
}

// plot_data.root already stores the selected SumE[1] and Delta_Ts[1]
// quantities as the scalar branches SumE [keV] and DeltaT_s [s].
TGraph *MakeScatterGraph(
    TTree *tree,
    const char *name,
    double xMinMeV,
    double xMaxMeV,
    double yMin,
    double yMax,
    bool applyFigure3Cut = false)
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
    int removedByFigure3Cut = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        if (tree->GetEntry(entry) <= 0)
            continue;

        const double energyMeV = energyKeV / 1000.0;
        if (energyMeV < xMinMeV || energyMeV > xMaxMeV)
            continue;
        if (deltaTSeconds < yMin || deltaTSeconds > yMax)
            continue;

        // Figure 3 only: remove the low-time rectangular region at
        // 10--20 MeV and DeltaT < 1 us.
        if (applyFigure3Cut &&
            energyMeV >= 10.0 && energyMeV <= 20.0 &&
            deltaTSeconds < 1.0e-6) {
            ++removedByFigure3Cut;
            continue;
        }

        graph->SetPoint(kept++, energyMeV, deltaTSeconds);
    }
    graph->Set(kept);

    tree->ResetBranchAddresses();
    tree->SetBranchStatus("*", true);

    std::cout << "[Info] " << name << ": kept " << kept << " / "
              << entries << " events (continuous Y values)";
    if (applyFigure3Cut)
        std::cout << "; Figure 3 cut removed " << removedByFigure3Cut;
    std::cout << std::endl;
    return graph;
}

void SetGraphStyle(TGraph *graph, Color_t color)
{
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetMarkerStyle(7);
    graph->SetMarkerSize(0.50);
}

TCanvas *DrawSingleGroup(
    TGraph *graph,
    const char *canvasName,
    Color_t color)
{
    TCanvas *canvas = new TCanvas(canvasName, "", 1100, 600);
    SetCanvasStyle(canvas);

    TH2D *frame = new TH2D(
        TString::Format("h_%s_frame", canvasName),
        "",
        125,
        kLowXMinMeV,
        kLowXMaxMeV,
        100,
        kYMinSeconds,
        kYMaxSeconds);
    frame->SetDirectory(nullptr);
    SetFrameStyle(frame);
    frame->Draw();

    SetGraphStyle(graph, color);
    graph->Draw("P SAME");

    canvas->Modified();
    canvas->Update();
    canvas->RedrawAxis();
    return canvas;
}

TCanvas *DrawBothGroups(TGraph *greenGraph, TGraph *redGraph)
{
    TCanvas *canvas = new TCanvas(
        "c_energy_time", "", 1100, 600);
    SetCanvasStyle(canvas);

    TH2D *frame = new TH2D(
        "h_energy_time_frame",
        "",
        500,
        kFullXMinMeV,
        kFullXMaxMeV,
        100,
        kYMinSeconds,
        kYMaxSeconds);
    frame->SetDirectory(nullptr);
    SetFrameStyle(frame);
    frame->Draw();

    SetGraphStyle(greenGraph, kGreen + 2);
    SetGraphStyle(redGraph, kRed + 1);

    // Draw the higher-statistics red group first so green remains visible.
    redGraph->Draw("P SAME");
    greenGraph->Draw("P SAME");

    canvas->Modified();
    canvas->Update();
    canvas->RedrawAxis();
    return canvas;
}

void SaveCanvas(TCanvas *canvas, const TString &baseName)
{
    for (const char *extension : {"pdf", "eps", "png"})
        canvas->SaveAs(baseName + "." + extension);
}

}  // namespace

void Draw_v4()
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
        std::cerr << "[Error] plot_data.root is missing one or more required "
                     "trees."
                  << std::endl;
        input->Close();
        delete input;
        return;
    }

    TGraph *greenLowGraph = MakeScatterGraph(
        lowGreen,
        "g_green_7_9p5MeV",
        kLowXMinMeV,
        kLowXMaxMeV,
        kYMinSeconds,
        kYMaxSeconds);
    TGraph *redLowGraph = MakeScatterGraph(
        lowRed,
        "g_red_7_9p5MeV",
        kLowXMinMeV,
        kLowXMaxMeV,
        kYMinSeconds,
        kYMaxSeconds);
    TGraph *greenFullGraph = MakeScatterGraph(
        fullGreen,
        "g_green_energy_time",
        kFullXMinMeV,
        kFullXMaxMeV,
        kYMinSeconds,
        kYMaxSeconds,
        true);
    TGraph *redFullGraph = MakeScatterGraph(
        fullRed,
        "g_red_energy_time",
        kFullXMinMeV,
        kFullXMaxMeV,
        kYMinSeconds,
        kYMaxSeconds,
        true);

    input->Close();
    delete input;

    if (!greenLowGraph || !redLowGraph ||
        !greenFullGraph || !redFullGraph) {
        std::cerr << "[Error] Failed to construct one or more plots."
                  << std::endl;
        return;
    }

    TCanvas *greenCanvas = DrawSingleGroup(
        greenLowGraph, "c_green_7_9p5MeV", kGreen + 2);
    TCanvas *redCanvas = DrawSingleGroup(
        redLowGraph, "c_red_7_9p5MeV", kRed + 1);
    TCanvas *combinedCanvas = DrawBothGroups(
        greenFullGraph, redFullGraph);

    SaveCanvas(
        greenCanvas, macroDirectory + "/Figure1_green_7_9p5MeV");
    SaveCanvas(
        redCanvas, macroDirectory + "/Figure2_red_7_9p5MeV");
    SaveCanvas(
        combinedCanvas, macroDirectory + "/Figure3_energy_time");

    std::cout << "\n======================================\n"
              << "Figure 1 points (green): " << greenLowGraph->GetN()
              << "\nFigure 2 points (red):   " << redLowGraph->GetN()
              << "\nFigure 3 points (green): " << greenFullGraph->GetN()
              << "\nFigure 3 points (red):   " << redFullGraph->GetN()
              << "\nSaved PDF/EPS/PNG files in:\n  " << macroDirectory
              << "\n======================================" << std::endl;
}
