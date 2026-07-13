#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TStyle.h"

void draw_DSSD_SSD()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameLineWidth(2);

    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadBottomMargin(0.14);
    gStyle->SetPadTopMargin(0.04);

    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.052, "XYZ");
    gStyle->SetTitleOffset(1.10, "X");
    gStyle->SetTitleOffset(1.25, "Y");

    gStyle->SetTickLength(0.025, "XY");
    gStyle->SetNdivisions(510, "XY");

    std::ifstream fin("data.dat");
    if (!fin.is_open()) {
        std::cerr << "Error: cannot open data.dat" << std::endl;
        return;
    }

    std::string line;
    std::getline(fin, line);  // skip header

    std::vector<double> x_zero, y_zero;
    std::vector<double> x_nonzero, y_nonzero;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::stringstream ss(line);

        int chain_no = 0;
        double Elab = 0;
        int runnum = 0;
        int x_strip = 0;
        int y_strip = 0;
        double ER_E = 0;
        double SF_E = 0;
        double DSSD_E1 = 0;
        double SSD_E1 = 0;
        double DeltaT_us = 0;

        ss >> chain_no
           >> Elab
           >> runnum
           >> x_strip
           >> y_strip
           >> ER_E
           >> SF_E
           >> DSSD_E1
           >> SSD_E1
           >> DeltaT_us;

        if (ss.fail()) continue;

        if (DSSD_E1 < 80.0 || DSSD_E1 > 300.0) continue;
        if (SSD_E1  < -5.0 || SSD_E1  > 80.0)  continue;

        if (SSD_E1 == 0.0) {
            x_zero.push_back(DSSD_E1);
            y_zero.push_back(SSD_E1);
        } else {
            x_nonzero.push_back(DSSD_E1);
            y_nonzero.push_back(SSD_E1);
        }
    }

    fin.close();

    TCanvas *c1 = new TCanvas("c1", "c1", 900, 750);
    c1->SetTicks(1, 1);

    TH1F *frame = c1->DrawFrame(80.0, -5.0, 300.0, 80.0);

    frame->GetXaxis()->SetTitle("DSSD Energy (MeV)");
    frame->GetYaxis()->SetTitle("SSD Energy (MeV)");

    frame->GetXaxis()->SetLabelOffset(0.008);
    frame->GetYaxis()->SetLabelOffset(0.008);

    frame->GetXaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetTitleOffset(1.25);

    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetNdivisions(510);

    TGraph *g_zero = new TGraph(x_zero.size());
    for (size_t i = 0; i < x_zero.size(); ++i) {
        g_zero->SetPoint(i, x_zero[i], y_zero[i]);
    }

    TGraph *g_nonzero = new TGraph(x_nonzero.size());
    for (size_t i = 0; i < x_nonzero.size(); ++i) {
        g_nonzero->SetPoint(i, x_nonzero[i], y_nonzero[i]);
    }

    // SSD_E1 == 0: solid square
    g_zero->SetMarkerStyle(21);
    g_zero->SetMarkerSize(1.25);
    g_zero->SetMarkerColor(kBlack);
    g_zero->SetLineColor(kBlack);

    // SSD_E1 != 0: solid circle
    g_nonzero->SetMarkerStyle(20);
    g_nonzero->SetMarkerSize(1.10);
    g_nonzero->SetMarkerColor(kBlack);
    g_nonzero->SetLineColor(kBlack);

    g_nonzero->Draw("P SAME");
    g_zero->Draw("P SAME");

    c1->RedrawAxis();

    c1->SaveAs("DSSD_SSD_correlation.pdf");
    c1->SaveAs("DSSD_SSD_correlation.png");

    std::cout << "SSD_E1 == 0 : " << x_zero.size() << std::endl;
    std::cout << "SSD_E1 != 0 : " << x_nonzero.size() << std::endl;
}