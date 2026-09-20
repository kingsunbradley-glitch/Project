{
    auto serv = new THttpServer("http:48080?loopback");

    // ============================================================
    // USER CONFIGURATION
    // ============================================================
    //
    // All appearance settings modified in this version are collected
    // here so that you only need to edit this section later.
    // ============================================================

    // ------------------------------------------------------------
    // Font
    // ------------------------------------------------------------
    // 133 = Times New Roman, regular, precision 3 (pixel-sized text).
    // A pixel size keeps all text physically identical across pads of
    // different heights and across the whole-canvas overlay.
    int fontCode = 133;
    double commonTextSize = 40.0;
    double panelLabelSize = 48.0; // Slightly larger publication-style (a), (b)

    // Panel boundary chosen so that the ACTUAL plotting frames in (a) and (b)
    // have the same height (and therefore the same aspect ratio).
    // With top/bottom margins used below, the exact value is about 0.5367.
    double panelBoundary = 0.5367;

    // Distance between tick labels and axis.
    // Larger value -> labels farther away from the axis.
    double xLabelOffset = 0.020;
    double yLabelOffset = 0.015;

    // Axis-title offsets.
    // X is slightly increased because the lower X ticks now point outward.
    double xTitleOffset = 1.20;
    double yTitleOffset = 1.05;

    // ------------------------------------------------------------
    // Line widths
    // ------------------------------------------------------------
    int frameLineWidth = 2;   // plot frame / border
    int histLineWidth  = 2;   // histogram lines
    int tickLineWidth  = 2;   // manually drawn tick marks
    int coverLineWidth = 2;   // black lines used to cover artificial baselines

    // ------------------------------------------------------------
    // Panel (a): decay-spectrum line
    // ------------------------------------------------------------
    int decayLineStyle = 1;   // 1=solid, 2=dashed, 3=dotted, ...

    // ------------------------------------------------------------
    // Manual tick lengths, in pad NDC units
    // ------------------------------------------------------------
    // Slightly longer than the previous version.
    double xTickMajor = 0.022;
    double xTickMinor = 0.014;

    // X-axis ticks: one tick every 10 MeV, with a major tick every 50 MeV.
    double xTickStep  = 10.0;
    double xTickMajorStep = 50.0;

    double yTickMajor = 0.013;
    double yTickMinor = 0.008;


    // ============================================================
    // Global style
    // ============================================================

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetTextFont(fontCode);
    gStyle->SetLabelFont(fontCode, "XYZ");
    gStyle->SetTitleFont(fontCode, "XYZ");

    gStyle->SetFrameLineWidth(frameLineWidth);


    // ============================================================
    // Open ROOT files
    // ============================================================

    TFile *fMap = TFile::Open(
        "Map_Total_Decay_Spectrum.root",
        "READ"
    );

    TFile *fDSSD = TFile::Open(
        "DSSDEnergy_run1_run2.root",
        "READ"
    );


    if (!fMap || fMap->IsZombie()) {
        cout
            << "[Error] Cannot open "
            << "Map_Total_Decay_Spectrum.root"
            << endl;

        return;
    }


    if (!fDSSD || fDSSD->IsZombie()) {
        cout
            << "[Error] Cannot open "
            << "DSSDEnergy_run1_run2.root"
            << endl;

        return;
    }


    // ============================================================
    // Get histograms
    //
    // Panel (a)
    //     htotal
    //     hdecay
    //
    // Panel (b)
    //     h_run1
    //     h_run2
    // ============================================================

    TH1D *htotal =
        (TH1D *)fMap->Get("htotal");

    TH1D *hdecay =
        (TH1D *)fMap->Get("hdecay");

    TH1D *h_run1 =
        (TH1D *)fDSSD->Get("h_run1");

    TH1D *h_run2 =
        (TH1D *)fDSSD->Get("h_run2");


    if (!htotal || !hdecay) {
        cout
            << "[Error] Cannot find htotal / hdecay"
            << endl;

        return;
    }


    if (!h_run1 || !h_run2) {
        cout
            << "[Error] Cannot find h_run1 / h_run2"
            << endl;

        return;
    }


    // ============================================================
    // Clone histograms
    // ============================================================

    TH1D *hTotal =
        (TH1D *)htotal->Clone(
            "hTotal_combined"
        );

    TH1D *hDecay =
        (TH1D *)hdecay->Clone(
            "hDecay_combined"
        );

    TH1D *hRun1 =
        (TH1D *)h_run1->Clone(
            "hRun1_combined"
        );

    TH1D *hRun2 =
        (TH1D *)h_run2->Clone(
            "hRun2_combined"
        );


    hTotal->SetDirectory(0);
    hDecay->SetDirectory(0);

    hRun1->SetDirectory(0);
    hRun2->SetDirectory(0);


    fMap->Close();
    fDSSD->Close();


    // ============================================================
    // Histogram styles
    // ============================================================

    // ------------------------------------------------------------
    // Panel (a)
    // ------------------------------------------------------------

    hTotal->SetLineColor(kBlack);
    hTotal->SetLineWidth(histLineWidth);

    hDecay->SetLineColor(kBlack);
    hDecay->SetLineWidth(histLineWidth);
    hDecay->SetLineStyle(decayLineStyle);


    // ------------------------------------------------------------
    // Panel (b)
    // ------------------------------------------------------------

    hRun1->SetLineColor(kGreen + 2);
    hRun1->SetLineWidth(histLineWidth);

    hRun2->SetLineColor(kRed + 1);
    hRun2->SetLineWidth(histLineWidth);


    // ============================================================
    // Y range : panel (a)
    // ============================================================

    double yminA = 0.2;

    // Fixed upper limit for panel (a).
    // 10e6 = 1.0e7
    double ymaxA = 1e6;


    hDecay->SetMinimum(yminA);
    hDecay->SetMaximum(ymaxA);


    // ============================================================
    // Y range : panel (b)
    //
    // Fixed:
    //
    //     0.1 - 10^3.15
    //
    // ============================================================

    double yminB = 0.1;

    double ymaxB =
        TMath::Power(
            10.0,
            3.15
        );


    hRun2->SetMinimum(yminB);
    hRun2->SetMaximum(ymaxB);


    // ============================================================
    // Canvas
    //
    // Each panel remains approximately equivalent to
    // original 1100 x 600 canvas.
    // ============================================================

    TCanvas *c1 = new TCanvas(
        "c1",
        "Combined spectra",
        1100,
        1200
    );


    // ============================================================
    // Pads
    //
    // Boundary selected from USER CONFIGURATION so that the
    // actual frame heights of panels (a) and (b) are equal.
    // ============================================================

    TPad *p1 = new TPad(
        "p1",
        "",
        0.0,
        panelBoundary,
        1.0,
        1.00
    );

    TPad *p2 = new TPad(
        "p2",
        "",
        0.0,
        0.00,
        1.0,
        panelBoundary
    );


    // ============================================================
    // Upper pad : panel (a)
    // ============================================================

    p1->SetLeftMargin(0.12);
    p1->SetRightMargin(0.04);

    p1->SetTopMargin(0.05);

    // No gap between two panels
    p1->SetBottomMargin(0.00);

    p1->SetLogy();

    p1->SetFrameLineWidth(frameLineWidth);

    p1->SetTickx(0);
    p1->SetTicky(0);


    // ============================================================
    // Lower pad : panel (b)
    // ============================================================

    p2->SetLeftMargin(0.12);
    p2->SetRightMargin(0.04);

    // No gap
    p2->SetTopMargin(0.00);

    p2->SetBottomMargin(0.18);

    p2->SetLogy();

    p2->SetFrameLineWidth(frameLineWidth);

    p2->SetTickx(0);
    p2->SetTicky(0);


    // ============================================================
    // Draw pads
    // ============================================================

    p1->Draw();
    p2->Draw();


    // ============================================================
    // Common histogram axis style
    // ============================================================

    auto SetAxisStyle =
        [&](TH1D *h)
        {
            // ----------------------------------------------------
            // Fonts
            // ----------------------------------------------------

            h->GetXaxis()->SetTitleFont(fontCode);
            h->GetYaxis()->SetTitleFont(fontCode);

            h->GetXaxis()->SetLabelFont(fontCode);
            h->GetYaxis()->SetLabelFont(fontCode);


            // ----------------------------------------------------
            // Font sizes
            // ----------------------------------------------------

            h->GetXaxis()->SetTitleSize(commonTextSize);
            h->GetYaxis()->SetTitleSize(commonTextSize);

            h->GetXaxis()->SetLabelSize(commonTextSize);
            h->GetYaxis()->SetLabelSize(commonTextSize);


            // ----------------------------------------------------
            // Title offset
            // ----------------------------------------------------

            h->GetXaxis()->SetTitleOffset(xTitleOffset);
            h->GetYaxis()->SetTitleOffset(yTitleOffset);

            h->GetXaxis()->SetLabelOffset(xLabelOffset);
            h->GetYaxis()->SetLabelOffset(yLabelOffset);


            // ----------------------------------------------------
            // Center titles
            // ----------------------------------------------------

            h->GetXaxis()->CenterTitle(true);
            h->GetYaxis()->CenterTitle(true);


            // ----------------------------------------------------
            // Divisions
            // ----------------------------------------------------

            // 5 major divisions over 0--250 MeV -> 50 MeV major spacing;
            // 5 subdivisions -> one tick every 10 MeV.
            h->GetXaxis()->SetNdivisions(505);
            h->GetYaxis()->SetNdivisions(510);


            // ----------------------------------------------------
            // ROOT ticks hidden
            //
            // Ticks will be drawn manually.
            // ----------------------------------------------------

            h->GetXaxis()->SetTickLength(0.0);
            h->GetYaxis()->SetTickLength(0.0);


            // ----------------------------------------------------
            // Axis color
            // ----------------------------------------------------

            h->GetXaxis()->SetAxisColor(kBlack);
            h->GetYaxis()->SetAxisColor(kBlack);
        };


    SetAxisStyle(hDecay);
    SetAxisStyle(hRun2);

    // ============================================================
    // Common X-axis definition
    // ============================================================

    double xMin = 0.0;
    double xMax = 250.0;


    // ============================================================
    // Position of panel labels
    //
    // IMPORTANT:
    //
    // These are FRACTIONS OF THE ACTUAL PLOT FRAME.
    //
    // Changing these two numbers moves BOTH (a) and (b).
    //
    // X = 0 -> left side of frame
    // X = 1 -> right side of frame
    //
    // Y = 0 -> bottom side of frame
    // Y = 1 -> top side of frame
    // ============================================================

    double labelFracX = 0.97;
    double labelFracY = 0.94;


    // ============================================================
    //
    //                       PANEL (a)
    //
    // ============================================================

    p1->cd();


    // ============================================================
    // No individual Y title
    //
    // A common Canvas-wide Y title is drawn later.
    // ============================================================

    hDecay->GetYaxis()->SetTitle("");


    // ============================================================
    // Upper panel does not show X labels/title
    // ============================================================

    hDecay->GetXaxis()->SetTitle("");

    hDecay->GetXaxis()->SetLabelSize(0.0);


    // ============================================================
    // Draw spectra
    // ============================================================

    // Only the decay spectrum is shown in panel (a).
    hDecay->Draw("hist");


    p1->Modified();
    p1->Update();


    // ============================================================
    // Actual frame boundaries : panel (a)
    // ============================================================

    double xLeftA =
        p1->GetLeftMargin();

    double xRightA =
        1.0 - p1->GetRightMargin();

    double yBottomA =
        p1->GetBottomMargin();

    double yTopA =
        1.0 - p1->GetTopMargin();


    // ============================================================
    // Cover artificial baseline from empty bins
    // ============================================================

    TLine *xCoverA = new TLine(
        xMin,
        yminA,
        xMax,
        yminA
    );

    xCoverA->SetLineColor(kBlack);
    xCoverA->SetLineWidth(coverLineWidth);

    xCoverA->Draw("same");


    // ============================================================
    // Cover left Y frame
    // ============================================================

    TLine *yCoverA = new TLine(
        xLeftA,
        yBottomA,
        xLeftA,
        yTopA
    );

    yCoverA->SetNDC();

    yCoverA->SetLineColor(kBlack);
    yCoverA->SetLineWidth(coverLineWidth);

    yCoverA->Draw();


    // ============================================================
    // Manual upper X-axis ticks : panel (a)
    //
    // Major = 50 MeV
    // Tick spacing = 10 MeV
    // Ticks point inward from the upper frame.
    // ============================================================

    for (
        double x = xMin;
        x <= xMax + 1e-9;
        x += xTickStep
    )
    {
        double frac =
            (x - xMin)
            /
            (xMax - xMin);


        double xNDC =
            xLeftA
            +
            frac
            *
            (xRightA - xLeftA);


        bool isMajor =
            (
                TMath::Abs(
                    x / xTickMajorStep
                    -
                    TMath::Nint(
                        x / xTickMajorStep
                    )
                )
                < 1e-6
            );


        double tickLength =
            isMajor
            ?
            xTickMajor
            :
            xTickMinor;


        TLine *tickTop =
            new TLine(
                xNDC,
                yTopA,
                xNDC,
                yTopA - tickLength
            );

        tickTop->SetNDC();
        tickTop->SetLineColor(kBlack);
        tickTop->SetLineWidth(tickLineWidth);
        tickTop->Draw();
    }


    // ============================================================
    // Manual logarithmic Y ticks : panel (a)
    // ============================================================

    double logMinA =
        TMath::Log10(yminA);

    double logMaxA =
        TMath::Log10(ymaxA);


    int decadeMinA =
        TMath::FloorNint(logMinA) - 1;

    int decadeMaxA =
        TMath::CeilNint(logMaxA) + 1;


    for (
        int decade = decadeMinA;
        decade <= decadeMaxA;
        ++decade
    )
    {
        double base =
            TMath::Power(
                10.0,
                decade
            );


        for (
            int m = 1;
            m <= 9;
            ++m
        )
        {
            double value =
                m * base;


            if (value < yminA)
                continue;

            if (value > ymaxA)
                continue;


            double frac =
                (
                    TMath::Log10(value)
                    -
                    logMinA
                )
                /
                (
                    logMaxA
                    -
                    logMinA
                );


            double yNDC =
                yBottomA
                +
                frac
                *
                (
                    yTopA
                    -
                    yBottomA
                );


            double tickLength =
                (m == 1)
                ?
                yTickMajor
                :
                yTickMinor;


            // ----------------------------------------------------
            // Left Y tick
            // ----------------------------------------------------

            TLine *tickLeft =
                new TLine(
                    xLeftA - tickLength,
                    yNDC,
                    xLeftA,
                    yNDC
                );

            tickLeft->SetNDC();

            tickLeft->SetLineColor(kBlack);
            tickLeft->SetLineWidth(tickLineWidth);

            tickLeft->Draw();


            // ----------------------------------------------------
            // Right Y tick
            // ----------------------------------------------------

            TLine *tickRight =
                new TLine(
                    xRightA,
                    yNDC,
                    xRightA - tickLength,
                    yNDC
                );

            tickRight->SetNDC();

            tickRight->SetLineColor(kBlack);
            tickRight->SetLineWidth(tickLineWidth);

            tickRight->Draw();
        }
    }


    // ============================================================
    // Spectrum label : panel (a)
    //
    // Put the label in the lower-middle part of the plotting frame.
    // The decay counts above the baseline are concentrated at the
    // left, so this position does not cover the data.
    // ============================================================
    double spectrumLabelFracX = 0.51;
    double spectrumLabelFracY = 0.34;

    TLatex spectrumLabelA;
    spectrumLabelA.SetNDC();
    spectrumLabelA.SetTextFont(fontCode);
    spectrumLabelA.SetTextSize(commonTextSize);
    spectrumLabelA.SetTextAlign(22);

    spectrumLabelA.DrawLatex(
        xLeftA
        + spectrumLabelFracX * (xRightA - xLeftA),
        yBottomA
        + spectrumLabelFracY * (yTopA - yBottomA),
        "Total decay spectrum"
    );


    // ============================================================
    // Label (a)
    //
    // Position calculated from ACTUAL frame size.
    // ============================================================

    TLatex latexA;

    latexA.SetNDC();
    latexA.SetTextFont(fontCode);
    latexA.SetTextSize(panelLabelSize);
    latexA.SetTextAlign(33);


    double labelXA =
        xLeftA
        +
        labelFracX
        *
        (xRightA - xLeftA);


    double labelYA =
        yBottomA
        +
        labelFracY
        *
        (yTopA - yBottomA);


    latexA.DrawLatex(
        labelXA,
        labelYA,
        "(a)"
    );


    // ============================================================
    //
    //                       PANEL (b)
    //
    // ============================================================

    p2->cd();


    // ============================================================
    // X title
    // ============================================================

    hRun2->GetXaxis()->SetTitle(
        "Energy (MeV)"
    );


    // ============================================================
    // No individual Y title
    // ============================================================

    hRun2->GetYaxis()->SetTitle("");


    // ============================================================
    // Draw spectra
    // ============================================================

    hRun2->Draw("hist");

    hRun1->Draw("hist same");


    p2->Modified();
    p2->Update();


    // ============================================================
    // Actual frame boundaries : panel (b)
    // ============================================================

    double xLeftB =
        p2->GetLeftMargin();

    double xRightB =
        1.0 - p2->GetRightMargin();

    double yBottomB =
        p2->GetBottomMargin();

    double yTopB =
        1.0 - p2->GetTopMargin();


    // ============================================================
    // Cover artificial baseline
    // ============================================================

    TLine *xCoverB = new TLine(
        xMin,
        yminB,
        xMax,
        yminB
    );

    xCoverB->SetLineColor(kBlack);
    xCoverB->SetLineWidth(coverLineWidth);

    xCoverB->Draw("same");


    // ============================================================
    // Cover left Y frame
    // ============================================================

    TLine *yCoverB = new TLine(
        xLeftB,
        yBottomB,
        xLeftB,
        yTopB
    );

    yCoverB->SetNDC();

    yCoverB->SetLineColor(kBlack);
    yCoverB->SetLineWidth(coverLineWidth);

    yCoverB->Draw();


    // ============================================================
    // Manual X ticks : panel (b)
    //
    // Major = 50 MeV
    // Tick spacing = 10 MeV
    // ============================================================

    for (
        double x = xMin;
        x <= xMax + 1e-9;
        x += xTickStep
    )
    {
        double frac =
            (x - xMin)
            /
            (xMax - xMin);


        double xNDC =
            xLeftB
            +
            frac
            *
            (xRightB - xLeftB);


        bool isMajor =
            (
                TMath::Abs(
                    x / xTickMajorStep
                    -
                    TMath::Nint(
                        x / xTickMajorStep
                    )
                )
                < 1e-6
            );


        double tickLength =
            isMajor
            ?
            xTickMajor
            :
            xTickMinor;


        // --------------------------------------------------------
        // Bottom X tick
        // --------------------------------------------------------

        TLine *tickBottom =
            new TLine(
                xNDC,
                yBottomB,
                xNDC,
                yBottomB - tickLength
            );

        tickBottom->SetNDC();

        tickBottom->SetLineColor(kBlack);
        tickBottom->SetLineWidth(tickLineWidth);

        tickBottom->Draw();


        // --------------------------------------------------------
        // Top X tick
        // --------------------------------------------------------

        TLine *tickTop =
            new TLine(
                xNDC,
                yTopB,
                xNDC,
                yTopB - tickLength
            );

        tickTop->SetNDC();

        tickTop->SetLineColor(kBlack);
        tickTop->SetLineWidth(tickLineWidth);

        tickTop->Draw();
    }


    // ============================================================
    // Manual logarithmic Y ticks : panel (b)
    // ============================================================

    double logMinB =
        TMath::Log10(yminB);

    double logMaxB =
        TMath::Log10(ymaxB);


    int decadeMinB =
        TMath::FloorNint(logMinB) - 1;

    int decadeMaxB =
        TMath::CeilNint(logMaxB) + 1;


    for (
        int decade = decadeMinB;
        decade <= decadeMaxB;
        ++decade
    )
    {
        double base =
            TMath::Power(
                10.0,
                decade
            );


        for (
            int m = 1;
            m <= 9;
            ++m
        )
        {
            double value =
                m * base;


            if (value < yminB)
                continue;

            if (value > ymaxB)
                continue;


            double frac =
                (
                    TMath::Log10(value)
                    -
                    logMinB
                )
                /
                (
                    logMaxB
                    -
                    logMinB
                );


            double yNDC =
                yBottomB
                +
                frac
                *
                (
                    yTopB
                    -
                    yBottomB
                );


            double tickLength =
                (m == 1)
                ?
                yTickMajor
                :
                yTickMinor;


            // ----------------------------------------------------
            // Left Y tick
            // ----------------------------------------------------

            TLine *tickLeft =
                new TLine(
                    xLeftB - tickLength,
                    yNDC,
                    xLeftB,
                    yNDC
                );

            tickLeft->SetNDC();

            tickLeft->SetLineColor(kBlack);
            tickLeft->SetLineWidth(tickLineWidth);

            tickLeft->Draw();


            // ----------------------------------------------------
            // Right Y tick
            // ----------------------------------------------------

            TLine *tickRight =
                new TLine(
                    xRightB,
                    yNDC,
                    xRightB - tickLength,
                    yNDC
                );

            tickRight->SetNDC();

            tickRight->SetLineColor(kBlack);
            tickRight->SetLineWidth(tickLineWidth);

            tickRight->Draw();
        }
    }


    // ============================================================
    // Legend : panel (b)
    // ============================================================

    TLegend *legB = new TLegend(
        0.15,
        0.80,
        0.37,
        0.94
    );

    legB->SetBorderSize(0);
    legB->SetFillStyle(0);

    legB->SetTextFont(fontCode);
    legB->SetTextSize(commonTextSize);


    legB->AddEntry(
        hRun1,
        "Run 1",
        "l"
    );

    legB->AddEntry(
        hRun2,
        "Run 2",
        "l"
    );


    legB->Draw();
    // ============================================================
    // Decay-time condition : panel (b)
    // ============================================================

    TLatex latexCut;

    latexCut.SetNDC();
    latexCut.SetTextFont(fontCode);
    latexCut.SetTextSize(commonTextSize);

    // Left/top anchoring for the cut text.
    latexCut.SetTextAlign(13);

    latexCut.DrawLatex(
        0.15,
        0.76,
        "#Delta#it{t}(ER-#alpha/SF)<10 ms"
    );

    // ============================================================
    // Spectrum label : panel (b)
    //
    // Position is specified directly in the panel's data coordinates.
    // ============================================================
    TLatex spectrumLabelB;
    spectrumLabelB.SetTextFont(fontCode);
    spectrumLabelB.SetTextSize(commonTextSize);
    spectrumLabelB.SetTextAlign(22);

    spectrumLabelB.DrawLatex(
        160.0,
        3.0,
        "SF of ^{242}Fm"
    );

    // ============================================================
    // Label (b)
    //
    // EXACT SAME relative position inside actual frame as (a)
    // ============================================================

    TLatex latexB;

    latexB.SetNDC();
    latexB.SetTextFont(fontCode);
    latexB.SetTextSize(panelLabelSize);
    latexB.SetTextAlign(33);


    double labelXB =
        xLeftB
        +
        labelFracX
        *
        (xRightB - xLeftB);


    double labelYB =
        yBottomB
        +
        labelFracY
        *
        (yTopB - yBottomB);


    latexB.DrawLatex(
        labelXB,
        labelYB,
        "(b)"
    );


    // ============================================================
    // Update physics pads
    // ============================================================

    p1->Modified();
    p1->Update();

    p2->Modified();
    p2->Update();


    // ============================================================
    // Transparent whole-canvas overlay
    //
    // Used for one common Y-axis title.
    // ============================================================

    c1->cd();


    TPad *pOverlay = new TPad(
        "pOverlay",
        "",
        0.0,
        0.0,
        1.0,
        1.0
    );


    pOverlay->SetFillStyle(4000);

    pOverlay->SetFrameFillStyle(0);
    pOverlay->SetFrameLineColor(0);

    pOverlay->SetBorderMode(0);
    pOverlay->SetBorderSize(0);

    pOverlay->SetLeftMargin(0.0);
    pOverlay->SetRightMargin(0.0);
    pOverlay->SetTopMargin(0.0);
    pOverlay->SetBottomMargin(0.0);


    pOverlay->Draw();

    pOverlay->cd();


    // ============================================================
    // Common Y-axis title
    //
    // Vertically centered on the WHOLE canvas.
    // ============================================================

    TLatex commonYTitle;

    commonYTitle.SetNDC();

    commonYTitle.SetTextFont(fontCode);

    commonYTitle.SetTextSize(commonTextSize);

    commonYTitle.SetTextAlign(22);

    commonYTitle.SetTextAngle(90);


    commonYTitle.DrawLatex(
        0.030,
        0.50,
        "Counts / 20 keV"
    );


    // ============================================================
    // Final update
    // ============================================================

    c1->cd();

    c1->Modified();
    c1->Update();


    // ============================================================
    // Save ROOT file
    // ============================================================

    TFile *fout = new TFile(
        "Combined_Spectrum.root",
        "RECREATE"
    );


    hTotal->Write();
    hDecay->Write();

    hRun1->Write();
    hRun2->Write();

    c1->Write();


    fout->Close();


    // ============================================================
    // Save PDF
    // ============================================================

    c1->SaveAs(
        "Combined_Spectrum.eps"
    );
        c1->SaveAs(
        "Combined_Spectrum.pdf"
    );

    // ============================================================
    // Output
    // ============================================================

    cout << endl;

    cout
        << "======================================"
        << endl;

    cout
        << "Panel (b) Y range = "
        << yminB
        << " - "
        << ymaxB
        << endl;

    cout
        << "Panel labels:"
        << endl;

    cout
        << "  labelFracX = "
        << labelFracX
        << endl;

    cout
        << "  labelFracY = "
        << labelFracY
        << endl;

    cout
        << "Saved:"
        << endl;

    cout
        << "  Combined_Spectrum.root"
        << endl;

    cout
        << "  Combined_Spectrum.pdf"
        << endl;

    cout
        << "======================================"
        << endl;
}
