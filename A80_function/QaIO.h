#pragma once
#include "A80Types.h"
#include <string>

class TFile;
class TCanvas;
class TH1D;
class TH2D;
class TLatex;

// Create directory recursively.
void EnsureDir(const std::string& dir);

// Write per-run objects into QA ROOT file (under run%05d directory).
void WriteRunObjects(
  TFile& fqa, int run,
  TH1D& hE,
  TH2D& hXY_all, TH2D& hXY1, TH2D& hXY2, TH2D& hXY3,
  TH1D& hX_all, TH1D& hY_all,
  TH1D& hX1, TH1D& hY1,
  TH1D& hX2, TH1D& hY2,
  TH1D& hX3, TH1D& hY3,
  TCanvas* cEnergy, TCanvas* cAll,
  TCanvas* cP1, TCanvas* cP2, TCanvas* cP3
);

// Draw a 4-panel page: XY heatmap + X/Y projections + text block (A80/X80/Y80 + ROI).
// NOTE: This function will allocate projection histograms and delete them before returning.
void DrawXYPage(
  TCanvas& c, TLatex& t,
  TH2D& hxy, const char* title,
  const PeakWin* pinfo,
  const AeqResult& Aall, const AeqResult& Xall, const AeqResult& Yall
);
