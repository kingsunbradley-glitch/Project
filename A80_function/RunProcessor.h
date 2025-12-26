#pragma once
#include "A80Types.h"
#include <string>

class TFile;

// Per-run processing: open input ROOT, fill histos under gates, fit peaks, compute A80/X80/Y80,
// write per-run QA objects into fqa, and produce per-run multi-page PDF.
class RunProcessor {
public:
  struct Config {
    double EgLo = 5000.0;
    double EgHi = 6000.0;
    double targetFrac = 0.80;
    double fitWindow = 80.0;

    int nx = 128;
    int ny = 48;
    double xlo = -0.5, xhi = 127.5;
    double ylo = -0.5, yhi = 47.5;

    PeakWin Ptpl[3];

    Config();
  };

  explicit RunProcessor(Config cfg = Config());

  // Returns true if the run was processed (input file existed and contained tr_map).
  // Returns false if skipped (missing file, open failure, missing tree).
  bool ProcessRun(const std::string& indir, int run, TFile& fqa, const std::string& pdfDir, SummaryRow& outRow) const;

private:
  Config cfg_;
};
