#pragma once
#include "Config.h"
#include <string>

class TFile;

// Per-run processing: open input ROOT, fill histos under gates, fit peaks, compute A80/X80/Y80,
// write per-run QA objects into fqa, and produce per-run multi-page PDF.
class RunProcessor {
public:
  using Config = A80Config;

  explicit RunProcessor(Config cfg = Config());

  // Returns true if the run was processed (input file existed and contained tr_map).
  // Returns false if skipped (missing file, open failure, missing tree).
  bool ProcessRun(const std::string& indir, int run, TFile& fqa, const std::string& pdfDir, SummaryRow& outRow) const;

private:
  Config cfg_;
};
