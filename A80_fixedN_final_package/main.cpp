// Modular version of process_runs_A80
// Build (one-liner):
//   g++ -O2 -std=c++17 src/*.cpp -Iinclude $(root-config --cflags --libs) -o process_runs_A80
//
// Run:
//   ./process_runs_A80 <indir> <runFirst> <runLast> <outSummary.root> <outQA.root> <pdfDir>

#include "RunProcessor.h"
#include "QaIO.h"
#include "A80Types.h"

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>

#include <cstdio>
#include <cstdlib>
#include <string>

static void BindSummaryBranches(TTree& t, SummaryRow& r) {
  t.Branch("runnum", &r.run, "runnum/I");
  t.Branch("nPeaks", &r.nPeaks, "nPeaks/I");

  t.Branch("N_all", &r.N_all, "N_all/L");

  t.Branch("A80_eq_all", &r.A80_eq_all, "A80_eq_all/D");
  t.Branch("A80_pix_all", &r.A80_pix_all, "A80_pix_all/L");
  t.Branch("A80_cov_all", &r.A80_cov_all, "A80_cov_all/D");

  t.Branch("X80_eq_all", &r.X80_eq_all, "X80_eq_all/D");
  t.Branch("X80_pix_all", &r.X80_pix_all, "X80_pix_all/L");
  t.Branch("X80_cov_all", &r.X80_cov_all, "X80_cov_all/D");

  t.Branch("Y80_eq_all", &r.Y80_eq_all, "Y80_eq_all/D");
  t.Branch("Y80_pix_all", &r.Y80_pix_all, "Y80_pix_all/L");
  t.Branch("Y80_cov_all", &r.Y80_cov_all, "Y80_cov_all/D");

  t.Branch("muE", r.muE, "muE[3]/D");
  t.Branch("sigE", r.sigE, "sigE[3]/D");

  t.Branch("Npk", r.Npk, "Npk[3]/L");
  t.Branch("mux", r.mux, "mux[3]/D");
  t.Branch("muy", r.muy, "muy[3]/D");

  t.Branch("A80eq", r.A80eq, "A80eq[3]/D");
  t.Branch("A80pix", r.A80pix, "A80pix[3]/L");
  t.Branch("A80cov", r.A80cov, "A80cov[3]/D");

  t.Branch("X80eq", r.X80eq, "X80eq[3]/D");
  t.Branch("X80pix", r.X80pix, "X80pix[3]/L");
  t.Branch("X80cov", r.X80cov, "X80cov[3]/D");

  t.Branch("Y80eq", r.Y80eq, "Y80eq[3]/D");
  t.Branch("Y80pix", r.Y80pix, "Y80pix[3]/L");
  t.Branch("Y80cov", r.Y80cov, "Y80cov[3]/D");

  t.Branch("Nfix", r.Nfix, "Nfix[3]/L");
  t.Branch("A80eq_fixN", r.A80eq_fixN, "A80eq_fixN[3]/D");
  t.Branch("A80pix_fixN", r.A80pix_fixN, "A80pix_fixN[3]/L");
  t.Branch("A80cov_fixN", r.A80cov_fixN, "A80cov_fixN[3]/D");
}

int main(int argc, char** argv) {
  if (argc < 7) {
    std::fprintf(stderr,
      "Usage:\n  %s <indir> <runFirst> <runLast> <outSummary.root> <outQA.root> <pdfDir>\n",
      argv[0]);
    return 2;
  }

  const std::string indir  = argv[1];
  const int runFirst       = std::atoi(argv[2]);
  const int runLast        = std::atoi(argv[3]);
  const std::string outSum = argv[4];
  const std::string outQA  = argv[5];
  const std::string pdfDir = argv[6];

  gStyle->SetOptStat(0);
  EnsureDir(pdfDir);

  TFile fqa(outQA.c_str(), "RECREATE");
  TFile fsum(outSum.c_str(), "RECREATE");
  TTree tsum("summary", "summary");

  SummaryRow row;
  row.reset(0);
  BindSummaryBranches(tsum, row);

  RunProcessor proc;

  for (int run = runFirst; run <= runLast; ++run) {
    if (!proc.ProcessRun(indir, run, fqa, pdfDir, row)) {
      continue;
    }
    fsum.cd();
    tsum.Fill();
  }

  fqa.cd(); fqa.Write(); fqa.Close();
  fsum.cd(); tsum.Write(); fsum.Close();

  std::printf("All done.\nOutput:\n  %s\n  %s\n  %s/\n",
              outSum.c_str(), outQA.c_str(), pdfDir.c_str());
  return 0;
}
