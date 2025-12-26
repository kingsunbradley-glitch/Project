#include "RunProcessor.h"

#include "PeakFinder.h"
#include "Metrics.h"
#include "QaIO.h"

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>

#include <cstdio>
#include <cstdlib>
#include <memory>

RunProcessor::Config::Config() {
  Ptpl[0].win_lo = 5100; Ptpl[0].win_hi = 5200;
  Ptpl[1].win_lo = 5435; Ptpl[1].win_hi = 5535;
  Ptpl[2].win_lo = 5754; Ptpl[2].win_hi = 5854;
}

RunProcessor::RunProcessor(Config cfg) : cfg_(std::move(cfg)) {}

static void ResetPeak(PeakWin& p, const PeakWin& tpl) {
  p = tpl;
  p.ok = false;
  p.mu = p.sig = 0.0;
  p.roi_lo = p.roi_hi = 0.0;
  p.N = 0;
  p.mux = p.muy = 0.0;
  p.A80_2D = {};
  p.X80_1D = {};
  p.Y80_1D = {};
}

bool RunProcessor::ProcessRun(const std::string& indir, int run, TFile& fqa, const std::string& pdfDir, SummaryRow& outRow) const {
  outRow.reset(run);

  PeakWin P[3];
  for (int i = 0; i < 3; ++i) ResetPeak(P[i], cfg_.Ptpl[i]);

  const TString fn = TString::Format("%s/run%05d_map.root", indir.c_str(), run);
  if (gSystem->AccessPathName(fn)) {
    return false; // missing file
  }

  TFile fin(fn, "READ");
  if (fin.IsZombie()) {
    std::printf("run%05d: cannot open, skip\n", run);
    return false;
  }

  TTree* tr = static_cast<TTree*>(fin.Get("tr_map"));
  if (!tr) {
    std::printf("run%05d: no tr_map, skip\n", run);
    return false;
  }

  // branches (only need first hit when multiplicity==1)
  UShort_t mx = 0, my = 0;
  Double_t Xch[256]{0}, Ych[256]{0};
  Double_t XE[256]{0};
  tr->SetBranchAddress("DSSDX_mul", &mx);
  tr->SetBranchAddress("DSSDY_mul", &my);
  tr->SetBranchAddress("DSSDX_Ch", Xch);
  tr->SetBranchAddress("DSSDY_Ch", Ych);
  tr->SetBranchAddress("DSSDX_E", XE);

  // Energy spectrum under global gate
  TH1D hE("hE", TString::Format("run%05d: DSSDX_E (keV);E_{x} (keV);Counts", run),
          1000, cfg_.EgLo, cfg_.EgHi);

  const Long64_t nent = tr->GetEntries();
  for (Long64_t i = 0; i < nent; ++i) {
    tr->GetEntry(i);
    if (mx != 1) continue;
    if (my != 1) continue;
    const double E = XE[0];
    if (E < cfg_.EgLo || E > cfg_.EgHi) continue;
    hE.Fill(E);
  }

  // Fit peaks independently
  for (int i = 0; i < 3; ++i) {
    if (FitOnePeakGaus(&hE, P[i], cfg_.fitWindow)) {
      outRow.nPeaks++;
      outRow.muE[i] = P[i].mu;
      outRow.sigE[i] = P[i].sig;
    }
  }

  // XY histograms (ROI union + per-peak)
  TH2D hXY_all("hXY_all", "All peaks (ROI union);DSSDX_Ch;DSSDY_Ch",
               cfg_.nx, cfg_.xlo, cfg_.xhi, cfg_.ny, cfg_.ylo, cfg_.yhi);
  TH2D hXY1("hXY1", "Peak1 XY;DSSDX_Ch;DSSDY_Ch", cfg_.nx, cfg_.xlo, cfg_.xhi, cfg_.ny, cfg_.ylo, cfg_.yhi);
  TH2D hXY2("hXY2", "Peak2 XY;DSSDX_Ch;DSSDY_Ch", cfg_.nx, cfg_.xlo, cfg_.xhi, cfg_.ny, cfg_.ylo, cfg_.yhi);
  TH2D hXY3("hXY3", "Peak3 XY;DSSDX_Ch;DSSDY_Ch", cfg_.nx, cfg_.xlo, cfg_.xhi, cfg_.ny, cfg_.ylo, cfg_.yhi);

  for (Long64_t i = 0; i < nent; ++i) {
    tr->GetEntry(i);
    if (mx != 1) continue;
    if (my != 1) continue;
    const double E = XE[0];
    if (E < cfg_.EgLo || E > cfg_.EgHi) continue;

    int pid = 0;
    for (int k = 0; k < 3; ++k) {
      if (P[k].ok && E >= P[k].roi_lo && E <= P[k].roi_hi) { pid = k + 1; break; }
    }
    if (pid == 0) continue;

    const double x = Xch[0];
    const double y = Ych[0];

    hXY_all.Fill(x, y);
    if (pid == 1) hXY1.Fill(x, y);
    if (pid == 2) hXY2.Fill(x, y);
    if (pid == 3) hXY3.Fill(x, y);
  }

  // Projections (IMPORTANT: keep ownership to avoid leaks)
  auto hX_all = std::unique_ptr<TH1D>(hXY_all.ProjectionX("hX_all"));
  auto hY_all = std::unique_ptr<TH1D>(hXY_all.ProjectionY("hY_all"));
  auto hX1    = std::unique_ptr<TH1D>(hXY1.ProjectionX("hX1"));
  auto hY1    = std::unique_ptr<TH1D>(hXY1.ProjectionY("hY1"));
  auto hX2    = std::unique_ptr<TH1D>(hXY2.ProjectionX("hX2"));
  auto hY2    = std::unique_ptr<TH1D>(hXY2.ProjectionY("hY2"));
  auto hX3    = std::unique_ptr<TH1D>(hXY3.ProjectionX("hX3"));
  auto hY3    = std::unique_ptr<TH1D>(hXY3.ProjectionY("hY3"));

  // Detach from any directory to keep lifetime controlled here
  for (auto* hh : {hX_all.get(), hY_all.get(), hX1.get(), hY1.get(), hX2.get(), hY2.get(), hX3.get(), hY3.get()}) {
    if (hh) hh->SetDirectory(nullptr);
  }

  // A80 / X80 / Y80
  const AeqResult Aall = CalcAeq2D(&hXY_all, cfg_.targetFrac);
  const AeqResult Xall = CalcAeq1D(hX_all.get(), cfg_.targetFrac);
  const AeqResult Yall = CalcAeq1D(hY_all.get(), cfg_.targetFrac);

  outRow.N_all = static_cast<long long>(hXY_all.GetEntries());

  outRow.A80_eq_all  = Aall.n_eq;
  outRow.A80_pix_all = Aall.n_pix;
  outRow.A80_cov_all = Aall.cov_at_npix;

  outRow.X80_eq_all  = Xall.n_eq;
  outRow.X80_pix_all = Xall.n_pix;
  outRow.X80_cov_all = Xall.cov_at_npix;

  outRow.Y80_eq_all  = Yall.n_eq;
  outRow.Y80_pix_all = Yall.n_pix;
  outRow.Y80_cov_all = Yall.cov_at_npix;

  // per-peak
  TH2D* hp2[3] = {&hXY1, &hXY2, &hXY3};
  TH1D* hpX[3] = {hX1.get(), hX2.get(), hX3.get()};
  TH1D* hpY[3] = {hY1.get(), hY2.get(), hY3.get()};

  for (int i = 0; i < 3; ++i) {
    if (!P[i].ok) continue;

    P[i].N   = static_cast<long long>(hp2[i]->GetEntries());
    P[i].mux = hp2[i]->GetMean(1);
    P[i].muy = hp2[i]->GetMean(2);

    P[i].A80_2D = CalcAeq2D(hp2[i], cfg_.targetFrac);
    P[i].X80_1D = CalcAeq1D(hpX[i], cfg_.targetFrac);
    P[i].Y80_1D = CalcAeq1D(hpY[i], cfg_.targetFrac);

    outRow.Npk[i] = P[i].N;
    outRow.mux[i] = P[i].mux;
    outRow.muy[i] = P[i].muy;

    outRow.A80eq[i] = P[i].A80_2D.n_eq; outRow.A80pix[i] = P[i].A80_2D.n_pix; outRow.A80cov[i] = P[i].A80_2D.cov_at_npix;
    outRow.X80eq[i] = P[i].X80_1D.n_eq; outRow.X80pix[i] = P[i].X80_1D.n_pix; outRow.X80cov[i] = P[i].X80_1D.cov_at_npix;
    outRow.Y80eq[i] = P[i].Y80_1D.n_eq; outRow.Y80pix[i] = P[i].Y80_1D.n_pix; outRow.Y80cov[i] = P[i].Y80_1D.cov_at_npix;
  }

  // ----- PDF output -----
  EnsureDir(pdfDir);
  const TString pdf = TString::Format("%s/run%05d_QA_A80.pdf", pdfDir.c_str(), run);

  TCanvas cEnergy("cEnergy", "cEnergy", 1100, 850);
  TCanvas cAll("cAll", "cAll", 1100, 850);
  TCanvas cP1("cP1", "cP1", 1100, 850);
  TCanvas cP2("cP2", "cP2", 1100, 850);
  TCanvas cP3("cP3", "cP3", 1100, 850);

  TLatex tx; tx.SetNDC(true);
  cEnergy.Print((pdf + "[").Data());

  // Page 1: energy + windows/ROIs (5000-6000 only)
  cEnergy.Clear();
  hE.Draw();

  TLine l; l.SetLineStyle(2);
  const double ymax = hE.GetMaximum();
  for (int i = 0; i < 3; ++i) {
    l.DrawLine(P[i].win_lo, 0, P[i].win_lo, ymax * 0.55);
    l.DrawLine(P[i].win_hi, 0, P[i].win_hi, ymax * 0.55);
    if (P[i].ok) {
      l.SetLineStyle(1);
      l.DrawLine(P[i].roi_lo, 0, P[i].roi_lo, ymax * 0.85);
      l.DrawLine(P[i].roi_hi, 0, P[i].roi_hi, ymax * 0.85);
      l.SetLineStyle(2);
    }
  }

  tx.SetTextSize(0.035);
  tx.DrawLatex(0.12, 0.92, "Cuts: DSSDX_mul==1 && DSSDY_mul==1 && DSSDX_E in [5000,6000] keV");
  for (int i = 0; i < 3; ++i) {
    const TString s = P[i].ok
      ? TString::Format("Peak%d: mu=%.1f keV, sigma=%.1f keV, ROI=[%.0f,%.0f]",
                        i + 1, P[i].mu, P[i].sig, P[i].roi_lo, P[i].roi_hi)
      : TString::Format("Peak%d: NOT FOUND", i + 1);
    tx.DrawLatex(0.12, 0.86 - 0.04 * i, s);
  }
  cEnergy.Print(pdf.Data());

  // Page 2: all peaks union
  DrawXYPage(cAll, tx, hXY_all, "All peaks (ROI union) XY", nullptr, Aall, Xall, Yall);
  cAll.Print(pdf.Data());

  // Page 3-5: each peak
  DrawXYPage(cP1, tx, hXY1, "Peak1 XY (ROI)", &P[0], Aall, Xall, Yall); cP1.Print(pdf.Data());
  DrawXYPage(cP2, tx, hXY2, "Peak2 XY (ROI)", &P[1], Aall, Xall, Yall); cP2.Print(pdf.Data());
  DrawXYPage(cP3, tx, hXY3, "Peak3 XY (ROI)", &P[2], Aall, Xall, Yall); cP3.Print(pdf.Data());

  cEnergy.Print((pdf + "]").Data());

  // write QA objects
  fqa.cd();
  WriteRunObjects(
    fqa, run, hE,
    hXY_all, hXY1, hXY2, hXY3,
    *hX_all, *hY_all,
    *hX1, *hY1,
    *hX2, *hY2,
    *hX3, *hY3,
    &cEnergy, &cAll, &cP1, &cP2, &cP3
  );

  // progress print
  std::printf("DONE run%05d: nPeaks=%d N_all=%lld A80_eq=%.3f X80_eq=%.3f Y80_eq=%.3f\n",
              run, outRow.nPeaks, outRow.N_all, outRow.A80_eq_all, outRow.X80_eq_all, outRow.Y80_eq_all);
  std::fflush(stdout);
  return true;
}
