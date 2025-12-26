#include "QaIO.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TString.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>

#include <memory>

void EnsureDir(const std::string& dir) {
  gSystem->mkdir(dir.c_str(), kTRUE);
}

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
) {
  const TString dname = TString::Format("run%05d", run);
  TDirectory* d = fqa.mkdir(dname);
  d->cd();

  hE.Write("hE_5000_6000");
  hXY_all.Write("hXY_all");
  hXY1.Write("hXY_p1");
  hXY2.Write("hXY_p2");
  hXY3.Write("hXY_p3");

  hX_all.Write("hX_all");
  hY_all.Write("hY_all");
  hX1.Write("hX_p1"); hY1.Write("hY_p1");
  hX2.Write("hX_p2"); hY2.Write("hY_p2");
  hX3.Write("hX_p3"); hY3.Write("hY_p3");

  if (cEnergy) cEnergy->Write("c_energy");
  if (cAll)    cAll->Write("c_all");
  if (cP1)     cP1->Write("c_p1");
  if (cP2)     cP2->Write("c_p2");
  if (cP3)     cP3->Write("c_p3");

  fqa.cd();
}

void DrawXYPage(
  TCanvas& c, TLatex& t,
  TH2D& hxy, const char* title,
  const PeakWin* pinfo,
  const AeqResult& Aall, const AeqResult& Xall, const AeqResult& Yall
) {
  c.Clear();
  c.Divide(2, 2);

  c.cd(1);
  gPad->SetRightMargin(0.14);
  hxy.SetTitle(title);
  hxy.Draw("COLZ");

  // RAII projections to avoid leaking histograms (important for long run batches)
  c.cd(2);
  std::unique_ptr<TH1D> hx(hxy.ProjectionX(Form("%s_px_tmp", hxy.GetName())));
  hx->SetDirectory(nullptr);
  hx->SetTitle("X projection;DSSDX_Ch;Counts");
  hx->Draw();

  c.cd(4);
  std::unique_ptr<TH1D> hy(hxy.ProjectionY(Form("%s_py_tmp", hxy.GetName())));
  hy->SetDirectory(nullptr);
  hy->SetTitle("Y projection;DSSDY_Ch;Counts");
  hy->Draw();

  c.cd(3);
  gPad->Clear();
  t.SetNDC(true);

  t.SetTextSize(0.045);
  t.DrawLatex(0.08, 0.85, title);

  t.SetTextSize(0.037);
  if (pinfo) {
    if (pinfo->ok) {
      t.DrawLatex(0.08, 0.72, TString::Format("N = %lld", pinfo->N));
      t.DrawLatex(0.08, 0.66, TString::Format("mux = %.3f, muy = %.3f", pinfo->mux, pinfo->muy));
      t.DrawLatex(0.08, 0.60, TString::Format("A80_eq = %.3f  (A80_pix=%lld, cov@pix=%.4f)",
                                               pinfo->A80_2D.n_eq, pinfo->A80_2D.n_pix, pinfo->A80_2D.cov_at_npix));
      t.DrawLatex(0.08, 0.54, TString::Format("X80_eq = %.3f (pix=%lld,cov=%.4f)",
                                               pinfo->X80_1D.n_eq, pinfo->X80_1D.n_pix, pinfo->X80_1D.cov_at_npix));
      t.DrawLatex(0.08, 0.48, TString::Format("Y80_eq = %.3f (pix=%lld,cov=%.4f)",
                                               pinfo->Y80_1D.n_eq, pinfo->Y80_1D.n_pix, pinfo->Y80_1D.cov_at_npix));
      t.DrawLatex(0.08, 0.42, TString::Format("ROI = [%.0f, %.0f] keV", pinfo->roi_lo, pinfo->roi_hi));
    } else {
      t.DrawLatex(0.08, 0.72, "Peak NOT FOUND");
    }
  } else {
    t.DrawLatex(0.08, 0.72, TString::Format("All peaks union: A80_eq=%.3f (pix=%lld,cov=%.4f)",
                                             Aall.n_eq, Aall.n_pix, Aall.cov_at_npix));
    t.DrawLatex(0.08, 0.66, TString::Format("X80_eq=%.3f (pix=%lld,cov=%.4f)",
                                             Xall.n_eq, Xall.n_pix, Xall.cov_at_npix));
    t.DrawLatex(0.08, 0.60, TString::Format("Y80_eq=%.3f (pix=%lld,cov=%.4f)",
                                             Yall.n_eq, Yall.n_pix, Yall.cov_at_npix));
  }
}
