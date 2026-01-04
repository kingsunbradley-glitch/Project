#include "PeakFinder.h"

#include <TH1D.h>
#include <TF1.h>

#include <algorithm>
#include <cmath>

bool FitOnePeakGaus(TH1D* hE, PeakWin& p, const A80Config& cfg) {
  if (!hE) return false;

  // seed: max bin in [win_lo, win_hi]
  const int bLo = hE->FindBin(p.win_lo);
  const int bHi = hE->FindBin(p.win_hi);
  int bMax = bLo;
  double yMax = -1.0;

  for (int b = bLo; b <= bHi; ++b) {
    const double y = hE->GetBinContent(b);
    if (y > yMax) { yMax = y; bMax = b; }
  }
  if (yMax <= 0.0) return false;
  if (yMax < cfg.minPeakMaxCount) return false;

  const double seed = hE->GetBinCenter(bMax);

  // local fit range: [seed-W, seed+W] clipped to search window
  const double W = cfg.fitWindow_keV;
  const double r1 = std::max(p.win_lo, seed - W);
  const double r2 = std::min(p.win_hi, seed + W);

  TF1 f1("f1", "gaus", r1, r2);
  f1.SetParameters(yMax, seed, cfg.seedSigma_keV);  // A, mu, sigma
  f1.SetParLimits(2, cfg.sigmaMin_keV, cfg.sigmaMax_keV);
  hE->Fit(&f1, "RQ0");

  double mu = f1.GetParameter(1);
  double sg = std::fabs(f1.GetParameter(2));
  if (!(sg > 0.0)) return false;

  // core refit in [mu-refitSigmaMult*sg, mu+refitSigmaMult*sg] clipped
  const double c1 = std::max(p.win_lo, mu - cfg.refitSigmaMult * sg);
  const double c2 = std::min(p.win_hi, mu + cfg.refitSigmaMult * sg);

  TF1 f2("f2", "gaus", c1, c2);
  f2.SetParameters(f1.GetParameter(0), mu, sg);
  f2.SetParLimits(2, cfg.sigmaMin_keV, cfg.sigmaMax_keV);
  hE->Fit(&f2, "RQ0");

  p.mu  = f2.GetParameter(1);
  p.sig = std::fabs(f2.GetParameter(2));
  if (!(p.sig > 0.0)) return false;

  // ROI according to config
  if (cfg.roiMode == A80Config::ROIMode::kSigma) {
    p.roi_lo = p.mu - cfg.roiSigmaMult * p.sig;
    p.roi_hi = p.mu + cfg.roiSigmaMult * p.sig;
  } else { // kFixed
    p.roi_lo = p.mu - cfg.roiHalfWidth_keV;
    p.roi_hi = p.mu + cfg.roiHalfWidth_keV;
  }

  p.ok = true;
  return true;
}
