#include "PeakFinder.h"

#include <TH1D.h>
#include <TF1.h>

#include <algorithm>
#include <cmath>

bool FitOnePeakGaus(TH1D* hE, PeakWin& p, double W) {
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

  const double seed = hE->GetBinCenter(bMax);

  // local fit range: [seed-W, seed+W] clipped to window
  const double r1 = std::max(p.win_lo, seed - W);
  const double r2 = std::min(p.win_hi, seed + W);

  TF1 f1("f1", "gaus", r1, r2);
  f1.SetParameters(yMax, seed, 20.0);  // A, mu, sigma
  f1.SetParLimits(2, 2.0, 200.0);
  hE->Fit(&f1, "RQ0");

  double mu = f1.GetParameter(1);
  double sg = std::fabs(f1.GetParameter(2));
  if (!(sg > 0.0)) return false;

  // core refit in [mu-2sg, mu+2sg] clipped to window
  const double c1 = std::max(p.win_lo, mu - 2.0 * sg);
  const double c2 = std::min(p.win_hi, mu + 2.0 * sg);

  TF1 f2("f2", "gaus", c1, c2);
  f2.SetParameters(f1.GetParameter(0), mu, sg);
  f2.SetParLimits(2, 2.0, 200.0);
  hE->Fit(&f2, "RQ0");

  p.mu  = f2.GetParameter(1);
  p.sig = std::fabs(f2.GetParameter(2));
  if (!(p.sig > 0.0)) return false;

  p.roi_lo = p.mu - 3.0 * p.sig;
  p.roi_hi = p.mu + 3.0 * p.sig;
  p.ok = true;
  return true;
}
