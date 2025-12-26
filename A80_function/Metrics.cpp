#include "Metrics.h"

#include <TH1D.h>
#include <TH2D.h>

#include <algorithm>

AeqResult CalcAeqFromCounts(std::vector<double> counts, double targetFrac) {
  AeqResult r;
  if (counts.empty()) return r;

  double tot = 0.0;
  for (double c : counts) tot += c;
  if (tot <= 0.0) return r;

  std::sort(counts.begin(), counts.end(), std::greater<double>());

  const double target = targetFrac * tot;

  double acc = 0.0, accPrev = 0.0;
  long long nFull = 0;

  for (double c : counts) {
    accPrev = acc;
    acc += c;
    ++nFull;
    if (acc >= target) {
      r.n_pix = nFull;
      r.cov_at_npix = acc / tot;

      const double need = target - accPrev; // (0..c]
      double f = (c > 0.0) ? (need / c) : 1.0;
      if (f < 0.0) f = 0.0;
      if (f > 1.0) f = 1.0;
      r.n_eq = (nFull - 1) + f;
      return r;
    }
  }

  // target not reached (shouldn't happen if tot>0)
  r.n_pix = static_cast<long long>(counts.size());
  r.n_eq = static_cast<double>(counts.size());
  r.cov_at_npix = 1.0;
  return r;
}

AeqResult CalcAeq2D(const TH2D* h, double targetFrac) {
  AeqResult r;
  if (!h) return r;

  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  std::vector<double> counts;
  counts.reserve(static_cast<size_t>(nx) * static_cast<size_t>(ny));

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      const double c = h->GetBinContent(ix, iy);
      if (c > 0.0) counts.push_back(c);
    }
  }
  return CalcAeqFromCounts(std::move(counts), targetFrac);
}

AeqResult CalcAeq1D(const TH1D* h, double targetFrac) {
  AeqResult r;
  if (!h) return r;

  const int nb = h->GetNbinsX();
  std::vector<double> counts;
  counts.reserve(static_cast<size_t>(nb));

  for (int b = 1; b <= nb; ++b) {
    const double c = h->GetBinContent(b);
    if (c > 0.0) counts.push_back(c);
  }
  return CalcAeqFromCounts(std::move(counts), targetFrac);
}
