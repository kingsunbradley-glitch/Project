#pragma once
#include "A80Types.h"

class TH1D;

// Fit one peak by Gaussian inside p.win_lo..p.win_hi.
// Seed = maximum bin inside the fixed window.
// Fit range: [seed-W, seed+W] clipped, then a core refit in [mu-2sigma, mu+2sigma] clipped.
// On success, fills p.mu/p.sig and ROI = mu ± 3sigma.
bool FitOnePeakGaus(TH1D* hE, PeakWin& p, double W = 80.0);
