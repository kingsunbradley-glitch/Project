#pragma once
#include "A80Types.h"
#include <vector>

class TH1D;
class TH2D;

// counts: per-bin counts (positive only recommended). Sorted descending internally.
AeqResult CalcAeqFromCounts(std::vector<double> counts, double targetFrac);

AeqResult CalcAeq2D(const TH2D* h, double targetFrac);
AeqResult CalcAeq1D(const TH1D* h, double targetFrac);
