#pragma once
#include "A80Types.h"
#include <vector>

class TH1D;
class TH2D;

// counts: per-bin counts (positive only recommended). Sorted descending internally.
AeqResult CalcAeqFromCounts(std::vector<double> counts, double targetFrac);

AeqResult CalcAeq2D(const TH2D* h, double targetFrac);
AeqResult CalcAeq1D(const TH1D* h, double targetFrac);


// 固定事件数（Fixed-N）抽样后的 2D Aeq 计算：
// 传入的是每个事件落入的像素 bin 下标 idx = iy*nx + ix（0-based）。
// 该函数会把这些 idx 统计成每像素计数，再调用 CalcAeqFromCounts。
AeqResult CalcAeq2DFromSampledBins(const std::vector<int>& binIdx, int nx, int ny, double targetFrac);
