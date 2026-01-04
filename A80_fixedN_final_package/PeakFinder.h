#pragma once
#include "Config.h"

class TH1D;

// 在 p.win_lo..p.win_hi 的搜索窗口内对单个峰做高斯拟合（两步：粗拟合 + 核心 refit）。
// - seed：窗口内最高 bin
// - 初次拟合： [seed - cfg.fitWindow_keV, seed + cfg.fitWindow_keV] 与 [win_lo,win_hi] 取交集
// - refit：    [mu - cfg.refitSigmaMult*sigma, mu + cfg.refitSigmaMult*sigma] 与 [win_lo,win_hi] 取交集
// - ROI：由 cfg.roiMode 决定：
//    * kSigma：mu ± cfg.roiSigmaMult*sigma
//    * kFixed：mu ± cfg.roiHalfWidth_keV
//
// 返回 true 表示拟合成功并写入 p.mu/p.sig/p.roi_lo/p.roi_hi/p.ok。
bool FitOnePeakGaus(TH1D* hE, PeakWin& p, const A80Config& cfg);
