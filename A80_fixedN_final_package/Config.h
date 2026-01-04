#pragma once
#include "A80Types.h"
#include <cstdint>

// A80 / X80 / Y80 QA 管线的“集中式参数配置”。
// 只需要在这里修改参数并重新编译即可，无需到各个 .cpp 里找魔法数字。
//
// 约定：
// - 能量单位：keV
// - 峰编号：0=Peak1，1=Peak2，2=Peak3
struct A80Config {
  // -------- 全局能量门（所有后续分析都在该门下进行） --------
  double EgLo = 5000.0;
  double EgHi = 6000.0;

  // -------- 能谱直方图分 bin（只影响找峰/拟合时的能谱形状） --------
  // 例如 Eg=[5000,6000]：
  // - nEbins=500  -> 2 keV/bin
  // - nEbins=1000 -> 1 keV/bin
  //
  // 注意：如果你用“峰窗口最大 bin 计数阈值”（minPeakMaxCount），
  // 改动 nEbins 会改变单 bin 计数大小；阈值可能需要相应调整。
  int nEbins = 500;

  // -------- A80 / X80 / Y80 的目标覆盖率 --------
  // A80eq 的定义是：覆盖 targetFrac(默认 80%) 计数所需的等效像素/等效 bin 数。
  double targetFrac = 0.80;

  // -------- 峰搜索与拟合相关参数 --------
  // 初次拟合窗口：以 seed 峰位为中心的局部窗口 [seed-W, seed+W]
  double fitWindow_keV = 80.0;

  // 二次拟合（refit）核心窗口： [mu - refitSigmaMult*sigma, mu + refitSigmaMult*sigma]
  // 典型取值：1.0（更“核心”、更抗背景）、2.0（更稳健）
  double refitSigmaMult = 1.0;

  // 峰搜索窗口内“最高 bin 的计数”小于该阈值时，认为统计不足，不进行拟合。
  // 这样可以避免低统计 run 在峰位/σ 上乱飘。
  double minPeakMaxCount = 20.0;

  // σ 的初值（用于拟合初始化）以及允许范围（keV）
  double seedSigma_keV = 20.0;
  double sigmaMin_keV = 2.0;
  double sigmaMax_keV = 200.0;

  // -------- ROI 定义（用于填 XY / 计数 Npk / 计算每峰 A80） --------
  enum class ROIMode { kSigma, kFixed };

  // 选择 ROI 模式：
  // - kSigma：ROI = mu ± roiSigmaMult * sigma   （随峰宽变化，统计口径会随 run 改变）
  // - kFixed：ROI = mu ± roiHalfWidth_keV      （固定宽度，便于不同 run 之间“统一口径”比较）
  ROIMode roiMode = ROIMode::kFixed;

  // roiMode = kSigma 时使用：
  double roiSigmaMult = 3.0;

  // roiMode = kFixed 时使用（例如统一成 ±100 keV）：
  double roiHalfWidth_keV = 100.0;

  // -------- Fixed-N A80（低统计偏差修正） --------
  // 背景：A80 是“覆盖 80% 计数所需的等效像素数/像素数”。当总计数很低时，
  //       事件天然只会落在少量像素上，A80 会“看起来更小”（假聚焦）。
  // 方法：对每个峰在 ROI 内的事件流做蓄水池抽样（reservoir sampling），
  //       均匀随机抽取固定数量 N0 个事件（N0=peakFixedN），只用这 N0 个事件计算 A80_fixN。
  // 规则：若该峰 ROI 内事件数 < N0，则不计算 A80_fixN（summary 里 Nfix=0、A80eq_fixN=-999）。
  bool enableFixedN = true;
  int  peakFixedN   = 5000;      // N0：固定抽样事件数
  uint64_t fixedNSeedBase = 12345ULL; // 基础随机种子（+runnum/峰号 生成每个 run 的可复现随机序列）

  // -------- XY 直方图分 bin（会直接影响 A80 的数值尺度） --------
  // 默认按像素化 index（X:128, Y:48）。如果改成更粗的 bin（例如 64x24），
  // A80 的绝对数值会随之变化（因为“一个 bin”代表更大的区域），比较前请统一标准。
  int nx = 128;
  int ny = 48;
  double xlo = -0.5, xhi = 127.5;
  double ylo = -0.5, yhi = 47.5;

  // -------- 三个峰的搜索窗口模板（keV） --------
  // 这里定义的是“搜索窗口”，不是 ROI；ROI 由拟合得到的 mu/sigma（或固定宽度）给出。
  PeakWin Ptpl[3];

  A80Config();
};

inline A80Config::A80Config() {
  Ptpl[0].win_lo = 5100; Ptpl[0].win_hi = 5200;
  Ptpl[1].win_lo = 5435; Ptpl[1].win_hi = 5535;
  Ptpl[2].win_lo = 5754; Ptpl[2].win_hi = 5854;
}
//search windows for the 3 peaks
