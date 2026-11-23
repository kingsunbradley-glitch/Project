#ifndef DSSD_NORM_CONFIG_H
#define DSSD_NORM_CONFIG_H

#include <setup.h>

namespace NormConfig
{
// 幅度窗口的单位与 tr_map 的 *_E 一致。
constexpr double MIN_X = 600.0;
constexpr double MIN_Y = 600.0;
constexpr double MIN_YH = 100.0;
constexpr double MAX_AMPLITUDE = 32100.0;
constexpr bool REJECT_PILEUP = true;
constexpr int MIN_PIXEL_EVENTS = 100;
constexpr double MIN_X_SPAN = 500.0; // 第90和第10百分位之差
constexpr unsigned MAX_PIXEL_SAMPLE = 3000; // 等概率 reservoir 抽样，0 表示全部
constexpr unsigned SEED = 20261115;
constexpr int SEED_TRIALS = 128;
constexpr double NOISE_FLOOR = 0.5; // X 幅度单位
// sigma_back / sigma_X，使用输入幅度单位的测量误差。
constexpr double SIGMA_Y_OVER_X = 1.0;
constexpr double SIGMA_YH_OVER_X = 1.0;
constexpr double MIN_INLIER_FRACTION = 0.6;
constexpr double MAX_REL_SLOPE_ERROR = 0.05;
constexpr int REFERENCE_X = -1; // 局部条号；-1 自动选择
constexpr double OFFSET_SCALE = 1000.0; // 全局拟合参数的数值缩放
}

#endif
