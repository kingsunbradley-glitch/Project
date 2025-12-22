#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include "TString.h"

// ===================================================================
//              全局配置文件 (Split Version: 分步运行版)
// ===================================================================

// --- 1. 运行配置 ---
inline constexpr int NUM_THREADS = 6; 

// --- 2. 文件路径配置 ---

// [输入] 原始 ROOT 文件 (TChain 通配符)
inline constexpr char INPUT_DIR_PATTERN[] = "/home/evalie2/Project/document/273Ds/nokb/SS032%05d_map.root";
inline constexpr int RUN_START = 2;
inline constexpr int RUN_END = 178;
inline constexpr char TREE_NAME[] = "tr_map";


// 预筛选输出的小文件（Step0 输出）
inline constexpr char DECAY_ONLY_FILE[]  = "SS032_decay_only.root";   // 衰变事件
inline constexpr char INJECT_ONLY_FILE[] = "SS032_inject_only.root";  // 注入事件

// [中间文件] 归一化程序(Step1)输出 -> 刻度程序(Step2)读取
// 解决报错: NORM_PARAM_FILE
inline constexpr char NORM_PARAM_FILE[] = "./SS032_Normalize_Params.txt";

// [诊断文件] 归一化过程的诊断图 (Step1 Output)
// 解决报错: NORM_DIAG_ROOT
inline constexpr char NORM_DIAG_ROOT[]  = "Diagnose_Normalize.root";

// [诊断文件] 刻度过程的诊断图 (Step2 Output)
// 解决报错: CALIB_DIAG_ROOT
inline constexpr char CALIB_DIAG_ROOT[] = "Diagnose_Calibration.root";

// [最终输出] 最终刻度系数 (Step2 Output)
inline constexpr char OUTPUT_DAT_FILE[] = "ener_Recal.dat";

// --- 3. 物理通道定义 ---
#define TOTAL_CH 224
inline constexpr int NUM_DSSDX_POS = 128; // 0-127
inline constexpr int NUM_DSSDY_POS = 48;  // 128-175
// YH 映射到 176-223

// --- 4. 归一化 (Step 1) 参数 ---
inline constexpr double NORM_MIN_E = 4000.0;
inline constexpr double NORM_MAX_E = 14000.0;
inline constexpr double NORM_E_DIFF = 200.0; // X-Y 能量差限制
inline constexpr int NORM_MIN_ENTRIES = 500; // 拟合最少计数
inline constexpr int NORM_REF_MIN_COUNTS = 2000;// 参考条至少计数


// 如果你想“手动指定”归一化参考条，把下面两行前面的 // 去掉：
// #define REF_STRIP_X 39   // 手动指定 X 面参考 strip
// #define REF_STRIP_Y 25   // 手动指定 Y 面参考 strip
// 如果这两行保持注释状态，则程序自动选择“统计量最多的条”作为参考。

// --- 5. 绝对刻度 (Step 2) 参数 ---
struct PeakDef {
    double energy;   // keV
    double window;   // keV (搜峰范围)
};

// 刻度用的 3 个参考峰
inline const std::vector<PeakDef> REF_PEAKS = {
    {6113, 50.0}, {7450, 50.0}, {8784, 50.0}
};

// 诊断用的 3 个峰 (算 FWHM)
inline const std::vector<PeakDef> DIAG_PEAKS = {
    {6113, 40.0}, {7129, 40.0}, {5307, 40.0}
};

// 绘图参数
inline constexpr int HIST_BINS = 2000;
inline constexpr double HIST_MIN = 4000.0;
inline constexpr double HIST_MAX = 14000.0;
inline constexpr bool ENABLE_BKG_SUB = true;
inline constexpr int BKG_ITERATIONS = 20;

#endif