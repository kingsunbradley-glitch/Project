#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include "TString.h"

// ===================================================================
//              全局配置文件 (Final Version)
// ===================================================================

// --- 1. 运行配置 ---
inline constexpr int NUM_THREADS = 10; 

// --- 2. 文件路径配置 ---

// [输入] 原始 ROOT 文件
inline constexpr char INPUT_DIR_PATTERN[] = "/home/evalie/Project/Document/After/SS032%05d_map.root";
inline constexpr int RUN_START = 2;
inline constexpr int RUN_END = 178;
inline constexpr char TREE_NAME[] = "tr_map";

// 预筛选输出的小文件（Step0 输出）
inline constexpr char DECAY_ONLY_FILE[]   = "SS032_decay_only.root";    // 衰变事件
// [修改] 注入事件文件名 Inject -> Implant
inline constexpr char IMPLANT_ONLY_FILE[] = "SS032_implant_only.root";  

// [中间文件] 归一化参数
inline constexpr char NORM_PARAM_FILE[] = "./SS032_Normalize_Params.txt";

// [诊断文件]
inline constexpr char NORM_DIAG_ROOT[]  = "Diagnose_Normalize.root";
inline constexpr char CALIB_DIAG_ROOT[] = "Diagnose_Calibration.root";

// [最终输出]
inline constexpr char OUTPUT_DAT_FILE[] = "ener_cal_Recal.dat";

// --- 3. 物理通道定义 ---
#define TOTAL_CH 224
inline constexpr int NUM_DSSDX_POS = 128; // 0-127
inline constexpr int NUM_DSSDY_POS = 48;  // 128-175
// YH 映射到 176-223

// --- 4. 归一化 (Step 1) 参数 ---
inline constexpr double NORM_MIN_E = 4000.0;
inline constexpr double NORM_MAX_E = 14000.0;
inline constexpr double NORM_E_DIFF = 200.0; 
inline constexpr int NORM_MIN_ENTRIES = 500; 
inline constexpr int NORM_REF_MIN_COUNTS = 2000;
// 手动指定归一化参考条：
#define REF_STRIP_X 10  //一定检查最热的条分辨如何,与Veto否决情况
#define REF_STRIP_Y 23

// --- 5. 绝对刻度 (Step 2) 参数 ---
struct PeakDef {
    double energy;   // keV
    double window;   // keV (搜峰范围)
};
double target_E = 6113.0; 
double fit_win  = 100.0; //仅仅搜索范围，double fit_radius = (win < 60.0) ? win/1.5 : 40.0;
// 刻度参考峰
inline const std::vector<PeakDef> REF_PEAKS = {
   // {8784, 30.0}, {7065, 30.0}, {7128, 30.0}
    {8783, 30.0}, {7065, 30.0}, {7128, 30.0}
};

// 诊断参考峰
inline const std::vector<PeakDef> DIAG_PEAKS = {
    {6113, 40.0}, {7686, 40.0}, {5304, 40.0}
};

// 绘图参数
inline constexpr int HIST_BINS = 2000;
inline constexpr double HIST_MIN = 4000.0;
inline constexpr double HIST_MAX = 14000.0;
inline constexpr bool ENABLE_BKG_SUB = false;
inline constexpr int BKG_ITERATIONS = 30;

#endif