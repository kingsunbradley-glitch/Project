#ifndef CONFIG_H
#define CONFIG_H

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "TString.h"

// ===================================================================
//              全局配置文件
// ===================================================================

// --- 0. 与 ANA_SHANS2 共享的 setup.h 对齐 ---
// 本程序原先自带一套道号常量，与 root2root / Cal_Auto_* 那套 setup.h 重复。
// 两份定义迟早会漂移（例如 DSSD_YH 的区间改了这边没改），所以这里直接吃共享
// 的 setup.h，再从 DSSD_*_CH_LOW/HIGH 推导本程序用到的量。
//
// 用 makefile 里的 ANA_SETUP_DIR (-I) 指过去；找不到就退回下面的内置兜底值，
// 两种情况都在启动日志里打印一行 [INFO]，不会静默。
#if defined(__has_include)
#  if __has_include("setup.h")
#    include "setup.h"
#    define HAVE_ANA_SETUP 1
#  endif
#endif

#ifndef HAVE_ANA_SETUP
#  define HAVE_ANA_SETUP 0
   // 内置兜底：与 2026-07 版 ANA setup.h 的 DSSD 段一致
#  define DSSD_X_CH_LOW    0
#  define DSSD_X_CH_HIGH   127
#  define DSSD_Y_CH_LOW    128
#  define DSSD_Y_CH_HIGH   175
#  define DSSD_YH_CH_LOW   176
#  define DSSD_YH_CH_HIGH  223
#  define RECONVERT_CALIBRATION_FILENAME "ener_cal_Recal.dat"
#endif

// --- 1. 运行配置 ---
inline constexpr int NUM_THREADS = 10;

// setup.h 快照相关的文件名（与 makefile 的 SETUP_LOCAL / SETUP_STAMP 保持一致）
inline constexpr char SETUP_LOCAL_FILE[] = "setup.h";
inline constexpr char SETUP_STAMP_FILE[] = "setup_stamp.txt"; 

// --- 2. 文件路径配置 ---
//
// 每台机器不同的量（数据路径、run 范围、归一化能量窗、参考条）建议写在 ./MyConfig.h 里：
// 升级程序时只覆盖 *.cpp / *.h，不用再改 Config.h，本地设置也不会被冲掉。
// 模板见 MyConfig.h.example（先 cp 一份再改）。
#if defined(__has_include)
#  if __has_include("MyConfig.h")
#    include "MyConfig.h"
#    define HAVE_MY_CONFIG 1
#  endif
#endif

#ifndef MY_INPUT_DIR_PATTERN
#  define MY_INPUT_DIR_PATTERN "/home/evalie/Project/Document/After/SS032%05d_map.root"
#endif
#ifndef MY_RUN_START
#  define MY_RUN_START 2
#endif
#ifndef MY_RUN_END
#  define MY_RUN_END 178
#endif

// [输入] 原始 ROOT 文件
inline constexpr char INPUT_DIR_PATTERN[] = MY_INPUT_DIR_PATTERN;
inline constexpr int RUN_START = MY_RUN_START;
inline constexpr int RUN_END   = MY_RUN_END;
inline constexpr char TREE_NAME[] = "tr_map";

// 预筛选输出的小文件（Step0 输出）
inline constexpr char DECAY_ONLY_FILE[]   = "SS032_decay_only.root";    // 衰变事件
inline constexpr char IMPLANT_ONLY_FILE[] = "SS032_implant_only.root";  // 注入事件

// [中间文件] 归一化参数
inline constexpr char NORM_PARAM_FILE[] = "./SS032_Normalize_Params.txt";

// [诊断文件]
inline constexpr char NORM_DIAG_ROOT[]  = "Diagnose_Normalize.root";
inline constexpr char CALIB_DIAG_ROOT[] = "Diagnose_Calibration.root";

// [最终输出] 名字来自共享 setup.h 的 RECONVERT_CALIBRATION_FILENAME，
// 因为读这个文件的是 root2root，两边必须一致。
inline constexpr char OUTPUT_DAT_FILE[] = RECONVERT_CALIBRATION_FILENAME;

// --- 3. 物理通道定义（全部由 setup.h 推导，不要在这里另写数字）---
inline constexpr int CH_BASE_X    = DSSD_X_CH_LOW;    //  0
inline constexpr int CH_BASE_Y    = DSSD_Y_CH_LOW;    //  128
inline constexpr int CH_BASE_YH   = DSSD_YH_CH_LOW;   //  176
inline constexpr int NUM_DSSDX_POS  = DSSD_X_CH_HIGH  - DSSD_X_CH_LOW  + 1;   // 128
inline constexpr int NUM_DSSDY_POS  = DSSD_Y_CH_HIGH  - DSSD_Y_CH_LOW  + 1;   //  48
inline constexpr int NUM_DSSDYH_POS = DSSD_YH_CH_HIGH - DSSD_YH_CH_LOW + 1;   //  48
// X/Y/YH 三个面的全局道上界（= 224）。注意不能叫 TOTAL_CH：共享 setup.h 里
// TOTAL_CH 是整个 DAQ 的道数（39*8 = 312），同名会冲突。
inline constexpr int TOTAL_DSSD_CH = std::max({DSSD_X_CH_HIGH, DSSD_Y_CH_HIGH, DSSD_YH_CH_HIGH}) + 1;

// --- 4. 归一化 (Step 1) 参数 ---
#ifndef MY_NORM_MIN_E
#  define MY_NORM_MIN_E 4000.0
#endif
#ifndef MY_NORM_MAX_E
#  define MY_NORM_MAX_E 14000.0
#endif
#ifndef MY_NORM_E_DIFF
#  define MY_NORM_E_DIFF 200.0
#endif
inline constexpr double NORM_MIN_E = MY_NORM_MIN_E;
inline constexpr double NORM_MAX_E = MY_NORM_MAX_E;
inline constexpr double NORM_E_DIFF = MY_NORM_E_DIFF;
inline constexpr int NORM_MIN_ENTRIES = 500; 
inline constexpr int NORM_REF_MIN_COUNTS = 2000;
// 手动指定归一化参考条（可以在 MyConfig.h 里改）：
#ifndef REF_STRIP_X
#  define REF_STRIP_X 10  //一定检查最热的条分辨如何,与Veto否决情况
#endif
#ifndef REF_STRIP_Y
#  define REF_STRIP_Y 23  //小心X:40，90的Veto缝隙
#endif

// --- 5. 绝对刻度 (Step 2) 参数 ---
struct PeakDef {
    double energy;   // keV
    double window;   // keV (搜峰范围)
};

// 参考峰表（运行期读取，换实验只改文件不用重编译）。
// 下面两个 _DEFAULT 只在 ref_peaks.txt 读不到或解析为空时兜底。
inline constexpr char REF_PEAK_FILE[] = "./ref_peaks.txt";

// 逐条一致性检查用的峰（输出表 FWHM / dE 两列），可被参考峰表的 [strip] 节覆盖
inline constexpr double target_E = 6113.0;
inline constexpr double fit_win  = 100.0; //仅仅搜索范围，double fit_radius = (win < 60.0) ? win/1.5 : 40.0;

inline const std::vector<PeakDef> REF_PEAKS_DEFAULT = {
    {8783, 30.0}, {7065, 30.0}, {7128, 30.0}
};

// 诊断参考峰（只用于刻度后验收，不参与拟合）
inline const std::vector<PeakDef> DIAG_PEAKS_DEFAULT = {
    {6113, 40.0}, {7686, 40.0}, {5304, 40.0}
};

// --- 5b. 绝对刻度的检查 ---
// 这里的量都不否决拟合结果，只决定日志怎么报告。原因：本程序也可能被用来处理
// 原始 ADC（不做外刻度），那种情况下 K/B 本来就该偏离 1/0，不能当成错误。
// 真正阻止使用拟合结果的只有两种情况：能定位的参考峰少于 2 个（直线定不出来），
// 或者拟合给出非物理的 K<=0 / 非有限值。
//
// 参考峰之间的最小间隔（以拟合出的 sigma 为单位）：小于该值的两点视为同一物理峰
inline constexpr double REF_DUP_SIGMA = 2.0;
// 定 (K,B) 所需的最少参考峰数（直线至少要 2 点）
inline constexpr int    REF_MIN_POINTS = 2;
// 有 3 点以上才有冗余，才能从残差看出模型对不对
inline constexpr int    REF_GOOD_POINTS = 3;
// 参考峰能量跨度低于该值时提醒"杠杆不足"（keV）
inline constexpr double REF_WARN_SPAN_KEV = 800.0;
// 单点残差超过该值时提醒（keV）
inline constexpr double REF_WARN_RESID = 8.0;
// 增益/零点偏离 1/0 超过以下范围时提醒（不否决）
inline constexpr double REF_WARN_K_TOL = 0.05;
inline constexpr double REF_WARN_B_TOL = 200.0;

// 绘图参数
inline constexpr int HIST_BINS = 2000;
inline constexpr double HIST_MIN = 4000.0;
inline constexpr double HIST_MAX = 14000.0;
inline constexpr bool ENABLE_BKG_SUB = false;
inline constexpr int BKG_ITERATIONS = 30;

// 启动时打印一次"用的哪份 setup.h / 哪个道号区间"，避免两份定义静默分叉
inline void Print_Channel_Map() {
    std::cout << "--> 本程序编译于 " << __DATE__ << " " << __TIME__ << std::endl;
#if HAVE_MY_CONFIG
    std::cout << "--> [MyConfig.h] 已加载本地配置" << std::endl;
#else
    std::cout << "--> [MyConfig.h] 未找到，使用 Config.h 默认路径" << std::endl;
#endif
    // 区分"本地快照"和"-I 共享"：sync-setup 会写下 setup_stamp.txt
    std::string stamp_src, stamp_when;
    {
        std::ifstream st(SETUP_STAMP_FILE);
        std::string line;
        while (std::getline(st, line)) {
            auto val = [&](const char* key) {
                std::string k(key);
                if (line.rfind(k, 0) == 0) {
                    auto p = line.find(':');
                    if (p != std::string::npos) {
                        std::string v = line.substr(p + 1);
                        auto a = v.find_first_not_of(" \t");
                        return (a == std::string::npos) ? std::string() : v.substr(a);
                    }
                }
                return std::string();
            };
            if (auto v = val("source"); !v.empty()) stamp_src = v;
            if (auto v = val("copied"); !v.empty()) stamp_when = v;
        }
    }

    std::cout << "--> [setup.h] ";
    if (!HAVE_ANA_SETUP) {
        std::cout << "未找到，使用 Config.h 内置兜底值";
    } else if (!stamp_src.empty()) {
        std::cout << "本地快照 ./" << SETUP_LOCAL_FILE;
        if (!stamp_when.empty()) std::cout << "（快照于 " << stamp_when << "）";
        std::cout << "，源: " << stamp_src;
    } else {
        std::cout << "已加载 ANA 共享 setup.h";
    }
    std::cout << "\n    DSSD X " << DSSD_X_CH_LOW << "-" << DSSD_X_CH_HIGH
              << ", Y "      << DSSD_Y_CH_LOW << "-" << DSSD_Y_CH_HIGH
              << ", YH "     << DSSD_YH_CH_LOW << "-" << DSSD_YH_CH_HIGH
              << " (共 " << TOTAL_DSSD_CH << " 道)" << std::endl;
}

#endif
