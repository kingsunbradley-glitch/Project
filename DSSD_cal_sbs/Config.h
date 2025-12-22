#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include "TString.h"

// ================= 用户配置区域 =================

// --- 多线程配置 ---
// 解决 'NUM_THREADS' was not declared 问题
static const int NUM_THREADS = 10; 

// --- 输入文件配置 ---
// 文件路径模板，%04d 会被替换为 0002, 0178 等
// 假设文件路径是: /home/.../SS03200002_map.root
static const char* INPUT_DIR_PATTERN = "/home/evalie2/Project/document/273Ds/inter_map/SS032%05d_map.root";
static const int RUN_START =2;   // 开始编号 (如 2)
static const int RUN_END = 178;   // 结束编号 (如 178)

// --- 输出文件 ---
static const char* OUTPUT_DAT_FILE = "ener_cal_Recal.dat";
static const char* OUTPUT_ROOT_FILE = "ener_cal_diagnostics.root";
static const char* TREE_NAME = "tr_map"; 

// --- 探测器通道定义 ---
#define DSSD_X_CH_LOW  0
#define DSSD_X_CH_HIGH 127
#define DSSD_Y_CH_LOW  128
#define DSSD_Y_CH_HIGH 175
#define DSSD_YH_CH_LOW 176
#define DSSD_YH_CH_HIGH 223
#define TOTAL_CH       224

// --- 峰值定义 ---
struct PeakDef {
    double energy;   
    double window;   
};

// 1. 参考峰 
static const std::vector<PeakDef> REF_PEAKS = {
    {7129, 20.0},  // 218Rn 
    {7312, 20.0},  // 219Fr 
    {8784, 20.0}   // 212Po 
};

// 2. 诊断峰 
static const std::vector<PeakDef> DIAG_PEAKS = {
    {6113, 15.0},    // 242Cm
    {6622, 15.0},    // 211Bi
    {7066, 15.0}     // 217At
};

// --- 其他参数 ---
static const bool ENABLE_BKG_SUB = false;
static const int BKG_ITERATIONS = 20;

static const int HIST_BINS = 2000;
static const double HIST_MIN = 4000.0;
static const double HIST_MAX = 14000.0;

#endif