#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

// ===================================================================
//                 项目全局配置文件 (Project Global Config)
// ===================================================================

// --- 1. 文件路径配置 ---
// 原始输入文件 (hadd合并后)
const char* input_root_file = "/home/evalie2/Project/document/273Ds/inter2/SS032_inter2_chain.root";
// 第1步程序输出、第2步程序输入的刻度参数文件 (更名以消除歧义)
const char* calibration_param_file = "/home/evalie2/Project/document/273Ds/inter2/SS032_SSD_calibration_params_pro.txt";
// 第1步程序输出的诊断图目录
const char* output_plot_dir = "/home/evalie2/Project/Function/273Ds/calibration_plots_pro";
// 第2步程序输出的最终ROOT文件
const char* final_calibrated_root_file = "/home/evalie2/Project/document/273Ds/inter2/SS032_inter2_chain_SSDcal_pro.root";


// --- 2. TTree 名称 ---
const char* tree_name = "tr_chain";


// --- 3. 刻度峰配置 ---
// 用于刻度的参考α粒子能量 (keV) 及其搜索窗口 (channel)
struct CalibrationPeak {
    double win_min;
    double win_max;
    double ref_E;
};

std::vector<CalibrationPeak> peaks_of_interest = {
    {6000, 6400, 6111.0},
    {7300, 7700, 7447.0},
    {5100, 5500, 5307.0}
};

// --- 3.1 背景扣除参数 ---
// 是否启用背景扣除功能 (true = 开启, false = 关闭)
const bool enable_background_subtraction = true; 
// TSpectrum 背景扣除的迭代次数
const int background_iterations = 20;

// --- 3.2 寻峰参数 ---
// TSpectrum 寻峰的高斯 sigma (单位: bin)
const double peak_search_sigma = 2.0;
// TSpectrum 寻峰的阈值 (0到1之间)
const double peak_search_threshold = 0.05;

// --- 4. 物理和程序参数 ---
// SSD总通道数
const int NUM_SSD_POS = 48;
// 每个通道最少的事例数，低于此值则不进行刻度
const int MIN_HIST_ENTRIES = 500;
// TTree中数组的最大长度 (用于防止越界)
const int MAX_MULT = 256; 
// 应用刻度时的最低能量阈值 (channel), 低于此值输出能量为0
const double ENERGY_THRESHOLD = 1.0;
// 输出ROOT文件的压缩级别 (0-9)
const int ROOT_COMPRESSION_LEVEL = 1;

// --- 5. 错误码定义 ---
// 定义统一的错误码，用于从main函数返回
enum ErrorCode {
    SUCCESS = 0,
    FILE_NOT_FOUND = 1,
    TREE_NOT_FOUND = 2,
    BRANCH_NOT_FOUND = 3,
    INSUFFICIENT_DATA = 4,
    CONFIG_ERROR = 5
};

#endif // CONFIG_H
