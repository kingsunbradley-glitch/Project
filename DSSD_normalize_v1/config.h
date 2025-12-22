#ifndef CONFIG_DSSD_H
#define CONFIG_DSSD_H

#include <vector>
#include <string>

// ===================================================================
//              项目全局配置文件 (DSSD Normalization Config)
// ===================================================================

// --- 1. 实验信息和文件路径配置 ---
// 实验编号，用于自动生成输出文件名和目录
const char* EXPERIMENT_ID = "SS032";

// [原始输入]: 指向包含原始 ADC/Channel 数据的 *_map.root 文件
// 这是 Sorter 程序 (第 0 步) 的输入
const char* input_root_file = "/home/evalie2/Project/document/273Ds/inter_map/SS0320001*_map.root";

// [中间文件]: 由 Sorter 程序生成的、经过筛选的 ROOT 文件
// 这是 Normalizer 程序 (第 1 步) 的输入
const char* sorted_root_file = "/home/evalie2/Project/Function/273Ds/DSSD_normalizev1/SS032_sorted_events.root";

// [最终输出]: 第 1 步程序输出的归一化参数文件
const char* normalize_param_file = "/home/evalie2/Project/Function/273Ds/DSSD_normalizev1/SS032_DSSD_normalize_params_pro2.txt";

// [最终输出]: 第 1 步程序输出的诊断图/ROOT文件目录
const char* output_plot_dir = "/home/evalie2/Project/Function/273Ds/DSSD_normalizev1/SS032_DSSD_normalize_pro2";

// TTree 名称 (与 Data_Map.h/cxx 中的一致)
const char* tree_name = "tr_map";


// --- 3. 物理和程序参数 ---  
// DSSD-X 总通道数
const int NUM_DSSDX_POS = 128; 
// DSSD-Y 和 DSSD-YH 总通道数
const int NUM_DSSDY_POS = 48;

// 每个通道进行归一化拟合时最少的事例数
const int MIN_FIT_ENTRIES = 100; 

// 最小能量阈值 (ADC Channel), 低于此值的数据不用于归一化拟合
const double MIN_ENERGY_THRESHOLD = 1000.0;
// 最大能量阈值 (ADC Channel)
const double MAX_ENERGY_THRESHOLD = 3200000.0;

// 归一化拟合的最小能量差要求 (可选，用于排除双重粒子)
const double ENERGY_DIFF_LIMIT = 200.0; 

// 是注入还是裂变 (MWPC_mul == 0 裂变 true ，MWPC_mul == 2 好注入 false)
const bool REQUIRE_MWPC_ANTI = true;

// 使用的线程数
const int num_threads = 12; // 默认使用 12 核

// --- 4. 错误码定义 ---
enum ErrorCode {
    SUCCESS = 0,
    FILE_NOT_FOUND = 1,
    TREE_NOT_FOUND = 2,
    INSUFFICIENT_DATA = 3,
    FIT_FAILED = 4,
    CONFIG_ERROR = 5
};

#endif // CONFIG_DSSD_H