#include <iostream>
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TSystem.h"
#include "config.h" // 引入配置文件

/**
 * @brief DSSD 事件筛选器 (Sorter)
 * 职责：
 * 1. 读取原始的、庞大的 *_map.root 文件。
 * 2. 应用全局的核心事件筛选条件。
 * 3. 将通过筛选的 "有效事件" 保存到一个新的、较小的 ROOT 文件中。
 * 这是整个分析流程的第 0 步，为后续的归一化分析做准备。
 */
ErrorCode DSSD_Sorter_Process() {
    ROOT::EnableImplicitMT(num_threads);
    std::cout << "\n--> Starting DSSD Event Sorter (Step 0)..." << std::endl;
    std::cout << "    Input files: " << input_root_file << std::endl;

    // 1. 从原始文件创建 RDataFrame
    ROOT::RDataFrame df(tree_name, input_root_file);
    auto n_initial = df.Count();
    
    if (*n_initial == 0) {
        std::cerr << "Error: No entries found in the input TTree." << std::endl;
        return INSUFFICIENT_DATA;
    }
    std::cout << "    Total initial entries: " << *n_initial << std::endl;

    // 2. 定义筛选条件 (与原 DSSD_normalize.cpp 中的核心过滤器一致)
    TString filter_str = TString::Format(
        "DSSDX_mul == 1 && DSSDY_mul == 1 " 
        "&& DSSDX_E[0] > %f && DSSDX_E[0] < %f "
        "&& DSSDY_E[0] > %f && DSSDY_E[0] < %f "
        "&& TMath::Abs(DSSDX_E[0] - DSSDY_E[0]) < %f",
        MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD,
        MIN_ENERGY_THRESHOLD, MAX_ENERGY_THRESHOLD,
        ENERGY_DIFF_LIMIT
    );

    if (REQUIRE_MWPC_ANTI) {
        filter_str += " && MWPC_mul == 0";
    }
    else {
        filter_str += " && MWPC_mul == 2"; // 不严格要求 MWPC
    }
    
    auto df_filtered = df.Filter(filter_str.Data(), "Core Event Filter");
    auto n_final = df_filtered.Count();

    // 3. 使用 Snapshot 将筛选后的数据保存到新的 ROOT 文件
    // Snapshot 会触发事件循环，并把结果写入新文件
    std::cout << "\n--> Filtering events and writing to new file..." << std::endl;
    std::cout << "    Output file: " << sorted_root_file << std::endl;
    
    // 检查输出目录是否存在，如果不存在则创建
    TString out_dir = gSystem->DirName(sorted_root_file);
    if (!gSystem->OpenDirectory(out_dir)) {
        gSystem->mkdir(out_dir, kTRUE);
    }

    df_filtered.Snapshot(tree_name, sorted_root_file);

    std::cout << "\n--> Sorter finished successfully!" << std::endl;
    if (*n_initial > 0) {
        std::cout << "    Entries after filtering: " << *n_final << " / " << *n_initial 
                  << " (" << 100.0 * (*n_final) / (*n_initial) << "%)" << std::endl;
    }

    return SUCCESS;
}

int main() {
    try {
        return DSSD_Sorter_Process();
    } catch (const std::exception& e) {
        std::cerr << "\nRuntime error in Sorter: " << e.what() << std::endl;
        return 2; 
    } catch (...) {
        std::cerr << "\nAn unknown error occurred in Sorter." << std::endl;
        return 3;
    }
}