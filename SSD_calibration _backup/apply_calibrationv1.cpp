// g++ apply_calibrationv1.cpp `root-config --cflags --libs` -o apply_calibration
// 该程序读取由 SSD_calibration.cpp 生成的参数文件，
// 并创建一个新的 ROOT 文件，其中 SSD_E 分支已被替换为刻度后的能量值。
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm> // For std::max
#include <memory>    // For std::unique_ptr

#include "config.h"

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

ErrorCode apply_calibration_replace() {
    // --- 1. 读取刻度参数文件 ---
    std::cout << "--> Reading calibration parameters from " << calibration_param_file << "..." << std::endl;
    std::ifstream param_in(calibration_param_file);
    if (!param_in.is_open()) {
        std::cerr << "Error: Cannot open parameter file: " << calibration_param_file << std::endl;
        return FILE_NOT_FOUND;
    }
    std::map<int, double> k_params, b_params;
    std::string line;
    while (std::getline(param_in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        int pos, conf;
        double k, b;
        if (ss >> pos >> k >> b >> conf) {
            k_params[pos] = k;
            b_params[pos] = b;
        }
    }
    param_in.close();
    std::cout << "--> Read " << k_params.size() << " parameter sets." << std::endl;

    // --- 2. 准备 ROOT 文件和 TTree ---
    std::cout << "\n--> Preparing ROOT files..." << std::endl;
    auto inFile = std::unique_ptr<TFile>(TFile::Open(input_root_file));
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Cannot open input ROOT file: " << input_root_file << std::endl;
        return FILE_NOT_FOUND;
    }

    TTree *inTree = (TTree*)inFile->Get(tree_name);
    if (!inTree) {
        std::cerr << "Error: Cannot find TTree with name '" << tree_name << "' in file: " << input_root_file << std::endl;
        return TREE_NOT_FOUND;
    }
    Long64_t n_entries = inTree->GetEntries();

    // --- 创建输出文件和克隆的TTree ---
    // 使用智能指针管理输出文件
    auto outFile = std::make_unique<TFile>(final_calibrated_root_file, "RECREATE", "", ROOT_COMPRESSION_LEVEL);
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Cannot create output ROOT file: " << final_calibrated_root_file << std::endl;
        return FILE_NOT_FOUND;
    }
    // CloneTree(0) 会完整复制inTree的所有分支结构
    TTree* outTree = inTree->CloneTree(0);
    outTree->SetDirectory(outFile.get());

    // --- 关键修改：重新链接输出树的 SSD_E 分支 ---
    UInt_t chain_mul = 0;
    Double_t SSD_E_raw[MAX_MULT], SSD_Pos[MAX_MULT]; // 用于读取原始值
    Double_t SSD_E_calibrated[MAX_MULT];            // 用于存放计算后的新值

    // 从输入树中读取原始值，并检查分支是否存在
    if (inTree->SetBranchAddress("chain_mul", &chain_mul) != 0) {
        std::cerr << "Error: Cannot find branch 'chain_mul' in input TTree." << std::endl;
        return BRANCH_NOT_FOUND;
    }
    if (inTree->SetBranchAddress("SSD_E", SSD_E_raw) != 0) {
        std::cerr << "Error: Cannot find branch 'SSD_E' in input TTree." << std::endl;
        return BRANCH_NOT_FOUND;
    }
    if (inTree->SetBranchAddress("SSD_Pos", SSD_Pos) != 0) {
        std::cerr << "Error: Cannot find branch 'SSD_Pos' in input TTree." << std::endl;
        return BRANCH_NOT_FOUND;
    }

    // 将输出树的 SSD_E 分支的地址指向我们存放新值的数组
    outTree->SetBranchAddress("SSD_E", SSD_E_calibrated);
    
    // --- 3. 遍历数据, 计算新值并填充新树 ---
    std::cout << "\n--> Processing events to create new calibrated file..." << std::endl;
    long calibrated_count = 0;
    long uncalibrated_count = 0;
    long total_hits = 0;

    for (Long64_t i = 0; i < n_entries; ++i) {
        // 改进进度显示，防止 n_entries 为 0 时崩溃
        if (n_entries > 0 && i > 0 && i % std::max(1LL, n_entries / 10) == 0) {
            std::cout << "    Processing event " << i << " of " << n_entries << " (" << int(100.0 * i / n_entries) << "%)" << std::endl;
        }
        inTree->GetEntry(i); // 读取当前事件的所有数据到inTree的缓冲区

        // 检查数组越界风险
        UInt_t safe_mult = (chain_mul > MAX_MULT) ? MAX_MULT : chain_mul;
        if (chain_mul > MAX_MULT) {
             std::cerr << "Warning: Event " << i << " has chain_mul (" << chain_mul 
                       << ") which exceeds MAX_MULT (" << MAX_MULT << "). Truncating to " << MAX_MULT << "." << std::endl;
        }

        for (UInt_t j = 0; j < safe_mult; ++j) {
            total_hits++;
            int pos_idx = static_cast<int>(TMath::Nint(SSD_Pos[j]));
            
            if (SSD_E_raw[j] < ENERGY_THRESHOLD) {
                SSD_E_calibrated[j] = 0.0;
                uncalibrated_count++; // 低于阈值的视为未刻度
            } else {
                if (k_params.count(pos_idx)) {
                    double calculated_E = SSD_E_raw[j] * k_params.at(pos_idx) + b_params.at(pos_idx);
                    SSD_E_calibrated[j] = (calculated_E > 0) ? calculated_E : 0.0;
                    calibrated_count++;
                } else {
                    SSD_E_calibrated[j] = SSD_E_raw[j]; // 没有刻度参数，使用原始值
                    uncalibrated_count++;
                }
            }
        }
        // 对于被截断的部分，用0填充
        for (UInt_t j = safe_mult; j < chain_mul; ++j) {
            SSD_E_calibrated[j] = 0.0;
        }

        outTree->Fill();
    }
    
    // --- 4. 写入并关闭文件 ---
    std::cout << "--> Writing and closing the output ROOT file..." << std::endl;
    outFile->cd();
    outTree->Write();
    outFile->Close();
    
    // --- 5. 打印总结信息 ---
    std::cout << "\n--> Process finished successfully!" << std::endl;
    std::cout << "--> New file created: " << final_calibrated_root_file << std::endl;
    std::cout << "\n--- Calibration Summary ---" << std::endl;
    std::cout << "Total hits processed: " << total_hits << std::endl;
    std::cout << "  - Calibrated hits:    " << calibrated_count << std::endl;
    std::cout << "  - Uncalibrated hits:  " << uncalibrated_count << std::endl;
    
    return SUCCESS;
}

int main() {
    return apply_calibration_replace();
}
