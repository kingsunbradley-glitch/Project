#include <TChain.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>    // <-- 修正：需要包含 TFile 的头文件！
#include <TTree.h>    // <-- 修正：需要包含 TTree 的头文件！
#include <TObjArray.h>
#include <TObjString.h>
#include <iostream>
#include <iomanip>

// 函数名与文件名保持一致，方便在ROOT中调用
void count_range() {
    // ----------------- 配置区 -----------------
    // 定义文件路径模式（TChain::Add可以直接接受通配符！）
    const char* filePattern = "/home/evalie2/Project/document/273Ds/inter2/SS03200*_chain.root";

    // 定义 TTree 名称
    const char* treeName = "tr_chain";
    
    // 定义选择条件（等同于Draw命令中的cut）
    const char* cutString = "SSD_E[1]==0 && Delta_Ts[1]>0 && DSSD_E[1] > 5000 && DSSD_E[1] < 10000 && log10(Delta_Ts[1]/1e9) > 1 && log10(Delta_Ts[1]/1e9) < 2";
    // ------------------------------------------

    std::cout << "============================================================" << std::endl;
    std::cout << "开始计数满足条件的事件..." << std::endl;
    std::cout << "文件模式: " << filePattern << std::endl;
    std::cout << "TTree 名称: " << treeName << std::endl;
    std::cout << "选择条件: " << cutString << std::endl;
    std::cout << "============================================================" << std::endl;

    TChain* chain = new TChain(treeName);
    
    // TChain::Add() 返回成功添加的文件数量
    Int_t nFiles = chain->Add(filePattern);

    if (nFiles == 0) {
        std::cerr << "错误：未找到或未能添加任何匹配的文件 (" << filePattern << ")。" << std::endl;
        delete chain;
        return;
    }
    
    std::cout << "成功加载 " << nFiles << " 个文件。" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    Long64_t totalCount = 0;
    
    // 遍历 TChain 内部的文件列表，获取每个文件的计数
    // 注意：TChain::GetNtrees() 返回的是包含的 TTree/TFile 数量
    for (int i = 0; i < chain->GetNtrees(); ++i) {
        // 关键步骤 1: 确保 TChain 加载了第 i 个 TTree 对应的文件
        // LoadTree(i) 返回的是 TTree 在 chain 里的全局 Entry 序号，这里不重要
        chain->LoadTree(chain->GetTreeOffset()[i]);
        
        // 关键步骤 2: 获取当前文件的 TFile 指针
        TFile* currentFile = chain->GetCurrentFile(); // GetCurrentFile() 比 GetFile() 更常用
        
        if (!currentFile) {
            std::cerr << "警告：无法获取第 " << i + 1 << " 个文件的指针，跳过。" << std::endl;
            continue;
        }

        // 获取文件名（不含路径）
        TString fileName = currentFile->GetName();
        TString baseName = gSystem->BaseName(fileName);
        
        // 关键步骤 3: 获取当前 TTree 的指针
        TTree* currentTree = chain->GetTree(); // GetTree() 返回当前加载的 TTree 指针
        
        if (!currentTree) {
            std::cerr << "警告：无法获取文件 " << baseName.Data() << " 中的 TTree，跳过。" << std::endl;
            continue;
        }

        // 使用 GetEntries(cut) 方法获取满足条件的事件数
        // 注意：这里是对当前 TTree 调用的 GetEntries()
        Long64_t count = currentTree->GetEntries(cutString);

        // 打印结果
        std::cout << "文件: " << std::left << std::setw(30) << baseName.Data() 
                  << " | 计数: " << std::right << std::setw(15) << count << std::endl;

        totalCount += count;
    }

    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "总计满足条件的事件数: " << totalCount << std::endl;
    std::cout << "============================================================" << std::endl;

    delete chain; // 清理
}