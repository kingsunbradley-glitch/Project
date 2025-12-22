#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include <vector>
#include <fstream>

// ROOT Classics
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TLine.h"
#include "TChain.h"  // <--- 解决 'TChain' was not declared

// ROOT RDataFrame
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RResultPtr.hxx"

class Calibrator {
public:
    Calibrator();
    ~Calibrator();

    void Run(); // <--- 解决 Run() 找不到匹配函数的问题 (现在不需要参数)

private:
    // 数据容器：主线程合并后的直方图
    std::vector<TH1F*> vec_h_raw;      
    std::vector<TH1F*> vec_h_cal;      
    TH1F* h_resolution;                

    // RDataFrame 填充逻辑
    void FillHistogramsRDF(TChain* chain);

    // 刻度处理逻辑
    void ProcessChannel(int chID, std::ofstream &ofs);
    
    // 算法辅助
    void InitialiseHistograms();
    TH1F* BackgroundSubtract(TH1F* h_in);
    bool FindSinglePeak(TH1F* h, double center, double window, double &mean, double &sigma);

    // 输出文件指针
    TFile* f_out_root;
};

#endif