#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include <vector>
#include <string>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"

struct NormParams { double k=1.0, b=0.0; };

// [修改] 增加 Res_RMS (Residual RMS, 刻度残差均方根)
struct PlaneResult { 
    double K=1.0; 
    double B=0.0; 
    double FWHM_Max=0.0; 
    double Res_RMS=0.0; // 新增：刻度后峰位与诊断峰的偏差
};

struct PeakInfo {
    double pos;    
    double height; 
};

class Calibrator {
public:
    Calibrator();
    ~Calibrator();
    void Run(); 

private:
    std::vector<NormParams> norm_params_; 
    TH1D *h_tot_X, *h_tot_Y, *h_tot_YH;
    TH1D *h_fwhm_summary; // 这就是逐条 FWHM 分布图

    void LoadNormParams(); 
    void FillSpectra(ROOT::RDataFrame& df); 
    
    PlaneResult CalibratePlane(TH1D* h_raw, const char* name, TFile* diag, std::ofstream& out); 
    
    void WriteOutput(const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH); 

    void FindAndListPeaks(TH1D* h, const char* title, std::ofstream& out);

    // 工具
    bool FindPeakGaussian(TH1D* h, double center, double win, double &pos, double &sigma);
    TH1D* SubtractBackground(TH1D* h);
    TH1D* ApplyCalibration(TH1D* h_in, double K, double B, const char* new_name);
};
#endif