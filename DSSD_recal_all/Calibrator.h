#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include <vector>
#include <string>
#include <fstream>
#include "TH1D.h"
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"
#include "Config.h"   // PeakDef / target_E / fit_win / 道号常量

struct NormParams { double k=1.0, b=0.0; };

struct PlaneResult { 
    double K=1.0; 
    double B=0.0; 
    double FWHM_Max=0.0; 
    double Res_RMS=0.0; 
    bool   valid=false;      // K,B 是否由本平面自己的参考峰拟合得到
    bool   inherited=false;  // 是否沿用别的平面的刻度（本平面参考峰不足/被否决）
    int    n_points=0;       // 实际参与拟合的参考峰数
    double span=0.0;         // 参与拟合的参考峰能量跨度 (keV)
    double fit_rms=0.0;      // 拟合残差 RMS (keV)
};

struct PeakInfo {
    double pos;    
    double height; 
    double sigma;   // 高斯宽度 (keV)；<0 表示二次拟合失败
};

class Calibrator { 
public:
    Calibrator();
    ~Calibrator();
    void Run(); 

private:
    std::vector<NormParams> norm_params_; 
    std::vector<double> strip_fwhms_; 
    std::vector<double> strip_dE_;      // 每条在 strip_peak_E_ 处的峰位偏差 (keV)，0 = 未测到

    // 参考峰表（运行期从 REF_PEAK_FILE 读，读不到时用 Config.h 的内置默认值）
    std::vector<PeakDef> ref_peaks_;
    std::vector<PeakDef> diag_peaks_;
    double strip_peak_E_   = target_E;   // 逐条一致性检查的峰
    double strip_peak_win_ = fit_win;
    
    TH1D *h_tot_X, *h_tot_Y, *h_tot_YH;
    TH1D *h_fwhm_summary; 

    void LoadNormParams(); 
    void LoadPeakTable(); 
    void FillSpectra(ROOT::RDataFrame& df); 
    
    // [修改] 增加 TFile* diag 参数
    void CalculateStripFWHM(ROOT::RDataFrame& df, const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH, TFile* diag);

    // fallback 非空时：本平面参考峰不足/拟合被否决 -> 沿用 fallback 的刻度
    PlaneResult CalibratePlane(TH1D* h_raw, const char* name, TFile* diag, std::ofstream& out,
                               const PlaneResult* fallback = nullptr); 
    
    void WriteOutput(const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH); 

    void FindAndListPeaks(TH1D* h, const char* title, std::ofstream& out);

    // 工具
    bool FindPeakGaussian(TH1D* h, double center, double win, double &pos, double &sigma);
    TH1D* SubtractBackground(TH1D* h);
    TH1D* ApplyCalibration(TH1D* h_in, double K, double B, const char* new_name);
};
#endif

