#ifndef CALIBRATOR_H
#define CALIBRATOR_H

#include "Config.h"

#include <ROOT/RDataFrame.hxx>
#include <TH1D.h>
#include <TFile.h>

#include <string>
#include <vector>

struct NormParam {
    double k = 1.0;
    double b = 0.0;
    int ok = 0;
    long n = 0;
};

struct PlaneCalib {
    bool ok = false;
    double K = FAIL_K; // Energy = K*ADC_norm + B
    double B = FAIL_B;
    double chi2ndf = FAIL_CHI2;
    double adc[3] = {0,0,0};
    double sig[3] = {0,0,0};
};

struct ChResult {
    double k = FAIL_K;
    double b = FAIL_B;
    double fwhm = FAIL_FWHM;
    double chi2ndf = FAIL_CHI2;
};

class Calibrator {
public:
    Calibrator(std::string inRoot,
               std::string normTxt,
               std::string outDat,
               std::string diagMain,
               std::string diagEd);

    int Run();

private:
    std::string inRoot_;
    std::string normTxt_;
    std::string outDat_;
    std::string diagMain_;
    std::string diagEd_;

    std::vector<NormParam> norm_;
    std::vector<ChResult> chRes_;

    PlaneCalib calX_, calY_, calYH_;

    bool LoadNorm();

    void BuildPlaneSummed(ROOT::RDataFrame& df, TH1D& hX, TH1D& hY, TH1D& hYH);
    PlaneCalib CalibratePlane(TH1D& h, const char* name, TFile& fout);

    static bool FindTop3PeaksTSpectrum(TH1D& h, std::vector<double>& peaks);
    static bool GaussianRefine(TH1D& h, double seed, double& mean, double& sigma, double& chi2ndf);

    void BuildFinalKB();
    void ComputePerChannelFWHM(ROOT::RDataFrame& df, TFile& outEd);
    bool FitFWHMOne(TH1D& h, double& fwhm, double& chi2ndf);

    void WriteEnerCalDat() const;
};

#endif
