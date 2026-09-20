#pragma once
#include "Waveform.h"
#include <limits>

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
struct AnalysisOptions {
    bool robust = false;
    int baselineSamples = 0, smooth = 1, minSeparation = 20;
    double threshold = 5, maxDelta = 2000, fixedT1 = NaN;
};
struct Pulse {
    int sample = -1;
    double x = NaN, height = NaN;
};
struct Processed {
    double baseline = NaN, noiseRMS = NaN, derivativeNoise = NaN;
    int polarity = 1, nBaseline = 0;
    std::vector<double> corrected, smoothed;
    std::vector<Pulse> pulses;
};
struct PulseTemplate {
    std::vector<double> x, y;
    std::string unit, source;
};
struct FitResult {
    bool valid = false;
    double t1 = NaN, t2 = NaN, deltaT = NaN, A1 = NaN, A2 = NaN;
    double chi1 = NaN, chi2 = NaN, aic1 = NaN, aic2 = NaN, bic1 = NaN, bic2 = NaN;
    int ndf1 = 0, ndf2 = 0;
    double fitSigma = NaN;
    std::string status = "not_requested";
    std::vector<double> single, dual, residual1, residual2, scanX, scanY, matchedPower;
};
Processed Process(const Waveform&, const AnalysisOptions&);
PulseTemplate MakeTemplate(const std::vector<const Waveform*>&, const AnalysisOptions&, const std::string&);
PulseTemplate LoadTemplate(const std::string&, const std::string&);
void SaveTemplate(const PulseTemplate&, const std::string&);
FitResult FitPileup(const Waveform&, const Processed&, const PulseTemplate&, const AnalysisOptions&);
