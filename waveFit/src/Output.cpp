#include "Output.h"
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TNamed.h>
#include <TTree.h>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace {
std::string Safe(std::string s) {
    for (auto& c : s) if (!std::isalnum(static_cast<unsigned char>(c)) && c != '-' && c != '_') c = '_';
    return s.empty() ? "unnamed" : s;
}
std::string Quote(const std::string& s) {
    std::string v = "\"";
    for (char c : s) { if (c == '"') v += '"'; v += c; }
    return v + '"';
}
void WriteGraph(const char* name, const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size() || x.empty()) return;
    TGraph g(static_cast<int>(x.size()), x.data(), y.data());
    g.SetName(name); g.Write();
}
}
std::string WaveID(const Waveform& w) {
    return Safe(w.canvasPath) + "_" + Safe(w.objectName) + "_w" + std::to_string(w.objectIndex);
}
void ExportWave(const Waveform& w, const std::filesystem::path& dir) {
    std::filesystem::create_directories(dir);
    std::ofstream out(dir / (WaveID(w) + ".csv"));
    if (!out) throw std::runtime_error("Cannot write waveform CSV");
    out << std::setprecision(17) << "sample,x,y\n";
    for (size_t i = 0; i < w.x.size(); ++i) out << i << ',' << w.x[i] << ',' << w.y[i] << '\n';
    if (!out) throw std::runtime_error("Error writing waveform CSV");
}
void DrawWave(const Waveform& w, const std::filesystem::path& dir, const std::string& unit) {
    if (w.x.empty()) return;
    std::filesystem::create_directories(dir);
    TCanvas c("extracted_canvas", "Extracted waveform", 1200, 600);
    TGraph g(static_cast<int>(w.x.size()), w.x.data(), w.y.data());
    g.SetName("extracted_waveform");
    g.SetTitle((w.canvasName + " / " + w.objectName + ";x (" + unit + ");ADC amplitude").c_str());
    g.SetLineWidth(2); g.Draw("AL"); c.SetGrid(); c.Update();
    c.SaveAs((dir / (WaveID(w) + ".png")).c_str());
    c.SaveAs((dir / (WaveID(w) + ".pdf")).c_str());
    TFile output((dir / (WaveID(w) + ".root")).c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot write drawing ROOT file");
    g.Write(); c.Write("redrawn_canvas");
}
void SaveAnalysis(const std::vector<AnalysisRow>& rows, const std::filesystem::path& dir,
                  const std::string& stem, const std::string& unit) {
    std::filesystem::create_directories(dir);
    TFile output((dir / (stem + "_analysis.root")).c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Cannot create analysis ROOT file");
    TNamed("x_unit", unit.c_str()).Write();
    TTree tree("WaveAnalysis", "Extracted waveform and pulse analysis");
    int chainNum, run, objectIndex, nSamples, nPulseCandidate, polarity, nBaseline, ndf1, ndf2;
    long long mapEntry;
    std::string canvasName, canvasTitle, canvasPath, padName, padPath, objectName, objectTitle, objectClass, fitStatus;
    double baseline, noiseRMS, derivativeNoise, t1, t2, deltaT, A1, A2, chi2_1pulse, chi2_2pulse;
    double deltaChi2, AIC1, AIC2, BIC1, BIC2, deltaBIC, fitSigma;
    std::vector<double> x, y, pulseX, pulseHeight;
    std::vector<int> pulseSample;
#define BRANCH(v) tree.Branch(#v, &v)
    BRANCH(chainNum); BRANCH(run); BRANCH(mapEntry); BRANCH(objectIndex); BRANCH(nSamples);
    BRANCH(canvasName); BRANCH(canvasTitle); BRANCH(canvasPath); BRANCH(padName); BRANCH(padPath);
    BRANCH(objectName); BRANCH(objectTitle); BRANCH(objectClass); BRANCH(fitStatus);
    BRANCH(baseline); BRANCH(noiseRMS); BRANCH(derivativeNoise); BRANCH(polarity); BRANCH(nBaseline);
    BRANCH(nPulseCandidate); BRANCH(t1); BRANCH(t2); BRANCH(deltaT); BRANCH(A1); BRANCH(A2);
    BRANCH(chi2_1pulse); BRANCH(chi2_2pulse); BRANCH(deltaChi2); BRANCH(AIC1); BRANCH(AIC2);
    BRANCH(BIC1); BRANCH(BIC2); BRANCH(deltaBIC); BRANCH(ndf1); BRANCH(ndf2); BRANCH(fitSigma);
    BRANCH(x); BRANCH(y); BRANCH(pulseX); BRANCH(pulseHeight); BRANCH(pulseSample);
#undef BRANCH
    std::ofstream csv(dir / (stem + "_analysis.csv"));
    if (!csv) throw std::runtime_error("Cannot create analysis CSV");
    csv << std::setprecision(17) << "chainNum,run,mapEntry,canvasName,padPath,objectName,objectIndex,nSamples,baseline,noiseRMS,polarity,nPulse,t1,t2,deltaT,A1,A2,chi2_1,chi2_2,deltaChi2,AIC1,AIC2,BIC1,BIC2,deltaBIC,ndf1,ndf2,fitSigma,fitStatus,xUnit\n";
    for (const auto& row : rows) {
        const auto& w = row.wave; const auto& p = row.processed; const auto& f = row.fit;
        chainNum = w.chainNum; run = w.run; mapEntry = w.mapEntry; objectIndex = w.objectIndex;
        canvasName = w.canvasName; canvasTitle = w.canvasTitle; canvasPath = w.canvasPath;
        padName = w.padName; padPath = w.padPath; objectName = w.objectName;
        objectTitle = w.objectTitle; objectClass = w.objectClass;
        nSamples = static_cast<int>(w.x.size()); x = w.x; y = w.y;
        baseline = p.baseline; noiseRMS = p.noiseRMS; derivativeNoise = p.derivativeNoise;
        polarity = p.polarity; nBaseline = p.nBaseline; nPulseCandidate = static_cast<int>(p.pulses.size());
        pulseSample.clear(); pulseX.clear(); pulseHeight.clear();
        for (const auto& pulse : p.pulses) { pulseSample.push_back(pulse.sample); pulseX.push_back(pulse.x); pulseHeight.push_back(pulse.height); }
        t1 = f.t1; t2 = f.t2; deltaT = f.deltaT; A1 = f.A1; A2 = f.A2;
        chi2_1pulse = f.chi1; chi2_2pulse = f.chi2; deltaChi2 = f.chi1 - f.chi2;
        AIC1 = f.aic1; AIC2 = f.aic2; BIC1 = f.bic1; BIC2 = f.bic2; deltaBIC = BIC1 - BIC2;
        ndf1 = f.ndf1; ndf2 = f.ndf2; fitSigma = f.fitSigma; fitStatus = f.status;
        tree.Fill();
        csv << chainNum << ',' << run << ',' << mapEntry << ',' << Quote(canvasName) << ',' << Quote(padPath)
            << ',' << Quote(objectName) << ',' << objectIndex << ',' << nSamples << ',' << baseline << ',' << noiseRMS
            << ',' << polarity << ',' << nPulseCandidate << ',' << t1 << ',' << t2 << ',' << deltaT << ',' << A1 << ',' << A2
            << ',' << chi2_1pulse << ',' << chi2_2pulse << ',' << deltaChi2 << ',' << AIC1 << ',' << AIC2 << ',' << BIC1
            << ',' << BIC2 << ',' << deltaBIC << ',' << ndf1 << ',' << ndf2 << ',' << fitSigma << ',' << Quote(fitStatus) << ',' << Quote(unit) << '\n';
        output.cd();
        auto* sub = output.mkdir(WaveID(w).c_str());
        if (!sub) throw std::runtime_error("Duplicate output waveform identity: " + WaveID(w));
        sub->cd();
        WriteGraph("raw_waveform", w.x, w.y);
        WriteGraph("baseline_subtracted", w.x, p.corrected);
        WriteGraph("best_single_fit", w.x, f.single); WriteGraph("best_double_fit", w.x, f.dual);
        WriteGraph("residual_single", w.x, f.residual1); WriteGraph("residual_double", w.x, f.residual2);
        WriteGraph("chi2_vs_deltaT", f.scanX, f.scanY);
        WriteGraph("matched_filter_power", f.scanX, f.matchedPower);
    }
    if (!csv) throw std::runtime_error("Error writing analysis CSV");
    output.cd(); tree.Write();
}
