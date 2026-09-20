#pragma once
#include "Analysis.h"
#include <filesystem>
#include <memory>
class TFile;
class TTree;
std::string WaveID(const Waveform&);
void ExportWave(const Waveform&, const std::filesystem::path&);
void DrawWave(const Waveform&, const std::filesystem::path&, const std::string&);
struct AnalysisRow {
    Waveform wave;
    Processed processed;
    FitResult fit;
};
void SaveAnalysis(const std::vector<AnalysisRow>&, const std::filesystem::path&, const std::string&, const std::string&);
