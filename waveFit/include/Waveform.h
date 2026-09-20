#pragma once
#include <string>
#include <vector>

struct EventInfo {
    int chainNum = -1, run = -1;
    long long mapEntry = -1;
};
struct Waveform : EventInfo {
    std::string canvasName, canvasTitle, canvasPath;
    std::string padName, padPath, objectName, objectTitle, objectClass;
    std::string xTitle, yTitle;
    int objectIndex = -1;
    bool isGraph = false;
    std::vector<double> x, y;
};
struct CanvasInfo : EventInfo {
    std::string name, title, path;
    std::vector<Waveform> waves;
};
struct ReadOptions {
    bool inspect = false, verbose = false, includeFrames = false;
};
std::vector<CanvasInfo> ReadCanvases(const std::string&, const ReadOptions&);
