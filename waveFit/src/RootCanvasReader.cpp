#include "Waveform.h"
#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TKey.h>
#include <TList.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <iostream>
#include <memory>
#include <regex>
#include <set>
#include <stdexcept>

namespace {
long long Parse(const std::string& s, const char* pattern) {
    std::smatch m;
    if (!std::regex_search(s, m, std::regex(pattern))) return -1;
    try { return std::stoll(m[1]); } catch (...) { return -1; }
}
int ParseInt(const std::string& s, const char* pattern) {
    const auto v = Parse(s, pattern);
    return v > 2147483647LL ? -1 : static_cast<int>(v);
}
void ScanPad(TVirtualPad* pad, CanvasInfo& c, const ReadOptions& opt,
             const std::string& path, int depth) {
    TIter next(pad->GetListOfPrimitives());
    while (auto* obj = next()) {
        auto* graph = dynamic_cast<TGraph*>(obj);
        auto* hist = dynamic_cast<TH1*>(obj);
        // ROOT DrawFrame creates an empty histogram called hframe. Retain all
        // other TH1 candidates, even when a graph is present in the same pad.
        bool frame = hist && hist->GetDimension() == 1 && std::string(hist->GetName()) == "hframe" &&
            hist->GetEntries() == 0;
        if (frame) {
            // GetMinimum/GetMaximum may return display limits set by DrawFrame,
            // so inspect bin contents instead of the axis range.
            for (int i = 1; i <= hist->GetNbinsX(); ++i)
                if (hist->GetBinContent(i) != 0) { frame = false; break; }
        }
        if (opt.inspect || opt.verbose) {
            std::cout << std::string(depth * 2, ' ') << '[' << obj->ClassName() << "] name="
                      << obj->GetName() << " title=" << obj->GetTitle() << " pad=" << path;
            if (graph) std::cout << " N=" << graph->GetN();
            if (hist) std::cout << " NbinsX=" << hist->GetNbinsX() << " dimension=" << hist->GetDimension();
            if (frame && !opt.includeFrames) std::cout << " [axis frame, skipped]";
            std::cout << '\n';
        }
        if (auto* child = dynamic_cast<TVirtualPad*>(obj)) {
            ScanPad(child, c, opt, path + "/" + child->GetName(), depth + 1);
            continue;
        }
        if (!graph && (!hist || hist->GetDimension() != 1 || (frame && !opt.includeFrames))) continue;
        Waveform w;
        static_cast<EventInfo&>(w) = static_cast<const EventInfo&>(c);
        w.canvasName = c.name; w.canvasTitle = c.title; w.canvasPath = c.path;
        w.padName = pad->GetName(); w.padPath = path;
        w.objectName = obj->GetName(); w.objectTitle = obj->GetTitle();
        w.objectClass = obj->ClassName(); w.isGraph = graph != nullptr; w.objectIndex = static_cast<int>(c.waves.size());
        if (graph) {
            for (int i = 0; i < graph->GetN(); ++i) {
                double x, y; graph->GetPoint(i, x, y); w.x.push_back(x); w.y.push_back(y);
            }
            if (graph->GetXaxis()) w.xTitle = graph->GetXaxis()->GetTitle();
            if (graph->GetYaxis()) w.yTitle = graph->GetYaxis()->GetTitle();
        } else {
            for (int i = 1; i <= hist->GetNbinsX(); ++i) {
                w.x.push_back(hist->GetBinCenter(i)); w.y.push_back(hist->GetBinContent(i));
            }
            w.xTitle = hist->GetXaxis()->GetTitle(); w.yTitle = hist->GetYaxis()->GetTitle();
        }
        c.waves.push_back(std::move(w));
    }
}
void ScanDirectory(TDirectory* dir, const std::string& path, const ReadOptions& opt,
                   std::vector<CanvasInfo>& out) {
    // Only the newest key cycle represents the current object.
    std::set<std::string> seen;
    TIter next(dir->GetListOfKeys());
    while (auto* key = dynamic_cast<TKey*>(next())) {
        if (!seen.insert(key->GetName()).second) continue;
        auto* cls = TClass::GetClass(key->GetClassName());
        const std::string full = path.empty() ? key->GetName() : path + "/" + key->GetName();
        if (opt.inspect) std::cout << "[Key] " << full << ';' << key->GetCycle() << " class=" << key->GetClassName() << '\n';
        if (!cls) continue;
        if (cls->InheritsFrom(TDirectory::Class())) {
            auto* sub = dir->GetDirectory(key->GetName());
            if (sub) ScanDirectory(sub, full, opt, out);
        } else if (cls->InheritsFrom(TCanvas::Class())) {
            std::unique_ptr<TObject> object(key->ReadObj());
            auto* canvas = dynamic_cast<TCanvas*>(object.get());
            if (!canvas) continue;
            CanvasInfo c; c.name = canvas->GetName(); c.title = canvas->GetTitle(); c.path = full;
            c.chainNum = ParseInt(c.name, "ChainNum_([0-9]+)");
            if (c.chainNum < 0) c.chainNum = ParseInt(c.title, "ChainNum_([0-9]+)");
            c.run = ParseInt(c.title, "\\brun\\s*=\\s*([0-9]+)");
            c.mapEntry = Parse(c.title, "\\bmap_entry\\s*=\\s*([0-9]+)");
            if (opt.inspect || opt.verbose)
                std::cout << "[Canvas] " << c.path << "\nTitle: " << c.title << "\nchainNum="
                          << c.chainNum << " run=" << c.run << " mapEntry=" << c.mapEntry << '\n';
            ScanPad(canvas, c, opt, c.path, 1);
            if (c.waves.empty()) std::cerr << "[WARNING] " << c.path << " contains no TGraph or TH1 waveform.\n";
            out.push_back(std::move(c));
        }
    }
}
}
std::vector<CanvasInfo> ReadCanvases(const std::string& name, const ReadOptions& opt) {
    std::unique_ptr<TFile> input(TFile::Open(name.c_str(), "READ"));
    if (!input || input->IsZombie()) throw std::runtime_error("Cannot open ROOT file: " + name);
    std::vector<CanvasInfo> result;
    ScanDirectory(input.get(), "", opt, result);
    return result;
}
