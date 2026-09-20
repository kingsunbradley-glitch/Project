#include "Analysis.h"
#include "Output.h"
#include "Waveform.h"
#include <TROOT.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <optional>
#include <regex>
#include <set>
#include <stdexcept>

namespace fs=std::filesystem;
namespace {
struct Options {
    std::string input, output, templateFile, pattern, unit="original_x";
    ReadOptions read;
    AnalysisOptions analysis;
    bool list=false, drawAll=false, analyzeAll=false, exportWaves=false;
    bool makeTemplate=false, overwrite=false;
    std::optional<int> draw, analyze;
    std::vector<int> templateChains;
    int graphIndex=-1;
    double samplePeriod=NaN;
};
void Help() {
    std::cout << R"HELP(ROOT Canvas waveform reader / pile-up analysis
Usage: ./wave_reader FILE.root [options]
  --inspect                  Print all keys and recursively inspect primitives
  --list                     List canvases and waveform candidates
  --export-waveforms         Export selected raw (x,y) points as CSV
  --draw N | --draw-all      Redraw extracted data to PNG, PDF and ROOT
  --analyze N | --analyze-all Baseline, polarity and rising-edge candidates
  --fit N                    Analyze N with a required --template FILE
  --make-template [N ...]    Average selected clean single-pulse waveforms
  --template FILE.root      Also run single/double template fits
  --object-pattern REGEX    Match object name/title or pad name/path
  --graph-index N           Select graph index (zero-based, per canvas)
  --include-frames          Include empty ROOT hframe histograms as candidates
  --baseline-samples N      Default min(200,max(2,Nsamples/10))
  --robust-baseline         Median + MAD instead of mean + sample deviation
  --smooth N                Odd moving-average window (default 1, detection only)
  --threshold K             Derivative-noise threshold (default 5)
  --min-separation N        Minimum edge separation in samples (default 20)
  --sample-period-ns VALUE  Convert sample-index x into ns: x_ns=x*VALUE
  --x-unit LABEL            Label unscaled native x (default original_x)
  --max-delta VALUE         Scan limit in current x units (default 2000)
  --t1 VALUE                Fix first edge time in current x units
  --output-dir DIR          Default output/<input_stem>
  --overwrite               Replace this program's existing output files
  --verbose                 Detailed reader diagnostics
  --help                    Show this help

No action defaults to --inspect. Input ROOT files are always opened READ.
Fits use positive amplitudes and a constant/linear residual baseline. The 1D
Delta_t scan fixes t1 to the first candidate (or --t1); model comparisons are
conditional on that reference and template quality. ADC heights are not energies.
)HELP";
}
double Number(const std::string& s) {
    size_t end=0; double v=std::stod(s,&end);
    if (end!=s.size() || !std::isfinite(v)) throw std::runtime_error("Invalid number: "+s);
    return v;
}
int Integer(const std::string& s,int min) {
    double v=Number(s);
    if (v<min || v>2147483647 || std::floor(v)!=v) throw std::runtime_error("Invalid integer: "+s);
    return static_cast<int>(v);
}
Options Parse(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;++i) {
        const std::string a=argv[i];
        auto value=[&]() { if (++i>=argc) throw std::runtime_error("Missing value for "+a);return std::string(argv[i]); };
        if (a=="--inspect") o.read.inspect=true;
        else if (a=="--list") o.list=true;
        else if (a=="--verbose") o.read.verbose=true;
        else if (a=="--include-frames") o.read.includeFrames=true;
        else if (a=="--draw") o.draw=Integer(value(),0);
        else if (a=="--draw-all") o.drawAll=true;
        else if (a=="--analyze" || a=="--fit") {o.analyze=Integer(value(),0);}
        else if (a=="--analyze-all") o.analyzeAll=true;
        else if (a=="--export-waveforms") o.exportWaves=true;
        else if (a=="--make-template") {
            o.makeTemplate=true;
            while(i+1<argc && std::string(argv[i+1]).rfind("--",0)!=0) o.templateChains.push_back(Integer(argv[++i],0));
        }
        else if (a=="--template") o.templateFile=value();
        else if (a=="--object-pattern") o.pattern=value();
        else if (a=="--graph-index") o.graphIndex=Integer(value(),0);
        else if (a=="--baseline-samples") o.analysis.baselineSamples=Integer(value(),2);
        else if (a=="--robust-baseline") o.analysis.robust=true;
        else if (a=="--smooth") o.analysis.smooth=Integer(value(),1);
        else if (a=="--threshold") o.analysis.threshold=Number(value());
        else if (a=="--min-separation") o.analysis.minSeparation=Integer(value(),1);
        else if (a=="--max-delta") o.analysis.maxDelta=Number(value());
        else if (a=="--t1") o.analysis.fixedT1=Number(value());
        else if (a=="--sample-period-ns") o.samplePeriod=Number(value());
        else if (a=="--x-unit") o.unit=value();
        else if (a=="--output-dir") o.output=value();
        else if (a=="--overwrite") o.overwrite=true;
        else if (a.rfind("-",0)==0) throw std::runtime_error("Unknown option: "+a);
        else if (o.input.empty()) o.input=a;
        else throw std::runtime_error("Unexpected argument: "+a);
    }
    if(o.input.empty()) throw std::runtime_error("Input ROOT file required; see --help");
    if(o.analysis.smooth%2==0 || o.analysis.threshold<=0 || o.analysis.maxDelta<0 ||
       (std::isfinite(o.samplePeriod) && o.samplePeriod<=0)) throw std::runtime_error("Invalid smoothing, threshold, scan limit or sample period");
    if(std::isfinite(o.samplePeriod)) o.unit="ns";
    if(o.unit.empty()) throw std::runtime_error("Empty x unit");
    if(o.output.empty()) o.output=(fs::path("output")/fs::path(o.input).stem()).string();
    bool fitRequested=false;
    for(int i=1;i<argc;++i) if(std::string(argv[i])=="--fit") fitRequested=true;
    if(fitRequested && o.templateFile.empty() && !o.makeTemplate) throw std::runtime_error("--fit requires --template FILE or --make-template");
    if(!o.read.inspect && !o.list && !o.draw && !o.drawAll && !o.analyze && !o.analyzeAll && !o.exportWaves && !o.makeTemplate) o.read.inspect=true;
    if(o.makeTemplate && !o.templateFile.empty()) throw std::runtime_error("Use either --make-template or --template, not both");
    return o;
}
void CheckOutput(const fs::path& path,const Options& o) {
    const auto dest=fs::weakly_canonical(path);
    if(dest==fs::weakly_canonical(o.input) || (!o.templateFile.empty() && dest==fs::weakly_canonical(o.templateFile)))
        throw std::runtime_error("Output would overwrite an input file: "+path.string());
    if(fs::exists(path) && !o.overwrite) throw std::runtime_error("Output exists: "+path.string()+"; choose a new --output-dir or use --overwrite");
}
}
int main(int argc,char** argv) {
    for(int i=1;i<argc;++i) if(std::string(argv[i])=="--help") {Help();return 0;}
    gROOT->SetBatch(true);
    try {
        const Options o=Parse(argc,argv);
        auto canvases=ReadCanvases(o.input,o.read);
        if(canvases.empty()) throw std::runtime_error("No TCanvas objects found");
        std::optional<std::regex> pattern;
        if(!o.pattern.empty()) pattern.emplace(o.pattern);
        std::vector<const Waveform*> selected;
        std::set<int> chains;
        for(auto& c:canvases) {
            chains.insert(c.chainNum);
            std::cout << "[Canvas] " << c.path << " title=" << c.title << " candidates=" << c.waves.size() << '\n';
            int graphIndex=0;
            for(auto& w:c.waves) {
                bool isGraph=w.isGraph;
                bool keep=o.graphIndex<0 || (isGraph && graphIndex==o.graphIndex);
                if(isGraph) ++graphIndex;
                if(pattern) keep=keep && (std::regex_search(w.objectName,*pattern) || std::regex_search(w.objectTitle,*pattern) ||
                                          std::regex_search(w.padName,*pattern) || std::regex_search(w.padPath,*pattern));
                std::cout << "  [" << w.objectIndex << "] " << w.objectClass << ' ' << w.objectName << " pad=" << w.padPath
                          << " N=" << w.x.size() << (keep ? " [selected]" : "") << '\n';
                if(std::isfinite(o.samplePeriod)) for(double& x:w.x) x*=o.samplePeriod;
                if(keep) selected.push_back(&w);
            }
        }
        auto requireChain=[&](int chain) {if(!chains.count(chain)) throw std::runtime_error("ChainNum_"+std::to_string(chain)+" not found");};
        if(o.draw) requireChain(*o.draw);
        if(o.analyze) requireChain(*o.analyze);
        for(int n:o.templateChains) requireChain(n);
        std::vector<const Waveform*> drawing,analyzing,templating;
        for(const auto* w:selected) {
            if(o.drawAll || (o.draw && w->chainNum==*o.draw)) drawing.push_back(w);
            if(o.analyzeAll || (o.analyze && w->chainNum==*o.analyze)) analyzing.push_back(w);
            if(o.makeTemplate && (o.templateChains.empty() || std::find(o.templateChains.begin(),o.templateChains.end(),w->chainNum)!=o.templateChains.end())) templating.push_back(w);
        }
        if((o.draw || o.drawAll) && drawing.empty()) throw std::runtime_error("No waveform matches drawing selection");
        if((o.analyze || o.analyzeAll) && analyzing.empty()) throw std::runtime_error("No waveform matches analysis selection");
        if(o.exportWaves && selected.empty()) throw std::runtime_error("No waveform matches export selection");
        if(o.makeTemplate && templating.empty()) throw std::runtime_error("No waveform matches template selection");
        if(o.makeTemplate && o.pattern.empty() && o.graphIndex<0) {
            for(const auto& c:canvases) if(c.waves.size()>1 &&
               std::any_of(templating.begin(),templating.end(),[&](const Waveform* w){return w->canvasPath==c.path;}))
                throw std::runtime_error("Multiple template candidates: select one detector using --object-pattern or --graph-index");
        }
        const fs::path out=o.output;
        const std::string stem=fs::path(o.input).stem().string();
        std::set<std::string> identities;
        for(const auto* w:selected) if(!identities.insert(WaveID(*w)).second) throw std::runtime_error("Colliding output identities: "+WaveID(*w));
        // Check every intended output before starting any writes.
        if(o.exportWaves) for(const auto* w:selected) CheckOutput(out/"waveforms"/(WaveID(*w)+".csv"),o);
        for(const auto* w:drawing) for(const auto& ext:{".png",".pdf",".root"}) CheckOutput(out/"drawings"/(WaveID(*w)+ext),o);
        if(o.makeTemplate) {CheckOutput(out/"pulse_template.root",o);CheckOutput(out/"pulse_template.csv",o);}
        if(!analyzing.empty()) {CheckOutput(out/(stem+"_analysis.root"),o);CheckOutput(out/(stem+"_analysis.csv"),o);}
        std::optional<PulseTemplate> pulseTemplate;
        if(!o.templateFile.empty()) pulseTemplate=LoadTemplate(o.templateFile,o.unit);
        if(o.makeTemplate) pulseTemplate=MakeTemplate(templating,o.analysis,o.unit);
        if(o.exportWaves) for(const auto* w:selected) ExportWave(*w,out/"waveforms");
        for(const auto* w:drawing) DrawWave(*w,out/"drawings",o.unit);
        if(o.makeTemplate) {fs::create_directories(out);SaveTemplate(*pulseTemplate,(out/"pulse_template.root").string());}
        std::vector<AnalysisRow> rows;
        bool failed=false;
        for(const auto* w:analyzing) {
            AnalysisRow row;row.wave=*w;
            try {
                row.processed=Process(*w,o.analysis);
                const auto& p=row.processed;
                std::cout << "[Analysis] " << WaveID(*w) << " Baseline=" << p.baseline << " Noise RMS=" << p.noiseRMS
                          << " polarity=" << p.polarity << " candidates=" << p.pulses.size() << '\n';
                for(size_t i=0;i<p.pulses.size();++i) {
                    const auto& pulse=p.pulses[i];
                    std::cout << "  Pulse candidate " << i+1 << ": sample=" << pulse.sample << " x=" << pulse.x
                              << ' ' << o.unit << " height=" << pulse.height << " ADC";
                    if(i) std::cout << " Delta_t=" << pulse.x-p.pulses[i-1].x << ' ' << o.unit;
                    std::cout << '\n';
                }
                if(pulseTemplate) {
                    row.fit=FitPileup(*w,p,*pulseTemplate,o.analysis);
                    const auto& f=row.fit;
                    std::cout << "  Fit: " << f.status << " deltaT=" << f.deltaT << ' ' << o.unit
                              << " chi2_1=" << f.chi1 << " chi2_2=" << f.chi2 << " reduced_chi2_2=" << f.chi2/f.ndf2 << " Delta BIC=" << f.bic1-f.bic2 << '\n';
                    if (f.valid && f.chi2/f.ndf2>5) std::cerr << "[WARNING] Large residual chi2: inspect template agreement and residuals before interpreting pile-up evidence.\n";
                }
            } catch(const std::exception& e) {
                failed=true;row.fit.status="analysis_error";
                std::cerr << "[ERROR] " << WaveID(*w) << ": " << e.what() << '\n';
            }
            rows.push_back(std::move(row));
        }
        if(!rows.empty()) SaveAnalysis(rows,out,stem,o.unit);
        std::cout << "[INFO] Canvases=" << canvases.size() << " selected_waveforms=" << selected.size()
                  << " analyzed=" << rows.size() << '\n';
        if(o.exportWaves || !drawing.empty() || !rows.empty() || o.makeTemplate) std::cout << "[INFO] Output: " << fs::absolute(out) << '\n';
        return failed ? 2 : 0;
    } catch(const std::exception& e) {std::cerr << "[ERROR] " << e.what() << '\n';return 1;}
}
