#include "Analysis.h"
#include <TDecompQRH.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMatrixD.h>
#include <TNamed.h>
#include <TVectorD.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>

namespace {
double Median(std::vector<double> v) {
    if (v.empty()) return 0;
    const size_t mid = v.size()/2;
    std::nth_element(v.begin(), v.begin()+mid, v.end());
    const double hi = v[mid];
    return v.size()%2 ? hi : (hi + *std::max_element(v.begin(), v.begin()+mid))/2;
}
double RMS(const std::vector<double>& v, double mean) {
    if (v.size() < 2) return 0;
    double sum = 0; for (double a : v) sum += (a-mean)*(a-mean);
    return std::sqrt(sum/(v.size()-1));
}
double MAD(const std::vector<double>& v) {
    const double m = Median(v);
    std::vector<double> d; for (double a : v) d.push_back(std::abs(a-m));
    return 1.482602218505602 * Median(d);
}
void Validate(const Waveform& w) {
    if (w.x.size() != w.y.size() || w.x.size() < 3) throw std::runtime_error("Need at least 3 waveform points");
    for (size_t i = 0; i < w.x.size(); ++i) {
        if (!std::isfinite(w.x[i]) || !std::isfinite(w.y[i])) throw std::runtime_error("Non-finite waveform point");
        if (i && w.x[i] <= w.x[i-1]) throw std::runtime_error("Analysis requires strictly increasing x; raw export remains available");
    }
}
double Step(const std::vector<double>& x) {
    std::vector<double> d;
    for (size_t i = 1; i < x.size(); ++i) d.push_back(x[i]-x[i-1]);
    return Median(d);
}
double Interpolate(const std::vector<double>& x, const std::vector<double>& y, double t) {
    if (t < x.front()) return 0;
    // A step-shaped preamplifier waveform need not have returned to baseline.
    // Hold the measured tail rather than inventing an abrupt falling edge.
    if (t >= x.back()) return y.back();
    auto it = std::upper_bound(x.begin(), x.end(), t);
    size_t j = static_cast<size_t>(it-x.begin()), i = j-1;
    return y[i]+(y[j]-y[i])*(t-x[i])/(x[j]-x[i]);
}
struct LinearFit {
    bool valid = false;
    double sse = std::numeric_limits<double>::infinity();
    double a1 = 0, a2 = 0;
    std::vector<double> model;
};
// Column-scaled Householder QR. Amplitudes are nonnegative; the two baseline
// coefficients remain unconstrained. No normal-equation inverse is formed.
LinearFit Project(const std::vector<double>& x, const std::vector<double>& y,
                  const PulseTemplate& h, double t1, double t2, bool dual) {
    const size_t n = x.size();
    std::vector<std::vector<double>> columns(4, std::vector<double>(n));
    const double center = (x.front()+x.back())/2, span = x.back()-x.front();
    for (size_t i = 0; i < n; ++i) {
        columns[0][i] = Interpolate(h.x, h.y, x[i]-t1);
        columns[1][i] = dual ? Interpolate(h.x, h.y, x[i]-t2) : 0;
        columns[2][i] = 1; columns[3][i] = (x[i]-center)/span;
    }
    LinearFit best;
    // The unconstrained solution normally suffices; enumerate boundary faces
    // only if a pulse amplitude is negative or columns are rank-deficient.
    const int fullMask = dual && std::abs(t2-t1) > 1e-10*Step(h.x) ? 3 : 1;
    const std::vector<int> masks = fullMask == 3 ? std::vector<int>{3,1,2,0} : std::vector<int>{1,0};
    for (int mask : masks) {
        std::vector<int> ids;
        if (mask & 1) ids.push_back(0);
        if (mask & 2) ids.push_back(1);
        ids.push_back(2); ids.push_back(3);
        const int k = static_cast<int>(ids.size());
        TMatrixD matrix(static_cast<int>(n), k); TVectorD rhs(static_cast<int>(n));
        std::vector<double> norm(k, 0);
        bool ok = true;
        for (int j = 0; j < k; ++j) {
            for (double a : columns[ids[j]]) norm[j] += a*a;
            norm[j] = std::sqrt(norm[j]);
            if (norm[j] < 1e-12) { ok = false; break; }
            for (size_t i = 0; i < n; ++i) matrix(static_cast<int>(i), j) = columns[ids[j]][i]/norm[j];
        }
        if (!ok) continue;
        for (size_t i = 0; i < n; ++i) rhs(static_cast<int>(i)) = y[i];
        TDecompQRH qr(matrix, 1e-10);
        Bool_t solved = false;
        TVectorD beta = qr.Solve(rhs, solved);
        if (!solved) continue;
        double coef[4] = {0,0,0,0};
        for (int j = 0; j < k; ++j) coef[ids[j]] = beta(j)/norm[j];
        if (!std::isfinite(coef[0]) || !std::isfinite(coef[1]) || coef[0] < 0 || coef[1] < 0) continue;
        LinearFit f; f.valid = true; f.a1 = coef[0]; f.a2 = coef[1]; f.sse = 0; f.model.resize(n);
        for (size_t i = 0; i < n; ++i) {
            double val = 0; for (int j = 0; j < 4; ++j) val += coef[j]*columns[j][i];
            f.model[i] = val; f.sse += (y[i]-val)*(y[i]-val);
        }
        if (f.sse < best.sse) best = std::move(f);
        if (mask == fullMask) break;
    }
    return best;
}
// Bounded Brent minimization (parabolic steps safeguarded by golden sections).
double Brent(const std::function<double(double)>& f, double a, double b, double tol) {
    constexpr double golden = 0.3819660112501051;
    double x = a+golden*(b-a), w=x, v=x, fx=f(x), fw=fx, fv=fx, d=0, e=0;
    for (int iter = 0; iter < 100; ++iter) {
        const double mid=(a+b)/2, tol1=tol+1e-12*std::abs(x), tol2=2*tol1;
        if (std::abs(x-mid) <= tol2-(b-a)/2) break;
        if (std::abs(e) > tol1) {
            double r=(x-w)*(fx-fv), q=(x-v)*(fx-fw), p=(x-v)*q-(x-w)*r;
            q=2*(q-r); if (q > 0) p=-p; q=std::abs(q);
            double old=e; e=d;
            if (q > 0 && std::abs(p) < std::abs(q*old/2) && p > q*(a-x) && p < q*(b-x)) {
                d=p/q; double u=x+d;
                if (u-a < tol2 || b-u < tol2) d=std::copysign(tol1, mid-x);
            } else { e=x<mid ? b-x : a-x; d=golden*e; }
        } else { e=x<mid ? b-x : a-x; d=golden*e; }
        const double u=x+(std::abs(d)>=tol1 ? d : std::copysign(tol1,d));
        const double fu=f(u);
        if (fu <= fx) {
            if (u < x) b=x; else a=x;
            v=w; fv=fw; w=x; fw=fx; x=u; fx=fu;
        } else {
            if (u < x) a=u; else b=u;
            if (fu <= fw || w == x) { v=w; fv=fw; w=u; fw=fu; }
            else if (fu <= fv || v == x || v == w) { v=u; fv=fu; }
        }
    }
    return x;
}
}
Processed Process(const Waveform& w, const AnalysisOptions& opt) {
    Validate(w);
    Processed p;
    const int n=static_cast<int>(w.x.size());
    if (opt.smooth>n) throw std::runtime_error("Smoothing window exceeds waveform length");
    p.nBaseline=opt.baselineSamples ? std::min(opt.baselineSamples,n) : std::min(200,std::max(2,n/10));
    const std::vector<double> base(w.y.begin(),w.y.begin()+p.nBaseline);
    p.baseline=opt.robust ? Median(base) : std::accumulate(base.begin(),base.end(),0.0)/base.size();
    p.noiseRMS=opt.robust ? MAD(base) : RMS(base,p.baseline);
    for (double y : w.y) p.corrected.push_back(y-p.baseline);
    const auto mm=std::minmax_element(p.corrected.begin(),p.corrected.end());
    p.polarity=std::abs(*mm.first)>*mm.second ? -1 : 1;
    for (double& y : p.corrected) y*=p.polarity;
    p.smoothed.resize(n);
    std::vector<double> prefix(n+1,0);
    for (int i=0;i<n;++i) prefix[i+1]=prefix[i]+p.corrected[i];
    const int half=opt.smooth/2;
    for (int i=0;i<n;++i) {
        const int lo=std::max(0,i-half),hi=std::min(n,i+half+1);
        p.smoothed[i]=(prefix[hi]-prefix[lo])/(hi-lo);
    }
    std::vector<double> d(n-1), baseD;
    for (int i=0;i<n-1;++i) {
        d[i]=p.smoothed[i+1]-p.smoothed[i];
        if (i>=half && i+1+half<p.nBaseline) baseD.push_back(d[i]);
    }
    const double meanD=baseD.empty() ? 0 : std::accumulate(baseD.begin(),baseD.end(),0.0)/baseD.size();
    p.derivativeNoise=opt.robust ? MAD(baseD) : RMS(baseD,meanD);
    const double maxY=*std::max_element(p.corrected.begin(),p.corrected.end());
    const double threshold=std::max(opt.threshold*p.derivativeNoise,1e-10*std::max(1.0,maxY));
    std::vector<int> candidates;
    // Collapse each connected above-threshold edge to its strongest derivative.
    for (int i=0;i<n-1;++i) {
        if (d[i]<=threshold) continue;
        int best=i;
        while (i+1<n-1 && d[i+1]>threshold) { ++i; if (d[i]>d[best]) best=i; }
        candidates.push_back(best);
    }
    std::sort(candidates.begin(),candidates.end(),[&](int a,int b){return d[a]>d[b];});
    std::vector<int> accepted;
    for (int i:candidates) {
        if (std::all_of(accepted.begin(),accepted.end(),[&](int j){return std::abs(i-j)>=opt.minSeparation;})) accepted.push_back(i);
    }
    std::sort(accepted.begin(),accepted.end());
    for (size_t k=0;k<accepted.size();++k) {
        const int i=accepted[k]; double offset=0;
        if (i>0 && i+1<n-1) {
            const double curvature=d[i-1]-2*d[i]+d[i+1];
            if (curvature<0) offset=std::clamp(0.5*(d[i-1]-d[i+1])/curvature,-0.5,0.5);
        }
        // Derivative d[i] is located at the midpoint of samples i and i+1.
        const double edge=(w.x[i]+w.x[i+1])/2+offset*(w.x[i+1]-w.x[i]);
        const int end=k+1<accepted.size() ? accepted[k+1]+1 : n;
        const double height=*std::max_element(p.corrected.begin()+i+1,p.corrected.begin()+end);
        p.pulses.push_back({i+1,edge,height});
    }
    if (!p.pulses.empty() && p.pulses.front().sample<p.nBaseline)
        std::cerr << "[WARNING] " << w.canvasName << "/" << w.objectName
                  << ": rising edge overlaps baseline window; use a shorter --baseline-samples.\n";
    return p;
}
PulseTemplate MakeTemplate(const std::vector<const Waveform*>& waves, const AnalysisOptions& opt,
                           const std::string& unit) {
    struct Clean { const Waveform* w; Processed p; double amplitude; };
    std::vector<Clean> clean;
    double lo=-std::numeric_limits<double>::infinity(),hi=std::numeric_limits<double>::infinity(),dt=0;
        for (const auto* w:waves) {
        auto p=Process(*w,opt);
        if (p.pulses.size()!=1) {
            std::cerr << "[WARNING] Skip template source " << w->canvasName << '/' << w->objectName
                      << ": " << p.pulses.size() << " pulse candidates\n"; continue;
        }
        const double amp=*std::max_element(p.corrected.begin(),p.corrected.end());
        if (!(amp>std::max(1e-12,opt.threshold*p.noiseRMS))) continue;
        const double spacing=Step(w->x);
        for (size_t i=1;i<w->x.size();++i)
            if (std::abs((w->x[i]-w->x[i-1])-spacing)>1e-5*spacing)
                throw std::runtime_error("Template sources must be uniformly sampled");
        if (dt && std::abs(spacing-dt)>1e-5*dt) throw std::runtime_error("Template source sample spacings differ");
        dt=spacing;
        lo=std::max(lo,w->x.front()-p.pulses[0].x);
        hi=std::min(hi,w->x.back()-p.pulses[0].x);
        clean.push_back({w,std::move(p),amp});
    }
    if (clean.empty() || hi-lo<3*dt) throw std::runtime_error("No usable clean single-pulse waveform for template");
    PulseTemplate h; h.unit=unit;
    const int start=static_cast<int>(std::ceil(lo/dt)), stop=static_cast<int>(std::floor(hi/dt));
    for (int j=start;j<=stop;++j) {
        const double t=j*dt; double val=0;
        for (const auto& c:clean) val+=Interpolate(c.w->x,c.p.corrected,t+c.p.pulses[0].x)/c.amplitude;
        h.x.push_back(t); h.y.push_back(val/clean.size());
    }
    // Peak-normalize the average so fitted A has an ADC-height interpretation.
    const double peak=*std::max_element(h.y.begin(),h.y.end());
    for (double& y:h.y) y/=peak;
    for (const auto& c:clean) h.source+=c.w->canvasPath+"/"+c.w->objectName+"\n";
    std::cout << "[INFO] Template: " << clean.size() << " clean waveforms, " << h.x.size() << " points\n";
    return h;
}
void SaveTemplate(const PulseTemplate& h,const std::string& path) {
    TFile file(path.c_str(),"RECREATE");
    if (file.IsZombie()) throw std::runtime_error("Cannot create template file");
    TGraph g(static_cast<int>(h.x.size()),h.x.data(),h.y.data()); g.SetName("pulse_template"); g.Write();
    TNamed("x_unit",h.unit.c_str()).Write(); TNamed("template_sources",h.source.c_str()).Write();
    const auto dot=path.find_last_of('.');
    std::ofstream csv(path.substr(0,dot)+".csv");
    if (!csv) throw std::runtime_error("Cannot create template CSV");
    csv << std::setprecision(17) << "x,y\n";
    for (size_t i=0;i<h.x.size();++i) csv << h.x[i] << ',' << h.y[i] << '\n';
    if (!csv) throw std::runtime_error("Error writing template CSV");
}
PulseTemplate LoadTemplate(const std::string& path,const std::string& unit) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(),"READ"));
    if (!file || file->IsZombie()) throw std::runtime_error("Cannot open template: "+path);
    auto* graph=dynamic_cast<TGraph*>(file->Get("pulse_template"));
    if (!graph || graph->GetN()<4) throw std::runtime_error("Template requires a TGraph named pulse_template with >=4 points");
    PulseTemplate h; h.unit=unit;
    if (auto* label=dynamic_cast<TNamed*>(file->Get("x_unit"))) {
        if (unit!=label->GetTitle()) throw std::runtime_error("Template x unit does not match waveform x unit");
    } else std::cerr << "[WARNING] Template has no x_unit metadata; assuming " << unit << '\n';
    Waveform check;
    for (int i=0;i<graph->GetN();++i) { double x,y;graph->GetPoint(i,x,y);check.x.push_back(x);check.y.push_back(y); }
    Validate(check); h.x=check.x; h.y=check.y;
    if (h.x.front()>0 || h.x.back()<0) throw std::runtime_error("Template x must include its rising-edge reference at zero");
    const double peak=*std::max_element(h.y.begin(),h.y.end());
    if (!(peak>0)) throw std::runtime_error("Template must have positive pulse polarity");
    for (double& y:h.y) y/=peak;
    return h;
}
FitResult FitPileup(const Waveform& w,const Processed& p,const PulseTemplate& h,const AnalysisOptions& opt) {
    FitResult result;
    if (w.x.size()<8) { result.status="too_few_samples";return result; }
    if (p.pulses.empty() && !std::isfinite(opt.fixedT1)) { result.status="no_pulse_candidate";return result; }
    const auto& x=w.x; const auto& y=p.corrected;
    const double dt=Step(x), sigma=p.noiseRMS>0 ? p.noiseRMS : 1.0;
    result.fitSigma=sigma;
    if (!p.noiseRMS) std::cerr << "[WARNING] Zero baseline noise; use 1 ADC for fit weights. Information criteria are diagnostic only.\n";
    // Triggered one-dimensional mode. The derivative reference matches the
    // template's alignment convention. An external trigger can override it.
    const double t1=std::isfinite(opt.fixedT1) ? opt.fixedT1 : p.pulses.front().x;
    if (t1<x.front() || t1>x.back()) { result.status="t1_out_of_range";return result; }
    const double maxDelta=std::min(opt.maxDelta,x.back()-t1);
    auto single=Project(x,y,h,t1,t1,false);
    if (!single.valid) { result.status="single_fit_failed";return result; }
    // Generalized matched-filter score: least-squares projection after removing
    // the constant/linear baseline. Time is fixed by the trigger in this mode.
    LinearFit best=single; double bestDelta=0;
    auto objective=[&](double delta) { return Project(x,y,h,t1,t1+delta,true).sse; };
    if (maxDelta/dt>100000) throw std::runtime_error("Delta scan exceeds 100000 steps; reduce --max-delta");
    const int steps=static_cast<int>(std::floor(maxDelta/dt));
    for (int i=0;i<=steps;++i) {
        const double delta=i*dt; auto fit=Project(x,y,h,t1,t1+delta,true);
        result.scanX.push_back(delta); result.scanY.push_back(fit.sse/(sigma*sigma));
        if (fit.sse<best.sse) { best=std::move(fit);bestDelta=delta; }
    }
    if (maxDelta>steps*dt) {
        auto fit=Project(x,y,h,t1,t1+maxDelta,true);
        result.scanX.push_back(maxDelta);result.scanY.push_back(fit.sse/(sigma*sigma));
        if (fit.sse<best.sse) {best=std::move(fit);bestDelta=maxDelta;}
    }
    // QR residualizes the delayed template against h1, 1 and t. The reduction
    // in weighted SSE equals the squared generalized matched-filter response
    // for the unconstrained interior solution; positivity clips boundary cases.
    for (double chi : result.scanY) result.matchedPower.push_back(std::max(0.0,single.sse/(sigma*sigma)-chi));
    if (maxDelta>0) {
        const double delta=Brent(objective,std::max(0.0,bestDelta-dt),std::min(maxDelta,bestDelta+dt),dt*1e-5);
        auto refined=Project(x,y,h,t1,t1+delta,true);
        if (refined.sse<best.sse) {best=std::move(refined);bestDelta=delta;}
    }
    result.valid=true;result.t1=t1;result.t2=t1+bestDelta;result.deltaT=bestDelta;
    result.A1=best.a1;result.A2=best.a2;
    result.chi1=single.sse/(sigma*sigma);result.chi2=best.sse/(sigma*sigma);
    // A,b0,b1 vs A1,A2,b0,b1,delta; count data-estimated t1 in both models.
    const int k1=std::isfinite(opt.fixedT1) ? 3 : 4, k2=k1+2;
    const double logN=std::log(x.size());
    result.ndf1=static_cast<int>(x.size())-k1;result.ndf2=static_cast<int>(x.size())-k2;
    result.aic1=result.chi1+2*k1;result.aic2=result.chi2+2*k2;
    result.bic1=result.chi1+k1*logN;result.bic2=result.chi2+k2*logN;
    result.status=bestDelta<=dt*1e-4 || best.a1==0 || best.a2==0 ? "single_or_unresolved" :
        (result.bic1-result.bic2>10 ? "pileup_candidate" : "no_strong_pileup_evidence");
    if (bestDelta>=maxDelta-dt*1e-4 && maxDelta>0) result.status="scan_boundary";
    result.single=std::move(single.model);result.dual=std::move(best.model);
    for (size_t i=0;i<y.size();++i) {result.residual1.push_back(y[i]-result.single[i]);result.residual2.push_back(y[i]-result.dual[i]);}
    return result;
}
