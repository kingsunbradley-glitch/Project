#include "Calibrator.h"

#include <TF1.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TSystem.h>
#include <TFitResult.h>
#include <TDirectory.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>

using ROOT::VecOps::RVec;

Calibrator::Calibrator(std::string inRoot,
                       std::string normTxt,
                       std::string outDat,
                       std::string diagMain,
                       std::string diagEd)
    : inRoot_(std::move(inRoot)),
      normTxt_(std::move(normTxt)),
      outDat_(std::move(outDat)),
      diagMain_(std::move(diagMain)),
      diagEd_(std::move(diagEd)) {
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);
    norm_.resize(TOTAL_CH);
    chRes_.resize(TOTAL_CH);
}

bool Calibrator::LoadNorm() {
    std::ifstream in(normTxt_);
    if (!in.is_open()) {
        std::cerr << "[WARNING] Cannot open norm file: " << normTxt_ << " (will use identity)\n";
        for (auto& p : norm_) { p.k = 1.0; p.b = 0.0; p.ok = 0; p.n = 0; }
        return false;
    }
    std::string line;
    int id; double k, b; int ok; long n;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        if (ss >> id >> k >> b >> ok >> n) {
            if (id >= 0 && id < TOTAL_CH) {
                norm_[id].k = k;
                norm_[id].b = b;
                norm_[id].ok = ok;
                norm_[id].n = n;
            }
        }
    }
    return true;
}

void Calibrator::BuildPlaneSummed(ROOT::RDataFrame& df, TH1D& hX, TH1D& hY, TH1D& hYH) {
    const unsigned int nSlots = df.GetNSlots();

    std::vector<std::unique_ptr<TH1D>> vX(nSlots), vY(nSlots), vYH(nSlots);
    for (unsigned int s = 0; s < nSlots; ++s) {
        vX[s] = std::unique_ptr<TH1D>((TH1D*)hX.Clone()); vX[s]->Reset(); vX[s]->SetDirectory(nullptr);
        vY[s] = std::unique_ptr<TH1D>((TH1D*)hY.Clone()); vY[s]->Reset(); vY[s]->SetDirectory(nullptr);
        vYH[s]= std::unique_ptr<TH1D>((TH1D*)hYH.Clone());vYH[s]->Reset();vYH[s]->SetDirectory(nullptr);
    }

    df.ForeachSlot(
        [&](unsigned int s,
            const RVec<double>& xE, const RVec<double>& xCh,
            const RVec<double>& yE, const RVec<double>& yCh,
            const RVec<double>& yhE, const RVec<double>& yhCh, unsigned short yhMul) {
            if (!xE.empty() && !xCh.empty()) {
                int xid = (int)xCh[0];
                if (xid >= 0 && xid < NUM_DSSDX_POS) {
                    double adc = norm_[xid].k * xE[0] + norm_[xid].b;
                    vX[s]->Fill(adc);
                }
            }
            if (!yE.empty() && !yCh.empty()) {
                int yid = NUM_DSSDX_POS + (int)yCh[0];
                if (yid >= NUM_DSSDX_POS && yid < NUM_DSSDX_POS + NUM_DSSDY_POS) {
                    double adc = norm_[yid].k * yE[0] + norm_[yid].b;
                    vY[s]->Fill(adc);
                }
            }
            if (yhMul == 1 && !yhE.empty() && !yhCh.empty()) {
                int yhid = NUM_DSSDX_POS + NUM_DSSDY_POS + (int)yhCh[0];
                if (yhid >= 0 && yhid < TOTAL_CH) {
                    double adc = norm_[yhid].k * yhE[0] + norm_[yhid].b;
                    vYH[s]->Fill(adc);
                }
            }
        },
        {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    for (unsigned int s = 0; s < nSlots; ++s) {
        hX.Add(vX[s].get());
        hY.Add(vY[s].get());
        hYH.Add(vYH[s].get());
    }
}

bool Calibrator::FindTop3PeaksTSpectrum(TH1D& h, std::vector<double>& peaks) {
    peaks.clear();
    TSpectrum spec(TS_MAX_PEAKS, 1.0);
    int nFound = spec.Search(&h, TS_SIGMA_BINS, "", TS_THRESHOLD);
    if (nFound <= 0) return false;

    // collect candidates: (height, pos)
    std::vector<std::pair<double,double>> cand;
    cand.reserve(nFound);
    double* xs = spec.GetPositionX();
    for (int i = 0; i < nFound; ++i) {
        double x = xs[i];
        int bin = h.GetXaxis()->FindBin(x);
        double y = h.GetBinContent(bin);
        cand.push_back({y, x});
    }
    std::sort(cand.begin(), cand.end(), [](auto& a, auto& b){return a.first > b.first;});

    // pick top3 with min separation
    for (auto& c : cand) {
        double x = c.second;
        bool tooClose = false;
        for (double p : peaks) {
            if (std::fabs(p - x) < PEAK_MIN_SEP) { tooClose = true; break; }
        }
        if (tooClose) continue;
        peaks.push_back(x);
        if (peaks.size() == 3) break;
    }
    if (peaks.size() != 3) return false;

    std::sort(peaks.begin(), peaks.end());
    return true;
}

bool Calibrator::GaussianRefine(TH1D& h, double seed, double& mean, double& sigma, double& chi2ndf) {
    chi2ndf = FAIL_CHI2;

    auto do_fit = [&](double m0, double win, double& m, double& s, double& c2) -> bool {
        double x1 = m0 - win;
        double x2 = m0 + win;
        TF1 f("f", "gaus(0)+pol0(3)", x1, x2);
        int b0 = h.GetXaxis()->FindBin(m0);
        double amp0 = std::max(1.0, h.GetBinContent(b0));
        f.SetParameters(amp0, m0, 10.0, 0.0);
        f.SetParLimits(1, x1, x2);
        f.SetParLimits(2, GAUS_SIGMA_MIN, GAUS_SIGMA_MAX);

        auto r = h.Fit(&f, "RQ0S");
        if ((int)r != 0) return false;
        m = f.GetParameter(1);
        s = std::fabs(f.GetParameter(2));
        double ndf = f.GetNDF();
        c2 = (ndf > 0) ? f.GetChisquare() / ndf : FAIL_CHI2;
        if (!std::isfinite(m) || !std::isfinite(s) || s <= 0) return false;
        return true;
    };

    // coarse
    double m1=0,s1=0,c21=FAIL_CHI2;
    if (!do_fit(seed, GAUS_COARSE_WIN, m1, s1, c21)) return false;

    // fine around +/- 1*sigma (but at least min win)
    double win2 = std::max(GAUS_FINE_MIN_WIN, s1);
    double m2=0,s2=0,c22=FAIL_CHI2;
    if (!do_fit(m1, win2, m2, s2, c22)) return false;

    mean = m2;
    sigma = s2;
    chi2ndf = c22;
    return true;
}

PlaneCalib Calibrator::CalibratePlane(TH1D& h, const char* name, TFile& fout) {
    PlaneCalib out;

    std::vector<double> seeds;
    if (!FindTop3PeaksTSpectrum(h, seeds)) {
        std::cerr << "[WARNING] " << name << ": TSpectrum did not yield 3 good peaks.\n";
        return out;
    }

    // refine each seed
    for (int i = 0; i < 3; ++i) {
        double m=0,s=0,c2=FAIL_CHI2;
        if (!GaussianRefine(h, seeds[i], m, s, c2)) {
            std::cerr << "[WARNING] " << name << ": Gaussian refine failed on peak " << i << "\n";
            return out;
        }
        out.adc[i] = m;
        out.sig[i] = s;
    }

    // linear fit Energy = K*ADC + B using 3 points
    TGraph g(3);
    for (int i = 0; i < 3; ++i) {
        g.SetPoint(i, out.adc[i], SOURCE_ENERGIES[i]);
    }
    TF1 fpol("fpol", "pol1", HIST_MIN, HIST_MAX);
    auto fr = g.Fit(&fpol, "Q0S");
    if ((int)fr != 0) {
        std::cerr << "[WARNING] " << name << ": pol1 fit failed\n";
        return out;
    }
    out.K = fpol.GetParameter(1);
    out.B = fpol.GetParameter(0);
    double ndf = fpol.GetNDF();
    out.chi2ndf = (ndf > 0) ? fpol.GetChisquare() / ndf : FAIL_CHI2;
    out.ok = std::isfinite(out.K) && std::isfinite(out.B) && out.K > 0;

    // write diagnostics
    fout.cd();
    h.SetTitle(Form("%s plane summed (normalized ADC);ADC_norm;Counts", name));
    h.Write(Form("h_%s_sum", name));
    g.Write(Form("g_%s_E_vs_ADC", name));

    return out;
}

void Calibrator::BuildFinalKB() {
    // fill per-channel final k,b based on plane
    for (int id = 0; id < TOTAL_CH; ++id) {
        const NormParam& np = norm_[id];
        ChResult& cr = chRes_[id];

        const PlaneCalib* pc = nullptr;
        if (id < NUM_DSSDX_POS) pc = &calX_;
        else if (id < NUM_DSSDX_POS + NUM_DSSDY_POS) pc = &calY_;
        else pc = &calYH_;

        if (!pc->ok) {
            cr.k = FAIL_K; cr.b = FAIL_B; continue;
        }
        // allow norm ok==0 but still apply identity; if fit failed, will use identity in file anyway
        cr.k = pc->K * np.k;
        cr.b = pc->K * np.b + pc->B;
        if (!std::isfinite(cr.k) || !std::isfinite(cr.b) || cr.k <= 0) {
            cr.k = FAIL_K; cr.b = FAIL_B;
        }
    }
}

bool Calibrator::FitFWHMOne(TH1D& h, double& fwhm, double& chi2ndf) {
    fwhm = FAIL_FWHM;
    chi2ndf = FAIL_CHI2;

    int b1 = h.FindBin(FWHM_REF_E - PEAK_FIT_WIN_E);
    int b2 = h.FindBin(FWHM_REF_E + PEAK_FIT_WIN_E);
    double sum = h.Integral(b1, b2);
    if (sum < MIN_PEAK_COUNTS_IN_WIN) return false;

    // initial mean from max bin in window
    int maxBin = b1;
    double maxY = -1;
    for (int b = b1; b <= b2; ++b) {
        double y = h.GetBinContent(b);
        if (y > maxY) { maxY = y; maxBin = b; }
    }
    double mu0 = h.GetXaxis()->GetBinCenter(maxBin);

    TF1 f("f", "gaus(0)+pol0(3)", mu0 - PEAK_FIT_WIN_E, mu0 + PEAK_FIT_WIN_E);
    f.SetParameters(std::max(1.0, maxY), mu0, 20.0, 0.0);
    f.SetParLimits(1, mu0 - PEAK_FIT_WIN_E, mu0 + PEAK_FIT_WIN_E);
    f.SetParLimits(2, 1.0, 200.0);

    auto r = h.Fit(&f, "RQ0S");
    if ((int)r != 0) return false;

    double sigma = std::fabs(f.GetParameter(2));
    if (!std::isfinite(sigma) || sigma <= 0) return false;

    fwhm = 2.355 * sigma;
    double ndf = f.GetNDF();
    chi2ndf = (ndf > 0) ? f.GetChisquare() / ndf : FAIL_CHI2;
    return true;
}

void Calibrator::ComputePerChannelFWHM(ROOT::RDataFrame& df, TFile& outEd) {
    const unsigned int nSlots = df.GetNSlots();

    // per slot, per channel hist
    std::vector<std::vector<std::unique_ptr<TH1D>>> hslot;
    hslot.resize(nSlots);
    for (unsigned int s = 0; s < nSlots; ++s) {
        hslot[s].resize(TOTAL_CH);
        for (int id = 0; id < TOTAL_CH; ++id) {
            hslot[s][id] = std::make_unique<TH1D>(Form("h_tmp_s%u_ch%d", s, id), "", E_HIST_BINS, E_HIST_MIN, E_HIST_MAX);
            hslot[s][id]->SetDirectory(nullptr);
        }
    }

    df.ForeachSlot(
        [&](unsigned int s,
            const RVec<double>& xE, const RVec<double>& xCh,
            const RVec<double>& yE, const RVec<double>& yCh,
            const RVec<double>& yhE, const RVec<double>& yhCh, unsigned short yhMul) {
            if (!xE.empty() && !xCh.empty()) {
                int id = (int)xCh[0];
                if (id >= 0 && id < NUM_DSSDX_POS) {
                    const ChResult& cr = chRes_[id];
                    double e = cr.k * xE[0] + cr.b;
                    hslot[s][id]->Fill(e);
                }
            }
            if (!yE.empty() && !yCh.empty()) {
                int id = NUM_DSSDX_POS + (int)yCh[0];
                if (id >= NUM_DSSDX_POS && id < NUM_DSSDX_POS + NUM_DSSDY_POS) {
                    const ChResult& cr = chRes_[id];
                    double e = cr.k * yE[0] + cr.b;
                    hslot[s][id]->Fill(e);
                }
            }
            if (yhMul == 1 && !yhE.empty() && !yhCh.empty()) {
                int id = NUM_DSSDX_POS + NUM_DSSDY_POS + (int)yhCh[0];
                if (id >= 0 && id < TOTAL_CH) {
                    const ChResult& cr = chRes_[id];
                    double e = cr.k * yhE[0] + cr.b;
                    hslot[s][id]->Fill(e);
                }
            }
        },
        {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    outEd.cd();
    auto* dir = outEd.mkdir("Strips");
    dir->cd();

    TH1D hFwhm("h_fwhm", "FWHM @5804.8 keV vs Channel;Channel;FWHM (keV)", TOTAL_CH, 0, TOTAL_CH);
    TH1D hChi2("h_chi2ndf", "Chi2/NDF of FWHM fit vs Channel;Channel;Chi2/NDF", TOTAL_CH, 0, TOTAL_CH);

    for (int id = 0; id < TOTAL_CH; ++id) {
        // merge slots into hslot[0][id]
        for (unsigned int s = 1; s < nSlots; ++s) {
            hslot[0][id]->Add(hslot[s][id].get());
        }
        TH1D* h = hslot[0][id].get();
        h->SetName(Form("h_ch%03d", id));
        h->SetTitle(Form("Channel %d calibrated;Energy (keV);Counts", id));

        double fwhm=FAIL_FWHM, chi2=FAIL_CHI2;
        bool ok = FitFWHMOne(*h, fwhm, chi2);
        if (!ok) {
            fwhm = FAIL_FWHM;
            chi2 = FAIL_CHI2;
        }
        chRes_[id].fwhm = fwhm;
        chRes_[id].chi2ndf = chi2;

        hFwhm.SetBinContent(id+1, fwhm);
        hChi2.SetBinContent(id+1, chi2);

        h->Write();
    }

    dir->cd();
    hFwhm.Write();
    hChi2.Write();
    outEd.cd();
}

void Calibrator::WriteEnerCalDat() const {
    std::ofstream out(outDat_);
    out << "# id  k  b  FWHM_keV  chi2ndf\n";
    out << std::fixed << std::setprecision(8);
    for (int id = 0; id < TOTAL_CH; ++id) {
        const auto& r = chRes_[id];
        out << id << " " << r.k << " " << r.b << " " << r.fwhm << " " << r.chi2ndf << "\n";
    }
}

int Calibrator::Run() {
    std::cout << "=== Step2: Cali ===\n";

    LoadNorm();

    if (gSystem->AccessPathName(inRoot_.c_str())) {
        std::cerr << "[ERROR] Cannot open input: " << inRoot_ << "\n";
        return 1;
    }

    ROOT::RDataFrame df(TREE_NAME, inRoot_);

    TH1D hX("h_tot_X", "X plane summed (normalized ADC);ADC_norm;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    TH1D hY("h_tot_Y", "Y plane summed (normalized ADC);ADC_norm;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    TH1D hYH("h_tot_YH", "YH plane summed (normalized ADC);ADC_norm;Counts", HIST_BINS, HIST_MIN, HIST_MAX);

    BuildPlaneSummed(df, hX, hY, hYH);

    TFile fMain(diagMain_.c_str(), "RECREATE");
    calX_ = CalibratePlane(hX, "X", fMain);
    calY_ = CalibratePlane(hY, "Y", fMain);
    calYH_= CalibratePlane(hYH,"YH", fMain);

    // If YH failed but looks usable, fall back to Y plane
    if (!calYH_.ok && calY_.ok) {
        calYH_ = calY_;
    }

    fMain.Close();

    std::cout << "Plane KB: X(" << (calX_.ok?"OK":"FAIL") << ") K=" << calX_.K << " B=" << calX_.B << "\n";
    std::cout << "          Y(" << (calY_.ok?"OK":"FAIL") << ") K=" << calY_.K << " B=" << calY_.B << "\n";
    std::cout << "         YH(" << (calYH_.ok?"OK":"FAIL") << ") K=" << calYH_.K << " B=" << calYH_.B << "\n";

    BuildFinalKB();

    TFile fEd(diagEd_.c_str(), "RECREATE");
    ComputePerChannelFWHM(df, fEd);
    fEd.Close();

    WriteEnerCalDat();

    std::cout << "Wrote: " << outDat_ << "\n";
    std::cout << "Wrote: " << diagMain_ << " and " << diagEd_ << "\n";
    std::cout << "=== Step2: Cali done ===\n";
    return 0;
}
