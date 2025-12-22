#include "Calibrator.h"
#include "Config.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TFitResult.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm> // 用于排序
#include <cmath>     // 用于 sqrt, pow

using namespace std;
using namespace ROOT;

Calibrator::Calibrator() {
    if(NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);
    norm_params_.resize(TOTAL_CH);
    h_tot_X  = new TH1D("h_tot_X", "Total X (Norm ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    h_tot_Y  = new TH1D("h_tot_Y", "Total Y (Norm ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    h_tot_YH = new TH1D("h_tot_YH", "Total YH (Raw ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    
    // FWHM Summary: X轴是StripID, Y轴是FWHM
    h_fwhm_summary = new TH1D("h_fwhm", "FWHM Distribution per Strip;Strip ID;FWHM (keV)", TOTAL_CH, 0, TOTAL_CH);
}

Calibrator::~Calibrator() {
    delete h_tot_X; delete h_tot_Y; delete h_tot_YH; delete h_fwhm_summary;
}

void Calibrator::LoadNormParams() {
    ifstream in(NORM_PARAM_FILE);
    if(!in.is_open()) return;
    int id; double k, b;
    string line;
    while(getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        if(ss >> id >> k >> b && id >= 0 && id < TOTAL_CH) {
            norm_params_[id] = {k, b};
        }
    }
}

void Calibrator::Run() {
    LoadNormParams(); 

    TChain chain(TREE_NAME);
    if (!gSystem->AccessPathName(DECAY_ONLY_FILE)) {
        chain.Add(DECAY_ONLY_FILE);
    } else {
        cerr << "[ERROR] Decay file not found: " << DECAY_ONLY_FILE << endl;
        return;
    }

    ROOT::RDataFrame df(chain);
    FillSpectra(df); 

    TFile* f_diag = new TFile(CALIB_DIAG_ROOT, "RECREATE");
    f_diag->cd(); 
    // 写入原始总谱
    h_tot_X->Write(); h_tot_Y->Write(); h_tot_YH->Write();

    ofstream out_txt("peaks_list.txt");
    out_txt << fixed << setprecision(2);
    out_txt << "=== DSSD Peak List Report (Fitted Precision) ===" << endl;
    
    // 3. 刻度 + 验证
    FindAndListPeaks(h_tot_X, "X-Plane Normalized (ADC)", out_txt);
    PlaneResult rX  = CalibratePlane(h_tot_X, "X", f_diag, out_txt);

    FindAndListPeaks(h_tot_Y, "Y-Plane Normalized (ADC)", out_txt);
    PlaneResult rY  = CalibratePlane(h_tot_Y, "Y", f_diag, out_txt);

    PlaneResult rYH = CalibratePlane(h_tot_YH, "YH", f_diag, out_txt); 

    if (rYH.K == 1.0 && rYH.B == 0.0) {
        rYH.K = rX.K;
        rYH.B = rX.B;
    }

    out_txt.close();
    cout << "--> Peaks saved to peaks_list.txt" << endl;

    // 4. 输出最终 dat 参数 (此时会填充 h_fwhm_summary)
    WriteOutput(rX, rY, rYH);

    // [新增] 将 FWHM 分布图写入 ROOT 文件
    f_diag->cd();
    h_fwhm_summary->SetLineWidth(2);
    h_fwhm_summary->SetLineColor(kBlue);
    h_fwhm_summary->Write();

    f_diag->Close(); delete f_diag;
}

void Calibrator::FillSpectra(ROOT::RDataFrame& df) {
    auto nSlots = df.GetNSlots();
    vector<TH1D*> vX(nSlots), vY(nSlots), vYH(nSlots);
    for(unsigned i=0; i<nSlots; ++i) {
        vX[i] = (TH1D*)h_tot_X->Clone(); vX[i]->Reset(); vX[i]->SetDirectory(0);
        vY[i] = (TH1D*)h_tot_Y->Clone(); vY[i]->Reset(); vY[i]->SetDirectory(0);
        vYH[i]= (TH1D*)h_tot_YH->Clone();vYH[i]->Reset();vYH[i]->SetDirectory(0);
    }

    df.ForeachSlot([&](unsigned int s, const RVec<double>& xe, const RVec<double>& xch,
                       const RVec<double>& ye, const RVec<double>& ych,
                       const RVec<double>& yhe, const RVec<double>& yhch, UShort_t yhm) {
        if(!xch.empty() && xch[0] < 128) {
            int c = (int)xch[0]; 
            vX[s]->Fill(xe[0]*norm_params_[c].k + norm_params_[c].b);
        }
        if(!ych.empty() && ych[0] < 48) {
            int c = 128 + (int)ych[0]; 
            vY[s]->Fill(ye[0]*norm_params_[c].k + norm_params_[c].b);
        }
        if(yhm==1 && !yhch.empty() && yhch[0] < 48) {
            int c = 176 + (int)yhch[0]; 
            vYH[s]->Fill(yhe[0]*norm_params_[c].k + norm_params_[c].b);
        }
    }, {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    for(auto h : vX) { h_tot_X->Add(h); delete h; }
    for(auto h : vY) { h_tot_Y->Add(h); delete h; }
    for(auto h : vYH){ h_tot_YH->Add(h); delete h; }
}

PlaneResult Calibrator::CalibratePlane(TH1D* h_raw, const char* name, TFile* diag, ofstream& out) {
    TH1D* h_sub = SubtractBackground(h_raw);
    
    vector<double> meas, ref;
    for(auto& p : REF_PEAKS) {
        double pos, sig;
        if(FindPeakGaussian(h_sub, p.energy, p.window * 2.0, pos, sig)) { 
            meas.push_back(pos); ref.push_back(p.energy);
        }
    }
    delete h_sub; 

    PlaneResult res;
    if(meas.size() >= 2) {
        TGraph g(meas.size(), &meas[0], &ref[0]);
        TF1 f("f","pol1"); g.Fit(&f, "Q");
        res.B = f.GetParameter(0); res.K = f.GetParameter(1);
        
        diag->cd();
        g.Write(TString::Format("gr_calib_%s", name));
    } else {
        return res; 
    }

    TH1D* h_calib = ApplyCalibration(h_raw, res.K, res.B, TString::Format("h_calib_%s", name));
    
    FindAndListPeaks(h_calib, TString::Format("%s-Plane Calibrated (keV)", name).Data(), out);

    // [修改] 计算 FWHM 和 Residual RMS (峰位偏差)
    TH1D* h_calib_sub = SubtractBackground(h_calib);
    
    double sum_sq_diff = 0.0;
    int count_diag = 0;

    for(auto& p : DIAG_PEAKS) {
        double pos, sig;
        if(FindPeakGaussian(h_calib_sub, p.energy, p.window, pos, sig)) {
            double fwhm = 2.355 * sig;
            if(fwhm > res.FWHM_Max) res.FWHM_Max = fwhm;

            // 计算偏差：(拟合值 - 理论值)
            double diff = pos - p.energy;
            sum_sq_diff += diff * diff;
            count_diag++;
        }
    }

    // 计算 RMS Error (均方根误差)
    if(count_diag > 0) {
        res.Res_RMS = std::sqrt(sum_sq_diff / count_diag);
    }

    diag->cd();
    h_calib->SetTitle(TString::Format("Calibrated %s (keV);Energy (keV);Counts", name));
    h_calib->Write();

    delete h_calib; delete h_calib_sub;
    return res;
}

void Calibrator::FindAndListPeaks(TH1D* h, const char* title, ofstream& out) {
    h->GetXaxis()->SetRange(0, 0); 
    TSpectrum s(200); 
    int nFound = s.Search(h, 4, "nobackground", 0.02); 

    Double_t *xpeaks = s.GetPositionX();
    Double_t *ypeaks = s.GetPositionY();

    vector<PeakInfo> peaks;
    for(int i=0; i<nFound; ++i) {
        double raw_pos = xpeaks[i];
        double final_pos = raw_pos;
        
        double fit_range = 30.0; 
        TF1 f("f_quick", "gaus", raw_pos - fit_range, raw_pos + fit_range);
        f.SetParameter(1, raw_pos);
        f.SetParLimits(1, raw_pos - fit_range, raw_pos - fit_range);
        
        TFitResultPtr r = h->Fit(&f, "RQN"); 
        if(int(r) == 0) {
            final_pos = f.GetParameter(1); 
        }

        peaks.push_back({final_pos, ypeaks[i]});
    }

    sort(peaks.begin(), peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.height > b.height;
    });

    out << "--- " << title << " ---" << endl;
    out << "Rank\tPos(Fitted)\tHeight(Raw)" << endl;
    for(size_t i=0; i<peaks.size(); ++i) {
        out << (i+1) << "\t" << peaks[i].pos << "\t" << peaks[i].height << endl;
    }
    out << "--------------------------------" << endl << endl;
}

TH1D* Calibrator::ApplyCalibration(TH1D* h_in, double K, double B, const char* new_name) {
    TH1D* h_out = new TH1D(new_name, new_name, HIST_BINS, HIST_MIN, HIST_MAX);
    h_out->SetDirectory(0);
    h_out->Sumw2(); 
    for(int i = 1; i <= h_in->GetNbinsX(); ++i) {
        double content = h_in->GetBinContent(i);
        if(content <= 0) continue; 
        double raw_center = h_in->GetBinCenter(i);
        double cal_energy = K * raw_center + B;
        h_out->Fill(cal_energy, content); 
    }
    h_out->SetEntries(h_in->GetEntries()); 
    return h_out;
}

bool Calibrator::FindPeakGaussian(TH1D* h, double center, double win, double &peakPos, double &sigma) {
    h->GetXaxis()->SetRange(0, 0); 
    
    TSpectrum s(100); 
    int nFound = s.Search(h, 4, "nobackground", 0.05);

    double targetPeak = -1;
    double maxHeight = -1.0;

    if (nFound > 0) {
        double* xpeaks = s.GetPositionX();
        for (int i = 0; i < nFound; ++i) {
            if (abs(xpeaks[i] - center) <= win) {
                int bin = h->FindBin(xpeaks[i]);
                double currentHeight = h->GetBinContent(bin);
                if (currentHeight > maxHeight) {
                    maxHeight = currentHeight;
                    targetPeak = xpeaks[i];
                }
            }
        }
    }

    if (targetPeak == -1) {
        h->GetXaxis()->SetRangeUser(center - win, center + win);
        int maxBin = h->GetMaximumBin();
        if (h->GetBinContent(maxBin) < 5) return false; 
        targetPeak = h->GetBinCenter(maxBin);
        if (abs(targetPeak - center) > win) return false;
    }

    peakPos = targetPeak; 
    
    double fit_radius = (win < 60.0) ? win/1.5 : 40.0; 
    double fit_min = peakPos - fit_radius;
    double fit_max = peakPos + fit_radius;
    
    h->GetXaxis()->SetRangeUser(fit_min, fit_max);

    TF1 func("f_comb", "gaus(0) + pol1(3)", fit_min, fit_max);
    
    int bin = h->FindBin(peakPos);
    double height = h->GetBinContent(bin);
    double bg_left  = h->GetBinContent(h->FindBin(fit_min));
    double bg_right = h->GetBinContent(h->FindBin(fit_max));
    double bg_const = (bg_left + bg_right) / 2.0;
    if (bg_const >= height) bg_const = 0;

    func.SetParameter(0, height - bg_const); 
    func.SetParameter(1, peakPos);           
    func.SetParameter(2, 12.0);          
    func.SetParLimits(2, 2.0, 100.0);

    func.SetParameter(3, bg_const);
    func.SetParameter(4, 0);

    TFitResultPtr r = h->Fit(&func, "RQNS");
    if (int(r) == 0 && r->IsValid()) {
        peakPos = func.GetParameter(1); 
        sigma   = func.GetParameter(2); 
        return true;
    }
    return false;
}

TH1D* Calibrator::SubtractBackground(TH1D* h) {
    if(!ENABLE_BKG_SUB) return (TH1D*)h->Clone();
    TH1D* h_out = (TH1D*)h->Clone(); h_out->SetDirectory(0);
    TSpectrum s; 
    TH1* bkg = s.Background(h_out, BKG_ITERATIONS);
    h_out->Add(bkg, -1);
    return h_out;
}

// [修改] 输出格式增加 RMS Error，并填充 h_fwhm_summary
void Calibrator::WriteOutput(const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH) {
    ofstream out(OUTPUT_DAT_FILE);
    out << fixed << setprecision(4);
    // 更新表头
    out << "Ch\tb\t\tk\t\tFWHM\tPeak_RMS_Err" << endl;
    
    auto w = [&](int start, int end, const PlaneResult& pr) {
        for(int i=start; i<end; ++i) {
            double fk = pr.K * norm_params_[i].k;
            double fb = pr.K * norm_params_[i].b + pr.B;
            
            // 输出格式: Ch, b, k, FWHM, RMS_Error
            out << i << "\t" << fb << "\t" << fk << "\t" << pr.FWHM_Max << "\t" << pr.Res_RMS << endl;
            
            // 填充 ROOT 里的 FWHM 分布图
            h_fwhm_summary->SetBinContent(i+1, pr.FWHM_Max);
        }
    };
    
    w(0, 128, rX);
    w(128, 176, rY);
    w(176, 224, rYH);
    out.close();
}