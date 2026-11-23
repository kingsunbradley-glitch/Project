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
#include <algorithm> 
#include <cmath>     

using namespace std;
using namespace ROOT;

Calibrator::Calibrator() {
    if(NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);
    norm_params_.resize(TOTAL_DSSD_CH);
    strip_fwhms_.resize(TOTAL_DSSD_CH, 0.0); 
    strip_dE_.resize(TOTAL_DSSD_CH, 0.0);

    h_tot_X  = new TH1D("h_tot_X", "Total X (Norm ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    h_tot_Y  = new TH1D("h_tot_Y", "Total Y (Norm ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    h_tot_YH = new TH1D("h_tot_YH", "Total YH (Norm ADC);ADC;Counts", HIST_BINS, HIST_MIN, HIST_MAX);
    
    h_fwhm_summary = new TH1D("h_fwhm", "FWHM Distribution per Strip;Strip ID;FWHM (keV)",
                              TOTAL_DSSD_CH, 0, TOTAL_DSSD_CH);
}

Calibrator::~Calibrator() {
    delete h_tot_X; delete h_tot_Y; delete h_tot_YH; delete h_fwhm_summary;
}

void Calibrator::LoadNormParams() {
    ifstream in(NORM_PARAM_FILE);
    if(!in.is_open()) {
        cerr << "[WARN] 归一化参数文件打不开：" << NORM_PARAM_FILE
             << "，全部通道按 (k,b)=(1,0) 处理" << endl;
        return;
    }
    int id; double k, b;
    string line;
    int n_x = 0, n_y = 0, n_yh = 0;
    while(getline(in, line)) {
        if(line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        if(ss >> id >> k >> b && id >= 0 && id < TOTAL_DSSD_CH) {
            norm_params_[id] = {k, b};
            if(id < CH_BASE_Y)        ++n_x;
            else if(id < CH_BASE_YH)  ++n_y;
            else                      ++n_yh;
        }
    }
    cout << "--> Normalization loaded: X=" << n_x << "/" << NUM_DSSDX_POS
         << "  Y=" << n_y << "/" << NUM_DSSDY_POS
         << "  YH=" << n_yh << "/" << NUM_DSSDYH_POS << endl;
    if(n_yh == 0) {
        cerr << "[WARN] 参数文件里没有 YH (" << CH_BASE_YH << "-" << (TOTAL_DSSD_CH-1)
             << ") 的系数，YH 面将保持原始幅度道。" << endl;
        cerr << "[WARN] 请重新运行 Step 1（本版本会输出 " << CH_BASE_YH << "-"
             << (TOTAL_DSSD_CH-1) << "）。" << endl;
    }
}

// 读参考峰表（ref_peaks.txt）。换实验只改文件，不用重编译。
// 任何一项解析失败都退回 Config.h 的默认值并打印原因，不静默。
void Calibrator::LoadPeakTable() {
    ref_peaks_        = REF_PEAKS_DEFAULT;
    diag_peaks_       = DIAG_PEAKS_DEFAULT;
    strip_peak_E_     = target_E;
    strip_peak_win_   = fit_win;

    ifstream in(REF_PEAK_FILE);
    if(!in.is_open()) {
        cout << "[WARN] 参考峰表 " << REF_PEAK_FILE << " 打不开，使用 Config.h 内置默认值" << endl;
        return;
    }

    vector<PeakDef> ref, diag;
    bool have_strip = false;
    double strip_E = 0.0, strip_win = 0.0;
    int section = 0;        // 1 = ref, 2 = diag, 3 = strip
    int n_bad = 0;
    string line;
    while(getline(in, line)) {
        size_t hash = line.find('#');
        if(hash != string::npos) line.erase(hash);
        // 去掉首尾空白
        size_t b = line.find_first_not_of(" \t\r\n");
        if(b == string::npos) continue;
        size_t e = line.find_last_not_of(" \t\r\n");
        line = line.substr(b, e - b + 1);
        if(line.empty()) continue;

        if(line[0] == '[') {
            if(line == "[ref]")        section = 1;
            else if(line == "[diag]")  section = 2;
            else if(line == "[strip]") section = 3;
            else {
                cerr << "[WARN] 参考峰表：" << REF_PEAK_FILE << " 里有未知的节 " << line
                     << "，该节的内容会被忽略" << endl;
                section = 0;
            }
            continue;
        }

        istringstream ss(line);
        double en = 0.0, win = 30.0;
        if(!(ss >> en)) { ++n_bad; continue; }
        if(!(ss >> win)) win = 30.0;      // 只给能量时用 30 keV 默认半窗
        if(en <= 0 || win <= 0 || win > 1000) { ++n_bad; continue; }

        if(section == 1)      ref.push_back({en, win});
        else if(section == 2) diag.push_back({en, win});
        else if(section == 3) {
            if(have_strip) continue;      // 只有第一行有效
            strip_E = en; strip_win = win; have_strip = true;
        } else {
            ++n_bad;                      // 没写节名就出现的数字行
        }
    }

    if(n_bad > 0)
        cerr << "[WARN] 参考峰表里有 " << n_bad << " 行无法解析（缺节名或数值非法），已跳过" << endl;

    if(ref.empty()) {
        cerr << "[WARN] 参考峰表里没有 [ref] 峰，沿用 Config.h 默认的 " << REF_PEAKS_DEFAULT.size()
             << " 个参考峰" << endl;
    } else {
        ref_peaks_ = ref;
    }
    if(diag.empty()) {
        cerr << "[WARN] 参考峰表里没有 [diag] 峰，沿用 Config.h 默认的 " << DIAG_PEAKS_DEFAULT.size()
             << " 个诊断峰" << endl;
    } else {
        diag_peaks_ = diag;
    }
    if(have_strip && strip_win > 0) {
        strip_peak_E_   = strip_E;
        strip_peak_win_ = strip_win;
    }

    cout << "--> 参考峰表 " << REF_PEAK_FILE << " (ref " << ref_peaks_.size()
         << " 个, diag " << diag_peaks_.size() << " 个, 逐条检查峰 "
         << strip_peak_E_ << " keV)" << endl;
    for(auto& p : ref_peaks_)
        cout << "      ref  " << p.energy << " keV (窗 ±" << p.window * 2.0 << " keV)" << endl;
    for(auto& p : diag_peaks_)
        cout << "      diag " << p.energy << " keV (窗 ±" << p.window << " keV)" << endl;
}

void Calibrator::Run() {
    Print_Channel_Map();
    LoadNormParams(); 
    LoadPeakTable();

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
    h_tot_X->Write(); h_tot_Y->Write(); h_tot_YH->Write();

    ofstream out_txt("peaks_list.txt");
    out_txt << fixed << setprecision(2);
    out_txt << "=== DSSD Peak List Report (Fitted Precision) ===" << endl;
    
    // 3. 全局刻度
    // Step 1 已把 X/Y/YH 三个面都桥接到 RefY 的原始刻度上，所以三者理论上共用同一个
    // 绝对刻度 (K,B)；各面单独再拟合一次只是为了吸收面间残余差异。参考峰不足或拟合
    // 被质量检查否决时（例如 YH 的两个候选峰被同一个物理峰吞掉），不独立定 (K,B)，
    // 直接沿用 X 面刻度，避免把退化拟合的结果写进表格。
    FindAndListPeaks(h_tot_X, "X-Plane Normalized (ADC)", out_txt);
    PlaneResult rX  = CalibratePlane(h_tot_X, "X", f_diag, out_txt);

    FindAndListPeaks(h_tot_Y, "Y-Plane Normalized (ADC)", out_txt);
    PlaneResult rY  = CalibratePlane(h_tot_Y, "Y", f_diag, out_txt);

    PlaneResult rYH = CalibratePlane(h_tot_YH, "YH", f_diag, out_txt, rX.valid ? &rX : nullptr);

    cout << "---- 绝对刻度汇总 ----" << endl;
    auto show = [](const char* n, const PlaneResult& r) {
        cout << "  " << n << ": K=" << r.K << " B=" << r.B;
        if (r.inherited)
            cout << "   [沿用 X 面刻度；本面只认到 " << r.n_points << " 个参考峰]";
        else if (r.valid)
            cout << "   [本面拟合: " << r.n_points << " 点, 跨度 " << r.span
                 << " keV, 残差 rms " << r.fit_rms << " keV]";
        else
            cout << "   [未刻度: 参考峰不足或被否决]";
        cout << endl;
    };
    show("X ", rX); show("Y ", rY); show("YH", rYH);
    // X 面是唯一的基准：它失败时另外两面也没有可继承的对象，整张表会退化成
    // k=1,b=0（等于"未刻度"）。这种情况必须喊出来，否则下游 root2root 会
    // 拿着看起来正常、其实是原始幅度道的表继续跑。
    if(!rX.valid && !rY.valid && !rYH.valid) {
        cerr << "[ERROR] 三个平面都没有得到独立刻度，输出表将全是 k=1, b=0（未刻度）。" << endl;
        cerr << "[ERROR] 请检查 " << REF_PEAK_FILE << " 的 [ref] 节：至少需要 "
             << REF_MIN_POINTS << " 个能定位的峰；跨度越大 K 越准（建议 >= "
             << REF_WARN_SPAN_KEV << " keV）。"
             << endl;
    }
    cout << "----------------------" << endl;
    
    out_txt.close();
    cout << "--> Peaks saved to peaks_list.txt" << endl;

    // 4. 逐条计算 FWHM 并保存每一条的图
    cout << "--> Calculating FWHM and Saving Individual Spectra..." << endl;
    CalculateStripFWHM(df, rX, rY, rYH, f_diag);

    // 5. 输出最终参数
    WriteOutput(rX, rY, rYH);

    f_diag->cd();
    h_fwhm_summary->SetLineWidth(2);
    h_fwhm_summary->SetLineColor(kBlue);
    h_fwhm_summary->Write();

    f_diag->Close(); delete f_diag;
}

// CalculateStripFWHM: 逐条 6113 keV 峰形。X/Y/YH 三个面都做——
// YH 现在有归一化和刻度了，如果跳过它，表里 176-223 的 FWHM/dE 永远是 0，
// 就没法发现坏 YH 条。
void Calibrator::CalculateStripFWHM(ROOT::RDataFrame& df, const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH, TFile* diag) {
    unsigned int nSlots = df.GetNSlots();

    int max_ch = TOTAL_DSSD_CH;   // X + Y + YH

    vector<vector<TH1D*>> h_slots(nSlots, vector<TH1D*>(max_ch, nullptr));
    
    // 1. 初始化
    for(unsigned int s=0; s<nSlots; ++s) {
        for(int i=0; i<max_ch; ++i) {
            h_slots[s][i] = new TH1D(Form("h_tmp_%d_%d", s, i), "", 500, 5000, 7500); 
            h_slots[s][i]->SetDirectory(0);
        }
    }

    // 2. 填充 (仅 X 和 Y)
    df.ForeachSlot([&](unsigned int s, const RVec<double>& xe, const RVec<double>& xch,
                       const RVec<double>& ye, const RVec<double>& ych,
                       const RVec<double>& yhe, const RVec<double>& yhch, UShort_t yhm) {
        // X Plane
        if(!xch.empty() && xch[0] < NUM_DSSDX_POS) {
            int c = CH_BASE_X + (int)xch[0];
            double e_cal = rX.K * (xe[0]*norm_params_[c].k + norm_params_[c].b) + rX.B;
            h_slots[s][c]->Fill(e_cal);
        }
        // Y Plane
        if(!ych.empty() && ych[0] < NUM_DSSDY_POS) {
            int c = CH_BASE_Y + (int)ych[0];
            double e_cal = rY.K * (ye[0]*norm_params_[c].k + norm_params_[c].b) + rY.B;
            h_slots[s][c]->Fill(e_cal);
        }
        // YH Plane（Y 面并行读出）
        if(yhm == 1 && !yhch.empty() && yhch[0] < NUM_DSSDYH_POS) {
            int c = CH_BASE_YH + (int)yhch[0];
            double e_cal = rYH.K * (yhe[0]*norm_params_[c].k + norm_params_[c].b) + rYH.B;
            h_slots[s][c]->Fill(e_cal);
        }
    }, {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    // 3. 处理 X / Y / YH
    TDirectory *dir = diag->mkdir("Strips_Calib");
    dir->cd();

   

    for(int i=0; i<max_ch; ++i) {
        TH1D* h_sum = h_slots[0][i];
        for(unsigned int s=1; s<nSlots; ++s) {
            h_sum->Add(h_slots[s][i]);
            delete h_slots[s][i];
        }

        h_sum->SetName(Form("h_strip_%d", i));
        const char* type_name = (i < CH_BASE_Y) ? "X" : (i < CH_BASE_YH ? "Y" : "YH");
        int local_i = (i < CH_BASE_Y) ? i : (i < CH_BASE_YH ? i - CH_BASE_Y : i - CH_BASE_YH);
        h_sum->SetTitle(Form("%s-Strip %d Calibrated;Energy (keV);Counts", type_name, local_i));

        if(h_sum->GetEntries() < 20) {
            strip_fwhms_[i] = 0.0;
            strip_dE_[i] = 0.0;
        } else {
            double pos, sig;
            if(FindPeakGaussian(h_sum, strip_peak_E_, strip_peak_win_, pos, sig)) {
                double fwhm = 2.355 * sig;
                // 6 MeV 重离子峰的 FWHM 不可能超过 200 keV（实测 X 中位 29、Y 36）。
                // 旧判据放到 300 keV，会把噪声拟合出来的 235 keV "峰"当成有效结果。
                if(fwhm > 0 && fwhm < 200.0) {
                    strip_fwhms_[i] = fwhm;
                    strip_dE_[i] = pos - strip_peak_E_;   // 该条自己还残留的峰位偏差
                } else {
                    strip_fwhms_[i] = 0.0;
                    strip_dE_[i] = 0.0;
                }
            } else {
                strip_fwhms_[i] = 0.0;
                strip_dE_[i] = 0.0;
            }
        }
        
        // [修改] 先 Unzoom 再保存
        h_sum->GetXaxis()->SetRange(0, 0); 
        h_sum->Write();
        delete h_sum;
    }
    
    diag->cd();
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
        if(!xch.empty() && xch[0] < NUM_DSSDX_POS) {
            int c = CH_BASE_X + (int)xch[0]; 
            vX[s]->Fill(xe[0]*norm_params_[c].k + norm_params_[c].b);
        }
        if(!ych.empty() && ych[0] < NUM_DSSDY_POS) {
            int c = CH_BASE_Y + (int)ych[0]; 
            vY[s]->Fill(ye[0]*norm_params_[c].k + norm_params_[c].b);
        }
        if(yhm==1 && !yhch.empty() && yhch[0] < NUM_DSSDYH_POS) {
            int c = CH_BASE_YH + (int)yhch[0]; 
            vYH[s]->Fill(yhe[0]*norm_params_[c].k + norm_params_[c].b);
        }
    }, {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    for(auto h : vX) { h_tot_X->Add(h); delete h; }
    for(auto h : vY) { h_tot_Y->Add(h); delete h; }
    for(auto h : vYH){ h_tot_YH->Add(h); delete h; }
}

PlaneResult Calibrator::CalibratePlane(TH1D* h_raw, const char* name, TFile* diag, ofstream& out,
                                       const PlaneResult* fallback) {
    PlaneResult res;
    TH1D* h_sub = SubtractBackground(h_raw);

    // ---- 1. 逐参考峰定位，剔除"同一物理峰被两个参考点认领"的重复 ----
    // 典型场景：7065 与 7128 只隔 63 keV，而 YH 面 FWHM 已有 40 keV 量级，两个
    // ±60 keV 搜索窗几乎完全重叠；FindPeakGaussian 在窗内取"最高峰"，于是两个
    // 参考点落到同一个峰上，三点拟合退化成两点，pol1 会给出完全离谱的增益。
    vector<double> meas, ref, sigs;
    for(auto& p : ref_peaks_) {
        double pos, sig;
        if(!FindPeakGaussian(h_sub, p.energy, p.window * 2.0, pos, sig)) {
            cout << "    [SKIP] " << name << ": 找不到 " << p.energy
                 << " keV 附近的峰 (窗 ±" << p.window*2.0 << " keV)" << endl;
            continue;
        }
        bool dup = false;
        for(size_t j=0; j<meas.size(); ++j) {
            if(std::abs(pos - meas[j]) < REF_DUP_SIGMA * std::max(sig, sigs[j])) {
                cout << "    [SKIP] " << name << ": " << p.energy << " keV 与 "
                     << ref[j] << " keV 认到同一个峰 (pos=" << pos << " keV, sigma="
                     << sig << ")，该点无独立信息，丢弃" << endl;
                dup = true; break;
            }
        }
        if(dup) continue;
        meas.push_back(pos); ref.push_back(p.energy); sigs.push_back(sig);
    }
    delete h_sub;

    res.n_points = (int)meas.size();
    if(!ref.empty()) {
        auto mm = std::minmax_element(ref.begin(), ref.end());
        res.span = *mm.second - *mm.first;
    }

    // 把被接受的参考点画出来，便于人工核对（旧版只在成功时才写，现在总是写）
    if(!meas.empty()) {
        diag->cd();
        TGraph g_pts((int)meas.size(), meas.data(), ref.data());
        g_pts.SetTitle(TString::Format("%s: accepted reference peaks;measured (keV);nominal (keV)", name));
        g_pts.SetMarkerStyle(20);
        g_pts.Write(TString::Format("gr_calib_%s", name));
    }

    // ---- 2. 拟合 (K,B)，然后报告（不否决）----
    // 只有"峰太少拟合不出来"和"K 非物理"才算失败；K/B 偏离 1/0 一律只警告，
    // 因为本程序也用于处理原始 ADC（不做外刻度），那时 K/B 本来就该偏。
    bool accepted = false;
    if(res.n_points < REF_MIN_POINTS) {
        cout << "    [FALLBACK] " << name << ": 可用参考峰 " << res.n_points
             << " < " << REF_MIN_POINTS << "，直线都定不出来，无法独立定 K/B" << endl;
    } else {
        TGraph g(res.n_points, meas.data(), ref.data());
        TF1 f("f","pol1"); g.Fit(&f, "Q");
        double K = f.GetParameter(1), B = f.GetParameter(0);

        double max_resid = 0.0, ss = 0.0;
        vector<size_t> idx(meas.size());
        for(size_t j=0; j<idx.size(); ++j) idx[j] = j;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){ return ref[a] < ref[b]; });
        for(size_t j : idx) {
            double d = K*meas[j] + B - ref[j];
            ss += d*d;
            max_resid = std::max(max_resid, std::abs(d));
        }
        res.fit_rms = std::sqrt(ss / meas.size());

        // 局域斜率检查：相邻两参考点测得的间距应接近它们的标称间距。
        // 这一项从不否决拟合（点太少时噪声大），但会明确指出误差被哪一对点放大。
        for(size_t j=1; j<idx.size(); ++j) {
            double d_ref  = ref[idx[j]] - ref[idx[j-1]];
            double d_meas = meas[idx[j]] - meas[idx[j-1]];
            if(d_ref <= 0) continue;
            double k_local = d_ref / d_meas;
            if(std::abs(k_local/K - 1.0) > 0.03)
                cout << "    [WARN] " << name << ": 相邻参考峰 " << ref[idx[j-1]] << "->"
                     << ref[idx[j]] << " keV 的局域增益 " << k_local << " 比全局 " << K
                     << " 偏 " << 100.0*(k_local/K - 1.0) << "%（该对点间距仅 "
                     << d_ref << " keV）" << endl;
        }

        if(!std::isfinite(K) || !std::isfinite(B) || K <= 0.0) {
            cout << "    [ERROR] " << name << ": 拟合给出非物理结果 K=" << K << " B=" << B
                 << "，弃用这套刻度" << endl;
        } else {
            res.K = K; res.B = B; res.valid = true;
            accepted = true;
            cout << "    [OK] " << name << ": K=" << K << " B=" << B << " keV (n="
                 << res.n_points << ", 跨度 " << res.span << " keV, rms="
                 << res.fit_rms << " keV)" << endl;

            // 以下都只是提示：原始 ADC 输入时 K/B 偏离 1/0 是正常的
            if(res.n_points < REF_GOOD_POINTS)
                cout << "    [WARN] " << name << ": 只有 " << res.n_points
                     << " 个参考峰，没有冗余，无法用残差检验模型" << endl;
            if(res.span < REF_WARN_SPAN_KEV)
                cout << "    [WARN] " << name << ": 参考峰跨度仅 " << res.span << " keV < "
                     << REF_WARN_SPAN_KEV << " keV，斜率杠杆不足，K 的误差会被放大" << endl;
            if(max_resid > REF_WARN_RESID)
                cout << "    [WARN] " << name << ": 单点最大残差 " << max_resid << " keV > "
                     << REF_WARN_RESID << " keV，参考峰可能有位移或被别的峰污染" << endl;
            if(std::abs(K - 1.0) > REF_WARN_K_TOL || std::abs(B) > REF_WARN_B_TOL)
                cout << "    [WARN] " << name << ": K=" << K << " B=" << B
                     << " 偏离 (1,0) 较远 —— 若输入是未做外刻度的原始 ADC，这是正常的"
                     << "（结果照用）；若不是，检查参考峰是否对错了线" << endl;
        }
    }

    if(!accepted) {
        if(fallback) {
            res.K = fallback->K; res.B = fallback->B;
            res.valid = true; res.inherited = true;
            cout << "    [INHERIT] " << name << ": 沿用 X 面刻度 K=" << res.K
                 << " B=" << res.B << " keV" << endl;
        } else {
            return res;   // K=1, B=0, valid=false：交给下游按"未刻度"处理
        }
    }

    TH1D* h_calib = ApplyCalibration(h_raw, res.K, res.B, TString::Format("h_calib_%s", name));
    
    FindAndListPeaks(h_calib, TString::Format("%s-Plane Calibrated (keV)", name).Data(), out);

    TH1D* h_calib_sub = SubtractBackground(h_calib);
    
    double sum_sq_diff = 0.0;
    int count_diag = 0;
    vector<pair<double,double>> diag_res;   // (标称能量, 实测-标称)

    for(auto& p : diag_peaks_) {
        double pos, sig;
        if(FindPeakGaussian(h_calib_sub, p.energy, p.window, pos, sig)) {
            double fwhm = 2.355 * sig;
            if(fwhm > res.FWHM_Max) res.FWHM_Max = fwhm;

            double diff = pos - p.energy;
            sum_sq_diff += diff * diff;
            count_diag++;
            diag_res.push_back({p.energy, diff});
        }
    }

    if(count_diag > 0) {
        res.Res_RMS = std::sqrt(sum_sq_diff / count_diag);
    }

    // 把本面的刻度结论写进报告，方便和上一版逐行比较
    out << TString::Format("--- %s Absolute Scale ---", name).Data() << endl;
    out << TString::Format("K=%.6f  B=%.4f keV  n_ref=%d  span=%.1f keV  fit_rms=%.2f keV%s",
                           res.K, res.B, res.n_points, res.span, res.fit_rms,
                           res.inherited ? "  (inherited from X)" : "").Data() << endl;
    for(auto& d : diag_res)
        out << TString::Format("  diag %7.1f keV : residual %+.1f keV", d.first, d.second).Data() << endl;
    out << TString::Format("  diag rms = %.2f keV", res.Res_RMS).Data() << endl;

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
        double final_sig = -1.0;

        // 二次拟合必须带本底：这些谱在大峰下面压着一条下降的连续谱，
        // 用无本底的纯 gaus 会把峰位系统性地往低能拽 1 keV 量级。
        double fit_range = 30.0;
        double hgt = h->GetBinContent(h->GetXaxis()->FindBin(raw_pos));
        TF1 f("f_quick", "gaus(0)+pol1(3)", raw_pos - fit_range, raw_pos + fit_range);
        f.SetParameter(0, hgt);
        f.SetParameter(1, raw_pos);
        f.SetParameter(2, 12.0);
        f.SetParLimits(2, 2.0, 100.0);
        f.SetParameter(3, 1.0);
        f.SetParameter(4, 0.0);

        // 注意：必须带 "S" 选项，否则 TFitResultPtr 是空的，IsValid() 会崩。
        TFitResultPtr r = h->Fit(&f, "RQNS");
        if(int(r) == 0 && r->IsValid()) {
            final_pos = f.GetParameter(1);
            final_sig = std::abs(f.GetParameter(2));
        }

        peaks.push_back({final_pos, ypeaks[i], final_sig});
    }

    sort(peaks.begin(), peaks.end(), [](const PeakInfo& a, const PeakInfo& b) {
        return a.height > b.height;
    });

    out << "--- " << title << " ---" << endl;
    out << "Rank\tPos(Fitted)\tHeight(Raw)\tSigma\tFWHM" << endl;
    for(size_t i=0; i<peaks.size(); ++i) {
        out << (i+1) << "\t" << peaks[i].pos << "\t" << peaks[i].height << "\t";
        if(peaks[i].sigma > 0)
            out << peaks[i].sigma << "\t" << 2.355 * peaks[i].sigma;
        else
            out << "-\t-";
        out << endl;
    }
    out << "--------------------------------" << endl << endl;
}

TH1D* Calibrator::ApplyCalibration(TH1D* h_in, double K, double B, const char* new_name) {
    // 精确移轴：直接把每个 bin 的边界映射为 K*edge+B，bin 宽度和内容原样保留。
    // 旧实现按 bin 中心重新 Fill，会把 5 keV 的栅格量化带进结果（峰位最多偏移
    // 2.5 keV），也让 peaks_list 里 "Normalized" 与 "Calibrated" 两段的差值
    // 呈现出整 30.00 keV 的假象。
    int nb = h_in->GetNbinsX();
    vector<double> edges(nb + 1);
    for(int i = 0; i <= nb; ++i)
        edges[i] = K * h_in->GetXaxis()->GetBinLowEdge(i + 1) + B;

    TH1D* h_out = new TH1D(new_name, new_name, nb, edges.data());
    h_out->SetDirectory(0);
    h_out->Sumw2();
    for(int i = 1; i <= nb; ++i) {
        h_out->SetBinContent(i, h_in->GetBinContent(i));
        h_out->SetBinError(i, h_in->GetBinError(i));
    }
    h_out->SetEntries(h_in->GetEntries());
    h_out->GetXaxis()->SetTitle(h_in->GetXaxis()->GetTitle());
    h_out->GetYaxis()->SetTitle(h_in->GetYaxis()->GetTitle());
    return h_out;
}

bool Calibrator::FindPeakGaussian(TH1D* h, double center, double win, double &peakPos, double &sigma) {
    h->GetXaxis()->SetRange(0, 0); 
    // SetRangeUser 只改显示范围、不改 Fit 范围，但会把调用者拿到的直方图留在
    // "放大"状态，后面所有按 bin 索引的操作都会被带偏。用 RAII 保证任何出口都还原。
    struct RangeReset {
        TH1D* hh;
        ~RangeReset() { hh->GetXaxis()->SetRange(0, 0); }
    } range_reset{h};
    
    
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

void Calibrator::WriteOutput(const PlaneResult& rX, const PlaneResult& rY, const PlaneResult& rYH) {
    ofstream out(OUTPUT_DAT_FILE);
    out << fixed << setprecision(4);
    // 文件格式：第一行起就是 TOTAL_DSSD_CH 行数据，表头放在文件末尾的注释块里。
    //
    // 为什么不把表头放最上面：root2root 的 Read_Par_From_File()
    //   while(fgets(...)) { sscanf(str,"%d\t%lf\t%lf\t",&ch,&par[0],&par[1]);
    //                       if(ch<TOTAL_CH) {...} }     // 这里的 TOTAL_CH 是 setup.h 的 DAQ 道数
    // 不检查 sscanf 的返回值。表头解析失败时 ch/par 保持"上一行"的值；第一行就失败
    // 时它们是未初始化的栈内容，可能写坏任意一个道。把数据放在最前面，后面所有
    // '#' 行解析失败时 ch 已经是最后一道、par 也是它自己的值，重复赋值无害。
    auto w = [&](int start, int end, const PlaneResult& pr) {
        for(int i=start; i<end; ++i) {
            double fk = pr.K * norm_params_[i].k;
            double fb = pr.K * norm_params_[i].b + pr.B;

            double fwhm = strip_fwhms_[i];
            double dE   = strip_dE_[i];

            // 前 3 列 (Ch 截距 斜率) 是 root2root 读的字段，顺序不能动
            out << i << "\t" << fb << "\t" << fk << "\t" << fwhm << "\t" << dE << endl;

            h_fwhm_summary->SetBinContent(i+1, fwhm);
        }
    };

    w(CH_BASE_X,  CH_BASE_Y,  rX);
    w(CH_BASE_Y,  CH_BASE_YH, rY);
    w(CH_BASE_YH, TOTAL_DSSD_CH, rYH);

    out << "# ------------------------------------------------------------" << endl;
    // 流里是 fixed/4 位，这个标签要 %g；注意 TString 必须活到函数结束，不能取 Data() 存指针
    TString e_tstr = TString::Format("%g", strip_peak_E_);
    const char* e_str = e_tstr.Data();
    out << "# 列: Ch  b(截距,keV)  k(斜率)  FWHM(keV)  dE_" << e_str << "(keV)" << endl;
    out << "# 前 3 列供 root2root Get_Calib_Energy() 使用: E = b + k*ADC (IntSour_Corr=1,0)" << endl;
    out << "# dE_" << e_str << " = 该条在 " << e_str
        << " keV 处的实测峰位减标称值；0 表示没测到峰" << endl;
    out << "# 平面级刻度结论（K = 增益, B = 零点, 单位 keV）" << endl;
    auto note = [&](const char* n, int lo, int hi, const PlaneResult& pr) {
        out << "# " << n << " (" << lo << "-" << (hi-1) << "): K=" << pr.K << "  B=" << pr.B;
        if(pr.inherited)       out << "  [沿用 X 面刻度]";
        else if(pr.valid)      out << "  [本面拟合: " << pr.n_points << " 点, 跨度 "
                                    << TString::Format("%g", pr.span) << " keV, 残差 rms "
                                    << TString::Format("%.2f", pr.fit_rms) << " keV]";
        else                   out << "  [未刻度: 参考峰不足或被否决, K=1 B=0]";
        out << "  诊断残差 rms=" << pr.Res_RMS << " keV, FWHM_max=" << pr.FWHM_Max << " keV";
        int nf = 0;
        for(int i=lo; i<hi; ++i) { if(strip_fwhms_[i] > 0) ++nf; }
        out << "  (" << nf << "/" << (hi-lo) << " 条测到 " << e_str << " 峰)" << endl;
    };
    note("X ", CH_BASE_X,  CH_BASE_Y,      rX);
    note("Y ", CH_BASE_Y,  CH_BASE_YH,     rY);
    note("YH", CH_BASE_YH, TOTAL_DSSD_CH,  rYH);
    out << "# ------------------------------------------------------------" << endl;
    out.close();

    cout << "--> 刻度系数已写入 " << OUTPUT_DAT_FILE << endl;
}