#include "Config.h"
#include "TFile.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TChain.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TFitResult.h"
#include "ROOT/RDataFrame.hxx"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace ROOT;

struct NormResult {
    double k = 1.0;
    double b = 0.0;
};

int main() {
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    cout << "=== Step 1: DSSD Normalization (With Bridge Correction) ===" << endl;
    Print_Channel_Map();

    // 1. 读取文件
    TChain chain(TREE_NAME);
    if (!gSystem->AccessPathName(IMPLANT_ONLY_FILE)) chain.Add(IMPLANT_ONLY_FILE);
    if (!gSystem->AccessPathName(DECAY_ONLY_FILE))   chain.Add(DECAY_ONLY_FILE);

    if (chain.GetListOfFiles()->GetEntries() == 0) {
        cerr << "[ERROR] No input files found (Check if Step 0 generated them)." << endl;
        return 1;
    }

    ROOT::RDataFrame df(chain);
    cout << "--> Total Entries for Normalization: " << df.Count().GetValue() << endl;

    // 2. 筛选 (X/Y均单条触发，无反符合，能量在范围内，且X-Y差值不大)
    // 说明：abs(DSSDX_E-DSSDY_E)<NORM_E_DIFF 这一刀切在"待测量的残差"上，严格说会
    // 引入截断偏置；实测该像素的残差分布 rms=28 keV，±200 keV 已是 ±7 sigma，
    // 实际几乎不切事件。保留它是为了挡住明显不匹配的本底（低能尾巴/噪声）。
    auto df_norm = df.Filter(
        Form("DSSDX_mul==1 && DSSDY_mul==1 && Veto_mul==0 && SSD_mul==0 && "
             "DSSDX_E[0]>%f && DSSDX_E[0]<%f && "
             "DSSDY_E[0]>%f && DSSDY_E[0]<%f && "
             "abs(DSSDX_E[0]-DSSDY_E[0])<%f",
             NORM_MIN_E, NORM_MAX_E,
             NORM_MIN_E, NORM_MAX_E,
             NORM_E_DIFF)
    ).Define("ChX", "(int)DSSDX_Ch[0]")
     .Define("EX",  "DSSDX_E[0]")
     .Define("ChY", "(int)DSSDY_Ch[0]")
     .Define("EY",  "DSSDY_E[0]");

    // 3. 寻找参考道 (Ref Strips) - 统计量最大的条
    cout << "--> Finding Reference Strips..." << endl;
    auto h_x_hits = df_norm.Histo1D({"h_x_hits", "X hits", NUM_DSSDX_POS, -0.5, NUM_DSSDX_POS-0.5}, "ChX");
    auto h_y_hits = df_norm.Histo1D({"h_y_hits", "Y hits", NUM_DSSDY_POS, -0.5, NUM_DSSDY_POS-0.5}, "ChY");

    // 如果 Config.h 里定义了手动参考条，可以在这里覆盖，否则自动寻找
    int ref_X = h_x_hits->GetMaximumBin() - 1;
    int ref_Y = h_y_hits->GetMaximumBin() - 1;

    #ifdef REF_STRIP_X
        ref_X = REF_STRIP_X;
    #endif
    #ifdef REF_STRIP_Y
        ref_Y = REF_STRIP_Y;
    #endif

    cout << "--> Ref X: " << ref_X << ", Ref Y: " << ref_Y << endl;

    // 4. 填充二维直方图 (多线程安全方式)
    cout << "--> Filling Histograms in Single Pass..." << endl;
    unsigned int nSlots = df.GetNSlots();
    vector<vector<TH2D*>> slots_X (nSlots, vector<TH2D*>(NUM_DSSDX_POS, nullptr));
    vector<vector<TH2D*>> slots_Y (nSlots, vector<TH2D*>(NUM_DSSDY_POS, nullptr));
    vector<vector<TH2D*>> slots_YH(nSlots, vector<TH2D*>(NUM_DSSDYH_POS, nullptr));

    for(unsigned int s=0; s<nSlots; ++s) {
        for(int i=0; i<NUM_DSSDX_POS; ++i) {
            slots_X[s][i] = new TH2D(Form("h2_X%d_s%d", i, s), "", 200, NORM_MIN_E, NORM_MAX_E, 200, NORM_MIN_E, NORM_MAX_E);
            slots_X[s][i]->SetDirectory(0);
        }
        for(int i=0; i<NUM_DSSDY_POS; ++i) {
            slots_Y[s][i] = new TH2D(Form("h2_Y%d_s%d", i, s), "", 200, NORM_MIN_E, NORM_MAX_E, 200, NORM_MIN_E, NORM_MAX_E);
            slots_Y[s][i]->SetDirectory(0);
        }
        // YH 面（Y 面并行读出，全局道 CH_BASE_YH 起）。旧版完全没有这一步，
        // 于是 Step 2 只能给 YH 用 (k,b)=(1,0)，YH 的谱其实是原始幅度道。
        for(int i=0; i<NUM_DSSDYH_POS; ++i) {
            slots_YH[s][i] = new TH2D(Form("h2_YH%d_s%d", i, s), "", 200, NORM_MIN_E, NORM_MAX_E, 200, NORM_MIN_E, NORM_MAX_E);
            slots_YH[s][i]->SetDirectory(0);
        }
    }

    // 填充逻辑：
    // X条校准：看它和 Ref_Y 的符合 (X_i vs Ref_Y)
    // Y条校准：看它和 Ref_X 的符合 (Y_i vs Ref_X)
    // YH条同上，也用 Ref_X 这一列
    df_norm.ForeachSlot([&](unsigned int slot, int cx, double ex, int cy, double ey,
                            const RVec<double>& yhe, const RVec<double>& yhch, UShort_t yhm){
        if(cy == ref_Y && cx >= 0 && cx < NUM_DSSDX_POS) slots_X[slot][cx]->Fill(ex, ey);
        if(cx == ref_X && cy >= 0 && cy < NUM_DSSDY_POS) slots_Y[slot][cy]->Fill(ey, ex);
        if(cx == ref_X && yhm == 1 && !yhch.empty()) {
            int cyh = (int)yhch[0];
            if(cyh >= 0 && cyh < NUM_DSSDYH_POS) slots_YH[slot][cyh]->Fill(yhe[0],  ex);
        }
    }, {"ChX", "EX", "ChY", "EY", "DSSDYH_E", "DSSDYH_Ch", "DSSDYH_mul"});

    // 5. 合并直方图并拟合
    cout << "--> Merging and Fitting..." << endl;
    TFile *f_diag = new TFile(NORM_DIAG_ROOT, "RECREATE");
    f_diag->cd();
    h_x_hits->Write(); h_y_hits->Write();

    vector<NormResult> x_res(NUM_DSSDX_POS), y_res(NUM_DSSDY_POS), yh_res(NUM_DSSDYH_POS);
    vector<int> bad_x, bad_y, bad_yh;

    // 合并各 slot 的同一条直方图
    auto merge_slots = [&](vector<vector<TH2D*>>& slots, int i) -> TH2D* {
        TH2D* h = slots[0][i];
        for(unsigned int s=1; s<nSlots; ++s) {
            h->Add(slots[s][i]);
            delete slots[s][i];
        }
        return h;
    };

    // 单条：ProfileX 后拟合 pol1。返回值即"该条是否可用"。
    // 增益合理性检查从 |k|>0.1 收紧到 0.5<k<2：一个 5% 的增益差才是物理上可能的，
    // 旧判据会让 k=0.3 这类明显失败的拟合也通过。
    auto fit_strip = [&](TH2D* h, const char* prof_name, int min_entries, NormResult& res) -> bool {
        if(h->GetEntries() < min_entries) return false;
        TProfile* p = h->ProfileX(prof_name);
        TFitResultPtr r = p->Fit("pol1", "QS0");
        bool ok = r->IsValid() &&
                  std::abs(r->Parameter(1)) > 0.5 && std::abs(r->Parameter(1)) < 2.0;
        if(ok) res = {r->Parameter(1), r->Parameter(0)};
        p->Write(); delete p;
        return ok;
    };

    // --- 处理 X 面 ---
    for(int i=0; i<NUM_DSSDX_POS; ++i) {
        TH2D* h_final = merge_slots(slots_X, i);
        h_final->SetName(Form("h2_X_%d", i));
        h_final->SetTitle(Form("X%d vs RefY%d;EX;RefEY", i, ref_Y));

        // 确保 Ref_X 本身能被拟合（放宽一点统计量要求）
        int min_entries = (i == ref_X) ? 100 : NORM_MIN_ENTRIES;
        if(!fit_strip(h_final, Form("pX_%d", i), min_entries, x_res[i])) bad_x.push_back(i);

        h_final->Write(); delete h_final;
    }

    // --- 处理 Y 面 ---
    for(int i=0; i<NUM_DSSDY_POS; ++i) {
        TH2D* h_final = merge_slots(slots_Y, i);
        h_final->SetName(Form("h2_Y_%d", i));
        h_final->SetTitle(Form("Y%d vs RefX%d;EY;RefEX", i, ref_X));

        if(!fit_strip(h_final, Form("pY_%d", i), NORM_MIN_ENTRIES, y_res[i])) bad_y.push_back(i);

        h_final->Write(); delete h_final;
    }

    // --- 处理 YH 面 ---
    for(int i=0; i<NUM_DSSDYH_POS; ++i) {
        TH2D* h_final = merge_slots(slots_YH, i);
        h_final->SetName(Form("h2_YH_%d", i));
        h_final->SetTitle(Form("YH%d vs RefX%d;YH_E;RefEX", i, ref_X));

        if(!fit_strip(h_final, Form("pYH_%d", i), NORM_MIN_ENTRIES, yh_res[i])) bad_yh.push_back(i);

        h_final->Write(); delete h_final;
    }

    f_diag->Close();
    delete f_diag;

    cout << "--> Fitted strips: X " << (NUM_DSSDX_POS-bad_x.size()) << "/" << NUM_DSSDX_POS
         << ", Y "  << (NUM_DSSDY_POS-bad_y.size()) << "/" << NUM_DSSDY_POS
         << ", YH " << (NUM_DSSDYH_POS-bad_yh.size()) << "/" << NUM_DSSDYH_POS << endl;

    // 6. 桥接修正 (Bridge Correction)
    // 目前：所有 X 条都对齐到了 Ref_Y，所有 Y/YH 条都对齐到了 Ref_X。
    // 目标：把 Y/YH 的参数转换一下，让它们也对齐到 Ref_Y (即 X 面的标准)。
    // 桥梁：Ref_X 和 Ref_Y 的交叉点。
    cout << "------------------------------------------------" << endl;
    cout << "--> Checking Bridge Correction (Ref_X <-> Ref_Y)..." << endl;

    double k_bridge = 1.0;
    double b_bridge = 0.0;
    bool bridge_ok = true;

    for(int bad : bad_x) {
        if(bad == ref_X) { bridge_ok = false; break; }
    }

    if(bridge_ok && x_res[ref_X].k > 0.5 && x_res[ref_X].k < 2.0) {
        // x_res[i] 的含义：把 Ch_i 的原始读数映射到参考轴
        //   X 面: E_refY_raw = k * E_Xi_raw + b
        //   Y 面: E_refX_raw = k_y * E_Yi_raw + b_y
        // 桥（i = Ref_X 那一行）：E_refY_raw = k_br * E_refX_raw + b_br
        // 联立得：E_refY_raw = (k_br*k_y)*E_Yi_raw + (k_br*b_y + b_br)
        k_bridge = x_res[ref_X].k;
        b_bridge = x_res[ref_X].b;

        cout << "    [OK] Bridge parameters (Ref_X -> Ref_Y): k=" << k_bridge << ", b=" << b_bridge << endl;

        auto apply_bridge = [&](vector<NormResult>& res, const vector<int>& bad) {
            for(int i=0; i<(int)res.size(); ++i) {
                bool is_bad = false;
                for(int b : bad) if(b == i) { is_bad = true; break; }
                if(is_bad) continue;
                double k_old = res[i].k, b_old = res[i].b;
                res[i].k = k_bridge * k_old;
                res[i].b = k_bridge * b_old + b_bridge;
            }
        };
        apply_bridge(y_res,  bad_y);
        apply_bridge(yh_res, bad_yh);
        cout << "--> Y-plane (and YH) parameters have been aligned to Ref_Y scale." << endl;

    } else {
        cout << "    [WARNING] Bridge calculation failed (Ref_X bad or k too small)." << endl;
        cout << "    Y parameters will NOT be corrected and may mismatch X scale." << endl;
    }
    cout << "------------------------------------------------" << endl;

    // 6b. 坏条处理
    // 旧版把坏条留在 (k,b)=(1,0)，这等于宣布"这条的增益和参考条一模一样"，
    // 而真实偏差是 0.5% 量级（几十 keV）。这里改成用左右最近好条按距离线性插值。
    auto interpolate = [&](vector<NormResult>& res, const vector<int>& bad, const char* face) {
        if(bad.empty()) return;
        vector<char> isbad(res.size(), 0);
        for(int i : bad) isbad[i] = 1;

        vector<int> good;
        for(int i=0; i<(int)res.size(); ++i) if(!isbad[i]) good.push_back(i);
        if(good.empty()) {
            cerr << "[WARN] " << face << " 面没有可用条，坏条无法插值" << endl;
            return;
        }

        for(int i : bad) {
            int lo = -1, hi = -1;
            for(int g : good) { if(g < i) lo = g; else { hi = g; break; } }
            if(lo < 0)       res[i] = res[hi];
            else if(hi < 0)  res[i] = res[lo];
            else {
                double w = double(hi - i) / double(hi - lo);   // lo 的权重
                res[i].k = w*res[lo].k + (1.0-w)*res[hi].k;
                res[i].b = w*res[lo].b + (1.0-w)*res[hi].b;
            }
        }
        cout << "--> " << face << ": " << bad.size() << " 条用邻道插值补齐" << endl;
    };
    interpolate(x_res,  bad_x,  "X");
    interpolate(y_res,  bad_y,  "Y");
    interpolate(yh_res, bad_yh, "YH");

    // 7. 输出结果
    ofstream out(NORM_PARAM_FILE, ios::trunc);
    if(!out.is_open()) {
        cerr << "[ERROR] Cannot open output file: " << NORM_PARAM_FILE << endl;
        return 1;
    }

    out << fixed << setprecision(6);
    out << "# DSSD normalization: E_norm = k * E_raw + b" << endl;
    out << "# " << CH_BASE_X << "-" << (CH_BASE_Y-1)  << " = DSSDX, "
        << CH_BASE_Y << "-" << (CH_BASE_YH-1) << " = DSSDY, "
        << CH_BASE_YH << "-" << (TOTAL_DSSD_CH-1) << " = DSSDYH (Y 面并行读出)" << endl;
    out << "# ID\tk\t\tb" << endl;
    for (int i : bad_x)  out << "# Bad X: "  << i << endl;
    for (int i : bad_y)  out << "# Bad Y: "  << i << endl;
    for (int i : bad_yh) out << "# Bad YH: " << i << endl;
    out << "# 坏条已用左右最近好条按距离线性插值，不是 (1,0)" << endl;

    for (int i = 0; i < NUM_DSSDX_POS;  ++i) out << (CH_BASE_X +i) << "\t" << x_res[i].k  << "\t" << x_res[i].b  << endl;
    for (int i = 0; i < NUM_DSSDY_POS;  ++i) out << (CH_BASE_Y +i) << "\t" << y_res[i].k  << "\t" << y_res[i].b  << endl;
    for (int i = 0; i < NUM_DSSDYH_POS; ++i) out << (CH_BASE_YH+i) << "\t" << yh_res[i].k << "\t" << yh_res[i].b << endl;

    out.close();
    cout << "--> Normalization Parameters saved to " << NORM_PARAM_FILE << endl;
    return 0;
}
