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

using namespace std;
using namespace ROOT;

struct NormResult {
    double k = 1.0;
    double b = 0.0;
};

int main() {
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    cout << "=== Step 1: DSSD Normalization (With Bridge Correction) ===" << endl;

    // 1. 读取文件
    TChain chain(TREE_NAME);
    // [修改] 这里改用 IMPLANT_ONLY_FILE，匹配 Config.h 的修改
    if (!gSystem->AccessPathName(IMPLANT_ONLY_FILE)) chain.Add(IMPLANT_ONLY_FILE);
    if (!gSystem->AccessPathName(DECAY_ONLY_FILE)) chain.Add(DECAY_ONLY_FILE);

    if (chain.GetListOfFiles()->GetEntries() == 0) {
        cerr << "[ERROR] No input files found (Check if Step 0 generated them)." << endl;
        return 1;
    }

    ROOT::RDataFrame df(chain);
    cout << "--> Total Entries for Normalization: " << df.Count().GetValue() << endl;

    // 2. 筛选 (X/Y均触发，无反符合，能量在范围内，且X-Y差值不大)
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
    vector<vector<TH2D*>> slots_X(nSlots, vector<TH2D*>(NUM_DSSDX_POS, nullptr));
    vector<vector<TH2D*>> slots_Y(nSlots, vector<TH2D*>(NUM_DSSDY_POS, nullptr));

    for(unsigned int s=0; s<nSlots; ++s) {
        for(int i=0; i<NUM_DSSDX_POS; ++i) {
            slots_X[s][i] = new TH2D(Form("h2_X%d_s%d", i, s), "", 200, NORM_MIN_E, NORM_MAX_E, 200, NORM_MIN_E, NORM_MAX_E);
            slots_X[s][i]->SetDirectory(0);
        }
        for(int i=0; i<NUM_DSSDY_POS; ++i) {
            slots_Y[s][i] = new TH2D(Form("h2_Y%d_s%d", i, s), "", 200, NORM_MIN_E, NORM_MAX_E, 200, NORM_MIN_E, NORM_MAX_E);
            slots_Y[s][i]->SetDirectory(0);
        }
    }

    // 填充逻辑：
    // X条校准：看它和 Ref_Y 的符合 (X_i vs Ref_Y)
    // Y条校准：看它和 Ref_X 的符合 (Y_i vs Ref_X)
    df_norm.ForeachSlot([&](unsigned int slot, int cx, double ex, int cy, double ey){
        if(cy == ref_Y && cx >= 0 && cx < NUM_DSSDX_POS) slots_X[slot][cx]->Fill(ex, ey);
        if(cx == ref_X && cy >= 0 && cy < NUM_DSSDY_POS) slots_Y[slot][cy]->Fill(ey, ex);
    }, {"ChX", "EX", "ChY", "EY"});

    // 5. 合并直方图并拟合
    cout << "--> Merging and Fitting..." << endl;
    TFile *f_diag = new TFile(NORM_DIAG_ROOT, "RECREATE");
    h_x_hits->Write(); h_y_hits->Write();

    vector<NormResult> x_res(NUM_DSSDX_POS), y_res(NUM_DSSDY_POS);
    vector<int> bad_x, bad_y;

    // --- 处理 X 面 ---
    for(int i=0; i<NUM_DSSDX_POS; ++i) {
        TH2D* h_final = slots_X[0][i];
        for(unsigned int s=1; s<nSlots; ++s) {
            h_final->Add(slots_X[s][i]);
            delete slots_X[s][i];
        }
        h_final->SetName(Form("h2_X_%d", i));
        h_final->SetTitle(Form("X%d vs RefY%d;EX;RefEY", i, ref_Y));
        
        // 确保 Ref_X 本身能被拟合（放宽一点统计量要求）
        int min_entries = (i == ref_X) ? 100 : NORM_MIN_ENTRIES;

        if(h_final->GetEntries() < min_entries) {
            bad_x.push_back(i);
        } else {
            // ProfileX 转为 1D 并拟合直线
            TProfile *p = h_final->ProfileX(Form("pX_%d", i));
            TFitResultPtr r = p->Fit("pol1", "QS0");
            if(r->IsValid() && abs(r->Parameter(1)) > 0.1) {
                // 存储参数: E_cal = k * E_raw + b
                x_res[i] = {r->Parameter(1), r->Parameter(0)};
            } else {
                bad_x.push_back(i);
            }
            p->Write(); delete p;
        }
        h_final->Write(); delete h_final;
    }

    // --- 处理 Y 面 ---
    for(int i=0; i<NUM_DSSDY_POS; ++i) {
        TH2D* h_final = slots_Y[0][i];
        for(unsigned int s=1; s<nSlots; ++s) {
            h_final->Add(slots_Y[s][i]);
            delete slots_Y[s][i];
        }
        h_final->SetName(Form("h2_Y_%d", i));
        h_final->SetTitle(Form("Y%d vs RefX%d;EY;RefEX", i, ref_X));

        if(h_final->GetEntries() < NORM_MIN_ENTRIES) {
            bad_y.push_back(i);
        } else {
            TProfile *p = h_final->ProfileX(Form("pY_%d", i));
            TFitResultPtr r = p->Fit("pol1", "QS0");
            if(r->IsValid() && abs(r->Parameter(1)) > 0.1) {
                y_res[i] = {r->Parameter(1), r->Parameter(0)};
            } else {
                bad_y.push_back(i);
            }
            p->Write(); delete p;
        }
        h_final->Write(); delete h_final;
    }

    f_diag->Close();
    delete f_diag;

    // 6. 桥接修正 (Bridge Correction)
    // 目前：所有 X 条都对齐到了 Ref_Y，所有 Y 条都对齐到了 Ref_X。
    // 目标：把所有 Y 条的参数转换一下，让它们也对齐到 Ref_Y (即 X面的标准)。
    // 桥梁：Ref_X 和 Ref_Y 的交叉点。
    
    cout << "------------------------------------------------" << endl;
    cout << "--> Checking Bridge Correction (Ref_X <-> Ref_Y)..." << endl;

    double k_bridge = 1.0;
    double b_bridge = 0.0;
    bool bridge_ok = true;

    // 检查 Ref_X 是否在坏条列表中
    for(int bad : bad_x) { 
        if(bad == ref_X) {
            bridge_ok = false;
            break;
        }
    }

    if(bridge_ok && abs(x_res[ref_X].k) > 0.01) {  
        // x_res[ref_X] 存的是: E_refX_cal = k * E_refX_raw + b = E_refY_raw
        // 这里的关系稍微有点绕，本质是利用 Ref_X 在 "X vs RefY" 图里的拟合参数
        // 实际上 x_res[ref_X] 描述了: Ref_X 读数 -> Ref_Y 读数 的关系
        
        // 我们需要把 Y面的 "对齐RefX" 转为 "对齐RefY"
        // 已知 Y[i]: E_Yi_cal = k_y * E_Yi_raw + b_y = E_refX_raw
        // 桥梁关系 (由 x_res[ref_X] 提供): E_refX_raw = (E_refY_raw - b_br) / k_br
        // ...这里采用直接数学倒推的逻辑:
        
        // 修正逻辑：
        // 1. 我们拟合得到 X_i vs Ref_Y => X 对齐到了 Ref_Y
        // 2. 我们拟合得到 Y_i vs Ref_X => Y 对齐到了 Ref_X
        // 3. 我们需要 Y_i vs Ref_Y
        // 4. 关键点：在 X_i vs Ref_Y 的循环中，当 i = Ref_X 时，
        //    我们得到了 Ref_X vs Ref_Y 的关系: Ref_X_val = k_ref * Ref_X_raw + b_ref ~ Ref_Y_raw
        //    Wait, 这里的 x_res[i] 含义是：将 Ch_i 的原始值 映射到 参考轴的值。
        //    对于 slots_X[ref_X], X轴是 Ref_X_raw, Y轴是 Ref_Y_raw
        //    所以 x_res[ref_X] 确实是: Ref_Y_raw = k * Ref_X_raw + b
        
        //    对于 slots_Y[i], X轴是 Y_i_raw, Y轴是 Ref_X_raw
        //    y_res[i] 是: Ref_X_raw = k_y * Y_i_raw + b_y
        
        //    联立：
        //    Ref_Y_raw = k_bridge * (k_y * Y_i_raw + b_y) + b_bridge
        //              = (k_bridge * k_y) * Y_i_raw + (k_bridge * b_y + b_bridge)
        
        //    所以修正系数很简单：
        k_bridge = x_res[ref_X].k;
        b_bridge = x_res[ref_X].b;
        
        cout << "    [OK] Bridge parameters (Ref_X -> Ref_Y): k=" << k_bridge << ", b=" << b_bridge << endl;

        for(int i=0; i<NUM_DSSDY_POS; ++i) {
            bool is_bad = false;
            for(int bad : bad_y) if(bad == i) { is_bad = true; break; }
            if(is_bad) continue;

            double k_old = y_res[i].k;
            double b_old = y_res[i].b;
            
            // 应用修正公式
            y_res[i].k = k_bridge * k_old;
            y_res[i].b = k_bridge * b_old + b_bridge;
        }
        cout << "--> Y-plane parameters have been aligned to Ref_Y scale." << endl;

    } else {
        cout << "    [WARNING] Bridge calculation failed (Ref_X bad or k too small)." << endl;
        cout << "    Y parameters will NOT be corrected and may mismatch X scale." << endl;
    }
    cout << "------------------------------------------------" << endl;

    // 7. 输出结果
    ofstream out(NORM_PARAM_FILE, ios::trunc); 
    if(!out.is_open()) {
        cerr << "[ERROR] Cannot open output file: " << NORM_PARAM_FILE << endl;
        return 1;
    }

    out << fixed << setprecision(6);
    out << "# ID\tk\t\tb" << endl;
    for (int i : bad_x) out << "# Bad X: " << i << endl;
    for (int i : bad_y) out << "# Bad Y: " << i << endl;
    
    // 输出
    for (int i = 0; i < NUM_DSSDX_POS; ++i) out << i << "\t" << x_res[i].k << "\t" << x_res[i].b << endl;
    for (int i = 0; i < NUM_DSSDY_POS; ++i) out << (128+i) << "\t" << y_res[i].k << "\t" << y_res[i].b << endl;
    
    out.close();
    cout << "--> Normalization Parameters saved to " << NORM_PARAM_FILE << endl;
    return 0;
}