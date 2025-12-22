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
    if (!gSystem->AccessPathName(INJECT_ONLY_FILE)) chain.Add(INJECT_ONLY_FILE);
    if (!gSystem->AccessPathName(DECAY_ONLY_FILE)) chain.Add(DECAY_ONLY_FILE);

    if (chain.GetListOfFiles()->GetEntries() == 0) {
        cerr << "[ERROR] No input files found." << endl;
        return 1;
    }

    ROOT::RDataFrame df(chain);
    cout << "--> Total Entries: " << df.Count().GetValue() << endl;

    // 2. 筛选
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

    // 3. 寻找参考道
    cout << "--> Finding Reference Strips..." << endl;
    auto h_x_hits = df_norm.Histo1D({"h_x_hits", "X hits", NUM_DSSDX_POS, -0.5, NUM_DSSDX_POS-0.5}, "ChX");
    auto h_y_hits = df_norm.Histo1D({"h_y_hits", "Y hits", NUM_DSSDY_POS, -0.5, NUM_DSSDY_POS-0.5}, "ChY");
    
    int ref_X = h_x_hits->GetMaximumBin() - 1;
    int ref_Y = h_y_hits->GetMaximumBin() - 1;
    cout << "--> Ref X: " << ref_X << ", Ref Y: " << ref_Y << endl;

    // 4. 填充直方图
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

    df_norm.ForeachSlot([&](unsigned int slot, int cx, double ex, int cy, double ey){
        if(cy == ref_Y && cx >= 0 && cx < NUM_DSSDX_POS) slots_X[slot][cx]->Fill(ex, ey);
        if(cx == ref_X && cy >= 0 && cy < NUM_DSSDY_POS) slots_Y[slot][cy]->Fill(ey, ex);
    }, {"ChX", "EX", "ChY", "EY"});

    // 5. 合并与拟合
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
        
        // 注意：这里我们放宽一点对 ref_X 的统计量要求，防止桥接计算被跳过
        int min_entries = (i == ref_X) ? 100 : NORM_MIN_ENTRIES;

        if(h_final->GetEntries() < min_entries) {
            bad_x.push_back(i);
        } else {
            TProfile *p = h_final->ProfileX(Form("pX_%d", i));
            TFitResultPtr r = p->Fit("pol1", "QS0");
            if(r->IsValid() && abs(r->Parameter(1)) > 0.1) {
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

if(bridge_ok && abs(x_res[ref_X].k) > 0.01) {  // 防止除零
    // x_res[ref_X] 存储的是: E_refX = k * E_refY + b
    // 我们需要反函数: E_refY = (E_refX - b) / k
    // 即: E_refY = (1/k) * E_refX + (-b/k)
    
    k_bridge = 1.0 / x_res[ref_X].k;
    b_bridge = -x_res[ref_X].b / x_res[ref_X].k;
    
    cout << "    [OK] Bridge calculated from intersection pixel (X" << ref_X << ", Y" << ref_Y << ")." << endl;
    cout << "    Original X[" << ref_X << "]: k=" << x_res[ref_X].k << ", b=" << x_res[ref_X].b << endl;
    cout << "    Inverted Bridge: k=" << k_bridge << ", b=" << b_bridge << endl;
    
    // 打印一个 Y 条修正前的样子,用于调试
    int sample_y = ref_Y;
    cout << "    [Debug] Y" << sample_y << " BEFORE: k=" << y_res[sample_y].k 
         << ", b=" << y_res[sample_y].b << endl;

    // 应用修正:将 Y 面参数从 "对齐 Ref_X" 转换到 "对齐 Ref_Y"
    // 原参数: E_Y = k_old * E_refX + b_old
    // 已知: E_refY = k_bridge * E_refX + b_bridge
    // 反解: E_refX = (E_refY - b_bridge) / k_bridge
    // 代入: E_Y = k_old * [(E_refY - b_bridge) / k_bridge] + b_old
    //          = (k_old / k_bridge) * E_refY + (b_old - k_old * b_bridge / k_bridge)
    
    for(int i=0; i<NUM_DSSDY_POS; ++i) {
        // 跳过坏条
        bool is_bad = false;
        for(int bad : bad_y) {
            if(bad == i) {
                is_bad = true;
                break;
            }
        }
        if(is_bad) continue;

        double k_old = y_res[i].k;
        double b_old = y_res[i].b;
        
        // 正确的转换公式
        y_res[i].k = k_old / k_bridge;
        y_res[i].b = b_old - k_old * b_bridge / k_bridge;
    }
    
    cout << "    [Debug] Y" << sample_y << " AFTER : k=" << y_res[sample_y].k 
         << ", b=" << y_res[sample_y].b << endl;
    cout << "--> Y-plane parameters have been aligned to Ref_Y scale." << endl;

} else {
    cout << "    [WARNING] Bridge calculation failed!" << endl;
    if(!bridge_ok) {
        cout << "    Reason: Ref_X (strip " << ref_X << ") was marked as BAD." << endl;
    } else {
        cout << "    Reason: k parameter too small (k=" << x_res[ref_X].k << ")." << endl;
    }
    cout << "    Y parameters will NOT be corrected and may mismatch X scale." << endl;
}
cout << "------------------------------------------------" << endl;

    // 6. 输出结果
    // 确保这里的 ofstream 会覆盖旧文件
    ofstream out(NORM_PARAM_FILE, ios::trunc); 
    if(!out.is_open()) {
        cerr << "[ERROR] Cannot open output file: " << NORM_PARAM_FILE << endl;
        return 1;
    }

    out << fixed << setprecision(6);
    out << "# ID\tk\t\tb" << endl;
    for (int i : bad_x) out << "# Bad X: " << i << endl;
    for (int i : bad_y) out << "# Bad Y: " << i << endl;
    
    // 输出修正后的参数
    for (int i = 0; i < NUM_DSSDX_POS; ++i) out << i << "\t" << x_res[i].k << "\t" << x_res[i].b << endl;
    for (int i = 0; i < NUM_DSSDY_POS; ++i) out << (128+i) << "\t" << y_res[i].k << "\t" << y_res[i].b << endl;
    
    out.close();
    cout << "--> Normalization Parameters saved to " << NORM_PARAM_FILE << endl;
    return 0;
}