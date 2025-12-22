#include "Config.h"
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT;

int main() {
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    cout << "=== Step 0: Preselect decay/implant events ===" << endl;

    TChain chain(TREE_NAME);
    for (int run = RUN_START; run <= RUN_END; ++run) {
        TString fname = TString::Format(INPUT_DIR_PATTERN, run);
        if (!gSystem->AccessPathName(fname)) {
            chain.Add(fname);
        }
    }

    if (chain.GetListOfFiles()->GetEntries() == 0) {
        cerr << "[ERROR] No input files found." << endl;
        return 1;
    }

    ROOT::RDataFrame df(chain);
    auto n_total = df.Count();
    cout << "--> Total entries in raw data: " << n_total.GetValue() << endl;

    vector<string> columns = {
        "DSSDX_Ch", "DSSDX_E", "DSSDX_mul",
        "DSSDY_Ch", "DSSDY_E", "DSSDY_mul",
        "DSSDYH_Ch","DSSDYH_E","DSSDYH_mul",
        "MWPC_mul", "Veto_mul", "SSD_mul"
    };

    // [修改] 注入/植入事件筛选 (Implant)
    auto df_implant = df.Filter(
        "MWPC_mul > 0 && DSSDX_mul == 1 && DSSDY_mul == 1 && "
        "Veto_mul == 0 && SSD_mul == 0"
    );

    // 衰变事件筛选
    auto df_decay = df.Filter(
        "MWPC_mul == 0 && DSSDX_mul == 1 && DSSDY_mul == 1 && "
        "Veto_mul == 0 && SSD_mul == 0"
    );

    // [修改] 输出文件
    cout << "--> Writing implant-only file: " << IMPLANT_ONLY_FILE << endl;
    df_implant.Snapshot(TREE_NAME, IMPLANT_ONLY_FILE, columns);

    cout << "--> Writing decay-only file:   " << DECAY_ONLY_FILE  << endl;
    df_decay.Snapshot(TREE_NAME, DECAY_ONLY_FILE,  columns);

    cout << "=== Preselect done ===" << endl;
    return 0;
}