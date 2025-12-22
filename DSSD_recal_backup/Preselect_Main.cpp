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
    // 如果想要并行，可以打开，但如果文件有问题容易报奇怪的错
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    cout << "=== Step 0: Preselect decay/inject events ===" << endl;

    // 1. 把所有 run 串成一个 TChain
    TChain chain(TREE_NAME);
    for (int run = RUN_START; run <= RUN_END; ++run) {
        TString fname = TString::Format(INPUT_DIR_PATTERN, run);
        // AccessPathName 返回 false 表示文件存在
        if (!gSystem->AccessPathName(fname)) {
            chain.Add(fname);
        } else {
            // 找不到文件属于正常情况（有些run可能没有），打个警告即可
            // cerr << "[WARN] File not found: " << fname << endl; 
        }
    }

    if (chain.GetListOfFiles()->GetEntries() == 0) {
        cerr << "[ERROR] No input files found with pattern: "
             << INPUT_DIR_PATTERN << endl;
        return 1;
    }

    // 2. 用 RDataFrame 包起来
    ROOT::RDataFrame df(chain);
    auto n_total = df.Count();
    cout << "--> Total entries in raw data: " << n_total.GetValue() << endl;

    // 3. 定义你后面可能用到的列
    // 【修改点】：删除了 "run" 和 "evt"，防止报错
    vector<string> columns = {
        "DSSDX_Ch", "DSSDX_E", "DSSDX_mul",
        "DSSDY_Ch", "DSSDY_E", "DSSDY_mul",
        "DSSDYH_Ch","DSSDYH_E","DSSDYH_mul",
        "MWPC_mul", "Veto_mul", "SSD_mul"
    };

    // 4. 按你给的条件筛注入 / 衰变事件
    // 注入事件：MWPC > 0 (且 DSSD X/Y 均为 1)
    auto df_inject = df.Filter(
        "MWPC_mul > 0 && DSSDX_mul == 1 && DSSDY_mul == 1 && "
        "Veto_mul == 0 && SSD_mul == 0"
    );

    // 衰变事件：MWPC == 0 (且 DSSD X/Y 均为 1)
    auto df_decay = df.Filter(
        "MWPC_mul == 0 && DSSDX_mul == 1 && DSSDY_mul == 1 && "
        "Veto_mul == 0 && SSD_mul == 0"
    );

    // 这一步 Count 会触发读取，可能会有点慢
    // 如果想快一点，可以把这两行注释掉，直接运行下面的 Snapshot
    // auto n_inject = df_inject.Count();
    // auto n_decay  = df_decay.Count();
    // cout << "--> Inject entries: "  << n_inject.GetValue() << endl;
    // cout << "--> Decay entries : "  << n_decay.GetValue()  << endl;

    // 5. Snapshot 写两个小 ROOT 文件
    // RDataFrame 会自动根据 columns 里的列名去原始 Tree 里找数据
    cout << "--> Writing inject-only file: " << INJECT_ONLY_FILE << endl;
    df_inject.Snapshot(TREE_NAME, INJECT_ONLY_FILE, columns);

    cout << "--> Writing decay-only file:  " << DECAY_ONLY_FILE  << endl;
    df_decay.Snapshot(TREE_NAME, DECAY_ONLY_FILE,  columns);

    cout << "=== Preselect done ===" << endl;
    return 0;
}