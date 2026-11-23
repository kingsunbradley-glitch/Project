#include "Config.h"
#include "ROOT/RDataFrame.hxx"
#include "TChain.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"
#include <iostream>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT;

int main() {
    Print_Channel_Map();
    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    cout << "=== Step 0: Preselect decay/implant events ===" << endl;

    TChain chain(TREE_NAME);
    vector<int> missing;
    for (int run = RUN_START; run <= RUN_END; ++run) {
        TString fname = TString::Format(INPUT_DIR_PATTERN, run);
        if (!gSystem->AccessPathName(fname)) {
            chain.Add(fname);
        } else {
            missing.push_back(run);
        }
    }

    const int n_expect = RUN_END - RUN_START + 1;
    const int n_found  = (int)chain.GetListOfFiles()->GetEntries();

    if (n_found == 0) {
        cerr << "[ERROR] No input files found." << endl;
        cerr << "        INPUT_DIR_PATTERN  = \"" << INPUT_DIR_PATTERN << "\"" << endl;
        cerr << "        RUN_START..RUN_END = " << RUN_START << " .. " << RUN_END
             << "  (" << n_expect << " 个运行号)" << endl;
        cerr << "        本程序想打开的名字（前 2 个）:" << endl;
        for (int r = RUN_START; r <= RUN_START + 1 && r <= RUN_END; ++r)
            cerr << "            " << TString::Format(INPUT_DIR_PATTERN, r) << endl;

        // 往下查一层：目录在不在、里面有什么
        TString pat(INPUT_DIR_PATTERN);
        int slash = (int)pat.Last('/');
        TString dir = (slash < 0) ? TString("./") : pat(0, slash + 1);

        cerr << "        当前工作目录 = " << gSystem->WorkingDirectory() << endl;
        cerr << "        模板所在目录 = " << dir << endl;

        // ROOT 的 TSystem 各版本 API 不一致，用 OpenDirectory 判断最稳
        auto is_dir = [](const TString& d) -> bool {
            void* h = gSystem->OpenDirectory(d);
            if (!h) return false;
            gSystem->FreeDirectory(h);
            return true;
        };

        if (!is_dir(dir)) {
            cerr << "        [诊断] 这个目录不存在或不可访问。" << endl;
            TString up = dir;
            while (up.Length() > 1 && !is_dir(up)) {
                int s2 = (int)up.Last('/');
                up = (s2 > 0) ? up(0, s2) : TString("/");
            }
            cerr << "        [诊断] 它最近一个存在的上级目录是: " << up << endl;
        } else {
            cerr << "        [诊断] 目录存在，下面是它里面真实的名字:" << endl;
            TSystemDirectory tdir(dir, dir);
            TList* lst = tdir.GetListOfFiles();
            int n_all = 0, n_map = 0;
            vector<TString> maps, others;
            if (lst) {
                for (int i = 0; i < lst->GetEntries(); ++i) {
                    TSystemFile* f = dynamic_cast<TSystemFile*>(lst->At(i));
                    if (!f || f->IsDirectory()) continue;
                    ++n_all;
                    TString n = f->GetName();
                    if (n.EndsWith("_map.root")) {
                        ++n_map;
                        if (maps.size() < 4) maps.push_back(n);
                    } else if (n.EndsWith(".root") && others.size() < 4) {
                        others.push_back(n);
                    }
                }
                delete lst;
            }
            cerr << "            共 " << n_all << " 个文件，其中 *_map.root 有 " << n_map << " 个" << endl;
            for (auto& n : maps)   cerr << "            [map]    " << n << endl;
            for (auto& n : others) cerr << "            [其它]   " << n << endl;
            if (n_map > 0)
                cerr << "        [诊断] 有 map 文件却匹配不上，请逐字对比上面"
                        "'想打开的名字'与'真实名字'的数字位数、前缀。" << endl;
            else
                cerr << "        [诊断] 这个目录里没有 *_map.root —— "
                        "注意 _tmp/_chain/_Rechain 都不是本程序的输入。" << endl;
        }

        cerr << "        若上面打印的模板不是你现在 Config.h 里的路径，" << endl;
        cerr << "        说明你跑的是旧二进制（本二进制的编译时间: "
             << __DATE__ << " " << __TIME__ << "），请重新 make。" << endl;
        return 1;
    }

    cout << "--> Input files: " << n_found << " / " << n_expect
         << " (run " << RUN_START << ".." << RUN_END << ")" << endl;
    if (!missing.empty()) {
        cout << "    [WARN] 缺 " << missing.size() << " 个运行号";
        cout << "（前 8 个:";
        for (size_t i = 0; i < missing.size() && i < 8; ++i) cout << " " << missing[i];
        cout << (missing.size() > 8 ? " ...）" : "）") << endl;
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