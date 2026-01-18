#include "Config.h"

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    std::string in = (argc >= 2) ? argv[1] : "";
    std::string out = (argc >= 3) ? argv[2] : PRESELECT_OUT_ROOT;
    if (in.empty()) {
        std::cerr << "Usage: " << argv[0] << " <runXXXXX_map.root> [preselected.root]\n";
        return 1;
    }

    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    std::cout << "=== Step0: Preselect ===\n";

    ROOT::RDataFrame df(TREE_NAME, in);
    auto nAll = *df.Count();
    auto dff = df.Filter("DSSDX_mul==1 && DSSDY_mul==1", "Xmul=1 & Ymul=1");
    auto nSel = *dff.Count();

    std::cout << "Input entries: " << nAll << "\n";
    std::cout << "Selected entries (Xmul=1 & Ymul=1): " << nSel << "\n";
    std::cout << "Writing: " << out << "\n";

    ROOT::RDF::RSnapshotOptions opt;
    opt.fMode = "RECREATE";
    dff.Snapshot(TREE_NAME, out, "", opt);

    std::cout << "=== Preselect done ===\n";
    return 0;
}
