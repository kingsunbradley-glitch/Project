#include "Calibrator.h"
#include "Config.h"

#include <TROOT.h>

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    std::string inRoot  = (argc >= 2) ? argv[1] : "";
    std::string normTxt = (argc >= 3) ? argv[2] : NORM_PARAM_FILE;
    std::string outDat  = (argc >= 4) ? argv[3] : OUTPUT_DAT_FILE;
    std::string diagMain= (argc >= 5) ? argv[4] : DIAG_MAIN_ROOT;
    std::string diagEd  = (argc >= 6) ? argv[5] : DIAG_ED_ROOT;

    if (inRoot.empty()) {
        std::cerr << "Usage: " << argv[0] << " <preselected.root> [SS032_Normalize_Params.txt] [ener_cal.dat] [diagMain.root] [diaged.root]\n";
        return 1;
    }

    // Batch mode: do not create GUI canvas
    gROOT->SetBatch(kTRUE);

    Calibrator cal(inRoot, normTxt, outDat, diagMain, diagEd);
    return cal.Run();
}
