#include "Config.h"

#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TH1D.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

struct Accum {
    long n = 0;
    double sx = 0.0;
    double sy = 0.0;
    double sxx = 0.0;
    double sxy = 0.0;
};

static void AddPoint(Accum& a, double x, double y) {
    a.n++;
    a.sx += x;
    a.sy += y;
    a.sxx += x*x;
    a.sxy += x*y;
}

static bool FitLine(const Accum& a, double& k, double& b) {
    if (a.n < NORM_MIN_ENTRIES) return false;
    double denom = (a.n * a.sxx - a.sx * a.sx);
    if (std::fabs(denom) < 1e-12) return false;
    k = (a.n * a.sxy - a.sx * a.sy) / denom;
    b = (a.sy - k * a.sx) / (double)a.n;
    if (!std::isfinite(k) || !std::isfinite(b)) return false;
    return true;
}

int main(int argc, char** argv) {
    std::string inRoot   = (argc >= 2) ? argv[1] : "";
    std::string outTxt   = (argc >= 3) ? argv[2] : NORM_PARAM_FILE;
    std::string outDiag  = (argc >= 4) ? argv[3] : NORM_DIAG_ROOT;

    if (inRoot.empty()) {
        std::cerr << "Usage: " << argv[0] << " <preselected.root> [SS032_Normalize_Params.txt] [Diagnose_Normalize.root]\n";
        return 1;
    }

    if (NUM_THREADS > 0) ROOT::EnableImplicitMT(NUM_THREADS);

    std::cout << "=== Step1: Normalize ===\n";

    ROOT::RDataFrame df(TREE_NAME, inRoot);
    auto nAll = *df.Count();
    std::cout << "Input entries: " << nAll << "\n";

    const unsigned int nSlots = df.GetNSlots();
    std::vector<std::vector<Accum>> slotAcc(nSlots, std::vector<Accum>(TOTAL_CH));

    // Use XY coincidence to define ref; independently regress each hit strip raw->ref
    df.ForeachSlot(
        [&](unsigned int s,
            const ROOT::VecOps::RVec<double>& xE,
            const ROOT::VecOps::RVec<double>& xCh,
            const ROOT::VecOps::RVec<double>& yE,
            const ROOT::VecOps::RVec<double>& yCh,
            const ROOT::VecOps::RVec<double>& yhE,
            const ROOT::VecOps::RVec<double>& yhCh,
            unsigned short yhMul) {
            if (xE.empty() || yE.empty() || xCh.empty() || yCh.empty()) return;

            const double xraw = xE[0];
            const double yraw = yE[0];
            if (xraw < NORM_RAW_MIN || xraw > NORM_RAW_MAX) return;
            if (yraw < NORM_RAW_MIN || yraw > NORM_RAW_MAX) return;

            const double ref = 0.5 * (xraw + yraw);

            const int xid = (int)xCh[0];
            const int yid = NUM_DSSDX_POS + (int)yCh[0];
            if (xid >= 0 && xid < NUM_DSSDX_POS) AddPoint(slotAcc[s][xid], xraw, ref);
            if (yid >= NUM_DSSDX_POS && yid < NUM_DSSDX_POS + NUM_DSSDY_POS) AddPoint(slotAcc[s][yid], yraw, ref);

            // YH: if present, also tie to the same ref
            if (yhMul == 1 && !yhE.empty() && !yhCh.empty()) {
                const double yhraw = yhE[0];
                if (yhraw >= NORM_RAW_MIN && yhraw <= NORM_RAW_MAX) {
                    const int yhid = NUM_DSSDX_POS + NUM_DSSDY_POS + (int)yhCh[0];
                    if (yhid >= 0 && yhid < TOTAL_CH) AddPoint(slotAcc[s][yhid], yhraw, ref);
                }
            }
        },
        {"DSSDX_E","DSSDX_Ch","DSSDY_E","DSSDY_Ch","DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

    // Merge slots
    std::vector<Accum> acc(TOTAL_CH);
    for (unsigned int s = 0; s < nSlots; ++s) {
        for (int id = 0; id < TOTAL_CH; ++id) {
            acc[id].n   += slotAcc[s][id].n;
            acc[id].sx  += slotAcc[s][id].sx;
            acc[id].sy  += slotAcc[s][id].sy;
            acc[id].sxx += slotAcc[s][id].sxx;
            acc[id].sxy += slotAcc[s][id].sxy;
        }
    }

    // Fit and write
    std::ofstream out(outTxt);
    out << "# id  k_norm  b_norm  ok  n\n";
    out << std::fixed << std::setprecision(8);

    std::vector<double> vK, vB;
    vK.reserve(TOTAL_CH);
    vB.reserve(TOTAL_CH);

    for (int id = 0; id < TOTAL_CH; ++id) {
        double k = 1.0, b = 0.0;
        bool ok = FitLine(acc[id], k, b);
        out << id << " " << k << " " << b << " " << (ok ? 1 : 0) << " " << acc[id].n << "\n";
        if (ok) { vK.push_back(k); vB.push_back(b); }
    }
    out.close();

    // Diagnostics
    TFile f(outDiag.c_str(), "RECREATE");
    TH1D hK("h_k_norm", "k_norm distribution;k_norm;Counts", 400, 0.0, 4.0);
    TH1D hB("h_b_norm", "b_norm distribution;b_norm (ADC);Counts", 400, -2000.0, 2000.0);
    TH1D hN("h_n_points", "N points per strip;N;Counts", 400, 0.0, 20000.0);

    for (int id = 0; id < TOTAL_CH; ++id) {
        double k = 1.0, b = 0.0;
        bool ok = FitLine(acc[id], k, b);
        if (ok) {
            hK.Fill(k);
            hB.Fill(b);
        }
        hN.Fill((double)acc[id].n);
    }

    hK.Write();
    hB.Write();
    hN.Write();
    f.Close();

    std::cout << "Saved norm params: " << outTxt << "\n";
    std::cout << "Saved diag ROOT:   " << outDiag << "\n";
    return 0;
}
