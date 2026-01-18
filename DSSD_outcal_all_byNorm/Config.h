#ifndef CONFIG_H
#define CONFIG_H

#include <vector>

// =========================
// Runtime
// =========================
inline constexpr int NUM_THREADS = 8;

// =========================
// Input tree
// =========================
inline constexpr char TREE_NAME[] = "tr_map";

// =========================
// Default outputs
// =========================
inline constexpr char PRESELECT_OUT_ROOT[] = "preselected.root";
inline constexpr char NORM_PARAM_FILE[]    = "./SS032_Normalize_Params.txt";
inline constexpr char NORM_DIAG_ROOT[]     = "Diagnose_Normalize.root";
inline constexpr char DIAG_MAIN_ROOT[]     = "diagMain.root";
inline constexpr char DIAG_ED_ROOT[]       = "diaged.root";
inline constexpr char OUTPUT_DAT_FILE[]    = "ener_cal.dat";

// =========================
// Channel mapping
// =========================
#define TOTAL_CH 224
inline constexpr int NUM_DSSDX_POS  = 128; // 0-127
inline constexpr int NUM_DSSDY_POS  = 48;  // 128-175
inline constexpr int NUM_DSSDYH_POS = 48;  // 176-223

// =========================
// Alpha source energies (keV)
// =========================
inline const std::vector<double> SOURCE_ENERGIES = {5156.6, 5485.6, 5804.8};

// =========================
// Step1 Normalize regression window (raw ADC)
// X/Y raw ~2000, YH raw ~500
// =========================
inline constexpr double NORM_RAW_MIN = 50.0;
inline constexpr double NORM_RAW_MAX = 3500.0;
inline constexpr int    NORM_MIN_ENTRIES = 200; // per-strip

// =========================
// Step2 Plane-summed spectra (normalized ADC)
// =========================
inline constexpr int    HIST_BINS = 2600;
inline constexpr double HIST_MIN  = 200.0;
inline constexpr double HIST_MAX  = 3200.0;

// =========================
// TSpectrum search (plane sum)
// =========================
inline constexpr int    TS_MAX_PEAKS = 80;
inline constexpr double TS_SIGMA_BINS = 3.0;
inline constexpr double TS_THRESHOLD  = 0.05; // higher to avoid noise peaks
inline constexpr double PEAK_MIN_SEP  = 30.0; // min separation in ADC_norm units

// =========================
// Gaussian refinement in ADC_norm
// =========================
inline constexpr double GAUS_COARSE_WIN = 60.0;
inline constexpr double GAUS_SIGMA_MIN  = 1.0;
inline constexpr double GAUS_SIGMA_MAX  = 120.0;
inline constexpr double GAUS_FINE_MIN_WIN = 15.0;

// =========================
// Per-channel calibrated spectra (Energy axis)
// =========================
inline constexpr int    E_HIST_BINS = 700;
inline constexpr double E_HIST_MIN  = 4800.0;
inline constexpr double E_HIST_MAX  = 6300.0;

// FWHM reference peak and fit window
inline constexpr double FWHM_REF_E     = 5804.8;
inline constexpr double PEAK_FIT_WIN_E = 120.0;
inline constexpr double MIN_PEAK_COUNTS_IN_WIN = 30.0;

// =========================
// Failure defaults
// =========================
inline constexpr double FAIL_K    = 1.0;
inline constexpr double FAIL_B    = 0.0;
inline constexpr double FAIL_FWHM = 0.0;
inline constexpr double FAIL_CHI2 = -1.0;

#endif
