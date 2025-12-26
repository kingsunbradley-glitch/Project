#pragma once
#include <cstdint>
#include <cstring>

// AeqResult: "equivalent pixels" for target coverage (e.g. 80%)
// n_pix: integer number of bins needed to reach target
// n_eq : fractional bins via linear interpolation in the last bin
// cov_at_npix: achieved coverage using n_pix bins
struct AeqResult {
  long long n_pix = 0;
  double    n_eq  = 0.0;
  double    cov_at_npix = 0.0;
};

// Peak window + fit result + ROI + per-peak metrics
struct PeakWin {
  double win_lo = 0.0;
  double win_hi = 0.0;

  bool   ok  = false;
  double mu  = 0.0;
  double sig = 0.0;

  double roi_lo = 0.0;
  double roi_hi = 0.0;

  long long N   = 0;
  double mux    = 0.0;
  double muy    = 0.0;

  AeqResult A80_2D;
  AeqResult X80_1D;
  AeqResult Y80_1D;
};

// One row in summary tree
struct SummaryRow {
  int run = 0;
  int nPeaks = 0;

  long long N_all = 0;

  double A80_eq_all = 0.0, X80_eq_all = 0.0, Y80_eq_all = 0.0;
  long long A80_pix_all = 0, X80_pix_all = 0, Y80_pix_all = 0;
  double A80_cov_all = 0.0, X80_cov_all = 0.0, Y80_cov_all = 0.0;

  double muE[3];
  double sigE[3];

  long long Npk[3];
  double mux[3];
  double muy[3];

  double A80eq[3];
  long long A80pix[3];
  double A80cov[3];

  double X80eq[3];
  long long X80pix[3];
  double X80cov[3];

  double Y80eq[3];
  long long Y80pix[3];
  double Y80cov[3];

  void reset(int run_) {
    run = run_;
    nPeaks = 0;
    N_all = 0;

    A80_eq_all = X80_eq_all = Y80_eq_all = 0.0;
    A80_pix_all = X80_pix_all = Y80_pix_all = 0;
    A80_cov_all = X80_cov_all = Y80_cov_all = 0.0;

    for (int i = 0; i < 3; ++i) {
      muE[i] = -999.0;
      sigE[i] = -999.0;

      Npk[i] = 0;
      mux[i] = -999.0;
      muy[i] = -999.0;

      A80eq[i] = X80eq[i] = Y80eq[i] = -999.0;
      A80pix[i] = X80pix[i] = Y80pix[i] = 0;
      A80cov[i] = X80cov[i] = Y80cov[i] = 0.0;
    }
  }
};
