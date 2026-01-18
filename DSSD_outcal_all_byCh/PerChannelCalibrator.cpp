// PerChannelCalibrator.cpp (robust per-channel 3-peak calibration)
// Build:
//   g++ -O2 -std=c++17 PerChannelCalibrator.cpp $(root-config --cflags --libs) -lSpectrum -o PerChannelCalibrator
//
// Run:
//   ./PerChannelCalibrator preselected.root tr_map originalSpectrum.root calibratedSpectrum.root ener_cal.dat

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>

#include <TFile.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TF1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using ROOT::VecOps::RVec;

// ------------------------
// User config
// ------------------------
static constexpr int TOTAL_CH = 224;

// raw spectrum binning (ADC)
static constexpr int    RAW_BINS = 4096;
static constexpr double RAW_MIN  = 0.0;
static constexpr double RAW_MAX  = 4096.0;

// energy spectrum binning (keV)
static constexpr int    E_BINS = 200;
static constexpr double E_MIN  = 5000.0;
static constexpr double E_MAX  = 6000.0;

// alpha energies (keV)
static const std::array<double,3> EREF = {5156.6, 5485.6, 5804.8};

// TSpectrum
static constexpr int    TS_MAX_PEAKS  = 200;
static constexpr double TS_SIGMA_BINS = 3.0;
static constexpr double TS_THRESHOLD  = 0.3; // fraction of max

// Gaussian fits (ADC axis)
static constexpr double COARSE_WIN = 50.0;     // +/- around seed for coarse fit
static constexpr double SIGMA_MIN  = 0.5;
static constexpr double SIGMA_MAX  = 200.0;

// peak selection constraints (ADC)
static constexpr double MIN_PEAK_SEP = 10.0;   // min separation between picked peaks
static constexpr double MAX_PEAK_SEP = 250.0;  // max separation (within dynamic window)
static constexpr double RATIO_TOL    = 0.6;    // require d2/d1 in [1-0.6, 1+0.6]

// min counts in fit window (integral)
static constexpr double MIN_WIN_COUNTS = 30.0;

// calibration sanity
// X/Y 典型 ~2-4 keV/ADC, YH 典型 ~8-15 keV/ADC；这里给宽一点
static constexpr double K_MIN = 0.5;
static constexpr double K_MAX = 30.0;

// 3-point line residual sanity (keV)
static constexpr double MAX_RESID_KEV = 150.0;

// FWHM rule
static constexpr double FWHM_MAX_OK = 300.0; // >300 -> -1

// failure defaults
static constexpr double FAIL_K    = 1.0;
static constexpr double FAIL_B    = 0.0;
static constexpr double FAIL_FWHM = 0.0;
static constexpr double FAIL_CHI2 = -1.0;

// channel mapping
static int chX(int xch)   { return xch; }
static int chY(int ych)   { return 128 + ych; }
static int chYH(int yhch) { return 176 + yhch; }

struct PeakFitResult {
  bool   ok = false;
  double mu = 0.0;
  double sigma = 0.0;
  double chi2ndf = 0.0;
};

static PeakFitResult FitPeakTwoStage(TH1D* h, double seed, double coarseWin) {
  PeakFitResult r;
  if (!h) return r;

  // 1) coarse fit in [seed-coarseWin, seed+coarseWin]
  double x1 = seed - coarseWin;
  double x2 = seed + coarseWin;
  if (x2 <= x1) return r;

  int b1 = h->FindBin(x1);
  int b2 = h->FindBin(x2);
  double sum = h->Integral(b1, b2);
  if (sum < MIN_WIN_COUNTS) return r;

  TF1 f1("f1", "gaus(0)", x1, x2);
  f1.SetParameters(h->GetMaximum(), seed, 8.0);
  f1.SetParLimits(1, seed - coarseWin, seed + coarseWin);
  f1.SetParLimits(2, SIGMA_MIN, SIGMA_MAX);

  if (h->Fit(&f1, "QNR") != 0) return r;

  double mu1 = f1.GetParameter(1);
  double s1  = std::fabs(f1.GetParameter(2));
  if (!std::isfinite(mu1) || !std::isfinite(s1) || s1 <= 0) return r;

  // 2) refine fit in [mu1-1sigma, mu1+1sigma] (minimum 8 ADC half-window)
  double w = std::max(8.0, 1.0 * s1);
  double y1 = mu1 - w;
  double y2 = mu1 + w;

  TF1 f2("f2", "gaus(0)", y1, y2);
  f2.SetParameters(f1.GetParameter(0), mu1, s1);
  f2.SetParLimits(1, mu1 - w, mu1 + w);
  f2.SetParLimits(2, SIGMA_MIN, SIGMA_MAX);

  if (h->Fit(&f2, "QNR") != 0) return r;

  double mu2 = f2.GetParameter(1);
  double s2  = std::fabs(f2.GetParameter(2));
  double chi2ndf = (f2.GetNDF() > 0) ? (f2.GetChisquare() / f2.GetNDF()) : 0.0;

  if (!std::isfinite(mu2) || !std::isfinite(s2) || s2 <= 0) return r;

  r.ok = true;
  r.mu = mu2;
  r.sigma = s2;
  r.chi2ndf = chi2ndf;
  return r;
}

struct CandPeak {
  double x = 0.0;
  double h = 0.0;
};

static bool ChooseBestTriplet(TH1D* hist,
                              const std::vector<CandPeak>& cands,
                              std::array<double,3>& seeds_out) {
  if (!hist) return false;
  if (cands.size() < 3) return false;

  // Use top M candidates by height
  std::vector<CandPeak> top = cands;
  std::sort(top.begin(), top.end(), [](const CandPeak& a, const CandPeak& b){
    return a.h > b.h;
  });
  const int M = std::min<int>((int)top.size(), 12);
  top.resize(M);

  double bestScore = -1e99;
  std::array<double,3> best = {0,0,0};

  for (int i=0;i<M;i++){
    for (int j=i+1;j<M;j++){
      for (int k=j+1;k<M;k++){
        // sort by x
        std::array<CandPeak,3> t = {top[i], top[j], top[k]};
        std::sort(t.begin(), t.end(), [](const CandPeak& a, const CandPeak& b){ return a.x < b.x; });

        double x1=t[0].x, x2=t[1].x, x3=t[2].x;
        double d1=x2-x1, d2=x3-x2;
        if (d1 < MIN_PEAK_SEP || d2 < MIN_PEAK_SEP) continue;
        if (d1 > MAX_PEAK_SEP || d2 > MAX_PEAK_SEP) continue;

        double ratio = d2/d1;
        if (ratio < (1.0 - RATIO_TOL) || ratio > (1.0 + RATIO_TOL)) continue;

        double sumH = t[0].h + t[1].h + t[2].h;
        // prefer ratio close to 1
        double penalty = std::fabs(ratio - 1.0);
        double score = sumH / (1.0 + 5.0*penalty);

        if (score > bestScore){
          bestScore = score;
          best = {x1,x2,x3};
        }
      }
    }
  }

  if (bestScore < -1e50) return false;
  seeds_out = best;
  return true;
}

static bool FindTop3PeaksRobust(TH1D* h, std::array<double,3>& seeds_out) {
  seeds_out = {0,0,0};
  if (!h) return false;
  if (h->GetEntries() < 200) return false;

  // Dynamic search window centered at maximum bin
  double xMax = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  double winL = std::max(RAW_MIN, xMax - 250.0);
  double winR = std::min(RAW_MAX, xMax + 250.0);
  if (winR <= winL) return false;

  TSpectrum sp(TS_MAX_PEAKS);
  int nfound = sp.Search(h, TS_SIGMA_BINS, "nobackground", TS_THRESHOLD);
  if (nfound <= 0) return false;

  std::vector<CandPeak> cands;
  cands.reserve(nfound);

  double* xpos = sp.GetPositionX();
  for (int i=0;i<nfound;i++){
    double x = xpos[i];
    if (x < winL || x > winR) continue;
    double hh = h->GetBinContent(h->FindBin(x));
    if (hh <= 0) continue;
    cands.push_back({x, hh});
  }

  if (cands.size() < 3) return false;

  // choose best triplet
  if (!ChooseBestTriplet(h, cands, seeds_out)) return false;

  // ensure ascending
  std::sort(seeds_out.begin(), seeds_out.end());
  return true;
}

static bool FitLine3pt(const std::array<double,3>& adc, double& k, double& b, double& maxAbsResid) {
  // analytic least squares for 3 points
  const double x1=adc[0], x2=adc[1], x3=adc[2];
  const double y1=EREF[0], y2=EREF[1], y3=EREF[2];

  const double Sx = x1+x2+x3;
  const double Sy = y1+y2+y3;
  const double Sxx = x1*x1 + x2*x2 + x3*x3;
  const double Sxy = x1*y1 + x2*y2 + x3*y3;

  const double n = 3.0;
  const double denom = (n*Sxx - Sx*Sx);
  if (std::fabs(denom) < 1e-12) return false;

  k = (n*Sxy - Sx*Sy) / denom;
  b = (Sy - k*Sx) / n;

  if (!std::isfinite(k) || !std::isfinite(b)) return false;

  // residual sanity
  double r1 = y1 - (k*x1 + b);
  double r2 = y2 - (k*x2 + b);
  double r3 = y3 - (k*x3 + b);
  maxAbsResid = std::max({std::fabs(r1), std::fabs(r2), std::fabs(r3)});
  return true;
}

struct CalibRow {
  int ch = -1;
  double k = FAIL_K;
  double b = FAIL_B;
  double fwhm = FAIL_FWHM;
  double chi2ndf = FAIL_CHI2; // use chi2/ndf of refined 5804.8 peak fit
  int ok = 0;
};

static void WriteFail(std::ofstream& out, CalibRow& row){
  out << row.ch << " " << FAIL_B << " " << FAIL_K << " " << FAIL_FWHM << " " << FAIL_CHI2 << "\n";
}

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cerr << "Usage: " << argv[0]
              << " <input.root> <tree> <originalSpectrum.root> <calibratedSpectrum.root> <ener_cal.dat>\n";
    return 1;
  }
  const std::string inRoot = argv[1];
  const std::string tree   = argv[2];
  const std::string outOrig= argv[3];
  const std::string outCal = argv[4];
  const std::string outDat = argv[5];

  ROOT::EnableImplicitMT();

  ROOT::RDataFrame df(tree, inRoot);

  // keep your clean selection (can remove if you want)
  auto dff = df.Filter("DSSDX_mul==1 && DSSDY_mul==1");

  const unsigned int nSlots = ROOT::GetThreadPoolSize() > 0 ? ROOT::GetThreadPoolSize() : 1;

  // slot-local hists
  std::vector<std::vector<std::unique_ptr<TH1D>>> hSlot(nSlots);
  for (unsigned int s=0;s<nSlots;s++){
    hSlot[s].resize(TOTAL_CH);
    for (int ch=0; ch<TOTAL_CH; ch++){
      hSlot[s][ch] = std::make_unique<TH1D>(
        Form("h_raw_s%u_ch%03d", s, ch),
        Form("raw ch %03d;ADC;Counts", ch),
        RAW_BINS, RAW_MIN, RAW_MAX
      );
      hSlot[s][ch]->SetDirectory(nullptr);
    }
  }

  // NOTE: your branches are: mul is unsigned short, Ch is double array
  auto fillRaw = [&](unsigned int slot,
                     const RVec<double>& xE, const RVec<double>& xCh, unsigned short xMul,
                     const RVec<double>& yE, const RVec<double>& yCh, unsigned short yMul,
                     const RVec<double>& yhE,const RVec<double>& yhCh,unsigned short yhMul){
    if (xMul==1 && xE.size()>0 && xCh.size()>0){
      int ch = chX((int)llround(xCh[0]));
      if (0<=ch && ch<TOTAL_CH) hSlot[slot][ch]->Fill(xE[0]);
    }
    if (yMul==1 && yE.size()>0 && yCh.size()>0){
      int ch = chY((int)llround(yCh[0]));
      if (0<=ch && ch<TOTAL_CH) hSlot[slot][ch]->Fill(yE[0]);
    }
    if (yhMul==1 && yhE.size()>0 && yhCh.size()>0){
      int ch = chYH((int)llround(yhCh[0]));
      if (0<=ch && ch<TOTAL_CH) hSlot[slot][ch]->Fill(yhE[0]);
    }
  };

  dff.ForeachSlot(fillRaw,
    {"DSSDX_E","DSSDX_Ch","DSSDX_mul",
     "DSSDY_E","DSSDY_Ch","DSSDY_mul",
     "DSSDYH_E","DSSDYH_Ch","DSSDYH_mul"});

  // merge slots
  std::vector<std::unique_ptr<TH1D>> hRaw(TOTAL_CH);
  for (int ch=0; ch<TOTAL_CH; ch++){
    hRaw[ch] = std::make_unique<TH1D>(
      Form("h_raw_ch%03d", ch),
      Form("raw ch %03d;ADC;Counts", ch),
      RAW_BINS, RAW_MIN, RAW_MAX
    );
    hRaw[ch]->SetDirectory(nullptr);
    for (unsigned int s=0;s<nSlots;s++){
      hRaw[ch]->Add(hSlot[s][ch].get());
    }
  }

  // write originalSpectrum.root
  {
    TFile f(outOrig.c_str(), "RECREATE");
    for (int ch=0; ch<TOTAL_CH; ch++) hRaw[ch]->Write();
    f.Close();
  }

  // output calibrated spectra + FWHM graph + tCalib
  TFile fCal(outCal.c_str(), "RECREATE");
  TH1D hFWHM("hFWHM_vsCh", "FWHM vs Channel;Ch;FWHM (keV)", TOTAL_CH, -0.5, TOTAL_CH-0.5);

  TTree t("tCalib","per-channel calibration");
  CalibRow row;
  t.Branch("ch",&row.ch,"ch/I");
  t.Branch("k",&row.k,"k/D");
  t.Branch("b",&row.b,"b/D");
  t.Branch("fwhm",&row.fwhm,"fwhm/D");
  t.Branch("chi2ndf",&row.chi2ndf,"chi2ndf/D");
  t.Branch("ok",&row.ok,"ok/I");

  std::ofstream out(outDat);
  out << "# ch  b  k  FWHM_keV  chi2ndf\n";

  for (int ch=0; ch<TOTAL_CH; ch++){
    row = CalibRow{};
    row.ch = ch;

    TH1D* h = hRaw[ch].get();
    if (!h || h->GetEntries() < 300) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // 1) robust peak seeds (windowed + best-triplet)
    std::array<double,3> seed{};
    if (!FindTop3PeaksRobust(h, seed)) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // 2) two-stage gaussian refine
    std::array<PeakFitResult,3> pf;
    bool okPeaks = true;
    for (int i=0;i<3;i++){
      pf[i] = FitPeakTwoStage(h, seed[i], COARSE_WIN);
      if (!pf[i].ok) { okPeaks = false; break; }
    }
    if (!okPeaks) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // 3) sort refined peaks by mu asc
    std::array<double,3> adc = {pf[0].mu, pf[1].mu, pf[2].mu};
    std::sort(adc.begin(), adc.end());

    // spacing sanity (avoid duplicate/shoulder)
    if ((adc[1]-adc[0] < MIN_PEAK_SEP) || (adc[2]-adc[1] < MIN_PEAK_SEP)) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // 4) line fit E = k*ADC + b + residual sanity
    double k=FAIL_K, b=FAIL_B, maxResid=1e9;
    if (!FitLine3pt(adc, k, b, maxResid)) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // k sanity
    if (!std::isfinite(k) || k < K_MIN || k > K_MAX) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }
    // residual sanity
    if (maxResid > MAX_RESID_KEV) {
      WriteFail(out, row);
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    row.k = k;
    row.b = b;

    // 5) choose the highest-ADC refined peak as 5804.8 reference
    int i_hi = 0;
    for (int i=1;i<3;i++){
      if (pf[i].mu > pf[i_hi].mu) i_hi = i;
    }
    double sigma_adc = pf[i_hi].sigma;
    row.chi2ndf = pf[i_hi].chi2ndf;

    // sigma sanity
    if (!std::isfinite(sigma_adc) || sigma_adc <= 0 || sigma_adc > SIGMA_MAX) {
      // treat as bad channel
      row.k = FAIL_K; row.b = FAIL_B; row.fwhm = 0.0; row.chi2ndf = FAIL_CHI2; row.ok = 0;
      out << ch << " " << row.b << " " << row.k << " " << row.fwhm << " " << row.chi2ndf << "\n";
      hFWHM.SetBinContent(ch+1, 0.0);
      t.Fill();
      continue;
    }

    // FWHM in keV
    double fwhm = 2.355 * (k * sigma_adc);

    if (!std::isfinite(fwhm) || fwhm <= 0) {
      fwhm = 0.0; // problematic -> 0
    } else if (fwhm > FWHM_MAX_OK) {
      fwhm = -1.0; // >300 -> -1
    }
    row.fwhm = fwhm;

    // 6) build calibrated histogram by bin-mapping
    auto hcal = std::make_unique<TH1D>(
      Form("h_cal_ch%03d", ch),
      Form("calibrated ch %03d;Energy (keV);Counts", ch),
      E_BINS, E_MIN, E_MAX
    );
    hcal->SetDirectory(&fCal);

    for (int ib=1; ib<=h->GetNbinsX(); ib++){
      double c = h->GetBinContent(ib);
      if (c <= 0) continue;
      double x = h->GetBinCenter(ib);
      double e = k*x + b;
      int be = hcal->FindBin(e);
      if (be>=1 && be<=hcal->GetNbinsX()) hcal->AddBinContent(be, c);
    }

    row.ok = 1;
    out << ch << " " << row.b << " " << row.k << " " << row.fwhm << " " << row.chi2ndf << "\n";

    // store FWHM vs ch (store 0 / -1 as you requested)
    hFWHM.SetBinContent(ch+1, row.fwhm);

    t.Fill();
    hcal->Write();
  }

  hFWHM.Write();
  t.Write();
  fCal.Close();
  out.close();

  std::cout << "Wrote: " << outOrig << "\n"
            << "Wrote: " << outCal << "\n"
            << "Wrote: " << outDat << "\n";
  return 0;
}
