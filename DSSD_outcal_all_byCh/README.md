# Per-channel 3-peak calibration (Ch=0..223)

This package calibrates each strip **independently** using 3 alpha peaks:
- 5156.6 keV
- 5485.6 keV
- 5804.8 keV

For each channel:
1. Build raw ADC histogram
2. TSpectrum peak search
3. Take the **top-3 highest peaks**
4. Gaussian fit (coarse) then refine within **±1σ** (fine)
5. Linear fit: `E = k*ADC + b`
6. FWHM at 5804.8 peak: `FWHM_keV = 2.355 * (k * sigma_ADC)`

Failure policy: `k=1, b=0, FWHM=0, chi2=-1`

FWHM rule (as requested):
- fit/statistics failure -> `FWHM = 0`
- computed `FWHM > 300 keV` -> `FWHM = -1`

## Build
```bash
make
```

## Run
```bash
make run IN=/path/to/preselected.root TREE=tr_map
```

Outputs:
- `originalSpectrum.root` : `h_raw_ch000` ... `h_raw_ch223`
- `calibratedSpectrum.root` : `h_cal_ch000` ... `h_cal_ch223`, plus `hFWHM_vsCh` and `tCalib`
- `ener_cal.dat` : text table per channel

## Notes on branch types
This version reads:
- `*_Ch` as `double` (and rounds to int)
- `*_mul` as `unsigned short`
so it matches your tree (the errors you saw were due to type mismatch).
