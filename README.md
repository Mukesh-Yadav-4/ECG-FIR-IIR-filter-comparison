# ECG FIR vs IIR Filter Comparison

A MATLAB-based comparative study of FIR and IIR bandpass filtering for ECG signal denoising, evaluated on real clinical data from the MIT-BIH Arrhythmia Database.

## Overview

ECG recordings are commonly corrupted by powerline interference (50 Hz), baseline wander (breathing-induced drift), and muscle artifact noise. This project simulates realistic noise on a real ECG recording, then compares two classic filtering approaches — FIR and IIR bandpass filters — on their ability to recover the clean signal, using both signal-quality metrics (SNR) and clinically meaningful metrics (R-peak detection accuracy).

## Method

1. Raw ECG signal loaded from MIT-BIH Arrhythmia Database, Record 100 (format 212, manually decoded — no external toolbox dependency)
2. Synthetic noise added: 50 Hz powerline interference, baseline wander, Gaussian muscle artifact noise
3. Two filters designed and applied via `filtfilt` (zero-phase):
   - 4th-order Butterworth IIR bandpass (0.5–40 Hz)
   - 100-tap equiripple FIR bandpass (0.5–40 Hz)
4. Evaluation: SNR (mean-centered to remove baseline offset bias) and R-peak detection accuracy (±50ms tolerance against reference peaks)

## Results

| Metric | Noisy | IIR Filtered | FIR Filtered |
|---|---|---|---|
| SNR (dB) | 3.95 | **10.12** | 7.00 |
| R-peak detection accuracy | — | 99.96% | **100%** |

**Finding**: IIR achieves superior raw noise suppression (higher SNR) for a much lower filter order, while FIR achieves marginally better R-peak preservation due to its linear-phase property — a classic and practically relevant trade-off for embedded/wearable ECG applications where compute budget and diagnostic fidelity must be balanced.

## Files

- `ECG_Filter_Analysis.m` — main reproducible script (data loading → noise simulation → filtering → evaluation)
- `fig1_raw_ecg.png` — raw ECG signal
- `fig3_clean_noisy_iir.png` — clean vs. noisy vs. IIR-filtered comparison
- `fig5_freq_response.png` — IIR filter frequency response
- `fig6_all_four_centered.png` — final 4-panel comparison (clean, noisy, IIR, FIR)

## Data

This project uses Record 100 from the MIT-BIH Arrhythmia Database (PhysioNet). Data files are not redistributed here — download `100.dat`, `100.hea`, `100.atr` directly from [physionet.org/content/mitdb](https://physionet.org/content/mitdb/1.0.0/) and place them in the same folder as the script to reproduce results.

**Citation**: Moody GB, Mark RG. The impact of the MIT-BIH Arrhythmia Database. IEEE Eng in Med and Biol 20(3):45-50 (May-June 2001).

## Author

Mukesh Yadav — B.Tech ECE, JSS Academy of Technical Education, Noida
