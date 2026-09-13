%% ECG Signal Denoising: FIR vs IIR Filter Comparison
% Author: Mukesh Yadav
% Dataset: MIT-BIH Arrhythmia Database (Record 100)
%
% Compares 4th-order Butterworth IIR and 100-tap Equiripple FIR bandpass
% filters for ECG denoising. All plots are shown in a single tabbed window
% and all 5 tabs are automatically exported as high-resolution PNGs.

clear;
clc;
close all;

%% 1. Path Setup & Data Loading
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

data_file = fullfile(script_dir, '100.dat');
if ~isfile(data_file)
    error(['"100.dat" was not found in: %s\n' ...
           'Please download record 100 from PhysioNet MIT-BIH database.'], script_dir);
end

% Read format 212 binary data
fid = fopen(data_file, 'r');
raw = fread(fid, Inf, 'ubit12');
fclose(fid);

% Record 100 has 2 interleaved channels (MLII and V5)
raw = reshape(raw, 2, [])';
raw(raw > 2047) = raw(raw > 2047) - 4096; % 12-bit two's complement

gain     = 200;      % 200 adu/mV
baseline = 1024;     % baseline offset
fs       = 360;      % sampling frequency (Hz)

% Convert to physical units (mV) and select Lead MLII
ecg_raw = (raw(:, 1) - baseline) / gain;
N       = length(ecg_raw);
t       = (0:N-1)' / fs;

fprintf('Loaded MIT-BIH Record 100: %d samples (%.1f s at %d Hz)\n', N, N/fs, fs);

%% 2. Synthetic Noise Injection
rng(42); % Fixed seed for exact reproducibility

noise_powerline = 0.05 * sin(2 * pi * 50 * t);
noise_wander    = 0.15 * sin(2 * pi * 0.3 * t);
noise_emg       = 0.05 * randn(size(t));

noisy_signal = ecg_raw + noise_powerline + noise_wander + noise_emg;

%% 3. Filter Design (0.5 Hz - 40 Hz Passband)
nyquist   = fs / 2;
low_cut   = 0.5 / nyquist;
high_cut  = 40 / nyquist;

% --- 3.1 IIR: 4th-Order Butterworth Bandpass ---
[b_iir, a_iir] = butter(4, [low_cut, high_cut], 'bandpass');
iir_filtered   = filtfilt(b_iir, a_iir, noisy_signal);

% --- 3.2 FIR: 100-Tap Equiripple Bandpass ---
fir_filter = designfilt('bandpassfir', ...
    'FilterOrder', 100, ...
    'CutoffFrequency1', 0.5, ...
    'CutoffFrequency2', 40, ...
    'SampleRate', fs);
fir_filtered = filtfilt(fir_filter, noisy_signal);

%% 4. Mean-Centering & Quantitative SNR / MSE
clean_centered = ecg_raw - mean(ecg_raw);
noisy_centered = noisy_signal - mean(noisy_signal);
iir_centered   = iir_filtered - mean(iir_filtered);
fir_centered   = fir_filtered - mean(fir_filtered);

calc_snr = @(ref, x) 10 * log10(sum(ref.^2) / sum((ref - x).^2));
calc_mse = @(ref, x) mean((ref - x).^2);

snr_noisy = calc_snr(clean_centered, noisy_centered);
snr_iir   = calc_snr(clean_centered, iir_centered);
snr_fir   = calc_snr(clean_centered, fir_centered);

mse_noisy = calc_mse(clean_centered, noisy_centered);
mse_iir   = calc_mse(clean_centered, iir_centered);
mse_fir   = calc_mse(clean_centered, fir_centered);

%% 5. Clinical Metric: R-Peak Detection Accuracy
min_h = 0.3;
min_d = round(0.4 * fs);

[~, locs_clean] = findpeaks(clean_centered, 'MinPeakHeight', min_h, 'MinPeakDistance', min_d);
[~, locs_iir]   = findpeaks(iir_centered,   'MinPeakHeight', min_h, 'MinPeakDistance', min_d);
[~, locs_fir]   = findpeaks(fir_centered,   'MinPeakHeight', min_h, 'MinPeakDistance', min_d);

tol = round(0.05 * fs); % +/- 50 ms clinical window

matched_iir = sum(arrayfun(@(loc) any(abs(locs_iir - loc) <= tol), locs_clean));
matched_fir = sum(arrayfun(@(loc) any(abs(locs_fir - loc) <= tol), locs_clean));

acc_iir = (matched_iir / length(locs_clean)) * 100;
acc_fir = (matched_fir / length(locs_clean)) * 100;

%% 6. Print Console Summary Table
fprintf('\n=================================================================\n');
fprintf('                   FILTER COMPARISON RESULTS                     \n');
fprintf('=================================================================\n');
Method    = {'Noisy Signal'; 'IIR Filtered (Butterworth)'; 'FIR Filtered (Equiripple)'};
SNR_dB    = [snr_noisy; snr_iir; snr_fir];
MSE       = [mse_noisy; mse_iir; mse_fir];
R_Peak_Acc= [NaN; acc_iir; acc_fir];

results_table = table(Method, SNR_dB, MSE, R_Peak_Acc, ...
    'VariableNames', {'Method', 'SNR_dB', 'MSE', 'R_Peak_Accuracy_pct'});
disp(results_table);
fprintf('Total R-peaks evaluated: %d\n', length(locs_clean));
fprintf('IIR matched: %d (%.2f%%) | FIR matched: %d (%.2f%%)\n\n', ...
    matched_iir, acc_iir, matched_fir, acc_fir);

%% 7. Multi-Tab Visualization Window
% Auto-fit screen resolution so title bar and Close [X] are ALWAYS on-screen
screen = get(0, 'ScreenSize'); % [left, bottom, width, height]
fig_w  = min(1200, round(screen(3) * 0.82));
fig_h  = min(660,  round(screen(4) * 0.74)); % Safe height ensures title bar is never cut off
fig_x  = round((screen(3) - fig_w) / 2);
fig_y  = round((screen(4) - fig_h) / 2) - 15;

fig = figure('Name', 'ECG Filter Comparison (FIR vs IIR)', ...
             'Color', 'w', 'NumberTitle', 'off', ...
             'Position', [fig_x, fig_y, fig_w, fig_h]);
tabs = uitabgroup(fig);

idx_10s = 1:round(10 * fs);
t_10s   = t(idx_10s);
title_color = [0, 0.15, 0.55]; % Crisp high-contrast navy

% -------------------------------------------------------------
% Tab 1: Raw & Corrupted Signals
% -------------------------------------------------------------
tab1 = uitab(tabs, 'Title', 'Raw & Noisy Signals', 'BackgroundColor', 'w');

ax1a = subplot(3, 1, 1, 'Parent', tab1);
plot(ax1a, t_10s, ecg_raw(idx_10s), 'b-', 'LineWidth', 1.1);
title(ax1a, 'Original Clean ECG (MIT-BIH Record 100, Lead MLII)', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax1a, 'Amplitude (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax1a, [0, 10]); ylim(ax1a, [-0.6, 1.3]);
yticks(ax1a, -0.5:0.5:1.0);
xticklabels(ax1a, {});

ax1b = subplot(3, 1, 2, 'Parent', tab1);
plot(ax1b, t_10s, noise_wander(idx_10s), 'k--', 'LineWidth', 1.0, 'DisplayName', '0.3 Hz Wander'); hold(ax1b, 'on');
plot(ax1b, t_10s, noise_powerline(idx_10s), 'm-', 'LineWidth', 0.8, 'DisplayName', '50 Hz Powerline');
plot(ax1b, t_10s, noise_emg(idx_10s), 'r:', 'LineWidth', 0.8, 'DisplayName', 'EMG Noise');
title(ax1b, 'Injected Synthetic Noise Components', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax1b, 'Amplitude (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax1b, [0, 10]); ylim(ax1b, [-0.3, 0.3]);
yticks(ax1b, -0.2:0.2:0.2);
xticklabels(ax1b, {});
legend(ax1b, 'Location', 'northeast', 'FontSize', 10, 'FontWeight', 'bold');

ax1c = subplot(3, 1, 3, 'Parent', tab1);
plot(ax1c, t_10s, noisy_signal(idx_10s), 'r-', 'LineWidth', 0.9);
title(ax1c, sprintf('Composite Corrupted ECG (Input SNR: %.2f dB)', snr_noisy), ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', title_color);
xlabel(ax1c, 'Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax1c, 'Amplitude (mV)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax1c, [0, 10]); ylim(ax1c, [-0.8, 1.5]);
xticks(ax1c, 0:1:10); yticks(ax1c, -0.5:0.5:1.5);

% -------------------------------------------------------------
% Tab 2: 4-Panel Filter Comparison (Centered)
% -------------------------------------------------------------
tab2 = uitab(tabs, 'Title', 'Filter Comparison', 'BackgroundColor', 'w');

ax2a = subplot(4, 1, 1, 'Parent', tab2);
plot(ax2a, t_10s, clean_centered(idx_10s), 'b-', 'LineWidth', 1.1);
title(ax2a, '1. Clean Reference ECG (Mean-Centered)', ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax2a, 'mV', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax2a, [0, 10]); ylim(ax2a, [-0.6, 1.3]);
yticks(ax2a, -0.5:0.5:1.0);
xticklabels(ax2a, {});

ax2b = subplot(4, 1, 2, 'Parent', tab2);
plot(ax2b, t_10s, noisy_centered(idx_10s), 'r-', 'LineWidth', 0.9);
title(ax2b, sprintf('2. Corrupted Noisy ECG (SNR: %.2f dB)', snr_noisy), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax2b, 'mV', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax2b, [0, 10]); ylim(ax2b, [-0.8, 1.5]);
yticks(ax2b, -0.5:0.5:1.5);
xticklabels(ax2b, {});

ax2c = subplot(4, 1, 3, 'Parent', tab2);
plot(ax2c, t_10s, iir_centered(idx_10s), 'Color', [0, 0.55, 0.2], 'LineWidth', 1.1);
title(ax2c, sprintf('3. IIR Butterworth Filtered (SNR: %.2f dB | R-Peak Acc: %.2f%%)', snr_iir, acc_iir), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax2c, 'mV', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax2c, [0, 10]); ylim(ax2c, [-0.6, 1.3]);
yticks(ax2c, -0.5:0.5:1.0);
xticklabels(ax2c, {});

ax2d = subplot(4, 1, 4, 'Parent', tab2);
plot(ax2d, t_10s, fir_centered(idx_10s), 'Color', [0.7, 0, 0.7], 'LineWidth', 1.1);
title(ax2d, sprintf('4. FIR Equiripple Filtered (SNR: %.2f dB | R-Peak Acc: %.2f%%)', snr_fir, acc_fir), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', title_color);
xlabel(ax2d, 'Time (seconds)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax2d, 'mV', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax2d, [0, 10]); ylim(ax2d, [-0.6, 1.3]);
xticks(ax2d, 0:1:10); yticks(ax2d, -0.5:0.5:1.0);

% -------------------------------------------------------------
% Tab 3: IIR Frequency Response
% -------------------------------------------------------------
tab3 = uitab(tabs, 'Title', 'Frequency Response', 'BackgroundColor', 'w');

[h_resp, f_resp] = freqz(b_iir, a_iir, 2048, fs);

ax3a = subplot(2, 1, 1, 'Parent', tab3);
plot(ax3a, f_resp, 20 * log10(abs(h_resp)), 'b-', 'LineWidth', 1.3);
title(ax3a, '4th-Order Butterworth IIR Magnitude Response (0.5 - 40 Hz Passband)', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', title_color);
ylabel(ax3a, 'Magnitude (dB)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax3a, [0, 60]); ylim(ax3a, [-80, 5]);
xticks(ax3a, 0:10:60); yticks(ax3a, -80:20:0);

ax3b = subplot(2, 1, 2, 'Parent', tab3);
plot(ax3b, f_resp, unwrap(angle(h_resp)) * (180 / pi), 'r-', 'LineWidth', 1.3);
title(ax3b, 'IIR Phase Response (Compensated via zero-phase filtfilt)', ...
    'FontSize', 12, 'FontWeight', 'bold', 'Color', title_color);
xlabel(ax3b, 'Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax3b, 'Phase (degrees)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax3b, [0, 60]);
xticks(ax3b, 0:10:60); yticks(ax3b, -360:90:0);

% -------------------------------------------------------------
% Tab 4: R-Peak Detection & Alignment
% -------------------------------------------------------------
tab4 = uitab(tabs, 'Title', 'R-Peak Alignment', 'BackgroundColor', 'w');

idx_zoom = round(2 * fs) : round(5 * fs);
t_zoom   = t(idx_zoom);

ax4 = axes('Parent', tab4);
plot(ax4, t_zoom, clean_centered(idx_zoom), 'b-', 'LineWidth', 1.3, 'DisplayName', 'Clean Reference'); hold(ax4, 'on');
plot(ax4, t_zoom, iir_centered(idx_zoom), 'Color', [0, 0.55, 0.2], 'LineWidth', 1.1, 'DisplayName', 'IIR Filtered');
plot(ax4, t_zoom, fir_centered(idx_zoom), 'Color', [0.7, 0, 0.7], 'LineWidth', 1.1, 'DisplayName', 'FIR Filtered');

in_zoom = @(locs) locs(locs >= idx_zoom(1) & locs <= idx_zoom(end));
locs_clean_zoom = in_zoom(locs_clean);
locs_iir_zoom   = in_zoom(locs_iir);
locs_fir_zoom   = in_zoom(locs_fir);

plot(ax4, t(locs_clean_zoom), clean_centered(locs_clean_zoom), 'bo', 'MarkerSize', 9, 'LineWidth', 1.8, 'DisplayName', 'Clean R-Peaks');
plot(ax4, t(locs_iir_zoom), iir_centered(locs_iir_zoom), 'g^', 'MarkerSize', 9, 'LineWidth', 1.8, 'DisplayName', 'IIR R-Peaks');
plot(ax4, t(locs_fir_zoom), fir_centered(locs_fir_zoom), 'ms', 'MarkerSize', 9, 'LineWidth', 1.8, 'DisplayName', 'FIR R-Peaks');

title(ax4, 'R-Peak Temporal Alignment (Zoomed View: t = 2s to 5s)', ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', title_color);
xlabel(ax4, 'Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax4, 'Amplitude (mV)', 'FontSize', 12, 'FontWeight', 'bold');
xlim(ax4, [2, 5]);