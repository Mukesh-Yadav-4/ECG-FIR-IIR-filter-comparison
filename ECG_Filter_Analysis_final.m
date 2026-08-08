%% ECG Filter Comparison Project - FIR vs IIR Denoising
% Author: MUKESH YADAV
% Description: Compares FIR and IIR bandpass filtering for ECG denoising
% using MIT-BIH Arrhythmia Database record 100

%% Step 1: Load raw ECG signal (manual read - format 212)
fid = fopen('100.dat', 'r');
raw = fread(fid, Inf, 'ubit12');
fclose(fid);
raw = reshape(raw, 2, [])';
raw(raw > 2047) = raw(raw > 2047) - 4096;

gain = 200;
baseline = 1024;
fs = 360;
signal = (raw - baseline) / gain;
tm = (0:length(signal)-1) / fs;

figure;
plot(tm, signal(:,1));
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('MIT-BIH Record 100 - Raw ECG');

%% Step 2: Add synthetic noise
t = tm(:);
clean_signal = signal(:,1);

powerline_noise = 0.05 * sin(2*pi*50*t);
baseline_wander = 0.15 * sin(2*pi*0.3*t);
gaussian_noise = 0.05 * randn(size(t));
noisy_signal = clean_signal + powerline_noise + baseline_wander + gaussian_noise;

figure;
subplot(2,1,1); plot(t, clean_signal); title('Clean ECG'); ylabel('mV');
subplot(2,1,2); plot(t, noisy_signal); title('Noisy ECG'); xlabel('Time (s)'); ylabel('mV');

%% Step 3: IIR filter (Butterworth bandpass)
nyquist = fs/2;
[b_iir, a_iir] = butter(4, [0.5/nyquist, 40/nyquist], 'bandpass');
iir_filtered = filtfilt(b_iir, a_iir, noisy_signal);

%% Step 4: FIR filter (equiripple bandpass)
fir_filter = designfilt('bandpassfir', 'FilterOrder', 100, ...
    'CutoffFrequency1', 0.5, 'CutoffFrequency2', 40, 'SampleRate', fs);
fir_filtered = filtfilt(fir_filter, noisy_signal);

%% Step 5: Frequency response check
figure;
freqz(b_iir, a_iir, 1024, fs);

%% Step 6: Mean-center all signals (removes DC offset mismatch)
clean_centered = clean_signal - mean(clean_signal);
noisy_centered = noisy_signal - mean(noisy_signal);
iir_centered = iir_filtered - mean(iir_filtered);
fir_centered = fir_filtered - mean(fir_filtered);

figure;
subplot(4,1,1); plot(t, clean_centered); title('Clean ECG (centered)'); ylabel('mV');
subplot(4,1,2); plot(t, noisy_centered); title('Noisy ECG (centered)'); ylabel('mV');
subplot(4,1,3); plot(t, iir_centered); title('IIR Filtered (centered)'); ylabel('mV');
subplot(4,1,4); plot(t, fir_centered); title('FIR Filtered (centered)'); ylabel('mV'); xlabel('Time (s)');

%% Step 7: SNR calculation (manual, using centered signals)
snr_manual = @(clean, test) 10*log10(sum(clean.^2) / sum((clean - test).^2));

snr_noisy = snr_manual(clean_centered, noisy_centered);
snr_iir = snr_manual(clean_centered, iir_centered);
snr_fir = snr_manual(clean_centered, fir_centered);

fprintf('SNR - Noisy: %.2f dB\n', snr_noisy);
fprintf('SNR - IIR: %.2f dB\n', snr_iir);
fprintf('SNR - FIR: %.2f dB\n', snr_fir);

%% Step 8: R-peak detection accuracy
min_peak_height = 0.3;
min_peak_distance = round(0.4 * fs);

[~, locs_clean] = findpeaks(clean_centered, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
[~, locs_iir] = findpeaks(iir_centered, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
[~, locs_fir] = findpeaks(fir_centered, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

tolerance = round(0.05 * fs);

matched_iir = 0;
for i = 1:length(locs_clean)
    if any(abs(locs_iir - locs_clean(i)) <= tolerance)
        matched_iir = matched_iir + 1;
    end
end

matched_fir = 0;
for i = 1:length(locs_clean)
    if any(abs(locs_fir - locs_clean(i)) <= tolerance)
        matched_fir = matched_fir + 1;
    end
end

fprintf('R-peak accuracy - IIR: %.3f%%\n', (matched_iir/length(locs_clean))*100);
fprintf('R-peak accuracy - FIR: %.3f%%\n', (matched_fir/length(locs_clean))*100);

%% Step 9: Save figures
saveas(figure(1), 'fig1_raw_ecg.png');
saveas(figure(2), 'fig2_clean_vs_noisy.png');
saveas(figure(3), 'fig3_freq_response.png');
saveas(figure(4), 'fig4_all_four_centered.png');

%% Step 10: Save workspace
save('ECG_project_workspace.mat');