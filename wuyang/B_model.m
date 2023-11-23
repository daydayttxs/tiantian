clear all;clc
% load E:\user-MATLAB\Test\实测\256合成孔径间隔1\data3\attenuation_model\tiqu_data\RFData_SA1.mat
% load E:\user-MATLAB\Test\实测\256合成孔径间隔1\测试体模能否使用数据\RFData_SA1_compare
load D:\实验室\Bmodule\实测脚本\5.15\data.mat
number_scan_lines=256;
length_input_signal=3200;
sample_frequency=11.9e6;
% sample_frequency=17.8571e6;
scan_lines = zeros(number_scan_lines, length_input_signal);
% scan_lines=RFData_SA1_wufei3';
scan_lines=RFData_SA1';
tone_burst_freq=3e6;
% tone_burst_freq=2.1e6;
c0=1500;
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(length_input_signal * 2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, 50 * 2), scan_line_win(1:length_input_signal - 50 * 2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);
% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(1, :) = scan_lines(end/2, :);

% -----------------------------
% Time Gain Compensation
% -----------------------------
% create radius variable assuming that t0 corresponds to the middle of the
% input signal
t0 = length_input_signal * 1/sample_frequency / 2;
r = c0 * ( (1:length_input_signal) * 1/sample_frequency - t0 ) / 2;    % [m]
% define absorption value and convert to correct units
tgc_alpha_db_cm = medium.alpha_coeff * (tone_burst_freq * 1e-6)^medium.alpha_power;
% tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 * 50;
tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 *40;%60
% create time gain compensation function based on attenuation value and
% round trip distance
tgc = exp(tgc_alpha_np_m * 3 * r);
% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);
% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(2, :) = scan_lines(end/2, :);

% -----------------------------
% Frequency Filtering
% -----------------------------
% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, sample_frequency, tone_burst_freq, 20, true);
% set(gca, 'XLim', [0, 6]);
scan_lines_harm = gaussianFilter(scan_lines, sample_frequency, 2 * tone_burst_freq, 10, true);
% set(gca, 'XLim', [0, 6]);
% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(3, :) = scan_lines_fund(end/2, :);

% Envelope Detection
% -----------------------------
% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);
% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(4, :) = scan_lines_fund(end/2, :);
% -----------------------------
% Log Compression
% -----------------------------
% normalised log compression
compression_ratio = 12;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);
% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(5, :) = scan_lines_fund(end/2, :);
% -----------------------------
% Scan Conversion
% -----------------------------
% upsample the image using linear interpolation
scale_factor = 2;
scan_lines_fund = interp2(1:length_input_signal, (1:number_scan_lines).', scan_lines_fund, 1:length_input_signal, (1:1/scale_factor:number_scan_lines).');
scan_lines_harm = interp2(1:length_input_signal, (1:number_scan_lines).', scan_lines_harm, 1:length_input_signal, (1:1/scale_factor:number_scan_lines).');
%%
subplot(2,3,1)
plot(scan_line_example(1,:));
title('原始数据')
subplot(2,3,2)
plot(scan_line_example(2,:));
title('时间补偿')
subplot(2,3,3)
plot(scan_line_example(3,:));
title('高斯滤波')
subplot(2,3,4)
plot(scan_line_example(4,:));
title('包络检测')
subplot(2,3,5) 
plot(scan_line_example(5,:));
title('对数压缩')
myfun 