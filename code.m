% 16-QAM调制：星座图、时域波形与频域分析
clc;
clear;
close all;

% 输入比特流
bit_stream = '111001001011110110100000111001001011110110101101111000101010110110100000111001011010000010111000111001001011111110101001111000101010000010111101';

% 参数设置
M = 16;             % 16-QAM调制
k = log2(M);        % 每个符号的比特数
sps = 100;          % 每符号采样点数
fs = 1;             % 符号率（Hz）
total_symbols = length(bit_stream)/k; % 总符号数

% 验证比特流长度
if mod(length(bit_stream), k) ~= 0
    error('比特流长度必须是4的倍数');
end

% 将比特流分组为4位一组
symbols_bin = reshape(bit_stream, k, total_symbols)';

% 16-QAM映射表（格雷编码）
mapping_table = [-3-3j, -3-1j, -3+3j, -3+1j, ...
                 -1-3j, -1-1j, -1+3j, -1+1j, ...
                  3-3j,  3-1j,  3+3j,  3+1j, ...
                  1-3j,  1-1j,  1+3j,  1+1j];

% 计算归一化因子（使平均功率为1）
norm_factor = sqrt(mean(abs(mapping_table).^2));
fprintf('归一化因子: %.4f\n', norm_factor);

% 映射过程
symbols = zeros(1, total_symbols);
for i = 1:total_symbols
    bin_str = symbols_bin(i, :);
    dec_index = bin2dec(bin_str);
    symbols(i) = mapping_table(dec_index + 1) / norm_factor;
end

% 生成I、Q分量
I = real(symbols);
Q = imag(symbols);

% 上采样生成波形
t_waveform = (0:total_symbols*sps-1)/sps; % 波形时间索引
I_waveform = repelem(I, sps);
Q_waveform = repelem(Q, sps);
complex_signal = I_waveform + 1i*Q_waveform;

% 计算信号频谱
N = length(complex_signal);           % 信号长度
frequencies = (-N/2:N/2-1)*(sps/N);  % 频率轴
signal_spectrum = fftshift(fft(complex_signal)) / N;

%% 信号解调
% 下采样获取符号点
received_symbols = complex_signal(1:sps:end);

% 对接收到的符号进行判决
decoded_symbols = zeros(1, total_symbols);
decoded_bits = '';

for i = 1:total_symbols
    % 寻找最近的星座点
    distances = abs(mapping_table/norm_factor - received_symbols(i));
    [~, idx] = min(distances);
    decoded_symbols(i) = mapping_table(idx)/norm_factor;
    
    % 转换为比特
    bin_str = dec2bin(idx-1, k);
    decoded_bits = [decoded_bits bin_str];
end

%% 绘图 - 完整分析界面
figure('Position', [100, 100, 1200, 900], 'Name', '16-QAM调制分析', 'NumberTitle', 'off');

% ================= 星座图 =================
subplot(3, 3, [1,2]);
plot(symbols, 'o', 'MarkerSize', 6, 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'k');
grid on; axis equal;
xlim([-1.8 1.8]); ylim([-1.8 1.8]);
title('16-QAM星座图 (格雷编码)', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('同相分量 (I)', 'FontWeight', 'bold');
ylabel('正交分量 (Q)', 'FontWeight', 'bold');
set(gca, 'FontSize', 10, 'GridAlpha', 0.3);

% 标记象限和典型点
text(1.2, 1.2, 'Q1', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
text(-1.5, 1.2, 'Q2', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
text(-1.5, -1.5, 'Q3', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
text(1.2, -1.5, 'Q4', 'FontSize', 14, 'Color', 'r', 'FontWeight', 'bold');
hold on;
plot([0,0], ylim, 'k--', 'LineWidth', 0.5);
plot(xlim, [0,0], 'k--', 'LineWidth', 0.5);
hold off;

% ================= 星座点映射说明 =================
subplot(3, 3, 3);
axis off;
text(0.05, 0.95, '16-QAM 比特映射表', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b');
text(0.05, 0.85, '第一象限 (Q1): 1010, 1011, 1110, 1111', 'FontSize', 10);
text(0.05, 0.75, '第二象限 (Q2): 0010, 0011, 0110, 0111', 'FontSize', 10);
text(0.05, 0.65, '第三象限 (Q3): 0000, 0001, 0100, 0101', 'FontSize', 10);
text(0.05, 0.55, '第四象限 (Q4): 1000, 1001, 1100, 1101', 'FontSize', 10);
text(0.05, 0.40, sprintf('归一化因子: %.4f', norm_factor), 'FontSize', 10, 'BackgroundColor', [0.9 0.95 1]);
text(0.05, 0.30, sprintf('符号总数: %d', total_symbols), 'FontSize', 10);
text(0.05, 0.20, sprintf('采样率: %d 样点/符号', sps), 'FontSize', 10);
text(0.05, 0.10, sprintf('总样点数: %d', N), 'FontSize', 10);

% ================= 时域波形 =================
num_show_symbols = 10;
show_samples = num_show_symbols * sps;

% I分量时域波形
subplot(3, 3, 4);
plot(t_waveform(1:show_samples), I_waveform(1:show_samples), 'b', 'LineWidth', 1.5);
title('同相分量 (I) 时域波形', 'FontSize', 11);
xlabel('时间 (符号周期)', 'FontSize', 9);
ylabel('幅度', 'FontSize', 9);
grid on;
xlim([0 num_show_symbols]);
ylim([-1.8 1.8]);
set(gca, 'FontSize', 8, 'XTick', 0:num_show_symbols);

% 添加符号标记
hold on;
for i = 1:num_show_symbols
    plot([i-1, i-1], [-2, 2], 'k:', 'LineWidth', 0.5);
    text(i-0.5, 1.9, sprintf('S%d', i), 'FontSize', 8, 'HorizontalAlignment', 'center');
end
hold off;

% Q分量时域波形
subplot(3, 3, 5);
plot(t_waveform(1:show_samples), Q_waveform(1:show_samples), 'r', 'LineWidth', 1.5);
title('正交分量 (Q) 时域波形', 'FontSize', 11);
xlabel('时间 (符号周期)', 'FontSize', 9);
ylabel('幅度', 'FontSize', 9);
grid on;
xlim([0 num_show_symbols]);
ylim([-1.8 1.8]);
set(gca, 'FontSize', 8, 'XTick', 0:num_show_symbols);

% 添加符号标记
hold on;
for i = 1:num_show_symbols
    plot([i-1, i-1], [-2, 2], 'k:', 'LineWidth', 0.5);
    text(i-0.5, 1.9, sprintf('S%d', i), 'FontSize', 8, 'HorizontalAlignment', 'center');
end
hold off;

% 复合信号包络
subplot(3, 3, 6);
plot(t_waveform(1:show_samples), abs(complex_signal(1:show_samples)), ...
    'Color', [0.6 0.2 0.8], 'LineWidth', 1.5);
title('复合信号包络', 'FontSize', 11);
xlabel('时间 (符号周期)', 'FontSize', 9);
ylabel('幅度', 'FontSize', 9);
grid on;
xlim([0 num_show_symbols]);
ylim([0 1.8]);
set(gca, 'FontSize', 8, 'XTick', 0:num_show_symbols);

% ================= 频域分析 =================
% 信号频谱幅度
subplot(3, 3, 7);
plot(frequencies, abs(signal_spectrum), 'LineWidth', 1.5, 'Color', [0 0.6 0.3]);
title('信号频谱幅度', 'FontSize', 11);
xlabel('频率 (Hz)', 'FontSize', 9);
ylabel('幅度', 'FontSize', 9);
grid on;
xlim([-1.5, 1.5]);
set(gca, 'FontSize', 8);

% 信号频谱相位
subplot(3, 3, 8);
plot(frequencies, angle(signal_spectrum)*180/pi, 'LineWidth', 1.5, 'Color', [0.9 0.5 0]);
title('信号频谱相位', 'FontSize', 11);
xlabel('频率 (Hz)', 'FontSize', 9);
ylabel('相位 (度)', 'FontSize', 9);
grid on;
xlim([-1.5, 1.5]);
set(gca, 'FontSize', 8);
yticks(-180:45:180);

% 信号功率谱密度
subplot(3, 3, 9);
window = hamming(512);
noverlap = 256;
nfft = 1024;
[Pxx, F] = pwelch(complex_signal, window, noverlap, nfft, sps);
plot(F, 10*log10(Pxx), 'LineWidth', 1.5, 'Color', [0.8 0.2 0.4]);
title('功率谱密度 (PSD)', 'FontSize', 11);
xlabel('频率 (Hz)', 'FontSize', 9);
ylabel('功率/频率 (dB/Hz)', 'FontSize', 9);
grid on;
xlim([0, 1.5]);
set(gca, 'FontSize', 8);

% 添加总标题
sgtitle('16-QAM调制系统全面分析：星座图、时域波形与频域特性', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.1 0.2 0.5]);

% 添加统计信息
info_str = sprintf('比特流长度: %d 比特 | 符号数: %d\n采样率: %d 样点/符号 | 总采样点数: %d', ...
                   length(bit_stream), total_symbols, sps, N);
annotation('textbox', [0.1, 0.01, 0.8, 0.03], 'String', info_str, ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10, ...
           'BackgroundColor', [0.95 0.95 0.95], 'HorizontalAlignment', 'center');

%% 添加比特流波形对比图
figure('Position', [150, 150, 1000, 400], 'Name', '比特流对比分析', 'NumberTitle', 'off');

% 创建时间轴
t_bits = (0:length(bit_stream)-1);

% 将字符串转换为数值数组
orig_bits = zeros(1, length(bit_stream));
for i = 1:length(bit_stream)
    orig_bits(i) = str2num(bit_stream(i));
end

decoded_bits_array = zeros(1, length(decoded_bits));
for i = 1:length(decoded_bits)
    decoded_bits_array(i) = str2num(decoded_bits(i));
end

% 原始比特流
subplot(2,1,1);
stem(t_bits, orig_bits, 'filled', 'LineWidth', 1.5, 'MarkerSize', 4);
title('原始比特流', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('比特索引', 'FontSize', 10);
ylabel('比特值', 'FontSize', 10);
grid on;
ylim([-0.2 1.2]);
set(gca, 'YTick', [0 1]);
set(gca, 'FontSize', 9);

% 解调后比特流
subplot(2,1,2);
stem(t_bits, decoded_bits_array, 'filled', 'LineWidth', 1.5, 'MarkerSize', 4, 'Color', [0.8 0.2 0.2]);
title('解调后比特流', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('比特索引', 'FontSize', 10);
ylabel('比特值', 'FontSize', 10);
grid on;
ylim([-0.2 1.2]);
set(gca, 'YTick', [0 1]);
set(gca, 'FontSize', 9);

% 添加比特错误统计
bit_errors = sum(orig_bits ~= decoded_bits_array);
error_rate = bit_errors / length(bit_stream);

sgtitle(sprintf('比特流对比分析 (错误率: %.2f%%)', error_rate*100), ...
    'FontSize', 14, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.6]);

% 添加统计信息
info_str = sprintf('总比特数: %d | 错误比特数: %d | 错误率: %.2f%%', ...
                   length(bit_stream), bit_errors, error_rate*100);
annotation('textbox', [0.1, 0.02, 0.8, 0.03], 'String', info_str, ...
           'FitBoxToText', 'on', 'EdgeColor', 'none', 'FontSize', 10, ...
           'BackgroundColor', [0.95 0.95 0.95], 'HorizontalAlignment', 'center');
