clc; clear; close all;
%% Initialization
M = 10;
d = 1;
D = 0: d: M-1;
fc = 150e6;
theta1 = 10; theta2= 20;
Theta = [theta1; theta2];
f1 = 20e3; f2 = 10e3;
fs = 1e6;
k = 2 * pi * fc / 3e8;
t = 0: 1 / fs: 1e-3 - 1 / fs;
s_1 = exp(1j * 2*pi * f1 * t);
s_2 = exp(1j * 2*pi * f2 * t);
A = exp(-1j * k * D .* sin(deg2rad(Theta))).';
S = [s_1; s_2];
Noise = 1 * randn(M, length(S));
X = A*S + Noise;

%% Part 1
[Q, G, V] = svd(X);

%% Part 2
theta = linspace(0, 90, 451);
a_theta = exp(-1j * k * D.' * sin(deg2rad(theta)));
Q_sig = Q(:, 1:2);
aQ_sig = a_theta' * Q_sig;
beamforming = sqrt(sum(aQ_sig .* conj(aQ_sig),2)).';
figure
findpeaks(beamforming, theta)
grid off
title('beamforming Method', Interpreter='latex')
xlabel("$\theta$", Interpreter="latex")
ylabel("$\left \| a^H Q_{sig} \right \|$", Interpreter="latex")
[values, beamforming_peaks] = findpeaks(beamforming, theta);
[~, ind] = sort(values);
beamforming_peaks = beamforming_peaks(ind(end-1: end));

%% Part 3
Q_null = Q(:, 3:end);
aQ_null = a_theta' * Q_null;
music = 1 ./ sqrt(sum(aQ_null .* conj(aQ_null), 2)).';
figure
findpeaks(music, theta)
grid off
title('music Method', Interpreter='latex')
xlabel("$\theta$", Interpreter="latex")
ylabel("$\frac{1}{\left \| a^H Q_{null} \right \|}$", Interpreter="latex")
[values, music_peaks] = findpeaks(music, theta);
[~, ind] = sort(values);
music_peaks = music_peaks(ind(end-1: end));

%% Part 4
f_var = 3e2: 1e1 :3e4;
s_f = exp(1j * 2*pi * f_var.' * t);
V_sig = V(:, 1:2);
sV_sig = s_f * V_sig;
beamforming_f = sqrt(sum(sV_sig .* conj(sV_sig), 2)).';
figure
findpeaks(beamforming_f, f_var)
grid off
title('beamforming Method', Interpreter='latex')
xlabel("$f$", Interpreter="latex")
ylabel("$\left \| s(f) V_{sig} \right \|$", Interpreter="latex")
[values, beamforming_f_peaks] = findpeaks(beamforming_f, f_var);
[~, ind] = sort(values);
beamforming_f_peaks = beamforming_f_peaks(ind(end-1: end));

%% Part 5
V_null = V(:, 3:end);
sV_null = s_f * V_null;
music_f = 1 ./ sqrt(sum(sV_null .* conj(sV_null), 2)).';
figure
findpeaks(music_f, f_var)
grid off
title('music Method', Interpreter='latex')
xlabel("$f$", Interpreter="latex")
ylabel("$\frac{1}{\left \| s(f) V_{null} \right \|}$", Interpreter="latex")
[values, music_f_peaks] = findpeaks(music_f, f_var);
[~, ind] = sort(values);
music_f_peaks = music_f_peaks(ind(end-1: end));
