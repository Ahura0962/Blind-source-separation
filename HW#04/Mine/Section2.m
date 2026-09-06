clc; clear; close all;

%% Initializzation
load hw4-X1.mat
load hw4-X2.mat
fs = 100;
t_plot = 0: 1/fs : length(X1)/fs - 1/fs;

%% Part 1
X1_shifted1 = circshift(X1, fs/4, 1);
X1_shifted2 = circshift(X1, fs/5, 1);
Rx_11 = X1 * X1_shifted1.';
Rx_12 = X1 * X1_shifted2.';
[B_T, Gamma] = eig(Rx_11, Rx_12);
[nuDiag , order] = sort(diag(Gamma), 'descend'); 
Gamma = diag(nuDiag);
B_T = B_T(:,order);
S1_hat = B_T.' * X1;
S1_hat = S1_hat ./ sqrt(sum(S1_hat.^2, 2));
subplot(2, 1, 1)
plot(t_plot, S1_hat(1, :))
title("First Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(t_plot, S1_hat(2, :))
title("Second Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_2$", Interpreter="latex",Rotation=0, FontSize=13)

%% Part 2
f = -fs/2 : fs / length(S1_hat): fs/2 - 1/length(S1_hat);
S1_hat_fft = [fftshift(fft(S1_hat(1,:))) ; fftshift(fft(S1_hat(2,:)))];
figure
subplot(2, 1, 1)
plot(f, abs(S1_hat_fft(1, :)))
title("Fourier of First Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$\hat{S}_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(f, abs(S1_hat_fft(2, :)))
title("Fourier of Second Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$\hat{S}_2$", Interpreter="latex",Rotation=0, FontSize=13)
X1_fft = [fftshift(fft(X1(1,:))) ; fftshift(fft(X1(2,:)))];
figure
subplot(2, 1, 1)
plot(f, abs(X1_fft(1, :)))
title("Fourier of First Observation Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$X_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(f, abs(X1_fft(2, :)))
title("Fourier of Second Observation Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$X_2$", Interpreter="latex",Rotation=0, FontSize=13)

%% Part 3
X2_shifted1 = circshift(X2, fs/4, 1);
X2_shifted2 = circshift(X2, fs/5, 1);
Rx_21 = X1 * X2_shifted1.';
Rx_22 = X1 * X2_shifted2.';
[B_T, Gamma] = eig(Rx_21, Rx_22);
[nuDiag , order] = sort(diag(Gamma), 'descend'); 
Gamma = diag(nuDiag);
B_T = B_T(:,order);
S2_hat = B_T.' * X2;
S2_hat = S2_hat ./ sqrt(sum(S2_hat.^2, 2));
figure
subplot(2, 1, 1)
plot(t_plot, S2_hat(1, :))
title("First Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(t_plot, S2_hat(2, :))
title("Second Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_2$", Interpreter="latex",Rotation=0, FontSize=13)

%% Part 4
S2_hat_fft = [fftshift(fft(S2_hat(1,:))) ; fftshift(fft(S2_hat(2,:)))];
figure
subplot(2, 1, 1)
plot(f, abs(S2_hat_fft(1, :)))
title("Fourier of First Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$\hat{S}_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(f, abs(S2_hat_fft(2, :)))
title("Fourier of Second Predicted Source Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$\hat{S}_2$", Interpreter="latex",Rotation=0, FontSize=13)
X1_fft = [fftshift(fft(X1(1,:))) ; fftshift(fft(X1(2,:)))];
figure
subplot(2, 1, 1)
plot(f, abs(X1_fft(1, :)))
title("Fourier of First Observation Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$X_1$", Interpreter="latex",Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(f, abs(X1_fft(2, :)))
title("Fourier of Second Observation Signal", Interpreter="latex", FontSize=12)
xlabel("$f$", Interpreter="latex", FontSize=13)
ylabel("$X_2$", Interpreter="latex",Rotation=0, FontSize=13)