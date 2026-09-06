clc; clear; close all;

%% Initialization
load hw7.mat
iterNum = 1e2;
L = 100;
K = 5;

%% Part 1
T = length(x1);
figure
plot(x1, "m")
title("x1 plot", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_1(t)$", Interpreter="latex", Rotation=0, FontSize=13)
[s_hat, comb, x1_hat] = Single_Channel__Blind_Deconvolution(x1, L, K, T, iterNum);
figure
subplot(3, 1, 1)
plot(s_hat, "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(comb, "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
xlim([1 T-L+1])
ylabel("$\Psi$", Interpreter="latex", FontSize=13)
subplot(3, 1, 3)
plot(1:T, x1, 'm', 1:T, x1_hat, 'k')
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$x_1$", "$\hat{x_1}$", Interpreter="latex", Box="off")

%% Part 2
T = length(x2);
figure
plot(x2, "m")
title("x2 plot", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_2(t)$", Interpreter="latex", Rotation=0, FontSize=13)
[s_hat, comb, x2_hat] = Single_Channel__Blind_Deconvolution(x2, L, K, T, iterNum);
figure
subplot(3, 1, 1)
plot(s_hat, "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(comb, "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
xlim([1 T-L+1])
ylabel("$\Psi$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(1:T, x2, 'm', 1:T, x2_hat, 'k')
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$x_2$", "$\hat{x_2}$", Interpreter="latex", Box="off")

%% Part 3
T = length(X);
N = size(X, 1);
figure
subplot(2, 1, 1)
plot(X(1, :), "m")
title("X1 plot", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$X_1(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(2, 1, 2)
plot(X(2, :), "m")
title("X2 plot", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$X_2(t)$", Interpreter="latex", Rotation=0, FontSize=13)
A_hat = randn(N, N);
A_hat = A_hat ./ sqrt(sum(A_hat .^ 2, 2));
s_hat = zeros(L, 2);
comb = zeros(2, T-L+1);
x_hat = zeros(2, T);
for iter = 1:iterNum
    % A: fixed
    B_hat = pinv(A_hat) * X;
    for i = 1:N
        b_hat = B_hat(i, :);
        [s_hat(:, i), comb(i, :), x_hat(i, :)] = Single_Channel__Blind_Deconvolution(b_hat, L, K, T, iterNum);
        B_hat(i, :) = x_hat(i, :);
    end
    % B: fixed
    A_hat = X * pinv(B_hat);
    A_hat = A_hat ./ sqrt(sum(A_hat .^ 2, 2));
end
figure
subplot(3, 1, 1)
plot(s_hat(:,1), "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_1$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(comb(1, :), "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
xlim([1 T-L+1])
ylabel("$\Psi_1$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(1:T, X(1, :), 'm', 1:T, x_hat(1, :), 'k')
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$X_1$", "$\hat{X_1}$", Interpreter="latex", Box="off")
figure
subplot(3, 1, 1)
plot(s_hat(:,2), "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$\hat{s}_2$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(comb(2, :), "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
xlim([1 T-L+1])
ylabel("$\Psi_2$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(1:T, X(2, :), 'm', 1:T, x_hat(2, :), 'k')
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$X_2$", "$\hat{X_2}$", Interpreter="latex", Box="off")

figure
X_hat = A_hat * B_hat;
subplot(2, 1, 1)
plot(1:T, X(1, :), "m", 1:T, X_hat(1, :), "k")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$X_1$", "$\hat{X_1}$", Interpreter="latex", Box="off")
subplot(2, 1, 2)
plot(1:T, X(2, :), "m", 1:T, X_hat(2, :), "k")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$X_2$", "$\hat{X_2}$", Interpreter="latex", Box="off")

%% Part 4
T = length(x1);
x1_f = fft(x1);
corr_max = zeros(1, K);
ind_corr_max = zeros(1, K);
S_hat_rand = rand(1, T);
S_hat_f = fft(S_hat_rand);
S_hat_f = S_hat_f ./ vecnorm(S_hat_f);
for iter = 1:iterNum
    % s: fixed
    comb_f = S_hat_f ./ x1_f;
    comb_t = real(ifft(comb_f));
    comb_t(comb_t < 0) = 0;
    for i = 1:K
        [corr_max(i), ind_corr_max(i)] = max(comb_t);
        min_L = max(ind_corr_max(i) - L, 1);
        max_L = min(ind_corr_max(i) + L, T);
        comb_t(min_L: max_L) = 0;
    end
    [tau_k, sort_ind] = sort(ind_corr_max, 'ascend');
    alpha = corr_max(sort_ind);
    alpha(alpha < 1e-5) = 1e2;
    comb_t = zeros(1, T);
    comb_t(tau_k) = alpha;
    comb_f = fft(comb_t);
    % comb: fixed
    S_hat_f = x1_f ./ comb_f;
    S_hat_f = S_hat_f ./ vecnorm(S_hat_f);
end
s_hat = ifft(ifftshift(S_hat_f));
s_hat = s_hat ./ sqrt(sum(s_hat.^2));
figure
% subplot(2, 1, 1)
% plot(1:2000, s_hat, "m")
% xlabel("$t$", Interpreter="latex", FontSize=13)
% ylabel("$\hat{s}$", Interpreter="latex", Rotation=0, FontSize=13)
% subplot(2, 1, 2)
plot(comb_t, "m")
xlabel("$t$", Interpreter="latex", FontSize=13)
xlim([1 T])
ylabel("$\Psi$", Interpreter="latex", FontSize=13)