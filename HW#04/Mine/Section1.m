clc; clear; close all;

%% Initialization
fs = 20;
c = [0.2 0.4 0.6 -0.1 -0.3];
d = [0.1 0.3 -0.2 0.5 -0.3];
k = 1:5;
t0 = 0:1/fs:1-1/fs;
t = (k.' - 1) + t0;

%% Part 1
A = [0.8 -0.6; 0.6 0.8];
S = [reshape(c .* sin(2*pi*t).', 1, []);
     reshape(d .* sin(4*pi*t).', 1, [])];
X = A * S;
figure 
t_plot = reshape(t.', 1, []);
subplot(2, 1, 1)
plot(t_plot, S(1, :))
title("First Source Signal", Interpreter="latex",FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$s_1$", Interpreter="latex",Rotation=0,FontSize=13)
subplot(2, 1, 2)
plot(t_plot, S(2, :))
title("Second Source Signal", Interpreter="latex",FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$s_2$", Interpreter="latex",Rotation=0,FontSize=13)

figure 
subplot(2, 1, 1)
plot(t_plot, X(1, :))
title("First Observation Signal", Interpreter="latex",FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_1$", Interpreter="latex",Rotation=0,FontSize=13)
subplot(2, 1, 2)
plot(t_plot, X(2, :))
title("Second Observation Signal", Interpreter="latex",FontSize=12)
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_2$", Interpreter="latex",Rotation=0,FontSize=13)

%% Part 2
window1_X = X(:, 1:1*length(t));
window1_R_x = window1_X * window1_X.';
window2_X = X(:, 1*length(t) + 1: 2*length(t));
window2_R_x = window2_X * window2_X.';
[B_T, Gamma] = eig(window1_R_x, window2_R_x);
[nuDiag , order] = sort(diag(Gamma), 'descend'); 
Gamma = diag(nuDiag);
B_T = B_T(:,order);
S_hat = B_T.' * X;
S_norm = S ./ sqrt(sum(S.^2, 2));
S_hat_norm = S_hat ./ sqrt(sum(S_hat.^2, 2));
R_s_sHat = S_norm * S_hat_norm.';
S_hat_norm = (2*(sum(R_s_sHat) > 0) - 1).' .* S_hat_norm;
figure
sgtitle("Part 2", Interpreter="latex", Color='m')
subplot(2, 1, 1)
plot(t_plot, S_norm(1,:), 'm', t_plot, S_hat_norm(1,:), '.k')
ylim([-0.3 0.3])
title("First Predicted Source Signal", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$s_1$", "$\hat{s}_1$", Interpreter="latex", FontSize=13, Location="northwest", Box="off")
subplot(2, 1, 2)
plot(t_plot, S_norm(2,:), 'm', t_plot, S_hat_norm(2,:), '.k')
ylim([-0.3 0.3])
title("Second Predicted Source Signal", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$s_2$", "$\hat{s}_2$", Interpreter="latex", FontSize=13, Location="northwest", Box="off")
E1 = sum((S_hat_norm - S_norm).^2, 2);

%% Part 3
observations_num = size(X, 1);
Bhat = orth(randn(observations_num));
[E2, S_hat_norm, func] = JD(Bhat, X, S_norm, 5);
figure
plot(func)
title("Error of Alternation Minimization", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
figure
sgtitle("Part 3", Interpreter="latex", Color='m')
subplot(2, 1, 1)
plot(t_plot, S_norm(1,:), 'm', t_plot, S_hat_norm(1,:), '.k')
ylim([-0.3 0.3])
title("First Predicted Source Signal", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$s_1$", "$\hat{s}_1$", Interpreter="latex", FontSize=13, Location="northwest", Box="off")
subplot(2, 1, 2)
plot(t_plot, S_norm(2,:), 'm', t_plot, S_hat_norm(2,:), '.k')
ylim([-0.3 0.3])
title("Second Predicted Source Signal", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$s_2$", "$\hat{s}_2$", Interpreter="latex", FontSize=13, Location="northwest", Box="off")

%% Part 4
SNR_dB = 20;
SNR_linear = 10 ^ (SNR_dB / 10);
power_X = sum(X.^2, 2);
sDeviation = sqrt(power_X / SNR_linear);
func = zeros(100, 100, 4);
E_list = zeros(2, 1, 4);
for runNum = 1:100
    W = randn(2, 100);
    W = W ./ sqrt(sum(W.^2, 2));
    Y = X + sDeviation .* W;
    [E_1, ~, func(runNum, :, 1)] = JD(Bhat, Y, S_norm, 2);
    [E_2, ~, func(runNum, :, 2)] = JD(Bhat, Y, S_norm, 3);
    [E_3, ~, func(runNum, :, 3)] = JD(Bhat, Y, S_norm, 4);
    [E_4, ~, func(runNum, :, 4)] = JD(Bhat, Y, S_norm, 5);
    E_list(:, :, 1) = E_list(:, :, 1) + E_1 ./ 100;
    E_list(:, :, 2) = E_list(:, :, 1) + E_2 ./ 100;
    E_list(:, :, 3) = E_list(:, :, 1) + E_3 ./ 100;
    E_list(:, :, 4) = E_list(:, :, 1) + E_4 ./ 100;
end
func = mean(func, 1);
figure
sgtitle("Part 4", Interpreter="latex", Color='m')
subplot(2, 2, 1)
plot(func(:, :, 1))
title("Error of A.M for $k = 2$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 2)
plot(func(:, :, 2))
title("Error of A.M for $k = 3$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 3)
plot(func(:, :, 3))
title("Error of A.M for $k = 4$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 4)
plot(func(:, :, 4))
title("Error of A.M for $k = 5$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)

%% Part 5
SNR_dB = 5:5:20;
SNR_linear = 10 .^ (SNR_dB ./ 10);
sDeviation = sqrt(power_X ./ SNR_linear);
func = zeros(100, 100, 4);
E_list = zeros(2, 1, 4);
for runNum = 1:100
    W = randn(2, 100);
    W = W ./ sqrt(sum(W.^2, 2));
    Y = X + reshape(sDeviation, 2, [], 4).* W;
    [E_1, ~, func(runNum, :, 1)] = JD(Bhat, Y(:, :, 1), S_norm, 5);
    [E_2, ~, func(runNum, :, 2)] = JD(Bhat, Y(:, :, 2), S_norm, 5);
    [E_3, ~, func(runNum, :, 3)] = JD(Bhat, Y(:, :, 3), S_norm, 5);
    [E_4, ~, func(runNum, :, 4)] = JD(Bhat, Y(:, :, 4), S_norm, 5);
    E_list(:, :, 1) = E_list(:, :, 1) + E_1 ./ 100;
    E_list(:, :, 2) = E_list(:, :, 1) + E_2 ./ 100;
    E_list(:, :, 3) = E_list(:, :, 1) + E_3 ./ 100;
    E_list(:, :, 4) = E_list(:, :, 1) + E_4 ./ 100;
end
func = mean(func, 1);
figure
sgtitle("Part 5", Interpreter="latex", Color='m')
subplot(2, 2, 1)
plot(func(:, :, 1))
title("Error of A.M for $SNR = 5 \left[dB\right]$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 2)
plot(func(:, :, 2))
title("Error of A.M for $SNR = 10 \left[dB\right]$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 3)
plot(func(:, :, 3))
title("Error of A.M for $SNR = 15 \left[dB\right]$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)
subplot(2, 2, 4)
plot(func(:, :, 4))
title("Error of A.M for $SNR = 20 \left[dB\right]$", Interpreter="latex")
xlabel("$Iteration$", Interpreter="latex", FontSize=13)
ylabel("$Error$", Interpreter="latex", FontSize=13)

%% Functions
function [E, S_hat_norm, func] = JD(Bhat, X, S_norm, k)
    X_paged = reshape(X, 2, [], 5);
    X_paged = X_paged(:, :, 1:k);
    R_x = pagemtimes(X_paged, pagetranspose(X_paged));
    observations_num = size(X, 1);
    
    iterNum = 100;
    func = zeros(1, iterNum);
    for itr = 1:iterNum
        f = 0;
        for i = 1:observations_num
            i_t = 1:observations_num;
            j = i_t ~= i;
            b_j = Bhat(j, :).';
            Rb = pagemtimes(R_x, b_j);
            R_i = sum(pagemtimes(Rb, pagetranspose(Rb)), 3);
            [U, ~] = eig(R_i);
            Bhat(i, :) = U(:, 1).';
            if i > 1
                BB = Bhat(1:i-1, :).';
                Bhat(i, :) = (eye(observations_num) - BB * BB.') * Bhat(i, :).';
                Bhat(i, :) = Bhat(i, :) / norm(Bhat(i, :));
            end
            f = f + Bhat(i, :) * R_i * Bhat(i, :).';
        end
        func(itr) = f;
    end
    S_hat = Bhat * X;
    S_hat_norm = S_hat ./ sqrt(sum(S_hat.^2, 2));
    R_s_sHat = S_norm * S_hat_norm.';
    if (abs(R_s_sHat(1,1)) < abs(R_s_sHat(1,2)))
        S_hat_norm = [S_hat_norm(2, :); S_hat_norm(1, :)];
    end
    S_hat_norm = (2*(sum(R_s_sHat) > 0) - 1).' .* S_hat_norm;
    E = sum((S_hat_norm - S_norm).^2, 2);
end