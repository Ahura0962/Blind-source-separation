clc; clear; close all;

%% Initialization
load hw6-part1.mat

%% Part 1
[M, N] = size(D);
fprintf("Part 1:\n")
D_1 = [D -D];
y = linprog(ones(1, 2*N), [], [], D_1, x, zeros(1, 2*N));
s = y(1: N) - y(N + 1: end);
S_hat = s(s ~= 0);
D_hat = x * pinv(S_hat);
err_func = x - D_hat * S_hat;
err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
disp("Error (Without Noise) = " + err)
fprintf("\n\n\n")

%% Part 2
fprintf("Part 2:\n")
err_thr = 1e-2;
x_r = x;
max_N_0 = M;
D_hat = zeros(M, N);
for i = 1:max_N_0
    [innerprodMax, ind] = max(x_r.' * D);
    D_hat(:, ind) = D(:, ind);
    x_r = x_r - innerprodMax * D(:, ind);
    S_hat = pinv(D_hat) * x;
    err_func = x - D_hat * S_hat;
    err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
    if err < err_thr
        break
    end
end
disp("Error (without Noise) = " + err)
fprintf("\n\n\n")

%% Part 3
fprintf("Part 3:\n")
D_1 = [D -D];
y = linprog(ones(1, size(D_1, 2)), [], [], D_1, x_noisy, zeros(1, size(D_1, 2)));
s = y(1: size(D, 2)) - y(size(D, 2) + 1: end);
S_hat = s(s ~= 0);
D_hat = x_noisy * pinv(S_hat);
err_func = x - D_hat * S_hat;
err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
disp("Error (with Noise) = " + err)
fprintf("\n\n\n")

%% Part 4
fprintf("Part 4:\n")
lambda = 0.125;
lambda_noisy = 4.2;
iterNum = 100;
s = randn(N, 1);
s_noisy = randn(N, 1);
for iter = 1:iterNum
    for i = 1:N
        D_filtered = [D(:, 1:i-1) D(:, i+1:end)];
        s_filtered = [s(1:i-1, :); s(i+1:end, :)];
        r = x - D_filtered * s_filtered;
        r_noisy = x_noisy - D_filtered * s_filtered;
        rho = r.' * D(:, i);
        rho_noisy = r_noisy.' * D(:, i);
        if rho > lambda / 2
            s(i) = rho - lambda / 2;
        elseif rho < -lambda / 2
            s(i) = rho + lambda / 2;
        else
            s(i) = 0;
        end
        if rho_noisy > lambda_noisy / 2
            s_noisy(i) = rho_noisy - lambda_noisy / 2;
        elseif rho_noisy < -lambda_noisy / 2
            s_noisy(i) = rho_noisy + lambda_noisy / 2;
        else
            s_noisy(i) = 0;
        end
    end
end
D_hat = x * pinv(s);
D_hat_noisy = x_noisy * pinv(s_noisy);
err_func = x - D_hat * s;
err_func_noisy = x - D_hat_noisy * s_noisy;
err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
err_noisy = sqrt(sum(err_func_noisy .^ 2)) / sqrt(sum(x .^ 2));
disp("Error (without Noise) = " + err)
disp("Error (with Noise) = " + err_noisy)
fprintf("\n\n\n")