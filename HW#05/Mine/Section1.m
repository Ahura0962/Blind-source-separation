%% Initialization
clc; clear; close all;
load hw5.mat

%% Part 1
N_0 = 3;
fprintf("Part 1:\n")
tic
err_thr = 1e-6;
combinations = nchoosek(1:size(D, 2), N_0);
D_hat = zeros(size(D, 1), N_0, size(combinations, 1));
[S_hat, Err] = deal([]);
for i = 1:size(combinations, 1)
    D_hat(:, :, i) = D(:, combinations(i, :));
    s_hat = pinv(D_hat(:, :, i)) * x;
    err_func = x - D_hat(:, :, i) * s_hat;
    err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
    if err < err_thr
        S_hat = cat(1, S_hat, s_hat);
        Err = cat(1, Err, err);
    end
end
toc
fprintf("S_hat = [%f; %f; %f]\n", S_hat)
disp("Error = " + Err)
fprintf("\n\n")

%% Part 2
s_hat = pinv(D) * x;
plot(1:length(s_hat), s_hat)
title("$\left \|s\right \|_2$ Using Output", Interpreter="latex")
xlabel("$t$", Interpreter="latex")
ylabel("$\hat{s}$", Interpreter="latex", Rotation=0, FontSize=16)

%% Part 3
fprintf("Part 3:\n")
tic
x_r = x;
D_hat = zeros(length(x), N_0);
for i = 1:N_0
    [innerprodMax, ind] = max(x_r.' * D);
    D_hat(:, i) = D(:, ind);
    x_r = x_r - innerprodMax*D(:, ind);
end
S_hat = pinv(D_hat) * x;
err_func = x - D_hat * S_hat;
Err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
toc
fprintf("S_hat = [%f; %f; %f]\n", S_hat)
disp("Error = " + Err)
fprintf("\n")
% Without Knowing N_0
tic
err_thr = 1e-1;
x_r = x;
max_N_0 = size(D, 1);
D_hat = zeros(length(x), max_N_0);
for i = 1:max_N_0
    [innerprodMax, ind] = max(x_r.' * D);
    D_hat(:, i) = D(:, ind);
    x_r = x_r - innerprodMax*D(:, ind);
    s_hat = pinv(D_hat) * x;
    err_func = x - D_hat * s_hat;
    err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
    if err < err_thr
        break
    end
end
toc
disp("N_0 = " + i)
disp("Error = " + err)
fprintf("\n\n")

 %% Part 4
fprintf("Part 4:\n")
tic
x_r = x;
D_hat = zeros(length(x), N_0);
for i = 1:N_0
    [innerprodMax, ind] = max(x_r.' * D);
    D_hat(:, i) = D(:, ind);
    if i > 1
        D_pr = D_hat(:, 1:i);
        s_pr = pinv(D_pr) * x_r;
        x_r = x_r - D_pr * s_pr;
    else
        x_r = x_r - innerprodMax*D(:, ind);
    end
end
S_hat = pinv(D_hat) * x;
err_func = x - D_hat * S_hat;
Err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
toc
fprintf("S_hat = [%f; %f; %f]\n", S_hat)
disp("Error = " + Err)
fprintf("\n")
% Without Knowing N_0
tic
err_thr = 1e-1;
x_r = x;
D_hat = zeros(length(x), max_N_0);
for i = 1:max_N_0
    [innerprodMax, ind] = max(x_r.' * D);
    D_hat(:, i) = D(:, ind);
    if i > 1
        D_pr = D_hat(:, 1:i);
        s_pr = pinv(D_pr) * x_r;
        x_r = x_r - D_pr * s_pr;
    else
        x_r = x_r - innerprodMax*D(:, ind);
    end
    s_hat = pinv(D_hat) * x;
    err_func = x - D_hat * s_hat;
    err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
    if err < err_thr
        break
    end
end
toc
disp("N_0 = " + i)
disp("Error = " + err)
fprintf("\n\n")

%% Part 5
fprintf("Part 5:")
tic
D_1 = [D -D];
y = linprog(ones(1, size(D_1, 2)), [], [], D_1, x, zeros(1, size(D_1, 2)));
s = y(1: size(D, 2)) - y(size(D, 2) + 1: end);
S_hat = s(s ~= 0);
D_hat = x * pinv(S_hat);
err_func = x - D_hat * S_hat;
err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
toc
fprintf("S_hat = [%f; %f; %f]\n", S_hat)
disp("Error = " + err)
fprintf("\n\n")

%% Part 6
fprintf("Part 6:\n")
tic
w = randn(size(D, 2), 1);
iterNum = 100;
for iter = 1:iterNum
    w = diag(w);
    lambda = 2 * (D / w * D.') \ x;
    S_hat = 0.5 * w \ D.' * lambda;
    w = 1 ./ abs(S_hat);
    w(abs(S_hat) > 1e3) = 1e-10;
    w(abs(S_hat) < 1e-3) = 1e10;
end
D_hat = x * pinv(S_hat);
err_func = x - D_hat * S_hat;
err = sqrt(sum(err_func .^ 2)) / sqrt(sum(x .^ 2));
inds = find(abs(S_hat) > 1e-3);
s_hat = S_hat(inds);
toc
fprintf("S_hat = [%f; %f; %f]\n", s_hat)
disp("Error = " + err)
fprintf("\n\n")