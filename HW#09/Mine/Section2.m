clc; clear; close all;

%% Initialization
load hw9.mat
warning("off")

%% Part 0
X_noNoise = A*S;
X = A*S + Noise;

%% Part 1
mu = .01;
B = randn(size(A));
B = B ./ vecnorm(B, 2, 2);
[M, T] = size(X);
K_Y = @(Y) [ones(1,T); Y; Y.^2; Y.^3; Y.^4; Y.^5];
dK_Y = @(Y) [zeros(1,T); ones(1,T); 2*Y; 3*Y.^2; 4*Y.^3; 5*Y.^4];
iterNum = 5e2;
score_func = zeros(1, iterNum);
for iter = 1:iterNum
    Y = B * X;
    Psi_hat = zeros(size(X));
    for i = 1:size(Y, 1)
        K_y = K_Y(Y(i, :));
        Theta_hat = (K_y * K_y.' / T) \ mean(dK_Y(Y(i,:)), 2);
        Psi_hat(i, :) = Theta_hat.' * K_y;
    end
    D = Psi_hat * Y.' ./ T;
    D = D - diag(diag(D));
    df_B = D * B;
    B = B - mu * df_B;
    B = B ./ vecnorm(B, 2, 2);
    S_hat = B * X;
    S_hat = S_hat ./ vecnorm(S_hat, 2, 2);
    S_hatt = S_hat;
    S_temp = S;
    for i = 1:size(Y, 1)
        corr = S_hatt(i,:) * S_temp.';
        [~, ind] = max(abs(corr));
        S_temp(ind, :) = 0;
        S_hat(ind, :) = sign(corr(ind)) * S_hatt(i, :);
        S_hat(ind, :) = S_hat(ind, :) * norm(S(i, :), "fro") ./ norm(S_hat(ind, :), "fro");
    end
    score_func(iter) = norm(df_B, "fro");
    E = norm(S_hat - S, "fro")^2 / norm(S, "fro")^2;
end
permutation_matrix = B * A;

%% Part 2
disp("Part 2:")
disp("Error (E) = " + E)
figure
subplot(3, 1, 1)
plot(1:T, S_hat(1, :), "m", 1:T, S(1,:), "k")
title("S(1, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$\hat{s}_1(t)$", "$s_1(t)$", Interpreter="latex", Box="off")
subplot(3, 1, 2)
plot(1:T, S_hat(2, :), "m", 1:T, S(2,:), "k")
title("S(2, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$\hat{s}_2(t)$", "$s_2(t)$", Interpreter="latex", Box="off")
subplot(3, 1, 3)
plot(1:T, S_hat(3, :), "m", 1:T, S(3,:), "k")
title("S(3, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
legend("$\hat{s}_3(t)$", "$s_3(t)$", Interpreter="latex", Box="off")

%% Part 3
figure
plot(1:iterNum, score_func)
title("Convergence Plot", Interpreter="latex")
xlabel("Iteration", Interpreter="latex", FontSize=13)
ylabel("Score Function", Interpreter="latex", FontSize=13)