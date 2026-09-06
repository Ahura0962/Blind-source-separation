clc; clear; close all;

%% Initialization
load hw10.mat
warning("off")
thr = 1e-15;

%% Part 0
X_noNoise = A*S;
X = A*S + Noise;

%% Part 1
R_x = X * X .';
[U, D] = eig(R_x);
[nuDiag , order] = sort(diag(D), 'descend'); 
D = diag(nuDiag);
U = U(:, order);
W = D^(-1/2) * U.';
Z = W * X;
B = randn(size(A));
B = B ./ vecnorm(B, 2, 2);
[M, T] = size(X);
iterNum = 5e4;
mu = 4;
score_func = zeros(M, iterNum);
for row = 1:size(B, 1)
    b = B(row, :).';
    y = b.' * Z;
    v = randn(1, T);
    for iter = 1:iterNum
        f = mean(-exp(-y.^2/2) - mean(exp(-v.^2/2)))^2;
        g = y .* exp(-y.^2/2);
        df = sqrt(f) * (Z * g.');
        b = b - mu * df;
        b = (eye(size(B)) - B(1:row-1, :).' * B(1:row-1, :)) * b;
        b = b ./ norm(b,"fro");
        score_func(row, iter) = f;
        y = b.' * Z;
    end
    B(row,:) = b.';
end
S_hat = B * X;
S_hat = S_hat ./ norm(S_hat, "fro");
S_hatt = S_hat;
S_temp = S;
for i = 1:size(B, 1)
    corr = S_hatt(i,:) * S_temp.';
    [~, ind] = max(abs(corr));
    S_temp(ind, :) = 0;
    S_hat(ind, :) = sign(corr(ind)) * S_hatt(i, :);
    S_hat(ind, :) = S_hat(ind, :) * norm(S(i, :), "fro") ./ norm(S_hat(ind, :), "fro");
end
E = norm(S_hat - S, "fro")^2 / norm(S, "fro")^2;
permutation_matrix = B * W * A;

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
subplot(3, 1, 1)
plot(1:iterNum, score_func(1, :))
title("$b_1$ Convergence Plot", Interpreter="latex")
subplot(3, 1, 2)
plot(1:iterNum, score_func(2, :))
title("$b_2$ Convergence Plot", Interpreter="latex")
ylabel("Score Function", Interpreter="latex", FontSize=13)
subplot(3, 1, 3)
plot(1:iterNum, score_func(3, :))
title("$b_3$ Convergence Plot", Interpreter="latex")
xlabel("Iteration", Interpreter="latex", FontSize=13)