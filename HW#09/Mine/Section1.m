clc; clear; close all;

%% Initialization
load hw9.mat
warning("off")

%% Part 0
X_noNoise = A*S;
X = A*S + Noise;
figure
subplot(3, 1, 1)
plot(S(1,:), "m")
title("S(1, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$s_1(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(S(2,:), "m")
title("S(2, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$s_2(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(S(3,:), "m")
title("S(3, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$s_3(t)$", Interpreter="latex", Rotation=0, FontSize=13)
figure
subplot(3, 1, 1)
plot(X_noNoise(1,:), "m")
title("X(noNoise)(1, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_1(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(X_noNoise(2,:), "m")
title("X(noNoise)(2, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_2(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(X_noNoise(3,:), "m")
title("X(noNoise)(3, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_3(t)$", Interpreter="latex", Rotation=0, FontSize=13)
figure
subplot(3, 1, 1)
plot(X(1,:), "m")
title("X(Noisy)(1, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_1(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 2)
plot(X(2,:), "m")
title("X(Noisy)(2, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_2(t)$", Interpreter="latex", Rotation=0, FontSize=13)
subplot(3, 1, 3)
plot(X(3,:), "m")
title("X(Noisy)(3, :)", Interpreter="latex")
xlabel("$t$", Interpreter="latex", FontSize=13)
ylabel("$x_3(t)$", Interpreter="latex", Rotation=0, FontSize=13)

%% Part 1
R_x = X * X .';
[U, D] = eig(R_x);
[nuDiag , order] = sort(diag(D), 'descend'); 
D = diag(nuDiag);
U = U(:,order);
Z = D^(-1/2) * U.' * X;
mu = .2;
B = randn(size(A));
B = orth(B);
[M, T] = size(X);
K_Y = @(Y) [ones(1,T); Y; Y.^2; Y.^3; Y.^4; Y.^5];
dK_Y = @(Y) [zeros(1,T); ones(1,T); 2*Y; 3*Y.^2; 4*Y.^3; 5*Y.^4];
iterNum = 5e2;
score_func = zeros(M, iterNum);
for row = 1:size(B, 1)
    b = B(row, :).';
    for iter = 1:iterNum
        Y = B * Z;
        y = Y(1, :);
        K_y = K_Y(y);
        Theta_hat = (K_y * K_y.' ./ T) \ mean(dK_Y(y), 2);
        Psi_hat = Theta_hat.' * K_y;
        df_B = Psi_hat * Z.' ./ T;
        b = b - mu * df_B.';
        b = (eye(size(B)) - B(1:row-1, :).' * B(1:row-1, :)) * b;
        b = b ./ norm(b, "fro");
        score_func(row, iter) = norm(df_B, "fro");
        B(row, :) = b.';
    end
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
score_func = sqrt(score_func(1, :).^2 + score_func(2, :).^2 + score_func(3, :).^2);
E = norm(S_hat - S, "fro")^2 / norm(S, "fro")^2;
permutation_matrix = B * (D^(-1/2) * U.') * A;

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