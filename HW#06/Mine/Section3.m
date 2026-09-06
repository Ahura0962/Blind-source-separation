clc; clear; close all;

%% Initialization
load hw6-part3.mat

%% Part 1
fprintf("Part 1:\n")
[mu_D, ind] = max(abs((D.' * D) - eye(size(D, 2))), [], 'all');
[row, col] = ind2sub(size(D,2), ind);
disp("mu_D = " + mu_D)
fprintf("\n\n\n")

%% Part 2:5
% Mod
fprintf("Part 2-5:\n")
[M, N] = size(D);
T = size(X, 2);
D_hat_init = randn(M, N);
tic
D_hat_MOD = D_hat_init;
D_hat_MOD = D_hat_MOD ./ sqrt(sum(D_hat_MOD.^ 2));
iterNum = 50;
RErr_MOD = zeros(1, iterNum);
for iter = 1:iterNum
    % D: fixed
    S_hat_MOD = zeros(N, T);
    for i = 1:T
        S_hat_MOD(:,i) = OMP(X(:,i), D_hat_MOD, N, N0);
    end
    % S: fixed
    D_hat_MOD = X * pinv(S_hat_MOD);
    D_hat_MOD = D_hat_MOD ./ sqrt(sum(D_hat_MOD.^2));
    RErr_MOD(iter) = trace((X - D_hat_MOD * S_hat_MOD) * (X - D_hat_MOD * S_hat_MOD).') / trace(X * X.');
end
toc
corr_MOD = D_hat_MOD.' * D;
SRR_MOD = sum(abs(corr_MOD) > 0.98, 'all') / N;
fprintf("MOD : Successful Recovery Rate: %.0f%%\n\n",SRR_MOD * 100);

% K-SVD
tic
D_hat_KSVD = D_hat_init;
D_hat_KSVD = D_hat_KSVD ./ sqrt(sum(D_hat_KSVD.^ 2));
RErr_KSVD = zeros(1, iterNum);
for iter = 1:iterNum
    % D: fixed
    S_hat_KSVD = zeros(N, T);
    for i = 1:T
        S_hat_KSVD(:,i) = OMP(X(:,i), D_hat_KSVD, N, N0);
    end
    % S: fixed
    for i = 1:N
        n = [1:(i-1) (i+1):N];
        Xr = X - D_hat_KSVD(:, n) * S_hat_KSVD(n, :);
        ind = (S_hat_KSVD(i,:) ~= 0);
        m_Xr = Xr(:,ind);
        [U, G, V] = svd(m_Xr);
        [m, n] = size(G);
        L = min(m, n);
        [~, index] = sort(diag(G(1:L, 1:L)));
        D_hat_KSVD(:,i) = U(:,index(end));
        D_hat_KSVD = D_hat_KSVD ./ sqrt(sum(D_hat_KSVD .^ 2));
        S_hat_KSVD(i, ind) = G(index(end), index(end)) * (V(:, index(end))).';
    end
    RErr_KSVD(iter) = trace((X - D_hat_KSVD * S_hat_KSVD) * (X - D_hat_KSVD * S_hat_KSVD).') / trace(X * X.');
end
toc
corr_KSVD = D_hat_KSVD.' * D;
SRR_KSVD = sum(abs(corr_KSVD) > 0.98, 'all') / N;
fprintf("K-SVD : Successful Recovery Rate: %.0f%%\n",SRR_KSVD * 100);
iterPlot = 1:iterNum;
figure
plot(iterPlot, RErr_MOD, iterPlot, RErr_KSVD)
title("MOD vs. K-SVD Algorithm", Interpreter="latex")
xlabel("Iteration", Interpreter="latex")
ylabel("Representation Error", Interpreter="latex")
legend("MOD", "K-SVD", Interpreter="latex")
fprintf("\n\n\n")

%% Part 6
fprintf("Part 6:\n")
[~ , ind_MOD] = max(abs(corr_MOD));
[~ , ind_KSVD] = max(abs(corr_KSVD));
S_sort_MOD = S_hat_MOD(ind_MOD, :);
S_sort_KSVD = S_hat_KSVD(ind_KSVD, :);
diff_MOD = abs(S_sort_MOD) - abs(S);
diff_KSVD = abs(S_sort_KSVD) - abs(S);
Err_MOD = trace(diff_MOD.' * diff_MOD) / trace(S.' * S);
Err_KSVD = trace(diff_KSVD.' * diff_KSVD) / trace(S.' * S);
disp("Error (MOD Algorithm) = " + Err_MOD)
fprintf("\n")
disp("Error (K-SVD Algorithm) = " + Err_KSVD)
fprintf("\n\n\n")