clc; clear; close all;

%% Part 1
N_max = 35;
mu_D = zeros(1, N_max - 1);
for N = 2:N_max
    step = pi / N;
    phi_N = 0: step: pi - step;
    vecs = [cos(phi_N); sin(phi_N)];
    vecs = vecs ./ sqrt(sum(vecs.^2));
    corrr = vecs.' * vecs;
    mu_D(N - 1) = max(abs(corrr - diag(diag(corrr))), [], 'all');
end
figure
plot(2:N_max, mu_D)
title("Mutual Coherence in 2D", Interpreter="latex")
xlabel("$N$", Interpreter="latex", FontSize=13)
xticks(2:3:N_max);
ylabel("$\mu(D)$", Interpreter="latex", Rotation=0, FontSize=13)

%% Part 2
N = 10;
vecs = randn(3, N);
vecs = vecs ./ sqrt(sum(vecs.^2));
iterNum = 5e4;
list_mu_D = zeros(1, iterNum);
corrr = vecs.' * vecs;
mu_D = max(abs(corrr - diag(diag(corrr))), [], "all");
for iter = 1:iterNum
    final_vecs = vecs;
    final_mu_D = mu_D;
    for col = 1:N
        test_vecs = vecs;
        test_vecs(:, col) = randn(3, 1);
        test_vecs = test_vecs ./ sqrt(sum(test_vecs.^2));
        corrr = test_vecs.' * test_vecs;
        test_mu_D = max(abs(corrr - diag(diag(corrr))), [], "all");
        if test_mu_D < final_mu_D
            final_vecs = test_vecs;
            final_mu_D = test_mu_D;
        end
    end
    vecs = final_vecs;
    mu_D = final_mu_D;
    list_mu_D(iter) = mu_D;
end
figure
plot(1:iterNum, list_mu_D)
title("Mutual Coherence in 3D", Interpreter="latex")
xlabel("Iteration", Interpreter="latex", FontSize=13)
ylabel("$\mu(D)$", Interpreter="latex", Rotation=0, FontSize=13)