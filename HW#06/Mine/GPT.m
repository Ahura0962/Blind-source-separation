clc; clear; close all;

%% Part 1
N_max = 35;
mu_D = zeros(1, N_max - 1);
for N = 2:N_max
    if mod(N, 2) ~= 0
        step = 2 * pi / N;
        phi_N = 0: step: 2*pi - step;
    else
        step = pi / N;
        phi_N = 0: step: pi - step;
    end
    vecs = [cos(phi_N); sin(phi_N)];
    corrr = vecs.' * vecs;
    mu_D(N - 1) = max(abs(corrr - diag(diag(corrr))), [], 'all');
end
figure
plot(2:N_max, mu_D)
title("Mutual Coherence in 2D", Interpreter="latex")
xlabel("$N$ (Iteration)", Interpreter="latex")
xticks(2:3:40);
ylabel("$\mu(D)$", Interpreter="latex", Rotation=0)

%%

function [vecs, minMu] = optimizeVectors(n, iterNum, candidates)
    % Initialize random vectors
    vecs = randn(3, n);
    vecs = vecs ./ sqrt(sum(vecs.^2));
    
    % Array to track the minimum mutual coherence
    list_mu_D = zeros(1, iterNum);
    
    % Initial mutual coherence
    corrr = vecs.' * vecs;
    mu_D = max(max(abs(corrr - diag(diag(corrr)))));
    
    for iter = 1:iterNum
        best_vecs = vecs;
        best_mu_D = mu_D;
        
        for col = 1:n
            % Generate multiple candidates for the current column
            new_vecs = randn(3, candidates);
            new_vecs = new_vecs ./ sqrt(sum(new_vecs.^2));
            

            test_vecs = vecs;
            test_vecs(:, col) = new_vecs(:, 1);
            test_vecs = test_vecs ./ sqrt(sum(test_vecs.^2));

            % Calculate the new mutual coherence
            corrr = test_vecs.' * test_vecs;
            test_mu_D = max(abs(corrr - diag(diag(corrr))), [], "all");

            % Update if the new mutual coherence is smaller
            if test_mu_D < best_mu_D
                best_vecs = test_vecs;
                best_mu_D = test_mu_D;
            end
        end
        
        vecs = best_vecs;
        mu_D = best_mu_D;
        list_mu_D(iter) = mu_D;
    end
    
    % Final minimum mutual coherence
    minMu = mu_D;
    
    % Plotting the convergence
    figure;
    plot(1:iterNum, list_mu_D);
    xlabel('Iteration');
    ylabel('Mutual Coherence');
    title('Convergence of Mutual Coherence');
    grid on;
end

% Example usage
n = 10;           % Number of vectors
iterNum = 1e5;    % Number of iterations
candidates = 1;  % Number of candidates per vector replacement
[optimizedVecs, minMu] = optimizeVectors(n, iterNum, candidates);
disp('Final Vectors:');
disp(optimizedVecs);
disp('Minimum Mutual Coherence:');
disp(minMu);
