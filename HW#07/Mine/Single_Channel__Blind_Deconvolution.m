function [s_hat, comb, x_hat] = Single_Channel__Blind_Deconvolution(x, L, K, T, iterNum)
    x1_resahpe = reshape(x, K, []);
    Y = x1_resahpe(:, 1:L);
    alpha = ones(K, 1);
    corr_max = zeros(1, K);
    ind_corr_max = zeros(1, K);
    for iter = 1:iterNum
        % comb function: fixed
        s_hat = (Y.' * alpha) / (alpha.' * alpha);
        s_hat = s_hat ./ sqrt(sum(s_hat .^ 2));
        % s: fixed
        Z = buffer(x, L, L-1, 'nodelay');
        b_T = pinv(s_hat) * Z;
        b_T(b_T < 0) = 0;
        for i = 1:K
            [corr_max(i), ind_corr_max(i)] = max(b_T);
            min_L = max(ind_corr_max(i) - L, 1);
            max_L = min(ind_corr_max(i) + L, T);
            b_T(min_L: max_L) = 0;
        end
        [tau_k, sort_ind] = sort(ind_corr_max, 'ascend');
        alpha = corr_max(sort_ind);
        alpha = reshape(alpha, K, 1);
        for i = 1:K
            Y(i, :) = x(tau_k(i): tau_k(i) + L - 1);
        end
    end
    comb = zeros(1, T-L+1);
    comb(ind_corr_max) = corr_max;
    x_hat = conv(s_hat, comb);
end