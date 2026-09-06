function s_hat = OMP(x, D, N, N0)
    x_r = x;
    ind = zeros(1, N0);
    s_hat = zeros(N, 1);
    for i = 1:N0
        innerprod = x_r.' * D;
        [~, ind(i)] = max(abs(innerprod));
        if i > 1
            D_pr = D(:, ind(1:i));
            s_hat(ind(1:i)) = pinv(D_pr) * x;
            x_r = x_r - D * s_hat;
        else
            x_r = x_r - innerprod(ind(i)) * D(:, ind(i));
        end
    end
end

