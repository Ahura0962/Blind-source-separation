clc; clear; close all;
%% Part 1
[x_1, x_2] = meshgrid(-10:0.02:10);
x = [x_1; x_2];
f = x_1.^2 + x_2.^2 - 4.*x_1 - 6.*x_2 + 13 + x_1.*x_2;
figure
subplot(1, 2, 1)
mesh(x_1, x_2, f)
title('Mesh plot in linear scale', Interpreter='latex')
subplot(1, 2, 2)
mesh(x_1, x_2, 10*log10(f))
title('Mesh plot in log scale', Interpreter='latex')

%% Part 2
figure
subplot(1, 2, 1)
contour(x_1, x_2, f, 'ShowText', 'on')
title('Contour plot in linear scale', Interpreter='latex')
subplot(1, 2, 2)
contour(x_1, x_2, 10*log10(f), 'ShowText', 'on')
title('Contour plot in log scale', Interpreter='latex')

%% Part 5
x_begin = [6; 6];
err = 1e-6;
x_sd_1 = steepest_descend(0.1, x_begin, err);
x_sd_2 = steepest_descend(0.01, x_begin, err);
figure
plot(0:length(x_sd_1) - 1, x_sd_1)
hold on
plot(0:length(x_sd_2) - 1, x_sd_2)
title("Steepest Descend algorithm", Interpreter="latex")
legend('$\mu = 0.1$', '$\mu = 0.01$', Interpreter="latex")

%% Part 6
x_nw = newton(x_begin, err);

%% Part 7
[x_am, x_x] = alternation_minimization(x_begin, err);
figure
contour(x_1, x_2, 10*log10(f), 'ShowText', 'on')
title('Contour plot in log scale', Interpreter='latex')
hold on
plot(x_x(1,:), x_x(2,:))

%% Part 8
x_sd_1 = gradient_projection(0.1, x_begin, err);

%% Functions
    % STEEPEST DESCEND
    function out = steepest_descend(mu, x, convergeLimit)

    %%% these lines should change based on the optimizing function
    grad = [2*x(1) + x(2) - 4; x(1) + 2*x(2) - 6];
    x_k = x; x_k_new = x_k - mu*grad;
    f = x_k(1).^2 + x_k(2).^2 - 4.*x_k(1) - 6.*x_k(2) + 13 + x_k(1).*x_k(2);
    f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
    %%% these lines should change based on the optimizing function

    out = [f f_new];
    while(abs(f_new - f) >= convergeLimit)
        grad = [2*x_k_new(1) + x_k_new(2) - 4; x_k_new(1) + 2*x_k_new(2) - 6];
        x_k = x_k_new;
        x_k_new = x_k - mu*grad;
        f = f_new;
        f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
        out = cat(2, out, f_new);
    end
    fprintf("S.D: final value of x for mu = %.2f is [%f, %f] ", mu, x_k_new)
    fprintf("in iteration: %i\n", length(out) - 1)
    end



    % NEWTON
    function out = newton(x, convergeLimit)

    %%% these lines should change based on the optimizing function
    grad = [2*x(1) + x(2) - 4; x(1) + 2*x(2) - 6];
    hessian = [2 1; 1 2];
    x_k = x; x_k_new = x_k - hessian\grad;
    f = x_k(1).^2 + x_k(2).^2 - 4.*x_k(1) - 6.*x_k(2) + 13 + x_k(1).*x_k(2);
    f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
    %%% these lines should change based on the optimizing function

    out = f;
    while(abs(f_new - f) >= convergeLimit)
        grad = [2*x_k_new(1) + x_k_new(2) - 4; x_k_new(1) + 2*x_k_new(2) - 6];
        x_k = x_k_new;
        x_k_new = x_k - hessian\grad;
        f = f_new;
        f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
        out = cat(2, out, f_new);
    end
    fprintf("\nNewton: final value of x for is [%f, %f] ", x_k_new)
    fprintf("in iteration: %i\n", length(out) - 1)
    end



    % ALTERNATION MINIMIZATION
    function [out, x_x] = alternation_minimization(x, convergeLimit)
    %%% these lines should change based on the optimizing function
    x_k = x; x_k_new = x;
    x_k_new(1) = 2 - x(2)/2;
    f = x_k(1).^2 + x_k(2).^2 - 4.*x_k(1) - 6.*x_k(2) + 13 + x_k(1).*x_k(2);
    f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
    %%% these lines should change based on the optimizing function

    out = [f f_new];
    x_x = [x_k x_k_new];
    while(abs(f_new - f) >= convergeLimit)
        if(mod(length(out), 2) == 0)
            x_k_new(2) = 3 - x_k_new(1)/2;
        else
            x_k_new(1) = 2 - x_k_new(2)/2;
        end
        f = f_new;
        f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
        out = cat(2, out, f_new);
        x_x = cat(2, x_x, x_k_new);
    end
    fprintf("\nAlternation Minimization: final value of x for is [%f, %f] ", x_k_new)
    fprintf("in iteration: %i\n", length(out) - 1)
    end



    % GRADIENT PROJECTION
    function out = gradient_projection(mu, x, convergeLimit)

    %%% these lines should change based on the optimizing function
    grad = [2*x(1) + x(2) - 4; x(1) + 2*x(2) - 6];
    x_k = x; x_k_new = x_k - mu*grad;
    x_k = x_k ./ sqrt(sum(x_k.^2)) ; x_k_new = x_k_new ./ sqrt(sum(x_k_new.^2));
    f = x_k(1).^2 + x_k(2).^2 - 4.*x_k(1) - 6.*x_k(2) + 13 + x_k(1).*x_k(2);
    f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
    %%% these lines should change based on the optimizing function

    out = [f f_new];
    while(abs(f_new - f) >= convergeLimit)
        grad = [2*x_k_new(1) + x_k_new(2) - 4; x_k_new(1) + 2*x_k_new(2) - 6];
        x_k = x_k_new;
        x_k_new = x_k - mu*grad;
        x_k_new = x_k_new ./ sqrt(sum(x_k_new.^2));
        f = f_new;
        f_new = x_k_new(1).^2 + x_k_new(2).^2 - 4.*x_k_new(1) - 6.*x_k_new(2) + 13 + x_k_new(1).*x_k_new(2);
        out = cat(2, out, f_new);
    end
    fprintf("\nG.P: final value of x for mu = %.2f is [%f, %f] ", mu, x_k_new)
    fprintf("in iteration: %i\n", length(out) - 1)
    end