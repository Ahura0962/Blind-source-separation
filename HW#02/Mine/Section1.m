clc; clear; close all;
%% Initialization
T = 1000;
s_1 = unifrnd(-3, 3, [1, T]);
s_2 = unifrnd(-2, 2, [1, T]);

A = [1 -2; 2 -1; 3 -2];
S = [s_1; s_2];
X = A*S;

%% Part 1
figure
scatter3(X(1, :), X(2, :), X(2, :), 'square')
title('Scatter plot of $X_{3\times T}$', Interpreter='latex')
xlabel("$x_1$", Interpreter="latex")
ylabel("$x_2$", Interpreter="latex")
zlabel("$x_3$", Interpreter="latex")

R_x = X * X .';
[U, D] = eig(R_x);
[nuDiag , order] = sort(diag(D), 'descend'); 
D = diag(nuDiag);
U = U(:,order);

%% Part 2
U_2d = U(:, 1:2);
D_2d = D(1:2, 1:2);
C = U_2d.' * A;
fprintf("Energy coverage ratio for 2-dimension is %.f%%\n", sum(D_2d(:)) / sum(D(:)) * 100)

%% Part 3
B_2d = D_2d^(-1/2) * U_2d.';
Z_2d = B_2d * X;
figure
scatter(Z_2d(1, :), Z_2d(2, :), 'square')
title('Scatter plot of $Z_{2\times T}$', Interpreter='latex')
xlabel("$z_1$", Interpreter="latex")
ylabel("$z_2$", Interpreter="latex")
figure
subplot(2, 1, 1)
plot(Z_2d(1, :));
title('Plot of $z_1(t)$ (2-D)', Interpreter='latex')
xlabel("$t$", Interpreter="latex")
ylabel("$z_1(t)$", Interpreter="latex")
subplot(2, 1, 2)
plot(Z_2d(2, :));
title('Plot of $z_2(t)$ (2-D)', Interpreter='latex')
xlabel("$t$", Interpreter="latex")
ylabel("$z_2(t)$", Interpreter="latex")

%% Part 4
[Q, G, V_T] = svd(X);

%% Part 5
F = S * Z_2d.';

%% Part 6
U_1d = U(:, 1);
D_1d = D(1, 1);
fprintf("Energy coverage ratio for 1-dimension is %.f%%\n", sum(D_1d(:)) / sum(D(:)) * 100)
B_1d = D_1d^(-1/2) * U_1d.';
Z_1d = B_1d * X;

figure
plot(Z_1d(:));
title('Plot of $z(t) (1-D)$', Interpreter='latex')
xlabel("$t$", Interpreter="latex")
ylabel("$z(t)$", Interpreter="latex")