clc; clear; close all;

%% Initaialization
load hw3-2.mat
fs = 250;
t = 0: 1/fs: 5 - 1/fs;

%% Part 1
[Xt1, Xt2, Xt3, Xt4, Xt5] = deal([]);
for i = 1:6 
    Xt1 = cat(1, Xt1, [sin(i*2*pi*freq(1)*t); cos(i*2*pi*freq(1)*t)]);
    Xt2 = cat(1, Xt2, [sin(i*2*pi*freq(2)*t); cos(i*2*pi*freq(2)*t)]);
    Xt3 = cat(1, Xt3, [sin(i*2*pi*freq(3)*t); cos(i*2*pi*freq(3)*t)]);
    Xt4 = cat(1, Xt4, [sin(i*2*pi*freq(4)*t); cos(i*2*pi*freq(4)*t)]);
    Xt5 = cat(1, Xt5, [sin(i*2*pi*freq(5)*t); cos(i*2*pi*freq(5)*t)]);
end
Xt = cat(3, Xt1, Xt2, Xt3, Xt4, Xt5);
rho = zeros(1, length(freq));
mylabel = zeros(1, size(data, 3));
R_X = pagemtimes(Xt, pagetranspose(Xt));
R_Y = pagemtimes(data, pagetranspose(data));
for i = 1:size(data, 3)
    R_XY = pagemtimes(Xt, pagetranspose(data(:, :, i)));
    R_YX = pagemtimes(data(:, :, i), pagetranspose(Xt));
    for j = 1:length(freq)
        R_x = R_X(:, :, j);
        R_y = R_Y(:, :, i);
        R_xy = R_XY(:, :, j);
        R_yx = R_YX(:, :, j);
        Sigma1 = R_x^(-1/2) * R_xy * R_y^(-1/2) * R_y^(-1/2) * R_yx * R_x^(-1/2);
        Sigma2 = R_y^(-1/2) * R_yx * R_x^(-1/2) * R_x^(-1/2) * R_xy * R_y^(-1/2);
        [c, D1] = eig(Sigma1, eye(12));
        [d, D2] = eig(Sigma2, eye(6));
        a = R_x^(-1/2) * c(:, 1);
        b = R_y^(-1/2) * d(:, 1);
        rho(j) = (a.' * Xt(:, :, j)) * (b.' * data(:, : ,i)).';
    end
    [~, ind] = max(abs(rho));
    mylabel(i) = freq(ind);
end
figure
stem(label, LineStyle="none", Marker="square", MarkerSize=8, Color="k");
hold on
stem(mylabel, Linestyle="none", Marker="x", MarkerSize=8, Color="m")
title("My Label vs. True Label", Interpreter="latex")
legend("True Label", "My Label", Interpreter="latex", Location = "southeast")
xlim([0 16])
ylim([6 12])
true_num = sum(label == mylabel);
fprintf("Prediction Rate is: %i%%\n\n", true_num / size(data, 3) * 100)