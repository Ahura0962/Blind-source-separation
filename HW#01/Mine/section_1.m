clc; clear; close all;
%% Initialization
T = 1000;
s_1 = unifrnd(-3, 3, [1, T]);
s_2 = unifrnd(-2, 2, [1, T]);

A = [0.6 0.8; 0.8 -0.6];
S = [s_1; s_2];
X = A*S;
x_1 = X(1, :);
x_2 = X(2, :);

%% Part 1
figure
% scatter(s_1, s_2, 'filled')
% xlim([-4 4])
% ylim([-4 4])
hold on
scatter(x_1, x_2, 'filled', 'd')
title('Scatter plot of X', Interpreter='latex')
xlabel('$x_1$', Interpreter='latex')
ylabel('$x_2$', Interpreter='latex')

%% Part 2
[max_x_1, ind_max_x_1] = max(x_1);
[min_x_1, ind_min_x_1] = min(x_1);
[max_x_2, ind_max_x_2] = max(x_2);
[min_x_2, ind_min_x_2] = min(x_2);

a_1_slope_pred_1 = (x_2(ind_max_x_1) - min_x_2) / (max_x_1 - x_1(ind_min_x_2));
a_1_slope_pred_2 = (max_x_2 - x_2(ind_min_x_1)) / (x_1(ind_max_x_2) - min_x_1);
a_2_slope_pred_1 = (x_2(ind_min_x_1) - min_x_2) / (min_x_1 - x_1(ind_min_x_2));
a_2_slope_pred_2 = (max_x_2 - x_2(ind_max_x_1)) / (x_1(ind_max_x_2) - max_x_1);

a_1_slope_pred = (a_1_slope_pred_1 + a_1_slope_pred_2) / 2;
a_2_slope_pred = (a_2_slope_pred_1 + a_2_slope_pred_2) / 2;

A_pred = ones(length(A));
A_pred(1, 2) = a_1_slope_pred; A_pred(2, 2) = a_2_slope_pred;

a_1_norm = sqrt(sum(A_pred(1, :).^2));
a_2_norm = sqrt(sum(A_pred(2, :).^2));

% we suppose that the norm of a1 and a2 is 1
A_pred(1, 1) = A_pred(1, 1) / a_1_norm;
A_pred(1, 2) = A_pred(1, 2) / a_1_norm;
A_pred(2, 1) = A_pred(2, 1) / a_2_norm;
A_pred(2, 2) = A_pred(2, 2) / a_2_norm;

% we have permutation ambiguity
disp('The result is:')
disp(A_pred)

%% Part 3
line_len = 100000;
dist_threshold = 0.8;

line_1_x = linspace(x_1(ind_min_x_2), max_x_1, line_len);
line_1_y = linspace(min_x_2, x_2(ind_max_x_1), line_len);
line_2_x = linspace(min_x_1, x_1(ind_max_x_2), line_len);
line_2_y = linspace(x_2(ind_min_x_1), max_x_2, line_len);
line_3_x = linspace(x_1(ind_max_x_2), max_x_1, line_len);
line_3_y = linspace(max_x_2, x_2(ind_max_x_1), line_len);
line_4_x = linspace(min_x_1, x_1(ind_min_x_2), line_len);
line_4_y = linspace(x_2(ind_min_x_1), min_x_2, line_len);

plot(line_1_x, line_1_y)
plot(line_2_x, line_2_y)
plot(line_3_x, line_3_y)
plot(line_4_x, line_4_y)

slope_line1 = slope_finder(x_1, x_2, line_len, dist_threshold, ...
                           x_1(ind_min_x_2), max_x_1, min_x_2, x_2(ind_max_x_1));
slope_line2 = slope_finder(x_1, x_2, line_len, dist_threshold, ...
                           min_x_1, x_1(ind_max_x_2), x_2(ind_min_x_1), max_x_2);
slope_line3 = slope_finder(x_1, x_2, line_len, dist_threshold, ...
                           x_1(ind_max_x_2), max_x_1, max_x_2, x_2(ind_max_x_1));
slope_line4 = slope_finder(x_1, x_2, line_len, dist_threshold, ...
                           min_x_1, x_1(ind_min_x_2), x_2(ind_min_x_1), min_x_2);

a_1_slope_pred = (slope_line1 + slope_line2) / 2;
a_2_slope_pred = (slope_line3 + slope_line4) / 2;

A_pred = ones(length(A));
A_pred(1, 2) = a_1_slope_pred; A_pred(2, 2) = a_2_slope_pred;

a_1_norm = sqrt(sum(A_pred(1, :).^2));
a_2_norm = sqrt(sum(A_pred(2, :).^2));

% we suppose that the norm of a1 and a2 is 1
A_pred(1, 1) = A_pred(1, 1) / a_1_norm;
A_pred(1, 2) = A_pred(1, 2) / a_1_norm;
A_pred(2, 1) = A_pred(2, 1) / a_2_norm;
A_pred(2, 2) = A_pred(2, 2) / a_2_norm;

% we have permutation ambiguity
disp('The result is:')
disp(A_pred)

%% Part 4 
bin = 20;
figure
histogram(x_1, bin)
title("Histogram of $x_1$", Interpreter="latex")

%% Part 5
figure
histogram(x_2, bin)
title("Histogram of $x_2$", Interpreter="latex")


%% Functions
function line_slope = slope_finder(x_1, x_2, line_len, dist_threshold, x_begin, x_end, y_begin, y_end)
    line_x = linspace(x_begin, x_end, line_len);
    line_y = linspace(y_begin, y_end, line_len);
    [dist, ind_line] = min(sqrt((line_x - x_1.').^2 + (line_y - x_2.').^2));
    ind_line = ind_line(dist <= dist_threshold);
    ind_line = unique(ind_line);
    line_slope = polyfit(x_1(ind_line), x_2(ind_line), 1);
    line_slope = line_slope(1);
end