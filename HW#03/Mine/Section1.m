clc; clear; close all;

%% Initialization
load hw3-1.mat
TrainData_class1 = TrainData_class1 - mean(TrainData_class1, [2 3]);
TrainData_class2 = TrainData_class2 - mean(TrainData_class2, [2 3]);
TestData = TestData - mean(TestData, [2 3]);

%% Part 1
trainData_size = size(TrainData_class1);
% Rx_1 = 0; Rx_2 = 0;
% for i = 1:trainData_size(3)
%     Rx_1 = Rx_1 + (TrainData_class1(:, :, i) * TrainData_class1(:, :, i).') ./ trainData_size(3);
%     Rx_2 = Rx_2 + (TrainData_class2(:, :, i) * TrainData_class2(:, :, i).') ./ trainData_size(3);
% end
Rx_1 = mean(pagemtimes(TrainData_class1, pagetranspose(TrainData_class1)), 3);
Rx_2 = mean(pagemtimes(TrainData_class2, pagetranspose(TrainData_class2)), 3);
[U, D] = eig(Rx_1, Rx_2);
[nuDiag , order] = sort(diag(D), 'descend'); 
D_csp = diag(nuDiag);
U = U(:,order);
W_csp = U ./ sqrt(sum(U.^2));
W1_exp49_r1 = W_csp(:, 1).' * TrainData_class1(:, :, 49); 
W2_exp49_r1 = W_csp(:, 1).' * TrainData_class2(:, :, 49); 
W1_exp49_r30 = W_csp(:, 30).' * TrainData_class1(:, :, 49); 
W2_exp49_r30 = W_csp(:, 30).' * TrainData_class2(:, :, 49);
figure
subplot(2, 1, 1)
plot(1:trainData_size(2) ,W1_exp49_r1, 'k', 1:trainData_size(2), W2_exp49_r1, 'm')
title("Filtered Signal for $W_{csp}(:, 1)$", Interpreter="latex")
xlabel("Samples", Interpreter="latex")
legend("$Var(Ch1)$", "$Var(Ch2)$", Interpreter="latex", Box="off")
xlim([1 trainData_size(2)])
subplot(2, 1, 2)
plot(1:trainData_size(2) ,W1_exp49_r30, 'k', 1:trainData_size(2), W2_exp49_r30, 'm')
title("Filtered Signal for $W_{csp}(:, 30)$", Interpreter="latex")
xlabel("Samples", Interpreter="latex")
legend("$Var(Ch1)$", "$Var(Ch2)$", Interpreter="latex", Box="off")
xlim([1 trainData_size(2)])
fprintf("Variances of Filtered Signal of 1st Filter is:\n\t%f for 1st Class\n\t%f for 2nd Class\n\n", ...
    var(W1_exp49_r1), var(W2_exp49_r1))
fprintf("Variances of Filtered Signal of 30th Filter is:\n\t%f for 1st Class\n\t%f for 2nd Class\n\n", ...
    var(W1_exp49_r30), var(W2_exp49_r30))

%% Part 2
figure
subplot(2, 1, 1)
stem(abs(W_csp(:, 1)), "m")
title("Abs of 1st Filter", Interpreter="latex")
xlabel("Channels", Interpreter="latex")
subplot(2, 1, 2)
stem(abs(W_csp(:, 30)), "k")
title("Abs of 30th Filter", Interpreter="latex")
xlabel("Channels", Interpreter="latex")

%% Part 3
important_filters = [W_csp(:, 1:7) W_csp(:, 24:30)];
num_important_filters = size(important_filters, 2);
% x_1_60_exps = zeros(num_important_filters, trainData_size(3));
% x_2_60_exps = zeros(num_important_filters, trainData_size(3));
% for i = 1:trainData_size(3)
%     X_1 = TrainData_class1(:, :, i);
%     X_2 = TrainData_class2(:, :, i);
%     W1_imp = important_filters.' * X_1;
%     W2_imp = important_filters.' * X_2;
%     x_1_60_exps(:, i) = var(W1_imp, 0, 2);
%     x_2_60_exps(:, i) = var(W2_imp, 0, 2);
% end
W1_imp = pagemtimes(important_filters.',TrainData_class1);
W2_imp = pagemtimes(important_filters.',TrainData_class2);
x_1_60_exps = var(W1_imp, 0, 2);
x_1_60_exps = reshape(x_1_60_exps, [num_important_filters, trainData_size(3)]);
x_2_60_exps = var(W2_imp, 0, 2);
x_2_60_exps = reshape(x_2_60_exps, [num_important_filters, trainData_size(3)]);
mu_1 = mean(x_1_60_exps, 2);
mu_2 = mean(x_2_60_exps, 2);
sigma_1 = ((x_1_60_exps - mu_1) * (x_1_60_exps - mu_1).') ./ trainData_size(3);
sigma_2 = ((x_2_60_exps - mu_2) * (x_2_60_exps - mu_2).') ./ trainData_size(3);
one = (mu_1 - mu_2) * (mu_1 - mu_2).';
two = sigma_1 + sigma_2;
[U, D] = eig(one, two);
[nuDiag , order] = sort(diag(D), 'descend'); 
D_lda = diag(nuDiag);
U = U(:,order);
W_lda = U ./ sqrt(sum(U.^2));
W_lda = W_lda(:, 1);
Mu1 = W_lda.' * mu_1;
Mu2 = W_lda.' * mu_2;
c = (Mu1 + Mu2) / 2;

%% Part 4
testData_size = size(TestData);
% var_num = zeros(num_important_filters, testData_size(3));
% for i = 1:testData_size(3)
%     var_signal = important_filters .' * TestData(:, :, i);
%     var_num(:, i) = var(var_signal, 0, 2);
% end
var_signal = pagemtimes(important_filters.', TestData);
var_num = var(var_signal, 0, 2);
var_num = reshape(var_num, [num_important_filters, testData_size(3)]);
MyLabel = W_lda.' * var_num;
MyLabel = 1 + (MyLabel < c);

%% Part 5
figure
stem(TestLabel, LineStyle="none", Marker="square", MarkerSize=8, Color="k");
hold on
stem(MyLabel, Linestyle="none", Marker="x", MarkerSize=8, Color="m")
title("My Label vs. True Label", Interpreter="latex")
legend("True Label", "My Label", Interpreter="latex")
xlim([0 41])
ylim([0 3])
true_num = sum(TestLabel == MyLabel);
fprintf("Prediction Rate is: %i%%\n\n", true_num / testData_size(3) * 100)