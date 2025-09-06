clc;
clear;
close all;
addpath(genpath('PR/'))
%-------------径向基代理模型-------------%
theta = [0.01 0.01];
lob = [1e-2 1e-2];
upb = [0.1 0.1];
variables_lob = [3 28];
variables_upb = [4 33];
range = variables_upb - variables_lob;
Decision_variables = importdata('x_1.txt');
testX = importdata('testX.txt');
testY1 = importdata('testY1.txt');
testY2 = importdata('testY2.txt');
Decision_variables =...
    (Decision_variables - variables_lob) ./ repmat(range, size(Decision_variables, 1), 1);

testX =(testX - variables_lob) ./ repmat(range, size(testX, 1), 1);

Objections_1 = importdata('y_1.txt');
Objections_2 = importdata('y_2.txt');
RBF_y1 = newrb(Decision_variables', Objections_1',0.1, 10);
RBF_y2 = newrb(Decision_variables', Objections_2',0.1, 10);
dmodel = [RBF_y1,RBF_y2];
% Plot models
X = gridsamp([0 0;1 1], 50);
X1 = reshape(X(: ,1), 50, 50);
X2 = reshape(X(: ,2), 50, 50);
YX = sim(RBF_y1,X');
YX = reshape(YX, size(X1));
figure(1)
mesh(X1, X2, YX)
hold on
plot3(Decision_variables(:, 1), Decision_variables(:, 2),Objections_1, '.k', 'markersize', 10)
hold on
plot3(testX(:, 1), testX(:, 2),testY1, '.r', 'markersize', 10)
hold off

X = gridsamp([0 0;1 1], 50);
X1 = reshape(X(: ,1), 50, 50);
X2 = reshape(X(: ,2), 50, 50);
% YX = sim(dmodel(2),X');
YX = sim(RBF_y2,X');
YX = reshape(YX, size(X1));
figure(2)
mesh(X1, X2, YX)
hold on
plot3(Decision_variables(:, 1), Decision_variables(:, 2), Objections_2, '.k', 'markersize', 10)
hold on
plot3(testX(:, 1), testX(:, 2),testY2, '.r', 'markersize', 10)


%--------------------误差验证-------------------%
y_true1 = testY1;
y_true2 = testY2;
y_pred1 = sim(RBF_y1,testX');
y_pred2 = sim(RBF_y2,testX');
[R2_1, RMSE_1, MARE_1] = evaluate_metrics(y_true1, y_pred1');
result_1 = [R2_1, RMSE_1, MARE_1];
[R2_2, RMSE_2, MARE_2] = evaluate_metrics(y_true2, y_pred2');
result_2 = [R2_2, RMSE_2, MARE_2];

