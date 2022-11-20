% add the solver path to the path
clear all; close all; clc;
%% Data Loading and Visulisation
train = importdata('train.mat');
X_train = train.X;
y_train = double(transpose(train.y));
d_train = double(transpose(train.z));
test = importdata('test.mat');
X_test = test.X;
y_test = double(transpose(test.y));
d_test = double(transpose(test.z));
clear train test;

figure();
histogram(y_train(d_train==1));
hold on;
histogram(y_train(d_train==0));
legend('Churned','Not Churned');
title('Histogram of Survival Times');

%% Linear Survival SVM
X_train = X_train(1:100,:);
y_train = y_train(1:100,:);
d_train = d_train(1:100,:);

observedX = X_train(d_train==1,:);
observedy = y_train(d_train==1);
y_comp = y_train(y_train>1);
X_comp = X_train(y_train>1,:);
X_comp_bar = [];
y_comp_bar = [];
for y=transpose(y_comp)
    [xa,ya] = get_comparable_data(y, observedX, observedy);
    X_comp_bar = [X_comp_bar; xa];
    y_comp_bar = [y_comp_bar; ya];
end
% YOUR CODE HERE:
% FORMULATE AND SOLVE DM
% FIND THE PREDICTION ON THE TRAINING SET
%prediction = ...

epsilon = size(X_comp_bar,1);
n = length(y_train);
minus = X_comp - X_comp_bar;
U = minus * minus';
W = X_train * minus';
V = X_train * X_train';
Z = diag(d_train);
r1 = 100 / epsilon;
r2 = 10000 / n;
U = U + 1e-9 * eye(size(U));
V = V + 1e-9 * eye(size(V));

alpha = sdpvar(epsilon,1);
beta = sdpvar(n,1);
gamma = sdpvar(n,1);

constraints = [alpha >= 0];
constraints = [constraints, alpha <= r1];
constraints = [constraints, beta >= 0];
constraints = [constraints, beta <= r2];
constraints = [constraints, gamma >= 0];
constraints = [constraints, gamma <= r2];
constraints = [constraints, abs(beta' * ones(length(beta), 1) - d_train' * gamma) <= 1e-12];

% Objective function
% obj = -1/2 * ([alpha; beta - Z * gamma])' * ...
%     ([U, transpose(W); W, V])' * [alpha; beta - Z * gamma] + ...
%     ([y_comp - y_comp_bar; y_train])' * [alpha; beta - Z * gamma];

% obj = -0.5 * alpha' * U * alpha - (beta - Z * gamma)' * W * alpha - ...
%     0.5 * (beta - Z * gamma)' * V * (beta - Z * gamma) + ...
%     (y_comp - y_comp_bar)' * alpha + y_train' * (beta - Z * gamma);

Q = [U, W'; W, V];
R = chol(Q);
x = [alpha; beta - Z * gamma];
h = sdpvar(length(x),1);
constraints = [constraints, h == R * x];
obj = -1/2 * h' * h + ...
    [y_comp - y_comp_bar; y_train]' * x;
obj = - obj;
% Specify solver settings, run solver, and retrieve optimal solution
ops = sdpsettings('solver', 'gurobi', 'verbose', 1);    % mosek	gurobi sdpt3
diagnosis = optimize(constraints, obj, ops);

w = alpha' * (X_comp - X_comp_bar) + (beta - gamma .* d_train)' * X_train;
b = y_train - X_train * w';


prediction = X_train * w' + 50;

%% Display Training Results
plot_results(y_train, prediction, d_train, 'Linear Survival (Training)')
disp('Linear Survival')
disp('Train Results')
compute_metrics(y_train, prediction, d_train)

%% Display Testing Results
% YOUR CODE HERE:
%prediction = ...
%FIND THE PREDICTION ON THE TEST SET





plot_results(y_test, prediction, d_test, 'Linear Survival (Testing)');
disp('Test Results')
compute_metrics(y_test, prediction, d_test)

%% Simple MAE Regression
% YOUR CODE HERE:
% SOLVE THE MAE REGRESSION PROBLEM
% FIND THE PREDICTION ON THE TRAINING SET
%prediction_reg = ...
%


%% Display Training Results
disp('Linear MAE')
disp('Train Results')
plot_results(y_train, prediction_reg, d_train, 'Linear MAE (Training)')
compute_metrics(y_train, prediction_reg, d_train)

%% Display Testing Results
disp('Test Results')
% YOUR CODE HERE
%prediction_reg = ...
% FIND THE PREDICTION ON THE TRAINING SET



plot_results(y_test, prediction_reg, d_test, 'Linear MAE (Testing)');
compute_metrics(y_test, prediction_reg, d_test)

%% Kernelized Version
% YOUR CODE HERE:
% FORMULATE AND SOLVE DM
% FIND THE PREDICTION ON THE TRAINING SET
%prediction = ...
%
    







%% Display Training Results
disp('Kernel Survival')
disp('Train Results')
plot_results(y_train, prediction, d_train, 'Kernel Survival (Training)')
compute_metrics(y_train, prediction, d_train)

%% Display Testing Results
disp('Test Results')
% YOUR CODE HERE:
%prediction = ...
%FIND THE PREDICTION ON THE TEST SET



plot_results(y_test, prediction, d_test, 'Kernel Survival (Testing)');
compute_metrics(y_test, prediction, d_test)