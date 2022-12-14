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

% %% Linear Survival SVM
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
% Set the parameters;
epsilon = size(X_comp_bar,1);
n = length(y_train);
r1 = 100 / epsilon;
r2 = 10000 / n;

minus = X_comp - X_comp_bar;
U = minus * minus';
W = X_train * minus';
V = X_train * X_train';
Z = diag(d_train);

% Guarantee positive definite
U = U + 1e-9 * eye(size(U));
V = V + 1e-9 * eye(size(V));

% Define decision variables
alpha = sdpvar(epsilon,1);
beta = sdpvar(n,1);
gamma = sdpvar(n,1);

% Add constraints
constraints = [alpha >= 0];
constraints = [constraints, alpha <= r1];
constraints = [constraints, beta >= 0];
constraints = [constraints, beta <= r2];
constraints = [constraints, gamma >= 0];
constraints = [constraints, gamma <= r2];
constraints = [constraints, abs(beta' * ones(length(beta), 1) - d_train' * gamma) <= 1e-12];

% Objective function
Q = [U, W'; W, V];
R = chol(Q);
h = sdpvar(length([alpha; beta - Z * gamma]),1);
constraints = [constraints, h == R * [alpha; beta - Z * gamma]];
obj = -1/2 * h' * h + ...
    [y_comp - y_comp_bar; y_train]' * [alpha; beta - Z * gamma];
obj = - obj;

% Specify solver settings, run solver, and retrieve optimal solution
ops = sdpsettings('solver', 'mosek', 'verbose', 0);    % mosek	gurobi sdpt3

% Run solver
diagnosis = optimize(constraints, obj, ops);

if diagnosis.problem == 0
 % Extract value
    alpha = value(alpha);
    beta = value(beta);
    gamma = value(gamma);
else
    display('Hmm, something went wrong!');
    diagnosis.info
    yalmiperror(diagnosis.problem)
end

% Calculate w and b based on KKT
w = alpha' * (X_comp - X_comp_bar) + (beta - gamma .* d_train)' * X_train;
b = [];
for i = 1:n
    if i == 60
        a=2;
    end
    if beta(i) > 0 + 1e-12 || ...
            (gamma(i) > 0 + 1e-12 && d_train(i) > 0)
        b = [b, y_train(i) - X_train(i,:) * w'];
    end
end
b_mean = mean(b);

% Fine the prediction on the training set
prediction = X_train * w' + b_mean;

%% Display Training Results
plot_results(y_train, prediction, d_train, 'Linear Survival (Training)')
disp('Linear Survival')
disp('Train Results')
compute_metrics(y_train, prediction, d_train)

%% Display Testing Results
% YOUR CODE HERE:
% Fine the prediction on the test set
prediction = X_test * w' + b_mean;

% Display test Results
plot_results(y_test, prediction, d_test, 'Linear Survival (Testing)');
disp('Test Results')
compute_metrics(y_test, prediction, d_test)

%% Simple MAE Regression
% YOUR CODE HERE:
% SOLVE THE MAE REGRESSION PROBLEM

% Define decision variables
w = sdpvar(15,1);
b = sdpvar(1,1);

% Objective function
obj = 0.5 * w' * w + r2 * sum(abs(y_train - X_train * w - b));

% Specify solver settings, run solver, and retrieve optimal solution
ops = sdpsettings('solver', 'mosek', 'verbose', 0);    % mosek	gurobi sdpt3

% Run solver
diagnosis = optimize([], obj, ops);

if diagnosis.problem == 0
 % Extract value
    w = value(w);
    b = value(b);
else
    display('Hmm, something went wrong!');
    diagnosis.info
    yalmiperror(diagnosis.problem)
end
b_mean = mean(b);

% Fine the prediction on the training set
prediction_reg = X_train * w + b_mean;

%% Display Training Results
disp('Linear MAE')
disp('Train Results')
plot_results(y_train, prediction_reg, d_train, 'Linear MAE (Training)')
compute_metrics(y_train, prediction_reg, d_train)

%% Display Testing Results
disp('Test Results')
% YOUR CODE HERE

% Fine the prediction on the test set
prediction_reg = X_test * w + b_mean;

% Display Testing Results
plot_results(y_test, prediction_reg, d_test, 'Linear MAE (Testing)');
compute_metrics(y_test, prediction_reg, d_test)

%% Kernelized Version
% YOUR CODE HERE:
% FORMULATE AND SOLVE DM

% Define kernel function    
kernel = @(x1,x2) exp(-0.5 * (norm(x1 - x2))^2 ./ 1.5);

% Set parameters
epsilon = size(X_comp_bar,1);
n = length(y_train);
r1 = 100 / epsilon;
r2 = 10000 / n;

% Calculate related inner products
U_1 = [];
for i = 1:epsilon
    for j = 1:epsilon
        U_1(i,j) = kernel(X_comp(i,:), X_comp(j,:));
    end
end
U_2 = [];
for i = 1:epsilon
    for j = 1:epsilon
        U_2(i,j) = kernel(X_comp(i,:), X_comp_bar(j,:));
    end
end
U_3 = [];
for i = 1:epsilon
    for j = 1:epsilon
        U_3(i,j) = kernel(X_comp_bar(i,:), X_comp_bar(j,:));
    end
end
U = U_1 - U_2 - U_2' + U_3;

W_1 = [];
for i = 1:n
    for j = 1:epsilon
        W_1(i,j) = kernel(X_train(i,:), X_comp(j,:));
    end
end
W_2 = [];
for i = 1:n
    for j = 1:epsilon
        W_2(i,j) = kernel(X_train(i,:), X_comp_bar(j,:));
    end
end
W = W_1 - W_2;

V = [];
for i = 1:n
    for j = 1:n
        V(i,j) = kernel(X_train(i,:), X_train(j,:));
    end
end

Z = diag(d_train);

% Guarantee positive definite
U = U + 1e-9 * eye(size(U));
V = V + 1e-9 * eye(size(V));

% Define decision variables
alpha = sdpvar(epsilon,1);
beta = sdpvar(n,1);
gamma = sdpvar(n,1);

% Add constraints
constraints = [alpha >= 0];
constraints = [constraints, alpha <= r1];
constraints = [constraints, beta >= 0];
constraints = [constraints, beta <= r2];
constraints = [constraints, gamma >= 0];
constraints = [constraints, gamma <= r2];
constraints = [constraints, abs(beta' * ones(length(beta), 1) - d_train' * gamma) <= 1e-12];

% Define objective function
Q = [U, W'; W, V];
R = chol(Q);
h = sdpvar(length([alpha; beta - Z * gamma]),1);
constraints = [constraints, h == R * [alpha; beta - Z * gamma]];
obj = -1/2 * h' * h + ...
    [y_comp - y_comp_bar; y_train]' * [alpha; beta - Z * gamma];
obj = - obj;

% Specify solver settings, run solver, and retrieve optimal solution
ops = sdpsettings('solver', 'mosek', 'verbose', 0);    % mosek	gurobi sdpt3

% Run solver
diagnosis = optimize(constraints, obj, ops);

if diagnosis.problem == 0
 % Extract value
    alpha = value(alpha);
    beta = value(beta);
    gamma = value(gamma);
else
    display('Hmm, something went wrong!');
    diagnosis.info
    yalmiperror(diagnosis.problem)
end

% Calculate b
b = [];
k_X_comp_X_train = [];
k_X_comp_bar_X_train = [];
k_X_train_X_train = [];
for k = 1:n
    for i = 1:epsilon
            k_X_comp_X_train(i,:) = kernel(X_comp(i,:), X_train(k,:));
    end
    for i = 1:epsilon
            k_X_comp_bar_X_train(i,:) = kernel(X_comp_bar(i,:), X_train(k,:));
    end
    for i = 1:n
            k_X_train_X_train(i,:) = kernel(X_train(i,:), X_train(k,:));
    end
    if beta(k) > 0 + 1e-12 || ...
            (gamma(k) > 0 + 1e-12 && d_train(k) > 0)
        b = [b, y_train(k) - alpha' * (k_X_comp_X_train - k_X_comp_bar_X_train) - (beta - Z * gamma)' * k_X_train_X_train];
    end
end
b_mean = mean(b);

% Fine the prediction on the training set
X_pre_1 = [];
for i = 1:epsilon
    for j = 1:n
        X_pre_1(i,j) = kernel(X_comp(i,:), X_train(j,:)) - kernel(X_comp_bar(i,:), X_train(j,:));
    end
end
X_pre_2 = [];
for i = 1:n
    for j = 1:n
        X_pre_2(i,j) = kernel(X_train(i,:), X_train(j,:));
    end
end
prediction = alpha' * X_pre_1 + (beta - Z * gamma)' * X_pre_2 + b_mean;
prediction = prediction';

%% Display Training Results
disp('Kernel Survival')
disp('Train Results')
plot_results(y_train, prediction, d_train, 'Kernel Survival (Training)')
compute_metrics(y_train, prediction, d_train)

%% Display Testing Results
disp('Test Results')
% YOUR CODE HERE:

% Fine the prediction on the test set
for i = 1:epsilon
    for j = 1:length(X_test)
        X_pre_1(i,j) = kernel(X_comp(i,:), X_test(j,:)) - kernel(X_comp_bar(i,:), X_test(j,:));
    end
end
for i = 1:n
    for j = 1:length(X_test)
        X_pre_2(i,j) = kernel(X_train(i,:), X_test(j,:));
    end
end
prediction = alpha' * X_pre_1 + (beta - Z * gamma)' * X_pre_2 + b_mean;
prediction = prediction';

% Display Testing Results
plot_results(y_test, prediction, d_test, 'Kernel Survival (Testing)');
compute_metrics(y_test, prediction, d_test);

%% save picture
for i = 1:13 
    saveas(i,['results_pic/fig',num2str(i),'.jpg']); 
end