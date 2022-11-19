
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-418: Convex Optimization | MOSEK Installation Test %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

% define decision variables
x = sdpvar(2,1);

% define objective function
objective = sum(x);

% define constraints
constraints = [-1 <= x(1); x(1) <= +1; -1 <= x(2); x(2) <= +1;];

% specify solver settings
opt_settings = sdpsettings('solver', 'mosek', 'verbose', 0);

% run solver
diagnosis = optimize(constraints, objective, opt_settings);

% display solver report
disp('solver report:');
disp(diagnosis);

% retireve and display optimal objective value
disp('optimal objective value:');
opt_objective = value(objective);
disp(opt_objective);

% retrieve and display optimal solution values
disp('optimal solution values:');
opt_x = value(x);
disp(opt_x);