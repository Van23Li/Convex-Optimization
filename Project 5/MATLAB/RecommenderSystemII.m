%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | Convex Optimization | Project 5 (graded), Solution (ii) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start completely fresh 
clear all; close all; clc;

%% Load the data

% Load the full (Rcompl) rating matrices (xxuyyi means #xx users #yy items)
load 100u10icompl;

% infer the size of the matrix (m x n)
[m,n] = size(Rcompl);

%% Question 5 (ii) (use Yalmip)
% Use your work from Question 5 (i)

% Generate several Rpartial, that is, matrices that
% contain p*100 percent of the data from Rcompl. 
% Compute the relative error (as given in the exercise) 
% for p=0.05, 0.1, 0.15, ...., 0.95. (steps of 0.05). 

Z = 19; 
Error = zeros(Z,1); 
for z=1:Z
    Rpartial = zeros(m,n);
    p = z*0.05
    for i=1:m
        for j=1:n
            if rand(1)>=1-p 
                Rpartial(i,j)=Rcompl(i,j); 
            end
        end
    end

    % Define the decision variables: 
    % Note! X is not symmetric. 
    % Without the 'full' constraint: wrong. 

    X = sdpvar(m, n, 'full');
    Lambda11 = sdpvar(n, n);
    Lambda22 = sdpvar(m, m);

    % Define the objective function:
    % Note! You cannot always group the Lambda matrices (dimension might
    % differ).
    obj = (trace(Lambda11) + trace(Lambda22)) / 2;


    % Define the constraints:
    % Constraints enforcing observed ratings (use Rpartial) 
    % that is, if Rij==0, the value is "not observed". 
    constr = [X(Rpartial~=0) == Rpartial(Rpartial~=0)];

    % Matrix version: 
    % k=find(Rpartial(:)); 
    % constr = [constr; X(k) == Rpartial(k)];

    % Define the constraints enforcing the range beween 1 and 5
    % constr = [constr, ismember(X, [1,2,3,4,5])];
    constr = [constr, X >= 1, X <= 5];

    % Define the semidefinite constraint
    constr = [constr, [Lambda11, -X'; -X, Lambda22] >= 0];
    
    % Solve the SDP

    % Make sure you have a SDP solver (like SDPT3) 
    % If you have bizzare memory errors, consider changing solvers or try to
    % run on a different operating system. 
    % 100u10i can take a few minutes (tic-toc measures the time) 
    % Set 'verbose' to 1 if you want to see the (any) progress. 
    tic
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);    % mosek	gurobi sdpt3
    diagnosis = optimize(constr, obj, ops);
    toc

    % Retrieve and round the optimal solution
    % Convice yourself of the validity of the rounding. 
    if diagnosis.problem == 0
         % Extract value
        Xopt_real = value(X); 
        Xopt_int = round(Xopt_real);
        Error(z,1) = norm((Xopt_int - Rcompl),'fro') / norm((Rcompl),'fro');
    else
        display('Hmm, something went wrong!');
        diagnosis.info
        yalmiperror(diagnosis.problem)
    end
end
%%  Plot the outcome
% Note! the outcome is evidently random,
% you can at most observe a trend. 

figp=figure;
p=5:5:95; 
plot(p,Error,'-ok'); 
grid on 
set(gca,'FontSize',12) 
xlabel('percentage observed','interpreter','latex','FontSize',20);
ylabel('relative error','interpreter','latex','FontSize',20);