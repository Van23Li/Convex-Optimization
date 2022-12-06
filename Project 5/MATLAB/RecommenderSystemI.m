%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | Convex Optimization | Project 5 (graded), Solution (i) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start completely fresh 
clear all; close all; clc;

%% Load the data

% Load the partially observed (Rpartial) rating matrices (xxuyyi means #xx users #yy items)
load 100u20ipartial;
% The data corresponds to user ratings songs from the Rolling Stone top 20.

% 1		Bob Dylan - Like a Rolling Stone
% 2		The Rolling Stones - Satisfaction
% 3		John Lennon - Imagine
% 4		Marvin Gaye - What's Going On
% 5		Aretha Franklin - Respect
% 6		The Beach Boys - Good Vibrations
% 7		Chuck Berry - Johnny B. Goode
% 8		The Beatles - Hey Jude
% 9		Nirvana - Smells Like Teen Spirit
% 10	Ray Charles - What'd I Say (Parts 1 And 2)
% 11	The Who - My Generation
% 12	Sam Cooke - A Change Is Gonna Come
% 13	The Beatles - Yesterday
% 14	Bob Dylan - Blowin' in the Wind
% 15	The Clash - London Calling
% 16	The Beatles - I Want to Hold Your Hand
% 17	Jimi Hendrix - Purple Haze
% 18	Chuck Berry - Maybellene
% 19	Elvis Presley - Hound Dog
% 20	The Beatles - Let It Be

% The rating is on a scale from 1 (very bad) to 5 (very good). A 0 means
% the entry is missing.

% infer the size of the matrix (m x n)
[m,n] = size(Rpartial);

% Eric likes Jimi Hendrix (song 17), what do you suggest for Eric? 
% Start by adding the user data of Eric to Rpartial
R_Eric = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0];
% Redefine Rpartial (and its dimension)
Rpartial = [Rpartial; R_Eric];
m = m + 1; 

%% Question 5 (i) (use Yalmip)

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

% Define the constraints enforcing the range beween 1 and 5
% constr = [constr, ismember(X, [1,2,3,4,5])];
constr = [constr, X >= 1, X <= 5];

% Define the semidefinite constraint
constr = [constr, [Lambda11, -X'; -X, Lambda22] >= 0];

%% Solve the SDP

% Make sure you have a SDP solver (like SDPT3) 
% If you have bizzare memory errors, consider changing solvers or try to
% run on a different operating system. 
% 100u20i can take a few minutes (tic-toc measures the time) 
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
else
    display('Hmm, something went wrong!');
    diagnosis.info
    yalmiperror(diagnosis.problem)
end

%% Suggestion to Eric
% What do you suggest to Eric? 
% Report the max indices 
Xopt_int(end,:)