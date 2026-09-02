%% Efficiency
% This script runs a numerical experiment to examine the runtime efficiency
% of solving Ax = b multiple times with a fixed matrix A and different
% vectors b using:
%       Method 1: Gaussian Elimination with Partial Pivoting (myGEpp)
%       Method 2: PLU Factorization (myPLU) + Forward/Backward
%       Substitution

clear; clc;

% Define Parameters
n = 100; % size of the linear system
m = 1000; % number of times we will solve the linear system

% Generate random linear systems for use in experiments
A = rand(n); % generate random nxn matrix A
x_true = rand(n,1);
b = A*x_true;
X = rand(n, m); % generate random matrix X whose columns are the true solutions
B = A*X; % this matrix stores the RHS vectors for multiple solves

%% Experiment 1: Timing a Single Solve for both Methods

fprintf('=====================================================\n');
fprintf('---     Running Experiment for Single Solve      --- \n')
fprintf('=====================================================\n');
fprintf('Matrix Size (n):         %d x %d\n', n, n);

% Method 1: GE with partial pivoting
tic
myGEpp(A,b);
time_GEpp = toc;

% Method 2: PLU Method 
tic
[P,L,U] = myPLU(A);
Pb = P*b;
y = myForwardSolve(L,Pb);
myBackSolve(U,y);
time_PLU = toc;

fprintf('-----------------------------------------------------\n');
fprintf('                  Results                            \n');
fprintf('-----------------------------------------------------\n');
fprintf('  Method 1 (Single GEpp):        %8.4f seconds\n', time_GEpp);
fprintf('  Method 2 (Single PLU + Single Solve):  %8.4f seconds\n', time_PLU);
fprintf('-----------------------------------------------------\n');
fprintf('\n');


%% Experiment 2: Timing Multiple Solves for both Methods

fprintf('=====================================================\n');
fprintf('---    Running Experiment for Multiple Solves    --- \n')
fprintf('=====================================================\n');
fprintf('Matrix Size (n):         %d x %d\n', n, n);
fprintf('Number of Solves (m):    %d\n\n', m);

% Solve the linear system m times using GEpp
tic
for k = 1:m
    myGEpp(A,B(:,k));
end
total_time_GEpp = toc;

% Solve the linear system m times using PLU decompositoin
tic
[P,L,U] = myPLU(A); % only need to calculate PLU decomp one time
for k = 1:m
    Pb = P*B(:,k);
    y = myForwardSolve(L,Pb);
    myBackSolve(U,y);
end
total_time_PLU = toc;

speedup = total_time_GEpp / total_time_PLU;

fprintf('-----------------------------------------------------\n');
fprintf('                  Results                            \n');
fprintf('-----------------------------------------------------\n');
fprintf('  Method 1 (Repeated GEpp):        %8.4f seconds\n', total_time_GEpp);
fprintf('  Method 2 (Single PLU + Multi. Solves):  %8.4f seconds\n', total_time_PLU);
fprintf('-----------------------------------------------------\n');
fprintf('  PLU Factorization Speedup:       %8.2fx faster\n', speedup);
fprintf('-----------------------------------------------------\n');