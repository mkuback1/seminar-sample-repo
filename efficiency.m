%% Efficiency 
% This script runs various numerical experiments to examine the runtime
% efficiency of solving Ax = b single and multiple times with a fixed
% matrix A and different vectors b using:
%       Method 1: Gaussian Elimination with Partial Pivoting (myGEpp)
%       Method 2: PLU Factorization (myPLU) + Forward/Backward
%       Substitution

clear; close; clc;

% Define Parameters
n_start = 10; % initial size of linear system
%n_total = 7; % will run experiments 1-2 for n_total matrices of increasing size
m = 50; % number of times we will solve the linear system in experiment 2


%% Experiment 1: Timing Single Solves for both Methods as n increases

fprintf('=====================================================\n');
fprintf('--- Experiment 1: Single Solve, for increasing n  --- \n')
fprintf('=====================================================\n');

% Pre-allocating time vectors
time1_GE = zeros(n_total,1);
time1_PLU = zeros(n_total,1);
n_size = zeros(n_total,1); % for storing matrix sizes

n_curr = n_start;

% Loop increasing matrix size and timing the solves
for k = 1:n_total
    n_size(k) = n_curr; % store current size
    A = rand(n_curr); % generate new random matrix A
    x_true = rand(n_curr,1); % generate true solution
    b = A*x_true; % generate RHS
    % Timing methods
    % Method 1
    tic
    myGEpp(A,b);
    time1_GE(k) = toc;
    % Method 2
    tic
    [P,L,U] = myPLU(A);
    Pb = P*b;
    y = myForwardSolve(L,Pb);
    myBackSolve(U,y);
    time1_PLU(k) = toc;
    % Update n
    %n_curr = n_curr + 10; % current matrix size
    n_curr = 2*n_curr;
end

%% Print table of results
Column_Names = {'n','Method 1 (GE)','Method 2 (PLU)'};
exp1_table = table(n_size,time1_GE,time1_PLU,'VariableNames',Column_Names);
disp(exp1_table);

%% Plot results on log scale so O(n^3) behavior appears as line w/slope 3
fig1 = figure(1);
hold on
loglog(n_size,time1_GE,'r','Linewidth',2)
loglog(n_size,time1_PLU,'b','LineWidth',2)
loglog(n_size(end-3:end),((n_size(end-3:end))./1000).^3,'--k','Linewidth',2)
hold off;
xlabel('Matrix size (n)')
ylabel('time (seconds)')
xscale('log'); yscale('log'); 
set(gca,'FontSize',14)
title('Time for Single Solve, for increasing n (log scale)')
legend({'Method 1 (GE)', 'Method 2 (PLU)', 'Line with Slope 3'},'Location','northwest')


%% Experiment 2: Timing Multiple Solves for both Methods, Increasing n
fprintf('=====================================================\n');
fprintf('--- Experiment 2: Multi Solves, for increasing n  --- \n')
fprintf('=====================================================\n');
fprintf('Number of Solves (m):    %d\n\n', m);

% Pre-allocating time vectors
time2_GE = zeros(n_total,1);
time2_PLU = zeros(n_total,1);
n_size = zeros(n_total,1); % for storing matrix sizes

n_curr = n_start;

% Loop over increasing matrix size and timing the solves
for k = 1:n_total
    n_size(k) = n_curr; % store current size
    A = rand(n_curr); % generate new random matrix A
    X = rand(n_curr,m); % generate m true solutions
    B = A*X; % RHS vectors for m true solutions
    % Timing methods
    % Method 1 (m solves)
    tic
    for j = 1:m
        myGEpp(A,B(:,j));
    end
    time2_GE(k) = toc;
    % Method 2 (m solves)
    tic
    [P,L,U] = myPLU(A);
    for j = 1:m
        Pb = P*B(:,j);
        y = myForwardSolve(L,Pb);
        myBackSolve(U,y);
    end
    time2_PLU(k) = toc;
    % Update n
    n_curr = 2*n_curr;
end

speedup2 = time2_GE ./ time2_PLU;

%% Print table of results
Column_Names = {'n','Method 1 (GE)','Method 2 (PLU)','Speedup'};
exp2_table = table(n_size,time2_GE,time2_PLU,speedup2,'VariableNames',Column_Names);
disp(exp2_table);

%% Plot results on log scale so O(n^3) behavior appears as line w/slope 3
fig2 = figure(2);
hold on
loglog(n_size,time2_GE,'r','Linewidth',2)
loglog(n_size,time2_PLU,'b','LineWidth',2)
loglog(n_size(end-3:end),m*((n_size(end-3:end))./1000).^3,'--k','Linewidth',2)
hold off;
xlabel('Matrix size (n)')
ylabel('time (seconds)')
xscale log; yscale log
set(gca,'FontSize',14)
title(['Time for ' num2str(m) '-Solves, for increasing n (log scale)'])
legend({'Method 1 (GE)', 'Method 2 (PLU)', 'Line with Slope 3'},'Location','northwest')
