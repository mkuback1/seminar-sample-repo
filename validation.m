%% Validation for GE and LU algorithms
% This script validates the algorithms myGE, myGEpp, myLU, and myPLU
% against an analytical solution derivated by hand. Using the built-in L2
% norm to evaluate errors.

clear; clc;
fprintf('======================================================\n')
fprintf('--- Validating against a known Analytical Solution ---\n')
fprintf('======================================================\n')

% Defining known test problem
A = [2 1 -1; 4 0 -1; -8 2 2];
x_true = [1; -3; 2];
b = A*x_true;

% Analytical LU decomposition
L_true = [1 0 0; 2 1 0; -4 -3 1];
U_true = [2 1 -1; 0 -2 1; 0 0 1];

% Analytical PLU decomposition
Lpp_true = [1 0 0; -1/4 1 0; -1/2 2/3 1];
Upp_true = [-8 2 2; 0 3/2 -1/2; 0 0 1/3];
P_true = [0 0 1; 1 0 0; 0 1 0];

%% GE and GEpp solver validation
fprintf('-------------------------------------\n')
fprintf('Testing Gaussian Elimination Solvers:\n')
fprintf('-------------------------------------\n')

xGE = myGE(A,b);
xGEpp = myGE(A,b);

error_GE = norm(x_true - xGE);
error_GEpp = norm(x_true - xGEpp);

fprintf(['Gaussian Elimination (myGE) error is ' num2str(error_GE) '\n']);
fprintf(['GE with Partial Pivoting (myGEpp) error is ' num2str(error_GEpp) '\n']);
fprintf('\n');
%% LU an PLU Decomposition Check
fprintf('-------------------------------------\n')
fprintf('Testing LU and PLU Decompositions:\n')
fprintf('-------------------------------------\n')

% LU decomposition
[L,U] = myLU(A);
error_L   = norm(L_true - L);
error_U   = norm(U_true - U);
error_LU  = norm(A - L*U);

fprintf('myLU Factor Breakdown: \n');
fprintf(['L error is ' num2str(error_L) '; U error is ' num2str(error_U) '\n']);
fprintf('\n');

% PLU decomposition
[P,Lpp,Upp] = myPLU(A);
error_P   = norm(P_true - P);
error_Lpp = norm(Lpp_true - Lpp);
error_Upp = norm(Upp_true - Upp);
error_PLU = norm(P*A - Lpp*Upp);

fprintf('myPLU Factor Breakdown: \n')
fprintf(['P error is ' num2str(error_P) '; Lpp error is ' num2str(error_Lpp) '; Upp error is ' num2str(error_Upp) '\n']);
fprintf('\n');

%% Verify solving with LU and PLU Decompositions
fprintf('-------------------------------------\n')
fprintf('Testing Solution using LU and PLU:\n')
fprintf('-------------------------------------\n')

% Solving with PLU
y = myForwardSolve(L,b);
xLU = myBackSolve(U,y);

error_solveLU = norm(x_true - xLU);
fprintf(['Error in solving with LU is ' num2str(error_solveLU) '\n']);

% Solving with PLU
y = myForwardSolve(Lpp,P*b);
xPLU = myBackSolve(Upp,y);

error_solvePLU =  norm(x_true - xPLU);
fprintf(['Error in solving with PLU is ' num2str(error_solvePLU) '\n']);



