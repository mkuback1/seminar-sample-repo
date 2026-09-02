function [xSoln] = myForwardSolve(L,b)
%
% [xSoln] = myForwardSolve(L,b)
%
% Performs forward substitution to solve lower triangular system Lx = b.
% Matrix L should be nxn lower triangular and b an nx1 vector.

n = length(b); % obtain size of system
xSoln = zeros(n,1); % initialize solution variable

xSoln(1) = b(1)/L(1,1); % calculate first entry of xSoln

for i = 2:n % loop forwards over entries of vector x
    sum = 0.0; % initialize sum 
    for j = 1:i-1 % loop over the entries in row i of U
        sum = sum + L(i,j)*xSoln(j); % update sum
    end
    xSoln(i) = (b(i)-sum)/L(i,i); % calculate entry i of x
end

end