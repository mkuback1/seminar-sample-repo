function [xSoln] = myBackSolve(U,b)
%
% [xSoln] = myBackSolve(U,b)
%
% Performs back substitution to solve an upper triangular system Ux = b.
% Matrix U should be nxn upper triangular and b an nx1 vector.

n = length(b); % obtain size of system
xSoln = zeros(n,1); % initialize solution variable

xSoln(n) = b(n)/U(n,n); % calculate last entry of xSoln

for i = n-1:-1:1 % loop backwards over entries of vector x
   xSoln(i) = (b(i)-U(i,i+1:n)*xSoln(i+1:n))/U(i,i);
end

% Same Routine, but written with nested for-loop
% for i = n-1:-1:1 % loop backwards over entries of vector x
%     sum = 0.0; % initialize sum 
%     for j = i+1:n % loop over the entries in row i of U
%         sum = sum + U(i,j)*xSoln(j); % update sum
%     end
%     xSoln(i) = (b(i)-sum)/U(i,i); % calculate entry i of x
% end

end