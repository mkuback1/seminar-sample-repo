function [L,U] = myLU(A)
%
% [L,U] = myLU(A)
%
% Performs LU decomposition without pivot to factor (if possible) a matrix
% A into A = LU.

% Compute dimensions of input
[m,n] = size(A);
if m~=n
    error('Matrix Must be Square')
end

% Create place holder for multipliers
L = eye(n,n);

% Elimination Phase (n-1 Steps)
for i = 1:n-1 % loop accross columns 1 through n-1 of augmented matrix

    if A(i,i) == 0 % stop algorithm if will divide by zero
        disp('Error: pivot is zero, cannot continue.')
        return % stop function routine
    end

    % compute multiplier vector for all entries below current pivot
    L(i+1:n,i) = A(i+1:n,i)/A(i,i);

    % perform row ops below current row across all columns
    A(i+1:n,:) = A(i+1:n,:)-L(i+1:n,i)*A(i,:); 
end 

U = A(:,1:n); % save U

end


