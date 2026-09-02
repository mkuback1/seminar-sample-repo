function [P,L,U] = myPLU(A)
%
% [P,L,U] = myPLU(A)
%
% Performs PLU decomposition with partial pivoting to factor (if possible)
% a matrix A into PA = LU.

% Compute dimensions of input
[m,n] = size(A);
if m~=n
    error('Matrix Must be Square')
end

P = eye(n,n);
L = eye(n,n);


% Elimination Phase (n-1 Steps)
for i = 1:n-1 % loop accross columns 1 through n-1 of augmented matrix
    [pivot,k] = max(abs(A(i:n,i))); 
    % k is the position of the pivot relative to row i
    % this means the pivot is in row i+k-1 in the overall matrix
    
    % swap row i (current pivot row) and the desired pivot row i+k-1
    if k>1 % if desired pivot row is below row i
        A([i,i+k-1],:) = A([i+k-1,i],:); % Swap rows i and i+k-1 of A
        P([i,i+k-1],:) = P([i+k-1,i],:); % Swap rows i and i+k-1 of P
        L([i,i+k-1],1:i-1) = L([i+k-1,i],1:i-1); % Swap rows up to col i-1
    end

    if A(i,i) == 0 % stop algorithm if will divide by zero
        disp('Error: pivot is zero, cannot continue.')
        return % stop function routine
    end

    % compute multipliers for all entries below current pivot
    L(i+1:n,i) = A(i+1:n,i)/A(i,i); 

    % perform row ops below current row across all columns
    A(i+1:n,:) = A(i+1:n,:)-L(i+1:n,i)*A(i,:); 
end 

U = A(:,1:n); % save U

end
