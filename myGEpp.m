function [x] = myGEpp(A,b)
%
% [x] = myGEpp(A,b)
%
% Performs Gaussian Elimination with partial pivoting to solve the linear
% system Ax = b.  A must be nxn and b nx1.

% Compute dimensions of input
[m,n] = size(A);
if m~=n
    error('Matrix Must be Square')
end

% Define augmented matrix
C = [A b]; % define augmented matrix

% Create place holder for multipliers
m = eye(n,n);

% Elimination Phase (n-1 Steps)
for i = 1:n-1 % loop accross columns 1 through n-1 of augmented matrix
    
    % find the largest entry and its location in column i
    % this will be our pivot at this step
    [pivot,k] = max(abs(C(i:n,i))); 
    % k is the position of the pivot relative to row i
    % this means the pivot is in row i+k-1 in the overall matrix
    
    % swap row i (current pivot row) and the desired pivot row i+k-1
    if k>1 % if desired pivot row is below row i
        C([i,i+k-1],:) = C([i+k-1,i],:); % Swap rows i and i+k-1
    end

    if C(i,i) == 0 % stop algorithm if will divide by zero
        disp(['Error: pivot is zero, cannot continue.'])
        return % stop function routine
    end

    % compute multiplier vector for all entries below current pivot
    m(i+1:n,i) = -C(i+1:n,i)/C(i,i);

    % perform row ops below current row across all columns
    C(i+1:n,:) = C(i+1:n,:)+m(i+1:n,i)*C(i,:); 
end 

if C(n,n) == 0 % if b(n) after elimination ends up being 0
    error('A is numerically singular') % stop algorithm before div by 0
end

% Backsubstitution Phase
x = zeros(n,1); % initialize solution vector x

x(n,1) = C(n,n+1)/C(n,n); % calculate x(n)

for i = n-1:-1:1 % Calculuate x(n-1) through x(1)
    x(i) = (C(i,n+1)-C(i,i+1:n)*x(i+1:n))/C(i,i);
end

end