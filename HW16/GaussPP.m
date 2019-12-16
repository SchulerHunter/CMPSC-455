function [x] = GaussPP(M,b)
%GAUSSPP A function to perform gaussian elimination
%   Utilizaes partial pivoting in gaussian elimination
%   This does not keep track of the r vector for performance
A = [M b];
n = size(A);
n = n(1);

% Diagonal iterator
for i = 1:(n-1)
    % Partial pivot to max value in the i column
    [val, row] = max(abs(A(i:end,i)));
    row = row + i - 1;
    temp = A(row,1:end);
    A(row, 1:end) = A(i, 1:end);
    A(i, 1:end) = temp;
    % Row iterator to gauss elimenate each row
    for j = (i+1):n
        mult = A(j,i)/A(i,i);
        if (mult ~= 0)
            % Column iterator for subtracting from the second column
            for k = i:n+1
                A(j,k) = A(j,k) - mult*A(i,k);
            end
        end
    end
end
% Backward substitution
x = zeros(n,1);
% Row iterator for backwards substitution
for i = n:-1:1
    x(i) = (A(i,n+1) - sum(A(i,i+1:n)*x(i+1:n))) ./ A(i,i);
end
end

