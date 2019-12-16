function [x] = Gauss(M, b)
%GAUSS Performs gaussian elimination on a matrix
%   Uses guassian elimination to solve
%   The system Mx=b
%   Finally uses back substitution to finish
A = [M b];
n = size(A);
n = n(1);
% Perform Gaussian elimination
for i = 1:(n-1)
    while(A(i,i) == 0)
        j = i + 1;
        if (A(j,i) ~= 0)
            temp = A(j,1:end);
            A(j, 1:end) = A(i, 1:end);
            A(i, 1:end) = temp;
        end
        if (j > n)
            fprintf('No Solution')
            return
        end
    end
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

