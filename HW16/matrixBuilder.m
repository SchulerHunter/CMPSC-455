function [M] = matrixBuilder(a,b, c, n)
%MATRIXBUILDER Builds a nxn matrix
%   Off diagonals are a and c
%   Diagonals are b
M = zeros(n);
for i = 1:n
    if (i-1 > 0)
        M(i, i-1) = a;
    end
    M(i,i) = b;
    if (i+1 < n+1)
        M(i, i+1) = c;
    end
end
end

