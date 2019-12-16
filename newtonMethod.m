function newtonMethod(f,fp,x0,tol, maxIter)
% Newton Method - A programatic implementation of the
%   newtons method of root finding
%   Accepts a function as f
%   Accepts the functions derivative as fp
%   Accepts an initial guess as x0
%   Accepts a tolerance as tol
%   Accepts a max iteration as maxIter
%   A sample of this would be:
%       syms x
%       newtonMethod(x - cos(x), 0, 10^-5, 100)
iter = 0;
xn = @(x) x - ( f(x) / fp(x) );

if (f(x0) == 0)
    fprintf('The initial guess: %10.6f is already the zero. 0 iterations\n', x0);
    return
end

fprintf('\tIteration\t|\t\t\tx\t\t|\t\tF(x)\t\n');
while (true)
    if (iter >= maxIter)
        fprintf("Max iteration reached. No root found.\n");
        return
    end
    if (abs(f(x0)) <= tol)
        fprintf('\t%d\t\t\t|\t%10.12f\t|\t%10.12f\t\n', [iter, x0, f(x0)]);
        fprintf("Root found within tolerance.\nRoot is %10.12f and took %d iterations.\nThe value of F(x) is %10.12f.\n", [x0, iter, f(x0)]);
        return
    end
    if (f(x0) == 0)
        fprintf('\t%d\t\t\t|\t%10.12f\t|\t%10.12f\t\n', [iter, x0, f(x0)]);
        fprintf("Root found!\nRoot is %10.12f and took %d iterations.\nThe value of F(x) is %10.12f.\n", [x0, iter, f(x0)]);
        return
    end
    fprintf('\t%d\t\t\t|\t%10.12f\t|\t%10.12f\t\n', [iter, x0, f(x0)]);
    x0 = xn(x0);
    iter = iter + 1;
end
