function [zero] = bisectionMethod(fx,bounds,tol)
%Bisection Method - A program to use bisection method to find zero
%   Accepts an anonymously defined function as fx
%   Accepts a list for bounds (only uses the first two points)
%   Accepts a tolerance for accuracy (minimum is 10^-16)
%   For question six in the homework, the function would run as:
%       bisectionMethod(@(x) sin(x)-cos(x),[0, 1], 10^-6)
a = bounds(1);
b = bounds(2);
if (fx(a) * fx(b) > 0)
    zero = "No zero can be determined according to IVT";
    return
end
while (abs(a-b) > tol)
    c = mean([a, b]);
    if (fx(c)*fx(a) < 0)
        b = c;
    else
        a = c;
    end
end
zero = sprintf('%10.6f',c);
end

