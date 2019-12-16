fx = @(x) sin(x)-cos(x);
tol = 10^-6;
a = 0;
b = 1;
if (fx(a) * fx(b) > 0)
    disp("No zero can be determined according to IVT");
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
fprintf('%10.6f\n',c);

