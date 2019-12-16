function MAIN()
%MAIN Main driver function for HW16
%   Builds all matrices and calls other functions
%   Outputs status code 0 upon completion

% Declare variables for root finding methods
maxIter = 100;
tol = 10^-12;

fprintf('\nMatrix 1 Part A, a=c=1, b=3\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t\n')
n = 10;
a = 1;
b = 3;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

fprintf('\nMatrix 1 Part B, a=c=1, b=3\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tJacobi\t\t|\tGauss-Seidel\t|\t\tJacobi\t\t|\tGauss-Seidel\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

fprintf('\nCalculating Optimal Omega for Matrix 1\n')
fprintf('\tSize\t|\tOmega\t|\tForward Infinity Error\t|\tOptimal Omega\t|\tRun Time\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;    
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end

fprintf('\nMatrix 2 Part A, a=c=-1, b=2\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t\n')
n = 10;
a = -1;
b = 2;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

fprintf('\nMatrix 2 Part B, a=c=-1, b=2\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tJacobi\t\t|\tGauss-Seidel\t|\t\tJacobi\t\t|\tGauss-Seidel\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

fprintf('\nCalculating Optimal Omega for Matrix 2\n')
fprintf('\tSize\t|\tOmega\t|\tForward Infinity Error\t|\tOptimal Omega\t|\tRun Time\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end

fprintf('\nMatrix 3 Part A, a=b=c=1\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t|\t\tGauss\t\t|\tGauss w/ PP\t\t\n')
n = 10;
a = 1;
b = 1;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x = m\y;
tic;
gX = Gauss(m, y);
gT = toc;
tic;
pX = GaussPP(m ,y);
pT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, gT, pT, norm(gX - x), norm(pX - x)])

fprintf('\nMatrix 3 Part B, a=b=c=1\n')
fprintf('\tSize\t\t|\t\t\t\tRun Time\t\t\t\t|\t\t\tForward Infinity Error\n')
fprintf('\t\t\t\t|\t\tJacobi\t\t|\tGauss-Seidel\t|\t\tJacobi\t\t|\tGauss-Seidel\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
tic;
jX = Jacobi(m, y, x0, maxIter, tol);
jT = toc;
tic;
gX = GaussS(m ,y, x0, maxIter, tol);
gT = toc;
fprintf('\t%d\t\t|\t%.6e\t|\t%.6e\t|\t%.6e\t|\t%.6e\n', [n, jT, gT, norm(jX - x), norm(gX - x)])

fprintf('\nCalculating Optimal Omega for Matrix 3\n')
fprintf('\tSize\t|\tOmega\t|\tForward Infinity Error\t|\tOptimal Omega\t|\tRun Time\n')
n = 10;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 100;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
n = 1000;
m = matrixBuilder(a,b,a,n);
y = ones([n,1]);
x0 = zeros([n,1]);
x = m\y;
minNorm = 10^32;
optOmega = 0;
for omega = -1.5:.1:1.5
    if (omega == 0)
        continue
    end
    tic;
    sX = SOR(m, y, omega, x0, maxIter, tol);
    sT = toc;
    if (norm(sX - x) < minNorm)
        minNorm = norm(sX - x);
        optOmega = omega;
    end
    fprintf('\t%d\t|\t%.2f\t|\t\t%.6e\t\t|\t\t%.2f\t\t|\t%.6e\n', [n, omega, norm(sX - x), optOmega, sT])
end
end