% Homework 15
% Author: Hunter Schuler
maxIter=1000;
% Question 1
fprintf('Question 1:\n')
x0=[0;0;0];
fprintf('\tStop Tol\t|\t\t\t# of Iterations\n')
fprintf('\t\t\t\t|\t\t\tJacobi\t\t|\t\tGauss-Seidel\n')
fprintf('\t\tx\t\t|\tx = 3\t|\tx = 2\t|\tx = 3\t|\tx = 2\n')
tol=10^-5;
fprintf('\t%.2e\t|\t%d\t\t|\t%d\t\t|\t%d\t\t|\t%d\n', [tol, jacobiMethod([3,-1,1;1,3,-1;-1,1,3],[3;3;3],x0,maxIter,tol), jacobiMethod([2,-1,1;1,2,-1;-1,1,2],[2;2;2],x0,maxIter,tol), gaussSeidelMethod([3,-1,1;1,3,-1;-1,1,3],[3;3;3],x0,maxIter,tol), gaussSeidelMethod([2,-1,1;1,2,-1;-1,1,2],[2;2;2],x0,maxIter,tol)]);
tol=10^-10;
fprintf('\t%.2e\t|\t%d\t\t|\t%d\t\t|\t%d\t\t|\t%d\n\n', [tol, jacobiMethod([3,-1,1;1,3,-1;-1,1,3],[3;3;3],x0,maxIter,tol), jacobiMethod([2,-1,1;1,2,-1;-1,1,2],[2;2;2],x0,maxIter,tol), gaussSeidelMethod([3,-1,1;1,3,-1;-1,1,3],[3;3;3],x0,maxIter,tol), gaussSeidelMethod([2,-1,1;1,2,-1;-1,1,2],[2;2;2],x0,maxIter,tol)]);

% Question 2
fprintf('Question 2:\n')
tol=10^-5;
b=[-1;-1;2];
fprintf('\tDelta\t\t|\t# of Iterations with Tol 10^-5\n')
fprintf('\t\t\t\t|\tJacobi\t\t|\tGauss-Seidel\n')
delta=1;
fprintf('\t%.2e\t|\t%d\t\t\t|\t%d\n', [delta, jacobiMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol), gaussSeidelMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol)]);
delta=10^-1;
fprintf('\t%.2e\t|\t%d\t\t\t|\t%d\n', [delta, jacobiMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol), gaussSeidelMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol)]);
delta=10^-2;
fprintf('\t%.2e\t|\t%d\t\t|\t%d\n', [delta, jacobiMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol), gaussSeidelMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol)]);
delta=0;
fprintf('\t%.2e\t|\t%d\t\t|\t%d\n\n', [delta, jacobiMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol), gaussSeidelMethod([1+delta,-1,0;-1,2+delta,-1;0,-1,1+delta],b,x0,maxIter,tol)]);

% Question 3
fprintf('Question 3:\n')
A=[10,1,1;1,-10,1;1,1,10];
b=[1;2;3];
fprintf('\tStop Tol\t|\t\t\t\t\t# of Iterations of SOR\n')
fprintf('\t\t\t\t|\tomega=1.1\t|\tomega=1.4\t|\tomega=1.6\t|\tomega=1.9\n')
tol=10^-5;
fprintf('\t%.2e\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\n', [tol, SOR(A,b,1.1,x0,maxIter,tol), SOR(A,b,1.4,x0,maxIter,tol), SOR(A,b,1.6,x0,maxIter,tol), SOR(A,b,1.9,x0,maxIter,tol)])
tol=10^-10;
fprintf('\t%.2e\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\n\n', [tol, SOR(A,b,1.1,x0,maxIter,tol), SOR(A,b,1.4,x0,maxIter,tol), SOR(A,b,1.6,x0,maxIter,tol), SOR(A,b,1.9,x0,maxIter,tol)])

% Question 4
fprintf('Question 4:\n')
fprintf('\tSize N\t|\t\t\t# of Iterations\n')
fprintf('\t\t\t\t|\t\t\tJacobi\t\t\t\t\t|\t\tGauss-Seidel\n')
fprintf('\t\t\t\t|\tlambda=100\t|\tlambda=10\t\t|\tlambda=100\t|\tlambda=10\n')
n=10;
R=rand(n);
I=eye(n);
b=rand(n,1);
x0=zeros(n,1);
fprintf('\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t\t|\t%d\t\t\t|\t%d\n', [n, jacobiMethod(R-100*I,b,x0,maxIter,tol), jacobiMethod(R-10*I,b,x0,maxIter,tol), gaussSeidelMethod(R-100*I,b,x0,maxIter,tol), gaussSeidelMethod(R-10*I,b,x0,maxIter,tol)])
n=100;
R=rand(n);
I=eye(n);
b=rand(n,1);
x0=zeros(n,1);
fprintf('\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\n', [n, jacobiMethod(R-100*I,b,x0,maxIter,tol), jacobiMethod(R-10*I,b,x0,maxIter,tol), gaussSeidelMethod(R-100*I,b,x0,maxIter,tol), gaussSeidelMethod(R-10*I,b,x0,maxIter,tol)])
n=1000;
R=rand(n);
I=eye(n);
b=rand(n,1);
x0=zeros(n,1);
fprintf('\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\t\t\t|\t%d\n', [n, jacobiMethod(R-100*I,b,x0,maxIter,tol), jacobiMethod(R-10*I,b,x0,maxIter,tol), gaussSeidelMethod(R-100*I,b,x0,maxIter,tol), gaussSeidelMethod(R-10*I,b,x0,maxIter,tol)])