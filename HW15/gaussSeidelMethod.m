function iter=gaussSeidelMethod(A,b,x0,maxIter,tol)
%A-- a nXn matrix
%b-- a nX1 vector
%x-- a solution of Ax=b
D=diag(A);
L=tril(A)-diag(D);
U=triu(A)-diag(D);
iter=0;
n=size(A,1);
x=x0;

%fprintf('\tIteration\t|\tResidual\t\n');
while (true)
    for i=1:n
        x(i)=(b(i)-L(i,:)*x-U(i,:)*x0)/D(i);
    end
    iter=iter+1;
    %fprintf('\t%d\t\t\t|\t%d\t\t\n', [iter, norm(A*x-b)]);
    if iter>maxIter
        %fprintf('Max iterations reached\n');
        break
    end
    if norm(x-x0)<tol
        %fprintf('Forward error less than tolerance\n');
        break
    end
    if norm(A*x-b)<tol
        %fprintf('Backward error less than tolerance\n');
        break
    end
    x0=x;
end