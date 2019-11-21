function iter=jacobiMethod(A,b,x0,maxIter,tol)
%A-- a nXn matrix
%b-- a nX1 vector
%x-- a solution of Ax=b
D=diag(A);
L=tril(A)-diag(D);
U=triu(A)-diag(D);
iter=0;

%fprintf('\tIteration\t|\tResidual\t\n');
while (true)
    x=-(L+U)*x0./D+b./D;
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