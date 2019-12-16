function x = sorPlot ( A, b, xold, W, TOL, Nmax )
            A=[7, -3, 0, 0, 0;
                -3, 9, 1, 0, 0;
                0, 1, 3, -1, 0;
                0, 0, -1, 10, -4;
                0, 0, 0, -4, 6];
            b=[4; -6; 3; 7; 2];
            TOL=5*10^(-6);
            xold=[0,0,0,0,0];
            Nmax=100;
            n = length ( b );
            for i=1:n
               for j=1:n
                 if i == j
                    D(i,j)=A(i,j);
                 else
                    D(i,j)=0;
                 end
               end
            end
            for i=1:n
               for j=1:n
                 if i<= j
                    L(i,j)=0;
                    U(i,j)=A(i,j);
                 else
                    L(i,j)=A(i,j);
                    U(i,j)=0;
                 end
               end
            end
            for i=1:n
               for j=1:n
                 if i>= j
                    U(i,j)=0;
                 else
                    U(i,j)=A(i,j);
                 end
               end
            end
            W=2/(1+sqrt(1-(max(abs(eig(inv(D)*(L+U)))))^2))
            [r c] = size ( A );
            if ( c ~= n ) 
               disp ( 'sor error: matrix dimensions and vector dimension not compatible' )
               return
            end;
            xnew = zeros ( 1, n );
            if ( nargout == 0 )
               s = sprintf ( '%3d \t %10f ', 0, xold(1) );
               for j = 2 : n 
                 s = sprintf ( '%s%10f ', s, xold(j) );
               end;
               disp ( s );
            end;
            for its = 1 : Nmax
                xnew(1) = ( 1 - W ) * xold(1) + W * ( b(1) - sum ( A(1,2:n) .* xold(2:n) ) ) / A(1,1);
              for i = 2 : n-1
                  xnew(i) = ( 1 - W ) * xold(i) + W * ( b(i) - sum ( A(i,1:i-1) .* xnew(1:i-1) ) - sum ( A(i,i+1:n) .* xold(i+1:n) ) ) / A(i,i);
              end;
              xnew(n) = ( 1 - W ) * xold(n) + W * ( b(n) - sum ( A(n,1:n-1) .* xnew(1:n-1) ) ) / A(n,n);
              
                if ( nargout == 0 )
                 s = sprintf ( '%3d \t %10f ', its, xnew(1) );
                 for j = 2 : n 
                     s = sprintf ( '%s%10f ', s, xnew(j) );
                 end;
                 disp ( s );
              end;
                conv = max ( abs ( xnew - xold ) );
              if ( conv < TOL )
                 x = xnew;
                 return
              else
                 xold = xnew;
              end;
            end;
            disp ( 'sor error: maximum number of iterations exceeded' );
            if ( nargout == 1 ) x = xnew; 
            end;