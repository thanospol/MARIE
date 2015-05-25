% Preconditioned GMRES internal iteration routine
%
% Solves A z = r for z, then returns x + z
% Generates Arnoldi relation A V(:,1:m) = V(:,1:m+1) H
%
% INPUT:  A      N-by-N matrix
%         x      current solution vector
%         r      N-by-1 preconditioned residual vector
%         m      number of GMRES iterations to perform
%         L1     first left preconditioner for A
%         L2     second left preconditioner for A
%         R1     first right preconditioner for A
%         R2     second right preconditioner for A
%         tol    specifies the tolerance of the method
% OUTPUT: x      updated solution vector
%         r      preconditioned residual vector
%         k      number of GMRES iterations actually performed
%         resvec vector containing norm of residual at each iteration of GMRES
%
%

function [x,r,k,resvec] = pgmres_iter(A,x,r,m,L1,L2,R1,R2,tol)

if(isempty(L1))
   existL1 = 0;
else
   existL1 = 1;
  [L1type,L1fun,L1fcnstr] = iterchk(L1);
end
if(isempty(L2))
   existL2 = 0;
else
   existL2 = 1;
   [L2type,L2fun,L2fcnstr] = iterchk(L2);
end
if(isempty(R1))
   existR1 = 0;
else
   existR1 = 1;
  [R1type,R1fun,R1fcnstr] = iterchk(R1);
end
if(isempty(R2))
   existR2 = 0;
else
   existR2 = 1;
   [R2type,R2fun,R2fcnstr] = iterchk(R2);
end

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);

% Preallocate and Initialize V
V = zeros(size(r,1),m+1);
V(:,1) = r / norm(r);

% Preallocate H and resvec
H = zeros(m+1,m);
resvec = zeros(m,1);

for k = 1:m
    
    w = V(:,k);
    
   % Find w using right preconditioning if available.
   % R1*R2 approx A, then w = inv(R1*R2)*v = R2\(R1\v)
   % w = R2 \ ( R1 \ V(:,k) );
   if(existR1) % first right preconditioning
       if strcmp(R1type,'matrix')
           w = R1 \ w;
       else
           w = iterapp('mtimes',R1fun,R1type,R1fcnstr,w);
       end
   end
   if(existR2) % second right preconditioning
       if strcmp(R2type,'matrix')
           w = R2 \ w;
       else
           w = iterapp('mtimes',R2fun,R2type,R2fcnstr,w);
       end
   end

   
   % Apply the system operator
   % w = A*w;
   w = iterapp('mtimes',afun,atype,afcnstr,w);
   
   % Apply the Left preconditioning if available.
   % w = L2 \ ( L1 \ w );
   if(existL1) % first left preconditioning
       if strcmp(L1type,'matrix')
           w = L1 \ w;
       else
           w = iterapp('mtimes',L1fun,L1type,L1fcnstr,w);
       end
   end
   if(existL2) % second left preconditioning
       if strcmp(L2type,'matrix')
           w = L2 \ w;
       else
           w = iterapp('mtimes',L2fun,L2type,L2fcnstr,w);
       end
   end
   
          
   % Create next column of V and H
   H(1:k,k) = V(:,1:k)' * w;
   w = w - V(:,1:k)*H(1:k,k);

   H(k+1,k) = norm(w);
   V(:,k+1) = w / H(k+1,k);

   % Initialize right hand side of least-squares system
   rhs = zeros(k+1,1);
   rhs(1) = norm(r);

   % Solve least squares system; Calculate residual norm
   Hc = H(1:k+1,1:k);
   y = Hc \ rhs;
   res = rhs - Hc * y;
   resvec(k) = norm(res);
   
   % check for early convergence
   if resvec(k) < tol
        break;
   end
   
end

% Calculate solution and residual
resvec = resvec(1:k);

% obtain update on solution
xupd = V(:,1:k)*y;
if(existR1) % first right preconditioning
    if strcmp(R1type,'matrix')
        xupd = R1 \ xupd;
    else
        xupd = iterapp('mtimes',R1fun,R1type,R1fcnstr,xupd);
    end
end
if(existR2) % second right preconditioning
    if strcmp(R2type,'matrix')
        xupd = R2 \ xupd;
    else
        xupd = iterapp('mtimes',R2fun,R2type,R2fcnstr, xupd);
    end
end   
% apply update x = x + V(:,1:k) * y;
x = x + xupd;
r = V(:,1:k+1) * res;
