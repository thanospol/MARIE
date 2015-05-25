
% modified GMRES with Left and Righ Preconditioners, and standard restart
%   Version by J. Fernandez Villena
%              Computational Prototyping Group, RLE at MIT
%   Can use double Left and/or Righ Preconditioners
%   Can use function handle for A and preconditioners
%
% INPUT:  A         N-by-N matrix or function
%         b         rhs
%         restart   inner iterations before restart
%         tol       relative tolerance for the residue
%         maxit     maximum number of outer iterations (total is restart*maxit)
%         L1        first left preconditioner for A, so that L1\A ~ EYE
%         L2        second left preconditioner for A, so that (L1*L2)\A ~ EYE
%         R1        first right preconditioner for A, so that A/R1 ~ EYE
%         R2        second right preconditioner for A, so that A/(R1*R2) ~ EYE
%         x0        initial guess
%
% OUTPUT: x         solution vector
%         flag      0 if converged, 1 if not
%         relres    final relative residue: norm(b - A*x)/norm(b) 
%         iter      vector with [current internal iterations, current external iterations]
%         resvec    vector containing norm of residual at each iteration of GMRES
%         
% EXAMPLE:  
%           A = sprand(1000,1000,0.05) +  spdiags(rand(1000,1),0,1000,1000);
%           b = rand(1000,1);
%           xtrue = A\b;
% 
%           [L,U] = luinc(A,0.01);
% 
%           % Approximation 1: no preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 1: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 2: left preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10,L,U);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 2: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 3: right preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10,[],[],L,U);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 3: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 4: left and right preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmres(A,b,30,1e-6,10,L,[],U,[]);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 4: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
%
%
%
% SUBFUNCTIONS: pgmres_iter.m, iterchk, iterapp
%
%

function [x,flag,relres,iter,resvec] = pgmres(A,b,restart,tol,maxit,L1,L2,R1,R2,x0)

% Initialize variables.
if(nargin < 2 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 3 || isempty(restart))
   restart = 50;
end
if(nargin < 4 || isempty(tol))
   tol = 1e-3;
end
if(nargin < 5 || isempty(maxit))
   maxit = 10;
end
if(nargin < 6 || isempty(L1))
   L1pre = 0; L1 = [];
else
   L1pre = 1;
   [L1type,L1fun,L1fcnstr] = iterchk(L1);
end
if(nargin < 7 || isempty(L2))
    L2pre = 0; L2 = [];
else
    L2pre = 1;
    [L2type,L2fun,L2fcnstr] = iterchk(L2);
end
if(nargin < 8 || isempty(R1))
   R1 = [];
end
if(nargin < 9 || isempty(R2))
   R2 = [];
end
if(nargin < 10 || isempty(x0))
   x0 = zeros(size(b));
end


% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);


% Calculate initial preconditioned residual.
% r = b - A*x0;
if (nnz(x0))
    r = b - iterapp('mtimes',afun,atype,afcnstr,x0);
else
    r = b;
end
% if there is Left preconditoning, the residual is L2\(L1\(B-A*X))
bp = b;
if L1pre
    if strcmp(L1type,'matrix')
        r = L1\r;
        bp = L1\bp;
    else
        r = iterapp('mtimes',L1fun,L1type,L1fcnstr,r);
        bp = iterapp('mtimes',L1fun,L1type,L1fcnstr,bp);
    end
end
if L2pre
    if strcmp(L2type,'matrix')
        r = L2\r;
        bp = L2\bp;
    else
        r = iterapp('mtimes',L2fun,L2type,L2fcnstr,r);
        bp = iterapp('mtimes',L2fun,L2type,L2fcnstr,bp);
    end
end

% Calculate rhs norm
bnorm = norm(bp);
clear bp;

% initialize the residue vector and other variables
resvec = zeros(restart*maxit,1);
it = 1;
outit = 0;
resvec(it) = norm(r)/bnorm;
x = x0;

% Restarted loop until convergence or maximum reached
while(resvec(it) > tol) && (outit < maxit)
  
    % increase external loop counter
    outit = outit+1;
    
   % Call the gmres to perform restart iterations
   [x,r,p,resvec_inner] = pgmres_iter(A,x,r,restart,L1,L2,R1,R2,tol*bnorm);
   resvec(it+1:it+p) = resvec_inner/bnorm;
   it = it + p;
   
    % % % Calculate relative residual.
    % % res = b - iterapp('mtimes',afun,atype,afcnstr,x);
    % % fprintf(1,'restart %d, ||r|| = %g ( resvec %g, tol %g, reltol %g)\n' ,outit, norm(res), resvec(it), tol, tol*bnorm);
    
end


% residual, and iteration count
resvec = resvec(1:it);
iter = [p outit];

% check connvergence
relres = norm(b - iterapp('mtimes',afun,atype,afcnstr,x))/norm(b);
if (relres <= tol ) % converged
    flag = 0;
else
    flag = 1;
end

