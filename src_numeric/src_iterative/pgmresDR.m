
% modified Preconditioned GMRES DR
%   Version by J. Fernandez Villena
%              Computational Prototyping Group, RLE at MIT
%   original code from ML Parks, Sandia Labs. 
%   Can use double Left and/or Righ Preconditioners
%   Can use function handle for A and preconditioners
%   Uses deflated restart
%
% INPUT:  A         N-by-N matrix or function
%         b         rhs
%         restart   inner iterations before restart
%         tol       relative tolerance for the residue
%         maxit     maximum number of outer iterations (total is restart*maxit)
%         k         size of the subspace kept at each restart for deflation
%         L1        first left preconditioner for A, so that L1\A ~ EYE
%         L2        second left preconditioner for A, so that L2\(L1\A) ~ EYE
%         R1        first right preconditioner for A, so that A/R1 ~ EYE
%         R2        second right preconditioner for A, so that (A/R1)/R2 ~ EYE
%         x0        initial guess
%
% OUTPUT: x         solution vector
%         flag      0 if converged, 1 if not
%         relres    final relative residue: norm(b - A*x)/norm(b) 
%         iter      vector with [current internal iterations, current external iterations]
%         resvec    vector containing norm of relative residual at each iteration of GMRES
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
%           [x1,flag,relres,iter,resvec] = pgmresdr(A,b,30,1e-6,10,15);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 1: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 2: left preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmresdr(A,b,30,1e-6,10,15,L,U);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 2: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 3: right preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmresdr(A,b,30,1e-6,10,15,[],[],L,U);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 3: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
% 
%           % Approximation 4: left and right preconditioning
%           tic_gmres    = tic;
%           [x1,flag,relres,iter,resvec] = pgmresdr(A,b,30,1e-6,10,15,L,[],U,[]);
%           Time_gmres = toc(tic_gmres);
%           r = b - A*x;
%           fprintf(1,'Approximation 4: Time gmres.  = %4.0d [sec], flag %d, relres %g, %3.0d iterations, residue %g \n' ,round(Time_gmres),flag,relres,length(resvec),norm(r));
%           fprintf(1,'                 Absolute Error %g, Relative Error %g \n\n' , max(abs(xtrue - x)), max(abs(xtrue -x)./abs(xtrue)));
%
%
%
% SUBFUNCTIONS: pgmresDR_iter.m, getHarmVecs1.m, getHarmVecs2.m, iterchk, iterapp
%
%

function [x,flag,relres,iter,resvec] = pgmresDR(A,b,restart,tol,maxit,ritz,L1,L2,R1,R2,x0)

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
if(nargin < 6 || isempty(ritz))
    ritz = min(10,round(restart/2));
else
    ritz = min(ritz,restart-1);
end
if(nargin < 7 || isempty(L1))
    existL1 = 0; L1 = [];
else
    existL1 = 1;
    [L1type,L1fun,L1fcnstr] = iterchk(L1);
end
if(nargin < 8 || isempty(L2))
    existL2 = 0; L2 = [];
else
    existL2 = 1;
    [L2type,L2fun,L2fcnstr] = iterchk(L2);
end
if(nargin < 9 || isempty(R1))
    existR1 = 0; R1 = [];
else
    existR1 = 1;
    [R1type,R1fun,R1fcnstr] = iterchk(R1);
end
if(nargin < 10 || isempty(R2))
    existR2 = 0; R2 = [];
else
    existR2 = 1;
    [R2type,R2fun,R2fcnstr] = iterchk(R2);
end
if(nargin < 11 || isempty(x0))
    x0 = zeros(size(b));
end


% -------------------------------------------------------------------------
% settings

% Determine whether A is a matrix or a function.
[atype,afun,afcnstr] = iterchk(A);

% Calculate initial preconditioned residual.
% r = b - A*x0;
r = b - iterapp('mtimes',afun,atype,afcnstr,x0);
% if there is Left preconditoning, the residual is L2\(L1\(B-A*X))
bp = b;
if existL1
    if strcmp(L1type,'matrix')
        r = L1\r;
        bp = L1\bp;
    else
        r = iterapp('mtimes',L1fun,L1type,L1fcnstr,r);
        bp = iterapp('mtimes',L1fun,L1type,L1fcnstr,bp);
    end
end
if existL2
    if strcmp(L2type,'matrix')
        r = L2\r;
        bp = L2\bp;
    else
        r = iterapp('mtimes',L2fun,L2type,L2fcnstr,r);
        bp = iterapp('mtimes',L1fun,L1type,L2fcnstr,bp);
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
C = [];


% -------------------------------------------------------------------------
% first call to an initial gmres to perform restart iterations

% Call the gmres to perform restart iterations
[V,H,~,p,resvec_inner] = pgmresDR_iter(A,r,restart,L1,L2,R1,R2,C,tol*bnorm);

% store vector with residues
resvec(it+1:it+p) = resvec_inner/bnorm;
it = it + p;

% obtain update on solution and residual
y = H\(V'*r);
x = V(:,1:p)*y;
r = r - V*(H*y);

% Generate the Ritz vectors for deflation ----------------------------
if (resvec(it) > tol)

    % Find the ritz smallest harmonic Ritz vectors.
    if ritz > p-1
        P = getHarmVecs1(p,p-1,H);
    else
        P = getHarmVecs1(p,ritz,H);
    end
    
    % Form the subspace to recycle: U
    U = V(:,1:p)*P;
    
    % Form orthonormalized C
    [C,R] = qr(H*P,0);
    C = V*C;
    
    % adjust U accordingly so that C = A*U
    U = U/R;
    
    % clean stuff
    clear V; clear H; clear P;
    
end


% -------------------------------------------------------------------------
% Restarted loop until convergence or maximum reached
while(resvec(it) > tol) && (outit < maxit)
  
    % increase external loop counter
    outit = outit+1;
        
    % Call the gmres to perform restart iterations ------------------------
    [V,H,B,p,resvec_inner] = pgmresDR_iter(A,r,restart,L1,L2,R1,R2,C,tol*bnorm);
        
    % store vector with residues
    resvec(it+1:it+p) = resvec_inner/bnorm;
    it = it + p;
           
   % Rescale U; 
   % Store inverses of the norms of columns of U in diagonal matrix D
   d = zeros(ritz,1);
   for i = 1:ritz
      d(i) = norm(U(:,i));
   end
   D = diag(1./d);
   U = U*D;

   % Form large H
   H2 = [D, B; zeros(size(H,1),size(D,2)), H];

   % obtain update on solution and residual
   y = H2\([C V]'*r);
   x = x + [U V(:,1:p)]*y;
   r = r - [C V]*(H2*y);

   
   if (resvec(it) > tol)
       % Calculate Harmonic Ritz vectors.
       P = getHarmVecs2(p+ritz,ritz,H2,V,U,C);
       
       % Form new U and C.
       U = [U V(:,1:p)]*P;
       
       % Form orthonormalized C and adjust U accordingly so that C = A*U
       [Q,R] = qr(H2*P,0);
       C = [C V] * Q;
       U = U/R;
   end
   
   
end

% Calculate final solution,
if(existR2) % second right preconditioning
    if strcmp(R2type,'matrix')
        x = R2\x;
    else
        x = iterapp('mtimes',R2fun,R2type,R2fcnstr, x);
    end
end
if(existR1) % first right preconditioning
    if strcmp(R1type,'matrix')
        x = R1 \ x;
    else
        x = iterapp('mtimes',R1fun,R1type,R1fcnstr,x);
    end
end
x = x0 + x;

% Calculate residual, and iteration count
resvec = resvec(1:it);
iter = [p outit];

% Calculate relative residual.
relres = norm(b - iterapp('mtimes',afun,atype,afcnstr,x))/norm(b);
if (relres <= tol ) % converged
    flag = 0;
else
    flag = 1;
end

