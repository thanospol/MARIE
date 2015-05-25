function [Jout,Sout,Eout,Bout,Gsar,Pabs,fN,fK] = Scat_Solver(Einc,Hinc,freq,r,epsilon_r,sigma_e,rho,tol,fN,fK,form,solver,precond,max_it,inner_it,outer_it,ritz,gpu_flag,logfile)
%%    Simple driver to completely solve EM problem
% _________________________________________________________________________
%
%       Completely solves the EM problem
%       Applies the JVIE formulation
%
% _________________________________________________________________________
%
%% INPUT
%   Einc       Incident Electric Field (LxMxNx3)
%   Hinc       Incident Magnetic Field (LxMxNx3)
%   freq       Frequency (in Hz)
%   r          4D (LxMxNx3) array with domain voxelized grid coordinates
%   epsilon_r  Relative epsilon (LxMxN)
%   sigma_e    Electric conductivity (LxMxN)
%
%
%% OPTIONAL INPUT
%   rho         Density (LxMxN)
%   tol         Relative tolerance for the method (DEFAULT 1e-3)
%   fN          N operator circulant
%   fK          K operator circulant
%   form        Choice of the formulation
%                   1 for JVIE I formulation
%                   2 for JVIE II formulation (DEFAULT)
%   solver      Choice of an Iterative Solver
%                   'G' for GMRES
%                   'D' for GMRES with Deflated Restart
%                   'B' for BiCGStab (DEFAULT)
%                   'F' for FFTW JVIEII BiCGStab (overrides form and gpu)
%                   'M' automatic selection
%   max_it      Number of maximum iterations (DEFAULT 1000)
%   inner_it    Number of internal iterations in restarted methods (DEFAULT 50)
%   outer_it    Number of external iterations in restarted methods (DEFAULT 200)
%   ritz        Number of vectors kept for deflated restart (DEFAULT 5)
%   precond     Use preconditioner
%                   'L' left preconditioning
%                   'R' right preconditioning
%                   otherwise (DEFAULT), no preconditioner
%   gpu_flag    if 0, forces to not use GPU (DEFAULT)
%               otherwise selects the GPU device with number gpu_flag (if possible)
%   logfile     file to print the data, if empty to the standard output (DEFAULT)
%
%
%% OUTPUT
%   Jout        Solution electric current (LxMxNx3)
%   Sout        Local SAR (LxMxN)
%   Eout        Solution electric field (LxMxNx3)
%   Bout        Solution magnetic field (LxMxNx3)
%   Gsar        Global SAR (number)
%   Pabs        Absorbed power (number)
%   fN          N operator circulant
%   fK          K operator circulant
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________



% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 6 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 7 || isempty(rho))
   rho = epsilon_r./epsilon_r; % all ones
end
if(nargin < 8 || isempty(tol))
   tol = 1e-3;
else
    if(tol < 1e-15) || (tol >= 1)
        tol = 1e-3;
    end
end
if(nargin < 9 || isempty(fN))
   fN = [];
end
if(nargin < 10 || isempty(fK))
   fK = [];
end
if(nargin < 11 || isempty(form))
   form = 2;
end
if(nargin < 12 || isempty(solver))
   solver = 'M';
else
    if ((solver ~= 'G') && (solver ~= 'D')) && ((solver ~= 'B') && (solver ~= 'F'))
        solver = 'M';
    end
end
if(nargin < 13 || isempty(precond))
   precond = 0;
   else
    if (precond ~= 'L') && (precond ~= 'R')
        precond = 0;
    end
end
if(nargin < 14 || isempty(max_it))
   max_it = 1000;
end
if(nargin < 15 || isempty(inner_it) || (inner_it < 1))
   inner_it = 50;
end
if(nargin < 16 || isempty(outer_it) || (outer_it < 1))
   outer_it = 200;
end
if(nargin < 17 || isempty(ritz) || (ritz < 1))
   ritz = 5;
end
if(nargin < 18 || isempty(gpu_flag))
   gpu_flag = 1;
end

if(nargin < 19 || isempty(logfile))
   fid = 1;
   logfile = [];
else
    fid = fopen(logfile, 'a');
    if (fid == -1)
        fid = 1; % standard output
    end
end

tini = tic;

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% properties
e_r = epsilon_r - 1j*sigma_e/(eo*omega);

% dimensions
dx = r(2,1,1,1) - r(1,1,1,1);
Gram = dx^3;

% -------------------------------------------------------------------------
%                  Generate circulants
% -------------------------------------------------------------------------

% compute the circulants
if (isempty(fN))
    [fN] = getOPERATORS(r,freq,'N');
end
if (isempty(fK))
    [fK] = getOPERATORS(r,freq,'K');
end


% -------------------------------------------------------------------------
%                 Solve the system 
% -------------------------------------------------------------------------

t1 = tic;

[Jout,~,~,~,~,~] = JVIE_Solver(Einc,e_r,r,freq,fN,tol,form,solver,precond,max_it,inner_it,outer_it,ritz,gpu_flag,logfile);

fprintf(1,'\n Volumetric currents generated:  \tElapsed time  = %.2f [sec]' , toc(t1));

% -------------------------------------------------------------------------
%                   Compute SAR
% -------------------------------------------------------------------------

t1 = tic;

[L,M,N,~] = size(r);
denom = 2*sigma_e.*rho;
idxsar = find(abs(denom(:))>1e-10);
idxzero = setdiff(1:L*M*N,idxsar);
Sout = Jout.*conj(Jout);
Sout = sum(Sout,4);
Sout(idxzero) = 0;
Sout(idxsar) = Sout(idxsar)./denom(idxsar);


% Sout= Jout.*conj(Jout);
% Sout = sum(Sout,4);
% idxsar = find(abs(sigma_e(:)>1e-10));
% Sout(idxsar) = Sout(idxsar)./(sigma_e(idxsar));
% idxsar = find(abs(rho(:)>1e-12));
% Sout(idxsar) = Sout(idxsar)./(2*rho(idxsar));

fprintf(1,'\n SAR computed:                   \tElapsed time  = %.2f [sec]' , toc(t1));

% -------------------------------------------------------------------------
%                   Compute E field
% -------------------------------------------------------------------------

t1 = tic;

[Eout] = E_field_Nop(Jout,fN,Gram,freq,Einc);

fprintf(1,'\n E fields computed:              \tElapsed time  = %.2f [sec]' , toc(t1));


% -------------------------------------------------------------------------
%                   Compute H field
% -------------------------------------------------------------------------

t1 = tic;

[Bout] = H_field_Kop(Jout,fK,Gram,Hinc);
Bout = 4*pi*1e-7*Bout;

fprintf(1,'\n B fields computed:              \tElapsed time  = %.2f [sec]' , toc(t1));

% -------------------------------------------------------------------------
%                   Compute Global SAR and Pabs
% -------------------------------------------------------------------------

t1 = tic;

Pabs = Eout.*conj(Jout);
Pabs = sum(Pabs,4);
Pabs = 0.5*real(sum(Pabs(:)))*Gram;
Gsar = sum(Sout(:))*Gram;

fprintf(1,'\n Power and GSAR computed:         \tElapsed time  = %.2f [sec]' , toc(t1));

% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

toverall = toc(tini);
fprintf(1,'\n --------------------------------------------------------------');
fprintf(1,'\n Overall Scattering Problem:      \tElapsed time  = %.2f [sec]', toverall);
fprintf(1,'\n --------------------------------------------------------------');
fprintf(1,'\n --------------------------------------------------------------\n');
