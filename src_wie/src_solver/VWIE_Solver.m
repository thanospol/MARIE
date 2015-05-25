function [ZP,Jcp,Jbp] = VWIE_Solver(RHBM,WCOIL,freq,tol,coup,gpu_flag,factorized,tolfactor,blockSize)
%%    Solver Routine for the VIE + WIE formulation
% _________________________________________________________________________
%
%   Def.: Fucntion to solve the VIE and WIE coupled system
%         Iterative approach with initial guess
%
%
%       Generates the port Z parameters of the coupled system
%       | Zmom   Zbc.' | | Jt |   | F |
%       |              | |    | = |   | Vin
%       | Zbc    Zbb   | | Jb |   | 0 |
%
%
%       Refs:
%           Villena, Polimeridis et al, MRM 2014
%
%       Subfunctions:
%           mv_AN
%           Assembly_WCOUP
%           Assembly_WIE
%           VWIE_op
%           pgmres, bicgstab
%
% _________________________________________________________________________
%
%% INPUT
%       RHBM  structure
%           r - mapping of the internal edge number to dof number
%           epsilon_r - voxel dielectric
%           sigma_e  - voxel conductivity
%           rho  - voxel density
%           ...
%       WCOIL structure
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           Dwire - diameter of wire
%           Rhocoil - resistivity of material
%           port - port definition
%       freq - frequency
%       tol - tolerance
%       coup - for the coupling, >0 explicitly form coupling matrix
%       gpu_flag - for GPU accelerated VIE operations
%
%
%% OUTPUT
%       ZP - Zparameters
%       Jcp - currents in the coil due to unitary port excitation
%       Jbp - volumetric currents in the body due to unitary port excitation
%
%       The currents in the coil due to Vp port voltages are Jc = Jcp*Vp
%       The currents in the body due to Vp port voltages are Jb = Jbp*Vps
%
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________

if(nargin < 3 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 4 || isempty(tol))
   tol = 1e-3;
else
    if(tol < 1e-15) || (tol >= 1)
        tol = 1e-3;
    end
end
if(nargin < 5 || isempty(coup))
   coup = 1;
end
if(nargin < 6 || isempty(gpu_flag))
   gpu_flag = 1;
end
if(nargin < 7 || isempty(factorized))
   factorized = 0;
end
if(nargin < 8 || isempty(tolfactor))
   tolfactor = 1e-5;
else
    if(tolfactor >= 1)
        tolfactor = 1e-5;
    end
end
if(nargin < 9 || isempty(blockSize))
   blockSize = 50;
else
    if(blockSize < 1)
        blockSize = 50;
    end
end

tini = tic;

% -------------------------------------------------------------------------
% Pre-proccess the RHBM and prepare for the VIE solve
% -------------------------------------------------------------------------

tic_i = tic;

[fVIEsolver,Gram,tau3,Scoord,dev] = fVIEsolver_Assembly(RHBM,freq,gpu_flag);

tVIE = toc(tic_i);


% -------------------------------------------------------------------------
% Generate the MoM Wire system
% -------------------------------------------------------------------------

tic_i = tic;

[Zcc,F] = Assembly_WIE(WCOIL,freq);

% coil coordinates are centers of segments
Ccoil = (WCOIL.Pcoil + WCOIL.Ncoil)/2; 

% segments 3D Cartesian components
Dcoil = WCOIL.Pcoil - WCOIL.Ncoil;

% get dimensions
Nc = size(Ccoil,1);
Nports = size(F,2);

tWIE = toc(tic_i);

% -------------------------------------------------------------------------
% Generate the MRGF perturbation
% -------------------------------------------------------------------------

tic_i = tic;

if ~isempty(RHBM.Dcoord)
    
    No = size(RHBM.Dcoord,1); % number of DEIM voxels
    
    % compute coupling of coil to DEIM points
    [Zbc] = Assembly_WCOUP(RHBM.Dcoord,Ccoil,Dcoil,freq);
    % recall: to transform into voxel fields we have to multiply by dV
    Zbc = reshape(Zbc,No*3,Nc); % reshape for MVP
    
    % apply the perturbation to the orginal matrix
    Zpp = Zcc - (Zbc).'*(Gram*RHBM.Um*(RHBM.Sm*(RHBM.Vm'*Zbc))); % use MRGF compressed matrices
    clear Zbc;
else
    Zpp = Zcc;    
end

tMRGF = toc(tic_i);



% -------------------------------------------------------------------------
% Factorize the coil matrix
% -------------------------------------------------------------------------


tic_i = tic;

[Lp,Up] = lu(Zpp); % keep the perturbation in factorized form
clear Zpp;

tLU = toc(tic_i);


% -------------------------------------------------------------------------
% Assemble the coupling (from coil to body)
% -------------------------------------------------------------------------


tic_i = tic;

No = size(Scoord,1);


% MAC does not allows to call memory
if ~ismac
    % check memory available and size of the coupling
    [~,MMv] = memory;
    % if memory available is smaller than 3 times the size of the coupling
    % use on the fly generation of coupling
    CoupMemEst = 1.5*3*No*Nc*16;
    if (MMv.PhysicalMemory.Available < CoupMemEst)
        factorized = 1;
        CoupMemEst = 3*No*16*Nc/1.5;
    end
    if (MMv.PhysicalMemory.Available < CoupMemEst)
        coup = -abs(coup);
        fprintf(1,'\n Available memory ~%.3gGB, coupling estimation ~%.3gGB:\n coupling will be computed on the fly,, procedure will be slower\n ', MMv.PhysicalMemory.Available/1e9, 3*No*Nc*16/1e9);
    else
        fprintf(1,'\n Available memory ~%.3gGB, coupling estimation ~%.3gGB:\n coupling will be pre-computed and stored\n ', MMv.PhysicalMemory.Available/1e9, 3*No*Nc*16/1e9);
    end
end

if (coup < 0)
    ZbcL = []; % use on-the-fly option inside MVP function
    ZbcR = [];
else
    
    try
        % generate the coupling matrix explicitly
        if (factorized) % use factorized approach
            
            blockSize = min(blockSize,Nc);
            [ZbcL,ZbcR] = Assembly_WCOUP_BlockFactor(Scoord,Ccoil,Dcoil,freq,tolfactor,blockSize);
            
        else
            
            [ZbcR] = Assembly_WCOUP(Scoord,Ccoil,Dcoil,freq);
            % recall: to transform into voxel fields we have to multiply by dV
            ZbcR = reshape(ZbcR,No*3,Nc); % reshape for MVP
            ZbcL = [];
        end  
    catch
       ZbcL = []; % use on-the-fly option inside MVP function
       ZbcR = []; 
    end
    
end

tCoup = toc(tic_i);

% MAC does not allows to call memory
if ~ismac
    infoCoup = whos('ZbcL');
    memCoup = infoCoup.bytes;
    infoCoup = whos('ZbcR');
    memCoup = memCoup + infoCoup.bytes;
    if (factorized)
        fprintf(1, '\n Coupling computed, Factorized size %d x %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, size(ZbcL,2), Nc, memCoup/(1024*1024*1000));
    else
        fprintf(1, '\n Coupling computed, Matrix size %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, Nc, memCoup/(1024*1024*1000));
    end
else
    if ~isempty(ZbcR)
        if (factorized)
            memCoup = (3*No*size(ZbcL,2) + size(ZbcL,2)*Nc)*16;
            fprintf(1, '\n Coupling computed, Factorized size %d x %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, size(ZbcL,2), Nc, memCoup/(1024*1024*1000));
        else
            memCoup = 3*No*Nc*16;
            fprintf(1, '\n Coupling computed, Matrix size %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, Nc, memCoup/(1024*1024*1000));
        end
    else
        fprintf(1,'\n Coupling will be computed on the fly, procedure will be slower\n ');
    end
end

    

% -------------------------------------------------------------------------
% Solve for the Y parameters of the system
% -------------------------------------------------------------------------
%
%       The coupled system is
%       | Zcc    Zbc.' | | Jt |   | F |
%       |              | |    | = |   | Vin
%       | Zbc    Zbb   | | Jb |   | 0 |
%
%       Needs to form
%       [YP] = F.' inv( Zcc  -  Zbc.' * inv(Zbb) * Zbc ) * F
%
%
%       We use an iterative method to solve
%       inv( Zcc  -  Zbc.' * inv(Zbb) * Zbc ) * F
%


% -------------------------------------------------------------------------
% generate the rhs and apply the batch of external solves
%       Jc = inv(Zcc - Zbc*inv(Zbb)*Zbc.')*F;

fprintf(1,'\n Applying the two-level nested iterative approach at freq %.3f MHz\n ', freq/1e6);

tic_i = tic;


% -------------------------------------------------------------------------
% Prepare the function with the internal MVP and iterative solver
fIB = @(Jc)VWIE_op(Jc,Zcc,fVIEsolver,ZbcL,ZbcR,freq,tau3,Gram,Scoord,Ccoil,Dcoil,RHBM,max(tol,1e-3),[],[]);
% [Vout] = VWIE_op(Jc,Zcc,fVIE,ZbcL,ZbcR,Scoord,Ccoord,Dcoil,tol,Jb,idxC)

% note, for the coupling use at most a tolerance of 1e-3, not need more 

% allocate solution
Jcp = zeros(Nc,Nports);
Jbp = zeros(3*No,Nports);

maxit = 1000;
restart = 50;

for ii = 1:Nports
    
    tic_iter = tic;
    
    % select current rhs
    rhs = F(:,ii);
    
    % we apply GMRES with restart and (right) preconditioning 
    [Jcp(:,ii),flag,~,~,resvec] = pgmres(@(J)fIB(J), rhs, restart, tol, maxit, [], [], Lp, Up);

    % other options for the solver
    % [Jcp(:,ii),flag,~,~,resvec] = pgmres(@(J)fIB(J), rhs, restart, tol, maxit);
    % [Jcp(:,ii),flag,~,~,resvec] = bicgstab(@(J)fIB(J), rhs, tol*1e-2, maxit, Lc, Uc);
    % [Jcp(:,ii),flag,~,~,resvec] = bicgstab(@(J)fIB(J), rhs, tol, maxit)
    
    
    % prepare right hand side for JVIE
    if ~isempty(ZbcR) % coupling is precomputed
        if isempty(ZbcL)
            vec_b = ZbcR*Jcp(:,ii);
        else
            vec_b = ZbcL*(ZbcR*Jcp(:,ii));
        end
    else % on the fly generation
        vec_b = WVIE_E_Coil2Scat_PM(Scoord,Ccoil,Dcoil,freq,Jcp(:,ii));
    end
    
    if isempty(RHBM.M)
        Jb = [];
    else
        % get the initial guess from MRGFs
        Jb = RHBM.M*(RHBM.X*(RHBM.P.'*vec_b)); % project field on DEIM points and scale basis
    end
    
    % get the RHS currents from fields
    vec_b = tau3.*(Gram*vec_b);
    
    % Solve the JVIE to get currents in the body
    Jbp(:,ii) = fVIEsolver(vec_b, tol, Jb);
    % [Jbp(:,ii),~,~,~,~] = bicgstab(@(J)fVIE(J,'notransp'), vec_b, tol, maxit, [], [], Jb);
    clear vec_b;
    
    % compute the residual for checking
    
    if ~isempty(ZbcR) % coupling is precomputed
        if isempty(ZbcL)
            residual_ = rhs - ( Zcc*Jcp(:,ii) - ZbcR.'*(Gram*Jbp(:,ii)) );
        else
            residual_ = rhs - ( Zcc*Jcp(:,ii) - ZbcR.'*(ZbcL.'*(Gram*Jbp(:,ii))) );
        end
    else % on the fly generation
        % compute the contribution to all dofs
        vec_c = WVIE_E_Scat2Coil_PM(Scoord,Ccoil,Dcoil,freq,Gram*Jbp(:,ii));
        % add the contribution to the internal variables
        residual_ = rhs - ( Zcc*Jcp(:,ii) - vec_c ); 
        clear vec_c;
    end
     
   fprintf(1,'\n Nested iterative method for RHS %d done\n  Elapsed time  = %.2f [sec], for %d external iterations (flag %d, residue %g, relative %g) \n' , ii,toc(tic_iter),length(resvec),flag,norm(residual_),norm(residual_)/norm(rhs));
   
    
end

tSolve = toc(tic_i);


% -------------------------------------------------------------------------
% form
%       YP = F.'*Jcp;
%       [ZP] = inv(YP);

% Generate the Z parameter matrix
ZP = inv(F.'*Jcp);
ZP = (ZP + ZP.')/2;

% -------------------------------------------------------------------------
% Clean and Report
% -------------------------------------------------------------------------

if (exist('dev', 'var') && (~isempty(dev))) % clean the GPU
    reset(dev); 
end

toverall = toc(tini);
fprintf(1,'\n VIE assembly:             \tElapsed time  = %.2f [sec]\n' , tVIE);
fprintf(1,'\n WIE assembly:             \tElapsed time  = %.2f [sec]\n' , tWIE);
fprintf(1,'\n MRGF assembly:            \tElapsed time  = %.2f [sec]' , tMRGF);
fprintf(1,'\n LU of Zcc:                \tElapsed time  = %.2f [sec]\n' , tLU);
fprintf(1,'\n Coupling assembly:        \tElapsed time  = %.2f [sec]\n' , tCoup);
fprintf(1,'\n Solve for %d rhs:         \tElapsed time  = %.2f [sec]\n', Nports, tSolve);
fprintf(1,'\n --------------------------------------------------------------');
fprintf(1,'\n Overall Iterative Solver: \tElapsed time  = %.2f [sec]\n', toverall);

