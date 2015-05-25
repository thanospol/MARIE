function [ZP,Jcp,Jbp] = VSIE_Solver(RHBM,SCOIL,freq,tol,coup,gpu_flag,factorized,tolfactor,blockSize)
%%    Solver Routine for the coupled VIE and SIE problem
% _________________________________________________________________________
%
%
%   Def.: Fucntion to solve the VIE and SIE coupled system
%         Iterative approach with preconditioner for the SIE vars
%         Can use a basis for initial guess of the VIE solution
%
%
%       Generates the port Z parameters of the coupled system
%       | Ztt    Zit.'  Zbt.' | | Jt |   | F |
%       | Zit    Zii    Zbi.' | | Ji | = | 0 | Vin
%       | Zbt    Zbi    Zbb   | | Jb |   | 0 |
%
% _________________________________________________________________________
%
%
%% INPUT
%       RHBM  structure
%           r - mapping of the internal edge number to dof number
%           epsilon_r - voxel dielectric
%           sigma_e  - voxel conductivity
%           rho  - voxel density
%           ...
%       SCOIL structure
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           ...
%       freq - frequency
%       tol - tolerance
%       coup - quadrature level for the coupling
%       gpu_flag - for GPU accelerated VIE operations
%
%
%% OUTPUT
%       ZP - Zparameters
%       Jcp - currents in the coil due to unitary port excitation
%       Jbp - volumetric currents in the body due to unitary port excitation
%
%       The currents in the coil due to Ip port currents are Jc = Jcp*Ip
%       The currents in the body due to Ip port currents are Jb = Jbp*Ip
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
   coup = 5;
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
% Generate the SIE system
% -------------------------------------------------------------------------

tic_i = tic;

[Zcc,F] = Assembly_SIE_par(SCOIL,freq);

tSIE = toc(tic_i);


% -------------------------------------------------------------------------
% Generate the MRGF perturbation
% -------------------------------------------------------------------------

tic_i = tic;

if ~(isempty(RHBM.Dcoord))
    
    No = size(RHBM.Dcoord,1); % number of DEIM voxels
    Nd = max(SCOIL.index); % number of dofs
    
    if (abs(coup) <= 1) % point matching explicit computation
        [Zbc] = Assembly_SCOUP_PM_par(RHBM.Dcoord,SCOIL.index,SCOIL.etod,SCOIL.Ct,SCOIL.Ln,SCOIL.Pn,freq);
    else
        % Quadrature        
        try % we try to use the PARALLEL MEXED Quadrature based function...
            % [Zbc] = Assembly_SCOUP_QMEX(RHBM.Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup))); % quadrature MEXED
            [Zbc] = Assembly_SCOUP_QMEX_par(RHBM.Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup))); % quadrature PARALLEL MEXED function
        catch % we used the PARALLEL NOT MEXED Quadrature based function...
            % [Zbc] = Assembly_SCOUP_Q(RHBM.Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)));  % standard Quadrature no mexed
            [Zbc] = Assembly_SCOUP_Q_par(RHBM.Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)));  % Quadrature PARALLEL no mexed
        end
    end
    Zbc = reshape(Zbc,No*3,Nd);
    
    % apply the perturbation to the orginal matrix
    Zpp = Zcc - (Zbc).'*(Gram*RHBM.Um*(RHBM.Sm*(RHBM.Vm'*Zbc))); % use MRGF compressed matrices
    clear Zbc;
else
    Zpp = Zcc;
end

tMRGF = toc(tic_i);

% -------------------------------------------------------------------------
% split according to the port information
% -------------------------------------------------------------------------

tic_i = tic;

idx = find(sum(abs(F),2));
Nports = size(F,2); % total number of external ports
Tports = length(idx); % number to terminals related to ports
Fp = F(1:Tports,1:Nports);
Fpinv = pinv(Fp); % pseudoinverse of the incidence matrix
Ztt = Zcc(1:Tports,1:Tports); % terminal subblock
Zii = Zcc(Tports+1:end,Tports+1:end); % internal subblock
Zit = Zcc(Tports+1:end,1:Tports);
clear Zcc;
Zpp = Zpp(Tports+1:end,Tports+1:end); % perturbed internal subblock
[Lpi,Upi] = lu(Zpp); % keep the perturbation in factorized form
Ni = size(Lpi,1);
clear Zpp;

tLU = toc(tic_i);

% -------------------------------------------------------------------------
% Assemble the coupling (from coil to head, left bottom block)
% -------------------------------------------------------------------------

tic_i = tic;

No = size(Scoord,1); % number of voxels
Nd = max(SCOIL.index); % number of dofs


% MAC does not allows to call memory
if ~ismac
    % check memory available and size of the coupling
    [~,MMv] = memory;
    % if memory available is smaller than 3 times the size of the coupling
    % use on the fly generation of coupling
    
    CoupMemEst = 3*No*16*Nd*1.5;
    if (MMv.PhysicalMemory.Available < CoupMemEst)
        factorized = 1;
        CoupMemEst = 3*No*16*Nd/1.5;
    end

    if (MMv.PhysicalMemory.Available < CoupMemEst)
        coup = -abs(coup);
        fprintf(1,'\n Available memory ~%.3gGB, coupling peak estimation ~%.3gGB:\n coupling will be computed on the fly, procedure will be slower\n ', MMv.PhysicalMemory.Available/1e9, CoupMemEst/1e9);
    else
        fprintf(1,'\n Available memory ~%.3gGB, coupling peak estimation ~%.3gGB:\n coupling will be pre-computed and stored\n ', MMv.PhysicalMemory.Available/1e9, CoupMemEst/1e9);
    end
end


if (coup < 0)
    ZbtR = [];
    ZbL = []; % on the fly approach
    ZbiR = [];
    idxT = 1:Tports; % set the index of the terminal currents
    idxC = Tports+1:Nd; % set the index of the internal currents
else
    
    try
        % Quadrature based function...
        
        if (factorized) % use factorized approach
            
            blockSize = min(blockSize,Nd);
            [ZbL,ZbR] = Assembly_SCOUP_Q_BlockFactor(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup),tolfactor,blockSize);
            
%             tic_ACA = tic;
%             aca_tol = tolfactor; aca_order = 1;
%             [ZbL,ZbR,~,~] = Assembly_SCOUP_Q_ACA(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup),aca_tol,aca_order);
%             tACA = toc(tic_ACA);
%             fprintf('\n Assembling with ACA');
%             fprintf('\n\tDimensions of U: %s', mat2str(size(ZbL)));
%             fprintf('\n\tDimensions of V: %s', mat2str(size(ZbR)));
%             fprintf('\n\tTime:            %g\n', tACA);
            
            % -------------------------------------------------------------------------
            % split coupling according to the port information
            ZbtR = ZbR(:,1:Tports); % body to terminals
            ZbiR = ZbR(:,Tports+1:end); % body to internal
            clear ZbR;
            % left factor is common to both cases
            
        else
            
            try % we used the PARALLEL MEXED Quadrature based function...
                % [Zbc] = Assembly_SCOUP_QMEX(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup)); % quadrature MEXED
                [Zbc] = Assembly_SCOUP_QMEX_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup)); % quadrature PARALLEL MEXED function
            catch
                if (coup <= 1) % point matching explicit computation
                    [Zbc] = Assembly_SCOUP_PM_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.Ct,SCOIL.Ln,SCOIL.Pn,freq);
                else
                    % we used the PARALLEL NOT MEXED Quadrature based function...
                    % [Zbc] = Assembly_SCOUP_Q(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup));  % standard Quadrature no mexed
                    [Zbc] = Assembly_SCOUP_Q_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,ceil(coup));  % Quadrature PARALLEL no mexed
                end
            end
            Zbc = reshape(Zbc,No*3,Nd);
            
            % -------------------------------------------------------------------------
            % split coupling according to the port information
            ZbtR = Zbc(:,1:Tports); % body to terminals
            % Zbi will be represented in factorized form Zbi = ZbiL*ZbiR
            ZbiR = Zbc(:,Tports+1:end); % body to internal
            ZbL = []; % not factorized
            clear Zbc;
            
        end
        
        % set idxT and idxC as empty: we only need it for on-the-fly
        idxT = [];
        idxC = [];
    catch
        ZbtR = [];
        ZbL = []; % on the fly approach
        ZbiR = [];
        idxT = 1:Tports; % set the index of the terminal currents
        idxC = Tports+1:Nd; % set the index of the internal currents
    end
    
end

tCoup = toc(tic_i);

% MAC does not allows to call memory
if ~ismac
    infoCoup = whos('ZbtR');
    memCoup = infoCoup.bytes;
    infoCoup = whos('ZbiR');
    memCoup = memCoup + infoCoup.bytes;
    infoCoup = whos('ZbL');
    memCoup = memCoup + infoCoup.bytes;
    if (factorized)
        fprintf(1, '\n Coupling computed, Factorized size %d x %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, size(ZbL,2), Nd, memCoup/(1024*1024*1000));
    else
        fprintf(1, '\n Coupling computed, Matrix size %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, Nd, memCoup/(1024*1024*1000));
    end
else
    if isempty(idxT)
        if (factorized)
            memCoup = (3*No*size(ZbL,2) + size(ZbL,2)*Nd)*16;
            fprintf(1, '\n Coupling computed, Factorized size %d x %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, size(ZbL,2), Nd, memCoup/(1024*1024*1000));
        else
            memCoup = 3*No*Nd*16;
            fprintf(1, '\n Coupling computed, Matrix size %d x %d,  Memory required %.3f GB (estimated) ' , 3*No, Nd, memCoup/(1024*1024*1000));
        end
    else
        fprintf(1,'\n Coupling will be computed on the fly, procedure will be slower\n ');
    end
end

% -------------------------------------------------------------------------
% Solve for the Z parameters of the system
% -------------------------------------------------------------------------
%
%       The coupled system is
%       | Ztt    Zit.'  Zbt.' | | Jt |   | F |
%       | Zit    Zii    Zbi.' | | Ji | = | 0 | Vin
%       | Zbt    Zbi    Zbb   | | Jb |   | 0 |
%
%       Needs to form
%       [ZP] = Fpinv * Ztt * Fpinv.' - Fpinv * Zta * inv(Zaa) * Zat * Fpinv.'
%
%       with
%       Fpinv = pinv(F)
%
%       Zaa = | Zii   Zbi.' |
%             | Zbi   Zbb   |
%
%       Zat = | Zit | = Zta.'
%             | Zbt |
%
%       We use a multistep iterative method
%       rhs_bt = Zbt*Fpinv.';
%       x_bt = inv(Zbb)*rhs_bt;
%       rhs_it = Zit*Fpinv.' - Zib*x_bt;
%       x_it = inv(Zii - Zib*inv(Zbb)*Zib.')*rhs_it;
%       Za = rhs_it.'*x_it;
%       [ZP] = Fpinv*Ztt*Fpinv.' -  Fpinv*Zbt.'x_bt - Za;
%


% -------------------------------------------------------------------------
% generate the rhs_bt and first batch of solves
%       rhs_bt = Zbt*Fpinv.';
%       x_bt = inv(Zbb)*rhs_bt;

fprintf(1,'\n Applying the two-level nested iterative approach at freq %.3f MHz\n ', freq/1e6);


tic_i = tic;

% generate the rhs
if isempty(idxT) % coupling is precomputed
    rhs_bt = ZbtR*Fpinv.';
    if ~isempty(ZbL)
        rhs_bt = ZbL*rhs_bt;
    end
    clear ZbtR;
else % on the fly generation
    rhs_bt = zeros(3*No,Nports);
    Jc = zeros(Nd,Nports);
    Jc(idxT,:) = Fpinv.';
    
    try % we used the PARALLEL MEXED Quadrature based function...
        for ii = 1:Nports
            rhs_bt(:,ii) = SVIE_E_Coil2Scat_QMEX_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc(:,ii)); % quadrature PARALLEL MEXED function
            % rhs_bt(:,ii) = SVIE_E_Coil2Scat_QMEX(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc(:,ii)); % quadrature MEXED
        end
    catch % we used the PARALLEL NOT MEXED Quadrature based function...
        for ii = 1:Nports
            rhs_bt(:,ii) = SVIE_E_Coil2Scat_Q_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc(:,ii)); % Quadrature PARALLEL no mexed
            % rhs_bt(:,ii) = SVIE_E_Coil2Scat_Q(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc(:,ii)); % standard Quadrature no mexed
        end
    end
    
    clear Jc;
end

% allocate solution
x_bt = zeros(3*No,Nports);

% maxit = 1000;
% restart = 50;

for ii = 1:Nports
    
    tic_iter = tic;
    
    % prepare rhs
    rhs = tau3.*(Gram*rhs_bt(:,ii));
    
    if isempty(RHBM.M)
        Jb = [];
    else
        % get the initial guess from MRGFs
        Jb = RHBM.M*(RHBM.X*(RHBM.P.'*rhs_bt(:,ii))); % project field on DEIM points and scale basis
    end
   
    % apply solver
    x_bt(:,ii) = fVIEsolver(rhs, tol, Jb);
    % [x_bt(:,ii),flag,~,~,resvec] = bicgstab(@(J)fVIE(J,'notransp'), rhs, tol, maxit, [], [], Jb);
    
    fprintf(1,'\n Body-Terminal system solved for RHS %d: Elapsed time  = %.2f [sec]' , ii,toc(tic_iter));
    
end

tBT = toc(tic_i);


% -------------------------------------------------------------------------
% generate the rhs_i and second batch of solves
%       rhs_it = Zit*Fpinv.' - Zib*x_bt;
%       x_it = inv(Zii - Zbi*inv(Zbb)*Zbi.')*rhs_it;

fprintf(1,'\n');

tic_i = tic;

% prepare right hand side
if isempty(idxC) % coupling is precomputed
    if isempty(ZbL) % not factorized
        rhs_it = Zit*Fpinv.' - ZbiR.'*(Gram*x_bt);
    else        
        rhs_it = Zit*Fpinv.' - ZbiR.'*(ZbL.'*(Gram*x_bt));
    end
    clear Zit;
else
    rhs_it = Zit*Fpinv.';
    
    try % we used the PARALLEL MEXED Quadrature based function...
        for ii = 1:Nports
            % compute the contribution to all dofs
            vec_ct = SVIE_E_Scat2Coil_QMEX_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bt(:,ii)); % quadrature PARALLEL MEXED function
            % vec_ct = SVIE_E_Scat2Coil_QMEX(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bt(:,ii)); % quadrature MEXED function
            % add the contribution to the internal variables
            rhs_it(:,ii) = rhs_it(:,ii) - vec_ct(idxC);
        end
    catch
        % we used the PARALLEL NOT MEXED Quadrature based function...
        for ii = 1:Nports
            % compute the contribution to all dofs
            vec_ct = SVIE_E_Scat2Coil_Q_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bt(:,ii)); % quadrature PARALLEL NOT MEXED function
            % vec_ct = SVIE_E_Scat2Coil_Q(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bt(:,ii)); % quadrature not mexed function
            % add the contribution to the internal variables
            rhs_it(:,ii) = rhs_it(:,ii) - vec_ct(idxC);
        end
    end
    
    clear vec_ct;
end
    


% -------------------------------------------------------------------------
% Prepare the function with the internal MVP and iterative solver
fIB = @(Jii)VSIE_op(Jii,Zii,fVIEsolver,ZbL,ZbiR,freq,tau3,Gram,Scoord,RHBM,SCOIL,max(tol,1e-3),[],idxC);

% note, for the coupling use at most a tolerance of 1e-3, not need more 

% allocate solution
x_it = zeros(Ni,Nports);
x_bi = zeros(3*No,Nports);

maxit = 1000;
restart = 50;

for ii = 1:Nports
    
    tic_iter = tic;
    rhs = rhs_it(:,ii);
    
     % we apply GMRES with restart and (right) preconditioning 
    [x_it(:,ii),flag,relres,~,resvec] = pgmres(@(J)fIB(J), rhs, restart, tol, maxit, [], [], Lpi, Upi);

    % other options for the solver
    % [x_it(:,ii),flag,relres,~,resvec] = pgmres(@(J)fIB(J), rhs, restart, tol, maxit);
    % [x_it(:,ii),flag,relres,~,resvec] = bicgstab(@(J)fIB(J), rhs, tol*1e-2, maxit, Li, Ui);
    % [x_it(:,ii),flag,relres,~,resvec] = bicgstab(@(J)fIB(J), rhs, tol, maxit)
    
    % prepare right hand side for JVIE
    if isempty(idxC) % coupling is precomputed
        if isempty(ZbL) % not factorized
            vec_bi = ZbiR*x_it(:,ii);
        else
            vec_bi = ZbL*(ZbiR*x_it(:,ii));
        end    
    else % on the fly generation
        Jc = zeros(Nd,1);
        Jc(idxC,1) = x_it(:,ii);
        
        try % we used the PARALLEL MEXED Quadrature based function...
            vec_bi = SVIE_E_Coil2Scat_QMEX_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc);% quadrature PARALLEL MEXED function
            % vec_bi = SVIE_E_Coil2Scat_QMEX(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc); % quadrature MEXED function
        catch % we used the PARALLEL NOT MEXED Quadrature based function...
            vec_bi = SVIE_E_Coil2Scat_Q_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc); % quadrature PARALLEL NOT MEXED function
            % vec_bi = SVIE_E_Coil2Scat_Q(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Jc); % quadrature not mexed function
        end
        
        clear Jc;
    end
    
    
    if isempty(RHBM.M)
        Jb = [];
    else
        % get the initial guess from MRGFs
        Jb = RHBM.M*(RHBM.X*(RHBM.P.'*vec_bi)); % project field on DEIM points and scale basis
    end
    
    % get the RHS currents from fields
    vec_bi = tau3.*(Gram*vec_bi);
    
    % solve the JVIE problem to generate the body solution
    x_bi(:,ii) = fVIEsolver(vec_bi, tol, Jb);
    % [x_bi(:,ii),~,~,~,~] = bicgstab(@(J)fVIE(J,'notransp'), vec_bi, tol, maxit, [], [], Jb);
    
    % compute the residual for checking
    
    if isempty(idxC) % coupling is precomputed
        if isempty(ZbL) % not factorized
            residual_ = rhs - ( Zii*x_it(:,ii) - ZbiR.'*(Gram*x_bi(:,ii)) );
        else
            residual_ = rhs - ( Zii*x_it(:,ii) - ZbiR.'*(ZbL.'*(Gram*x_bi(:,ii))) );
        end        
    else % on the fly generation
        
        % compute the contribution to all dofs
        try % we used the PARALLEL MEXED Quadrature based function...
            vec_ci = SVIE_E_Scat2Coil_QMEX_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bi(:,ii)); % quadrature PARALLEL MEXED function
            % vec_ci = SVIE_E_Scat2Coil_QMEX(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bi(:,ii)); % quadrature MEXED function
        catch % we used the PARALLEL NOT MEXED Quadrature based function...
            vec_ci = SVIE_E_Scat2Coil_Q_par(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bi(:,ii)); % quadrature PARALLEL NOT MEXED function
            % vec_ci = SVIE_E_Scat2Coil_Q(Scoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,abs(ceil(coup)),Gram*x_bi(:,ii)); % quadrature not mexed function
        end
        
        % add the contribution to the internal variables
        residual_ = rhs - ( Zii*x_it(:,ii) - vec_ci(idxC) );
        clear vec_ci;
    end
    
    fprintf(1,'\n Nested iterative method for RHS %d done\n  Elapsed time  = %.2f [sec], for %d external iterations (flag %d, residue %g, relative %g, gmres %g) \n' , ii,toc(tic_iter),length(resvec),flag,norm(residual_),norm(residual_)/norm(rhs),relres);
    
end

tIT = toc(tic_i);


% -------------------------------------------------------------------------
% form
%       Za = rhs_it.'*x_it;
%       [ZP] = Fpinv*Ztt*Fpinv.' -  rhs_bt.'x_bt - Za;

% Generate the Z parameter matrix
ZP = Fpinv*Ztt*Fpinv.' - rhs_bt.'*(Gram*x_bt) - rhs_it.'*x_it;
ZP = (ZP + ZP.')/2;

% -------------------------------------------------------------------------
% return currents in i edges and body due to port currents

Jbp = x_bt + x_bi; % body currents, combination of the effect of t and i
Jcp = [Fpinv.'; x_it]; % terminal currents due to Ip and internal edges currents due to Ip


% -------------------------------------------------------------------------
% Clean and Report
% -------------------------------------------------------------------------

if (exist('dev', 'var') && (~isempty(dev))) % clean the GPU
    reset(dev); 
end

toverall = toc(tini);

fprintf(1,'\n VIE assembly:             \tElapsed time  = %.2f [sec]' , tVIE);
fprintf(1,'\n SIE assembly:             \tElapsed time  = %.2f [sec]' , tSIE);
fprintf(1,'\n MRGF assembly:            \tElapsed time  = %.2f [sec]' , tMRGF);
fprintf(1,'\n LU of Zcc:                \tElapsed time  = %.2f [sec]' , tLU);
fprintf(1,'\n Coupling assembly:        \tElapsed time  = %.2f [sec]' , tCoup);
fprintf(1,'\n Terminal Solve (%d rhs):  \tElapsed time  = %.2f [sec]', Nports, tBT);
fprintf(1,'\n Internal Solve (%d rhs):  \tElapsed time  = %.2f [sec]', Nports, tIT);
fprintf(1,'\n --------------------------------------------------------------');
fprintf(1,'\n Overall Iterative Solver: \tElapsed time  = %.2f [sec]\n', toverall);

