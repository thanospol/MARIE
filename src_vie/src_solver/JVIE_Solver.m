function [Je,flag,relres,iter,resvec,solvert] = JVIE_Solver(Einc,e_r,r,freq,fN,tol,form,solver,precond,max_it,inner_it,outer_it,ritz,gpu_flag,logfile)
%%    Solver Routine for the Electric JVIE formulation
% _________________________________________________________________________
%
%       Solves the VIE problem, by applying operator N
%       Generates the Electric currents J given an incident E field
%
% _________________________________________________________________________
%
%% INPUT
%   Einc        Incident Electric Field (LxMxNx3)
%   e_r         Relative electric permittivity (LxMxN)
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%   freq        Frequency (in Hz)
%   fN          N operator: FFT of circulant matrix
%
%% OPTIONAL INPUT
%   tol         Relative tolerance for the method (DEFAULT 1e-3)
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
%   gpu_flag    if 0, forces to not use GPU
%               otherwise selects the GPU device with number gpu_flag (if possible)
%               if empty DEFAULT is 1
%   logfile     file to print the data, if empty to the standard output (DEFAULT)
%
%
%% OUTPUT
%   Je          Solution electric current (LxMxNx3)
%   flag        0 if converged, 1 if not
%   relres      final relative residue: norm(b - A*x)/norm(b) 
%   iter        vector with [current internal iterations, current external iterations]
%   resvec      vector containing norm of relative residual at each iteration of GMRES
%   solvert     elapsed time by the solver, in seconds
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


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 5 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 6 || isempty(tol))
   tol = 1e-3;
else
    if(tol < 1e-15) || (tol >= 1)
        tol = 1e-3;
    end
end
if(nargin < 7 || isempty(form))
   form = 2;
end
if(nargin < 8 || isempty(solver))
   solver = 'B';
else
    if ((solver ~= 'G') && (solver ~= 'D'))
        solver = 'B';
    end
end
if(nargin < 9 || isempty(precond))
   precond = 0;
   else
    if (precond ~= 'L') && (precond ~= 'R')
        precond = 0;
    end
end
if(nargin < 10 || isempty(max_it))
   max_it = 1000;
end
if(nargin < 11 || isempty(inner_it) || (inner_it < 1))
   inner_it = 50;
end
if(nargin < 12 || isempty(outer_it) || (outer_it < 1))
   outer_it = 200;
end
if(nargin < 13 || isempty(ritz) || (ritz < 1))
   ritz = 5;
end
if(nargin < 14 || isempty(gpu_flag))
   gpu_flag = 1;
end

if(nargin < 15 || isempty(logfile))
   fid = 1;
else
    fid = fopen(logfile, 'a');
    if (fid == -1)
        fid = 1; % standard output
    end
end

tini = tic;
fprintf(1,'\n');

% -------------------------------------------------------------------------
% Prepare the domain and EM variables
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% Compute the relative permittivity and suceptibility
Mr = e_r;
Mc = e_r - 1.0;

% voxel volume as Gram matrix
dx = r(2,1,1,1) - r(1,1,1,1);
Gram = dx^3;

% Domain dimensions
[L,M,N,~] = size(r);
nD = L*M*N; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(Mc(:)) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions for 3 Cartesian components

% set the multiplier to convert E fields to the solver rhs (currents)
tau = 1j*omega*eo*Mc;
tau3 = [tau(:); tau(:); tau(:)]; % 3 Cartesian components in vector form

% get the RHS for the solver
% notice that only the non-air positions in Cartesian components are used
Vrhs = Gram.*tau3(idxS3).*Einc(idxS3);
% *** IMPORTANT *** the fields are transformed into voxel density

fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n Domain:                %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n # DOFS:                %d', 3*nD);
fprintf(fid, '\n # DOFS in Scatterer:   %d', 3*length(idxS));
fprintf(fid, '\n Operating Frequency:   %.2f MHz',freq/1e6);
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');



fAGPU = [];
fPGPU = [];
% -------------------------------------------------------------------------
% Select the formulation and Prepare the EM variables accordingly
% -------------------------------------------------------------------------

if (form == 1) 
    
    % ---------------------------------------------------------------------
    % JVIE_I formulation
    % ---------------------------------------------------------------------
    
    % Check the possibility of using a GPU card
    infofN = whos('fN');
    memestimated = 4*infofN.bytes;
    
    fprintf(fid,'\n Estimated MVP Peak Memory %.3f MB.\n' , memestimated/(1024*1024));
    
    if (gpuDeviceCount) && (gpu_flag)  % check GPU existence and user opt
        
        if (ceil(gpu_flag) <= gpuDeviceCount)  % Try to select gpu_flag device
            dev = gpuDevice(ceil(gpu_flag));
        else
            dev = gpuDevice(1); % get the device one by default
        end
        
        % reset the device
        reset(dev);
        
        if (dev.TotalMemory > memestimated ) % verify memory of GPU and problem size
            
            try
                % -------------------------------------------------------------
                % Apply the GPU based solver

                % send data to GPU memory
                fN_gpu = gpuArray(fN);
                Mr_gpu = gpuArray(Mr);
                Mc_gpu = gpuArray(Mc);
                
                % define the main function for the matrix vector product
                fAGPU   = @(J,transp_flag)mv_AN(J, fN_gpu, Mr_gpu, Mc_gpu, Gram, transp_flag, idxS3, 1);
                % define the preconditioner
                fPGPU   = @(J,transp_flag)mv_AN(J, fN_gpu, 1./Mr_gpu, (1./Mr_gpu - 1.0), Gram, transp_flag, idxS3, 1);
                
            catch
                
                fAGPU = [];
                fPGPU = [];
                
            end
            
        end
               
    end
        
    % -----------------------------------------------------------------
    % option for CPU based solver
    
    % define the main function for the matrix vector product
    fACPU   = @(J,transp_flag)mv_AN(J, fN, Mr, Mc, Gram, transp_flag, idxS3, 0);
    % define the preconditioner
    fPCPU   = @(J,transp_flag)mv_AN(J, fN, 1./Mr, (1./Mr - 1.0), Gram, transp_flag, idxS3, 0);

        
else
    
    % ---------------------------------------------------------------------
    % JVIE_II formulation (DEFAULT)
    % ---------------------------------------------------------------------
    
    % Modify the Mr, Mc, and RHS according to the formulation
    Mc(:) = Mc(:)./Mr(:);
    Vrhs(:) = Vrhs(:)./[Mr(idxS); Mr(idxS); Mr(idxS)]; % notice that Vrhs is only defined for non-air components
    Mr(:) = Mr(:)./Mr(:);

        
    % Check the possibility of using a GPU card
    infofN = whos('fN');
    memestimated = 4*infofN.bytes;
    
    fprintf(fid,'\n Estimated MVP Peak Memory %.3f MB.\n' , memestimated/(1024*1024));
    
    if (gpuDeviceCount) && (gpu_flag)  % check GPU existence and user opt
        
        if (ceil(gpu_flag) <= gpuDeviceCount)  % Try to select gpu_flag device
            dev = gpuDevice(ceil(gpu_flag));
        else
            dev = gpuDevice(1); % get the device one by default
        end
        
        % reset the device
        reset(dev);
        
        if (dev.TotalMemory > memestimated ) % verify memory of GPU and problem size
                        
            try
                % -------------------------------------------------------------
                % Apply the GPU based solver
                
                % send data to GPU memory
                fN_gpu = gpuArray(fN);
                Mr_gpu = gpuArray(Mr);
                Mc_gpu = gpuArray(Mc);
                
                % define the main function for the matrix vector product
                fAGPU   = @(J,transp_flag)mv_AN(J, fN_gpu, Mr_gpu, Mc_gpu, Gram, transp_flag, idxS3, 1);
                % define the preconditioner
                fPGPU   = @(J,transp_flag)mv_AN(J, fN_gpu, 1./Mr_gpu, (1./Mr_gpu - 1.0), Gram, transp_flag, idxS3, 1);
                
            catch
                
                fAGPU = [];
                fPGPU = [];
                
            end
            
        end
        
    end
            
    % -------------------------------------------------------------
    % for CPU based solver
    
    % define the main function for the matrix vector product
    fACPU   = @(J,transp_flag)mv_AN(J, fN, Mr, Mc, Gram, transp_flag, idxS3, 0);
    % define the preconditioner
    fPCPU   = @(J,transp_flag)mv_AN(J, fN, 1./Mr, (1./Mr - 1.0), Gram, transp_flag, idxS3, 0);
    
end
  

% -------------------------------------------------------------------------
% Select the iterative method and solve
% -------------------------------------------------------------------------
   
switch solver
    

    case 'G'
        
        % ---------------------------------------------------------------------
        % GMRES
        % ---------------------------------------------------------------------
        
        switch precond
            
            case 'R'
               
                % -------------------------------------------------------------
                % GMRES: w/ Right preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, [], [], @(J)fPCPU(J,'notransp'));
                    fprintf(fid,'\n CPU GMRES solver (%d restart), w/ Right preconditioner', inner_it);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, [], [], @(J)fPGPU(J,'notransp'));
                        fprintf(fid,'\n GPU GMRES solver (%d restart), w/ Right preconditioner', inner_it);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, [], [], @(J)fPCPU(J,'notransp'));
                        fprintf(fid,'\n CPU GMRES solver (%d restart), w/ Right preconditioner', inner_it);
                    end
                end
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
                
            case 'L'
                
                % -------------------------------------------------------------
                % GMRES: w/ Left preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, @(J)fPCPU(J,'notransp'));
                    fprintf(fid,'\n CPU GMRES solver (%d restart), w/ Left preconditioner', inner_it);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, @(J)fPGPU(J,'notransp'));
                        fprintf(fid,'\n GPU GMRES solver (%d restart), w/ Left preconditioner', inner_it);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, @(J)fPCPU(J,'notransp'));
                        fprintf(fid,'\n CPU GMRES solver (%d restart), w/ Left preconditioner', inner_it);
                    end
                end
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
                
            otherwise
            
                % -----------------------------------------------------------------
                %  GMRES: no preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it);
                    fprintf(fid,'\n CPU GMRES solver (%d restart), wo/ preconditioner', inner_it);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it);
                        fprintf(fid,'\n GPU GMRES solver (%d restart), wo/ preconditioner', inner_it);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it);
                        fprintf(fid,'\n CPU GMRES solver (%d restart), wo/ preconditioner', inner_it);
                    end
                end
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
            
        end
        
        
    case 'D'
        
        % ---------------------------------------------------------------------
        % GMRES with Deflated Restart
        % ---------------------------------------------------------------------
        
        switch precond
                
            case 'R'
                
                % -------------------------------------------------------------
                % GMRES DR: w/ Right preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, [], [], @(J)fPCPU(J,'notransp'));
                    fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), w/ Right preconditioner \n', inner_it, ritz);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, [], [], @(J)fPGPU(J,'notransp'));
                        fprintf(fid,'\n GPU GMRES solver (%d restart, deflated %d), w/ Right preconditioner \n', inner_it, ritz);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, [], [], @(J)fPCPU(J,'notransp'));
                        fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), w/ Right preconditioner \n', inner_it, ritz);
                    end
                end                
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
                
            case 'L'
                
                % -------------------------------------------------------------
                % GMRES DR: w/ Left preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, @(J)fPCPU(J,'notransp'));
                    fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), w/ Left preconditioner \n', inner_it, ritz);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, @(J)fPGPU(J,'notransp'));
                        fprintf(fid,'\n GPU GMRES solver (%d restart, deflated %d), w/ Left preconditioner \n', inner_it, ritz);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz, @(J)fPCPU(J,'notransp'));
                        fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), w/ Left preconditioner \n', inner_it, ritz);
                    end
                end
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
                
            otherwise
                  
                % -----------------------------------------------------------------
                %  GMRES DR: no preconditioning
                ta = tic;
                if isempty(fAGPU)
                    [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz);
                    fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), wo/ preconditioner \n', inner_it, ritz);
                else
                    try
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fAGPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz);
                        fprintf(fid,'\n GPU GMRES solver (%d restart, deflated %d), wo/ preconditioner \n', inner_it, ritz);
                    catch
                        [vsol, flag, relres, iter, resvec] = pgmresDR(@(J)fACPU(J,'notransp'), Vrhs, inner_it, tol, outer_it, ritz);
                        fprintf(fid,'\n CPU GMRES solver (%d restart, deflated %d), wo/ preconditioner \n', inner_it, ritz);
                    end
                end
                solvert = toc(ta);
                fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
        end
        
        
        
    otherwise
        
        % ---------------------------------------------------------------------
        % BICGSTAB
        % ---------------------------------------------------------------------
        
        if (precond)
            
            % -----------------------------------------------------------------
            % MATLAB BICGSTAB: w/ preconditioning (LEFT)
            ta = tic;
            if isempty(fAGPU)
                [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fACPU(J,'notransp'), Vrhs, tol, max_it, @(J)fPCPU(J,'notransp'));
                fprintf(fid,'\n CPU BICGSTAB solver, w/ preconditioner');
            else
                try
                    [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fAGPU(J,'notransp'), Vrhs, tol, max_it, @(J)fPGPU(J,'notransp'));
                    fprintf(fid,'\n GPU BICGSTAB solver, w/ preconditioner');
                catch
                    [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fACPU(J,'notransp'), Vrhs, tol, max_it, @(J)fPCPU(J,'notransp'));
                    fprintf(fid,'\n CPU BICGSTAB solver, w/ preconditioner');
                end
            end
            solvert = toc(ta);
            fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
            
        else
            
            % -----------------------------------------------------------------
            % MATLAB BICGSTAB: no preconditioning
            ta = tic;
            if isempty(fAGPU)
                [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fACPU(J,'notransp'), Vrhs, tol, max_it);
                fprintf(fid,'\n CPU BICGSTAB solver, wo/ preconditioner');
            else
                try
                    [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fAGPU(J,'notransp'), Vrhs, tol, max_it);
                    fprintf(fid,'\n GPU BICGSTAB solver, wo/ preconditioner');
                catch
                    [vsol, flag, relres, iter, resvec] = bicgstab(@(J)fACPU(J,'notransp'), Vrhs, tol, max_it);
                    fprintf(fid,'\n CPU BICGSTAB solver, wo/ preconditioner');
                end
            end
            solvert = toc(ta);
            fprintf(fid,' Time  = %.2f [sec] for %d iterations, flag %d, residue %g, relative residue %g \n' ,solvert,length(resvec),flag,resvec(end),relres);
            
        end
        

    
end


% -------------------------------------------------------------------------
% Assign veriables for return and clean whatever created... in GPU
% -------------------------------------------------------------------------


if (exist('dev', 'var')) % clean the GPU
    reset(dev); 
end

Je = zeros(L,M,N,3);
Je(idxS3) = vsol ; % return to global variables

tend = toc(tini);
fprintf(fid,' VIE Solve finished. Elapsed Time  = %.2f [sec] \n' ,tend);

if (fid ~= 1)
    fclose(fid);
end
