function [fVIE,Gram,tau3,Scoord,dev] = fVIE_Assembly(RHBM,freq,gpu_flag)
%%    Assembly Routine to generate the fVIE handle
% _________________________________________________________________________
%
%   Prepares the RHBM data for solving
%   Checks if there is already circulant and dimensions
% _________________________________________________________________________
%
%% INPUT
%       RHBM - structure with data
%       freq - frequency
%       gpu_flag - for GPU accelerated VIE operations
%
%
%% OUTPUT
%       fVIE - VIE MVP handle function
%       Gram - voxel volume
%       tau3 - multiplier to get VIE rhs
%       Scoord - coordinates of non-air voxels
%       dev - device if GPU is used, empty otherwise
%
%
% -------------------------------------------------------------------------
%
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

fid = 1;

% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*freq;

% -------------------------------------------------------------------------
% Pre-proccess the RHBM and sample points
% -------------------------------------------------------------------------

% properties
e_r = RHBM.epsilon_r - 1j*RHBM.sigma_e/(eo*omega);

% voxel volume as Gram matrix
dx = RHBM.r(2,1,1,1) - RHBM.r(1,1,1,1);
Gram = dx^3;

% Domain dimensions
[L,M,N,~] = size(e_r);
ND = L*M*N; % number of variables in the system

% get the positions of the non-air voxels in the 3D grid
idxS = find(abs(e_r(:)-1) > 1e-12); % these are the indexes of the non-air voxel positions
idxS3 = [idxS; ND+idxS; 2*ND+idxS]; % the vector of non-air positions for 3 Cartesian components

% x,y,z coordinates of the domain
xd = RHBM.r(:,:,:,1);
yd = RHBM.r(:,:,:,2);
zd = RHBM.r(:,:,:,3);

%   S are the coordinates of the RHBM points
Scoord = [xd(idxS), yd(idxS), zd(idxS)];
clear xd; clear yd; clear zd;

% -------------------------------------------------------------------------
%      Check circulant sizes of interesting subdomain
% -------------------------------------------------------------------------

% define subdomain of interest
xd = RHBM.r(:,:,:,1);
xmin = min(xd(idxS))-dx; xmax = max(xd(idxS))+dx;
yd = RHBM.r(:,:,:,2);
ymin = min(yd(idxS))-dx; ymax = max(yd(idxS))+dx;
zd = RHBM.r(:,:,:,3);
zmin = min(zd(idxS))-dx; zmax = max(zd(idxS))+dx;

% cut grid to subdomain dimensions
idxX = find( (RHBM.r(:,1,1,1) >= xmin) & (RHBM.r(:,1,1,1) <= xmax));
idxY = find( (RHBM.r(1,:,1,2) >= ymin) & (RHBM.r(1,:,1,2) <= ymax));
idxZ = find( (RHBM.r(1,1,:,3) >= zmin) & (RHBM.r(1,1,:,3) <= zmax));

% -------------------------------------------------------------------------
%                  Create new grid
% -------------------------------------------------------------------------

% truncate e_r to new reduced grid
e_rnew = e_r(idxX,idxY,idxZ);

% get dimensions
[Lf,Mf,Nf] = size(e_rnew);

% -------------------------------------------------------------------------
%      Check if frequency is the same, or original domain is too large   
% -------------------------------------------------------------------------

if ((Lf*Mf*Nf < 0.7*ND) || isempty(RHBM.fN) || (RHBM.freqfN ~= freq))
    % we need to compute a new circulant and set new data
    
    % crop domain
    rnew = RHBM.r(idxX,idxY,idxZ,:);
    
    % compute the circulant
    [fN] = getOPERATORS(rnew,freq,'N');
    
    % re-assign data with cropped domain data
    e_r = e_rnew;
    [L,M,N,~] = size(e_r);
    ND = L*M*N; % number of variables in the system
    
    % prepare domain
    idxS = find(abs(e_r(:)-1)> 1e-12);
    idxS3 = [idxS; ND+idxS; 2*ND+idxS]; % the vector of non-air positions for 3 Cartesian components

else
    % use the one available
    fN = RHBM.fN;    

end

% -------------------------------------------------------------------------
% Prepare the vars and functions for the VIE solve (JVIE II form)
% -------------------------------------------------------------------------

% Compute the relative permittivity and suceptibility for second
% formulation
Mr = e_r./e_r;
Mc = (e_r - 1.0)./e_r;

% multiplier for transforming the fields into rhs (JVIE II form)
tau = 1j*omega*eo*Mc;
tau3 = [tau(idxS); tau(idxS); tau(idxS)]; % 3 Cartesian components in vector form

% operator
transp_flag = 'notransp';

% Check the possibility of using a GPU card
infofN = whos('fN');
memestimated = 6*infofN.bytes;

fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n VIE Domain:            %dx%dx%d voxels',L,M,N);
fprintf(fid, '\n Resolution:            %.2fmm',dx*1000);
fprintf(fid, '\n # DOFS:                %d', 3*ND);
fprintf(fid, '\n # DOFS in Scatterer:   %d', length(idxS3));
fprintf(fid, '\n Operating Frequency:   %.2f MHz',freq/1e6);
fprintf(fid, '\n MVP Peak Memory:       %.3f MB (estimated) ' , memestimated/(1024*1024));
fprintf(fid, '\n');
fprintf(fid, '\n');

if (gpuDeviceCount) && (gpu_flag)  % check GPU existence and user opt
    
    if (ceil(gpu_flag) <= gpuDeviceCount)  % Try to select gpu_flag device
        dev = gpuDevice(ceil(gpu_flag));
    else
        dev = gpuDevice(1); % get the device one by default
    end
    
    % reset the device
    reset(dev);
    
    if (dev.FreeMemory > memestimated ) % verify memory of GPU and problem size
        
        % -------------------------------------------------------------
        % Apply the GPU based solver
        
        % send data to GPU memory
        fN_gpu = gpuArray(fN);
        Mr_gpu = gpuArray(Mr);
        Mc_gpu = gpuArray(Mc);
        
        % define the main function for the matrix vector product
        fVIE = @(J)mv_AN(J, fN_gpu, Mr_gpu, Mc_gpu, Gram, transp_flag, idxS3, 1);
        fprintf(fid, '\n GPU VIE op handle generated (GPU memory %.3f MB)', dev.FreeMemory/(1024*1024));
        
    end
    
end
        
% -----------------------------------------------------------------
% Not available GPU 
% Not enough GPU memory for the problem: go to CPU based solver
% or user opted for CPU based solver
if isempty(fVIE)
    
    % -------------------------------------------------------------
    % Not enough GPU memory for the problem: go to CPU based solver
    
    % define the main function for the matrix vector product
    fVIE   = @(J,transp_flag)mv_AN(J, fN, Mr, Mc, Gram, transp_flag, idxS3, 0);
    fprintf(fid, '\n CPU VIE op handle generated (No GPU available)');
    
    dev = [];

end
       
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------');

