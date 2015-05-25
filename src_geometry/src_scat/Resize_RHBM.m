function [RHBM] = Resize_RHBM(RHBM,resolution,xlimits,ylimits,zlimits,trans) 
%%    Resize and rediscretize a RHBM
% _________________________________________________________________________
%
%   Gets the original RHBM, resizes and rediscretizes the domain
%
% _________________________________________________________________________
%
%
%% INPUT
%       RHBM - initial RHBM
%       resolution - desired resolution for the discretization (min if empty)
%       limits - limits of domain (automatic if left empty)
%       trans = translation in x, y and z
%
%
%% OUTPUT
%       RHBM -  Realistic Human Body Model structure with
%           name        - name/description of model
%           r           - mapping of the internal edge number to dof number
%           epsilon_r   - voxel dielectric
%           sigma_e     - voxel conductivity
%           rho         - voxel density
%           freq        - frequency for which the data below is computed
%           fN          - N operator circulant in fft domain
%           fK          - K operator circulant in fft domain
%           Scoord      - coordinates of the DEIM points (Npx3)
%           P           - incidence matrix of the DEIM points
%           Um          - MRGF for scattering  (Um*Sm*Vm')
%           Sm          - MRGF for scattering
%           Vm          - MRGF for scattering
%           Xin         - MRGF deim interpolation matrix
%           M           - MRGF J basis
%
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


if(nargin < 2)
   resolution = [];
end
if(nargin < 3)
   xlimits = [];
end
if(nargin < 4)
   ylimits = [];
end
if(nargin < 5)
   zlimits = [];
end
if(nargin < 6)
   trans = [];
end

% -------------------------------------------------------------------------
% Open file and read subfiles
% -------------------------------------------------------------------------

% get the coordinates
x = unique(RHBM.r(:,:,:,1));
y = unique(RHBM.r(:,:,:,2));
z = unique(RHBM.r(:,:,:,3));
        
% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------

rediscretize = 0; % flag for the rediscretization

% limits of domain
if isempty(xlimits)
    xmin = min(x); xmax = max(x);
else
    xmin = min(xlimits);
    xmax = max(xlimits);
    rediscretize = 1;
end

if isempty(ylimits)
    ymin = min(y); ymax = max(y);
else
    ymin = min(ylimits);
    ymax = max(ylimits);
    rediscretize = 1;
end

if isempty(zlimits)
    zmin = min(z); zmax = max(z);
else
    zmin = min(zlimits);
    zmax = max(zlimits);
    rediscretize = 1;
end

% define discretization
if isempty(resolution)
    dx = min(x(2:end) - x(1:end-1)); % min discretization
    dy = min(y(2:end) - y(1:end-1)); % min discretization
    dz = min(z(2:end) - z(1:end-1)); % min discretization
    dx = min([dx, dy, dz]); % get the minimum
else % force resolution
    dx = resolution;
    rediscretize = 1;
end

if (rediscretize)
    
    % number of voxels in each direction
    nX = ceil((xmax - xmin)/dx);
    nY = ceil((ymax - ymin)/dx);
    nZ = ceil((zmax - zmin)/dx);
    
    % make sure we cover the whole domain (extra voxel in case non exact)
    xnew = xmin:dx:xmin+nX*dx;
    ynew = ymin:dx:ymin+nY*dx;
    znew = zmin:dx:zmin+nZ*dx;
    
    % generate the 3D grid
    r = grid3d(xnew,ynew,znew);
    
    % interpolate the nearest
    epsilon_r = interp3(y,x,z,RHBM.epsilon_r-1,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
    idxnan = isnan(epsilon_r(:)); epsilon_r(idxnan) = 0;
    epsilon_r = epsilon_r + 1;
    sigma_e = interp3(y,x,z,RHBM.sigma_e,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
    idxnan = isnan(sigma_e(:)); sigma_e(idxnan) = 0;
    rho = interp3(y,x,z,RHBM.rho,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
    idxnan = isnan(rho(:));rho(idxnan) = 0;
    
    % mu_r = interp3(y,x,z,mu_r-1,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
    % mu_r = mu_r + 1;
    % sigma_h = interp3(y,x,z,sigma_h,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
    
end

idxS = find( abs(epsilon_r(:) - 1 + 1j*sigma_e(:)) > 1e-20);

% translate domain
if ~isempty(trans)
    trans = squeeze(trans);
    r(:,:,:,1) = r(:,:,:,1) - trans(1);
    r(:,:,:,2) = r(:,:,:,2) - trans(2);
    r(:,:,:,3) = r(:,:,:,3) - trans(3);
end

% -------------------------------------------------------------------------
% Assign data to structure
% -------------------------------------------------------------------------

RHBM = struct('name', RHBM.name,...
              'r',  r,...
              'epsilon_r', epsilon_r, ...
              'sigma_e', sigma_e, ...
              'rho', rho, ...
              'idxS', idxS, ...
              'freqfN', [], ...
              'freqfK', [], ...
              'fN', [], ...
              'fK', [], ...
              'freqM', [], ...
              'Dcoord', [], ...
              'P', [], ...
              'Um', [], ...
              'Sm', [], ...
              'Vm', [], ...
              'X', [], ...
              'M', []);
           
           
