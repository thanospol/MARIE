function [RHBM] = Import_RHBM(vsmfile,resolution,xlimits,ylimits,zlimits) 
%%    Loads a RHBM structure from .vmm file
% _________________________________________________________________________
%
%   Reads the voxelized marie model (.vmm) file
%   and generates the data for the RHBM
%
% _________________________________________________________________________
%
%
%% INPUT
%       vsmfile - name of the .vmm file, with path
%       resolution - desired resolution for the discretization (min if empty)
%       limits - limits of domain (automatic if left empty)
%
%
%% OUTPUT
%       RHBM -  Realistic Human Body Model structure with
%           name        - name/description of model
%           r           - mapping of the internal edge number to dof number
%           epsilon_r   - voxel dielectric
%           sigma_e     - voxel conductivity
%           rho         - voxel density
%           idxS        - positions of non-air voxels
%           freqfN      - frequency for which the data below is computed
%           freqfK      - frequency for which the data below is computed
%           fN          - N operator circulant in fft domain
%           fK          - K operator circulant in fft domain
%           freqM       - frequency of the MRGFs
%           Dcoord      - coordinates of the DEIM points (Npx3)
%           P           - incidence matrix of the DEIM points
%           Um          - MRGF for scattering  (Um*Sm*Vm')
%           Sm          - MRGF for scattering
%           Vm          - MRGF for scattering
%           X           - MRGF deim interpolation matrix
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
%           Copyright � 2014
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


% -------------------------------------------------------------------------
% Open file and read subfiles
% -------------------------------------------------------------------------

fid = fopen(vsmfile, 'r');
if (fid < 0)
    fprintf(1, ' \n WARNING!: Unable to open file: %s\n', vsmfile);
    RHBM = [];
    return;
end


modelname = fgetl(fid); % read first line with name
line = fgetl(fid); % read line
nchar = length(line);
extension = line(nchar-3:nchar);

if strcmp(extension, '.dat') % dat file
    
    voxelfile = line; % voxel dat file
    xfile = fgetl(fid); % x coord dat file
    yfile = fgetl(fid); % x coord dat file
    zfile = fgetl(fid); % x coord dat file
        
    % -------------------------------------------------------------------------
    %   Parse
    % -------------------------------------------------------------------------
	voxelfile = which(voxelfile);
	xfile = which(xfile);
	yfile = which(yfile);
	zfile = which(zfile);
    
    [x,y,z,epsilon_r,sigma_e,rho,~,~,~] = RHBMDat_ParseMEX(voxelfile,xfile,yfile,zfile);

else
    
    if strcmp(extension, '.mat') % mat saved file
    
        load(line); % load the data
        % it should have r, epsilon_r, sigma_e and rho of the model
        
        % get the coordinates in case of truncation
        x = unique(r(:,:,:,1));
        y = unique(r(:,:,:,2));
        z = unique(r(:,:,:,3));
            
    else
        
        fprintf(1, ' \n WARNING!: Invalid file format\n');
        RHBM = [];
        return;
    
    end
    
end
        
% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------

fprintf(1, '\n Pre-processing the model to fit regular grid\n');

% limits of domain
if isempty(xlimits)
    xmin = min(x); xmax = max(x);
else
    xmin = min(xlimits);
    xmax = max(xlimits);
end

if isempty(ylimits)
    ymin = min(y); ymax = max(y);
else
    ymin = min(ylimits);
    ymax = max(ylimits);
end

if isempty(zlimits)
    zmin = min(z); zmax = max(z);
else
    zmin = min(zlimits);
    zmax = max(zlimits);
end

% define discretization
if isempty(resolution)
    dx = min(x(2:end) - x(1:end-1)); % min discretization
    dy = min(y(2:end) - y(1:end-1)); % min discretization
    dz = min(z(2:end) - z(1:end-1)); % min discretization
    dx = min([dx, dy, dz]); % get the minimum
else % force resolution
    dx = resolution;
end

% number of voxels in each direction
nX = ceil((xmax - xmin)/dx);
nY = ceil((ymax - ymin)/dx);
nZ = ceil((zmax - zmin)/dx);

fprintf(1, ' domain %dx%dx%d, resolution %g\n', nX,nY,nZ,dx);

% make sure we cover the whole domain (extra voxel in case non exact)
xnew = xmin:dx:xmin+nX*dx;
ynew = ymin:dx:ymin+nY*dx;
znew = zmin:dx:zmin+nZ*dx;

% generate the 3D grid
r = grid3d(xnew,ynew,znew);

% interpolate the nearest
epsilon_r = interp3(y,x,z,epsilon_r-1,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
idxnan = isnan(epsilon_r(:)); epsilon_r(idxnan) = 0;
epsilon_r = epsilon_r + 1;
sigma_e = interp3(y,x,z,sigma_e,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
idxnan = isnan(sigma_e(:)); sigma_e(idxnan) = 0;
rho = interp3(y,x,z,rho,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
idxnan = isnan(rho(:));rho(idxnan) = 0;

% clear possible approximation errors in rho
idxD = find( abs(rho(:)) > 1e-20);
idxS = find( abs(epsilon_r(:) - 1 + 1j*sigma_e(:)) > 1e-20);
idxA = setdiff(idxD,idxS); % air indexes
rho(idxA) = 0;

% mu_r = interp3(y,x,z,mu_r-1,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');
% mu_r = mu_r + 1;
% sigma_h = interp3(y,x,z,sigma_h,r(:,:,:,2),r(:,:,:,1),r(:,:,:,3), 'nearest');

fprintf(1, ' ... done!\n');

% -------------------------------------------------------------------------
% Assign data to structure
% -------------------------------------------------------------------------

RHBM = struct('name', modelname,...
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
           
           
% -------------------------------------------------------------------------
% If there is MRGF
% -------------------------------------------------------------------------

mrgfflag = 0;
if exist('Dcoord', 'var')
    RHBM.Dcoord = Dcoord;
    mrgfflag = 1;    
end

if exist('freqM', 'var')
    RHBM.freqM = freqM;
else
    mrgfflag = 0;
end

if exist('Sm', 'var')
    RHBM.Sm = Sm;
else
    mrgfflag = 0;
end

if exist('Um', 'var')
    RHBM.Um = Um;
else
    mrgfflag = 0;
end

if exist('Vm', 'var')
    RHBM.Vm = Vm;
else
    mrgfflag = 0;
end


if (mrgfflag == 0)
    RHBM.freqM = [];
    RHBM.Dcoord = [];
    RHBM.Sm = [];
    RHBM.Um = [];
    RHBM.Vm = [];
else
    if exist('M', 'var')
        RHBM.M = M;
    else
        mrgfflag = 0;
    end
    
    if exist('P', 'var')
        RHBM.P = P;
    else
        mrgfflag = 0;
    end
    
    if exist('X', 'var')
        RHBM.X = X;
    else
        mrgfflag = 0;
    end
    
    if (mrgfflag == 0)
        RHBM.X = [];
        RHBM.P = [];
        RHBM.M = [];
    end
    
end


          
