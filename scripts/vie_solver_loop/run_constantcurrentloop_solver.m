%%    Script that solves the VIE for a Constant Current Loop excitation
% _________________________________________________________________________
%
%       Illustrates how to generate a constant current loop
%       Assumes Jc current amplitude
%       Solves the scattering problem and 
%       Generates the E and H fields due to the constant current loop
%
%       Also illustrates how to save the same geometry in a .wmm model
%       Solves the exact same problem using the MARIE GUI wrappers
%       and the saved .wmm
%
%       Also solves the same problem with a realistic wire coil model 
%       using the MR_Solve routine of MARIE
%       i.e. solving the WIE + VIE problem (WVIE solver)
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


clear
close all
clc

% -------------------------------------------------------------------------
% Add the corresponding path 
% -------------------------------------------------------------------------

% find the current folder
currentfolder = pwd;

% find the MARIE folder (two folders up)
if ispc
    idx = find(currentfolder == '\');
else
    idx = find(currentfolder == '/');
end
mariefolder = currentfolder(1:idx(end-1)-1);

% obtain the string with the recursive paths
p = genpath(mariefolder);

% add the path tree to the current path
addpath(p);

% -------------------------------------------------------------------------
% Set the frequency and resolution
% -------------------------------------------------------------------------

freq = 298.2e6; % 7T
resolution = 3e-3; % 3mm

% -------------------------------------------------------------------------
%                   initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% initialize figure counter
figidx = 1;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Load the corresponding body model
% -------------------------------------------------------------------------

% we use the same way of loading as standard GUI

if ispc
    RHBMfile = '.\body_model\RHBM_HT_5mm.vmm';
else
    RHBMfile = './body_model/RHBM_HT_5mm.vmm';
end

% we just set some dimensions to limit the size of the model
[RHBM] = Import_RHBM(RHBMfile,resolution, [-0.153, 0.153], [-0.153, 0.153], [-0.153, 0.153]);

% find scatterer 
[L,M,N,~] = size(RHBM.r);
idxS = find(abs(RHBM.epsilon_r(:)-1 + RHBM.sigma_e(:)) > 1e-12); % non-air positions
nS = length(idxS);
nD = L*M*N; % number of voxels in the complete domain
idxSS = [idxS; nD+idxS; 2*nD+idxS]; % the vector of non-air positions in 3D grid

% coordinates of domain
xd = RHBM.r(:,:,:,1);
yd = RHBM.r(:,:,:,2);
zd = RHBM.r(:,:,:,3);
Dcoord = [xd(:), yd(:), zd(:)];

% scatterer
xs = xd(idxS);
ys = yd(idxS);
zs = zd(idxS);
Scoord = [xs(:), ys(:), zs(:)];

figure(figidx);
plot3(xs, ys, zs ,'bs');
hold off;
axis equal;

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

% properties at the given frequency
e_r = RHBM.epsilon_r - 1j*RHBM.sigma_e/(eo*omega);

% to simplify
r = RHBM.r; 

% get the voxel side
dx = r(2,1,1,1) - r(1,1,1,1);
% Gram matrix (value)
Gram = dx^3; % volume of the voxel


% -------------------------------------------------------------------------
%                  Generate circulants
% -------------------------------------------------------------------------

% compute the circulants
[fN] = getOPERATORS(r,freq,'N',[],'DEMCEM');
[fK] = getOPERATORS(r,freq,'K',[],'DEMCEM');




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% FIRST WAY OF DOING IT... 
% generate the loop by hand and use low level functions to solve
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Generate a Loop geometry
% -------------------------------------------------------------------------

RAD = 0.040; % 40mm radius of the coil
NP = 100; % 100 segments
Pcoil = zeros(NP,3); % positive edges of the segments
Ncoil = zeros(NP,3); % negative edges of the segments
ALPHA = 2*pi/NP;
for ii = 1:NP-1
    
    xp = 0.12; % position in the x axis
    yp = RAD*cos(ALPHA*ii); % just a vertical loop
    zp = RAD*sin(ALPHA*ii);
    
    Pcoil(ii,:) = [xp, yp, zp];
    Ncoil(ii+1,:) = [xp, yp, zp];
    
end

% do not forget to close the loop
Ncoil(1,:) = [0.12,RAD,0];
Pcoil(NP,:) = [0.12,RAD,0];

nC = size(Pcoil,1); % number of elements in the 'coil'


% plot the geometry
orange = [1 0.5 0.2];
for ii = 1:nC
    
    clear x; clear y; clear z;
    x(1) = Ncoil(ii,1); x(2) = Pcoil(ii,1);
    y(1) = Ncoil(ii,2); y(2) = Pcoil(ii,2);
    z(1) = Ncoil(ii,3); z(2) = Pcoil(ii,3);

    figure(figidx);
    hold on;
    plot3(x,y,z,'Color', orange, 'LineWidth', 2.0);
    
end
hold off;
axis equal;
grid on;
view(3);


% -------------------------------------------------------------------------
% Solve directly
% -------------------------------------------------------------------------

tini = tic;
fprintf(1, '\n  -- Constant Current Loop starting');

Jc = 1; % constant current in the loop

% -------------------------------------------------------------------------
% Scale the loop constant current Jc by the length of each segment
% -------------------------------------------------------------------------

Ccoil = (Pcoil + Ncoil)/2; % center of coil segments
Dcoil = Pcoil - Ncoil; % length of each segment

Jcc = Dcoil; % scale the current by the length of each segment
for ii = 1:3
    Jcc(:,ii) = Jcc(:,ii).*Jc; % for each component
end


% -------------------------------------------------------------------------
% Evaluate the DGF at the domain positions
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%  Compute incident fields due to Jc

% [Einc] = E_field_DGF(Jcc,Ccoil,Dcoord,ko);
% [Hinc] = H_field_DGF(Jcc,Ccoil,Dcoord,ko);
[Einc,Hinc] = eval_DGF(Jcc,Ccoil,Dcoord,ko);

% in the whole domain, now reshape them to LMN3 size
Einc = reshape(Einc,L,M,N,3);
Hinc = reshape(Hinc,L,M,N,3);


% -------------------------------------------------------------------------
%                 Solve the system 
% -------------------------------------------------------------------------

% options for the solver
%
%
%   tol         Relative tolerance for the method (DEFAULT 1e-3)
%   form        Choice of the formulation
%                   1 for JVIE I formulation
%                   2 for JVIE II formulation (DEFAULT)
%   solver      Choice of an Iterative Solver
%                   'G' for GMRES
%                   'D' for GMRES with Deflated Restart
%                   'B' for BiCGStab (DEFAULT)
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
%
%
% for more technical details see:
%
%   A. G. Polimeridis, J. Fernández Villena, L. Daniel and J. K. White.
%   "Stable FFT-JVIE Solvers for Fast Analysis of Highly Inhomogeneous Dielectric Objects"
%   Journal of Computational Physics, 269:280-296, 2014.
%

tol = 1e-3;
form = 2;
solver = 'B';
precond = 'N';
[Jout1,flag,relres,iter,resvec,solvert] = JVIE_Solver(Einc,e_r,r,freq,fN,tol,form,solver,precond);


% -------------------------------------------------------------------------
%                   Compute total E field
% -------------------------------------------------------------------------

[Eout1] = E_field_Nop(Jout1,fN,Gram,freq,Einc);

% -------------------------------------------------------------------------
%                   Compute total H field
% -------------------------------------------------------------------------

[Hout1] = H_field_Kop(Jout1,fK,Gram,Hinc);


fprintf(1, '\n          Done (explicit functions)');
fprintf(1, '\n          Elapsed time  = %g [sec]', toc(tini));
fprintf(1, '\n');




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% SECOND WAY OF DOING IT... 
% save the same geometry as a .wsm file
% and use the MARIE constant current loop functionalities
% to solve the same problem
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% save the same geometry as a .wsm file
% -------------------------------------------------------------------------

% Constant Current Loop file name
name = 'verticalloop';
description = 'Vertical loop, 40mm radius, distance in x 120mm';
Dwire = 0.500e-3; % diameter of the wire
Rhocoil = 1.2579e-08; % resistivity of copper

% save the data in a .wmm file
if ispc
    fid = fopen( strcat('.\loop_model\', name, '.wmm'), 'W'); % the discretized filename
else
    fid = fopen( strcat('./loop_model/', name, '.wmm'), 'W'); % the discretized filename
end

% write the data
fprintf(fid, '%s\n', description);
fprintf(fid, '%g %f\n', Rhocoil, Dwire);
fprintf(fid, '%s.wsd', name);
fclose(fid);

% create the discretized file name
if ispc
    cclfile = strcat('.\loop_model\', name, '.wsd'); % the discretized filename
else
    cclfile = strcat('./loop_model/', name, '.wsd'); % the discretized filename
end

% set the ports: here are not used, but are needed for the format
port = [1]; % 1 port, on the first segment
Bcoil = zeros(NP,1);
Bcoil(port) = 1; % assign on the Bcoil the entry of the first segment to 1

% prepare the data of the coil to save
DATAcoil = [Pcoil, Ncoil, Bcoil];

% export the data in ascii mode
save(cclfile, 'DATAcoil', '-ascii');


% clear loop related data
clear Pcoil; clear Ncoil; clear Bcoil;
clear DATAcoil; clear port;
clear description; clear Rhocoil; clear Dwire;


% -------------------------------------------------------------------------
% Solve using the MARIE GUI wrappers
% -------------------------------------------------------------------------

tini = tic;
fprintf(1, '\n  -- Constant Current Loop starting');


% -------------------------------------------------------------------------
% Load loop and create excitation with MARIE GUI function
% -------------------------------------------------------------------------

% this is pretty much the same as before, but reading from the .wsd file

% Exctype for constant current loop is 'L'
% vec1 is the name of the loop file: the .wsd (in cclfile)
% vec2 is the amplitude of the constant current
Exctype = 'L';
[Einc,Hinc] = Generate_Excitation(RHBM.r,freq,Exctype,cclfile,Jc);


% -------------------------------------------------------------------------
%  Solve using the MARIE GUI solver 
% -------------------------------------------------------------------------

% add the circulant to the RHBM struct
RHBM.fN = fN; RHBM.freqfK = freq;
RHBM.fK = fK; RHBM.freqfN = freq;

% call the solver (tolerance is the same as before)
[Jout,Sout,Eout,Bout,Gsar,Pabs,~,~] = Scat_Solver(Einc,Hinc,freq,RHBM.r,RHBM.epsilon_r,RHBM.sigma_e,RHBM.rho,tol,RHBM.fN,RHBM.fK);


fprintf(1, '\n          Done (MARIE GUI wrappers)');
fprintf(1, '\n          Elapsed time  = %g [sec]', toc(tini));
fprintf(1, '\n');


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% NOW, solve the same geometry as it was a realistic wire coil... 
% save the same geometry as a .wsm file
% and use the MARIE WVIE functionalities
% to solve the same problem
% Note that we use an unitary excitation at the port
% and the current is no constant along the loop
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Generate the Wire Coil geometry for the wmm that we created before
% -------------------------------------------------------------------------

if ispc
    WCOILfile = strcat('.\loop_model\', name, '.wmm');
else
    WCOILfile = strcat('./loop_model/', name, '.wmm');
end

WCOIL = Import_COIL(WCOILfile);

% -------------------------------------------------------------------------
% Solve with the WIRE approach (WVIE)
% -------------------------------------------------------------------------

% generate the system matrix
fprintf(1, '\n  -- Applying VWIE solver');
tini = tic;

% assign the coil structure
COIL = struct('name', name,...
    'type', 'W',...
    'Rhocoil', WCOIL.Rhocoil,...
    'Thickness', [],...
    'index',  [],...
    'etod', [], ...
    'node', [], ...
    'edge', [], ...
    'elem', [], ...
    'index_elem', [], ...
    'Ct', [], ...
    'Ln', [], ...
    'Pn', [], ...
    'Pcoil', WCOIL.Pcoil,...
    'Ncoil', WCOIL.Ncoil,...
    'Dwire', WCOIL.Dwire,...
    'port', WCOIL.port);


% call the MARIE solver with the coil
[ZPw,Jcw,Jbw,Sbw,Ebw,Bbw,Gsarw,Pabsw] = MR_Solver(RHBM,COIL,freq,tol);

fprintf(1, '\n          Elapsed time  = %g [sec]', toc(tini));
fprintf(1, '\n');


% -------------------------------------------------------------------------
% Done
% -------------------------------------------------------------------------

maximumval = 0; % if 0: autoscale in the plot
minimumval = 0;
mapvec = maximumval*ones(L,M); % to set maximum
mapvec(L,M) = minimumval; % to set minimum

% for the comparison we just weight the Jout so that the current is approx
% the same as in the wire case
Wamp = -max(Jcw(:));
Jout = Jout*Wamp;

figidx = figidx + 1;
figure(figidx); vec = sum(conj(Jbw).*Jbw,4);
imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(:,:,round(N/2))))); hold off;
colorbar; colormap('hot'); title('RMS(J) for VWIE');% wire
figidx = figidx + 1;
figure(figidx); vec = sum(conj(Jout).*Jout,4);
imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(:,:,round(N/2))))); hold off;
colorbar; colormap('hot'); title('RMS(J) for CCL');% wire % ccl
figidx = figidx + 1;
figure(figidx); vec = sum(conj(Jout-Jbw).*(Jout-Jbw),4);
imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(:,:,round(N/2))))); hold off;
colorbar; colormap('hot'); title('RMSE(J) CCL minus VWIE');% wire % ccl% error


