%%    SPHERE test vs MIE series
% _________________________________________________________________________
%
%       Script to test the VIE solver with a sphere
%       Compare with MIE series
%
%       It explains some advanced options for the solver
%
% _________________________________________________________________________
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

close all
clear
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
%            Define the Model to Simulate
% -------------------------------------------------------------------------

modelname = 'Sphere';

% frequency
freq = 298.2e6;

% sphere dimensions and properties
Rad = 0.10; % 10cm of radius
Cnt = [0 0 0]; % center at origin
er = 52;  % value of epsilon_r within sphere
se = 0.55; % value of sigma_e within sphere
dens = 1000; % value of density within sphere

% Resolution
Res = 5e-3; % 5mm

% -------------------------------------------------------------------------
%                   initialize stuff
% -------------------------------------------------------------------------

% generate EM constants
EMconstants;

% initialize figure counter
figidx = 1;

% -------------------------------------------------------------------------
%                   Define domain and sphere
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generatedomain(Res,Rad);

% generate sphere as homogeneous object
[idx,epsilon_r,sigma_e,rho] = homogensphere(r,Cnt,Rad,er,se,dens);

% plot geometry
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

figure(figidx);
plot3(xd(idx), yd(idx), zd(idx), 's');
axis equal;
grid on;

% -------------------------------------------------------------------------
%                   Define Excitation
% -------------------------------------------------------------------------

% Define excitation - a plane wave
Eo = [1,0,0]; % polarization of E-field
k = ko * [0,0,1];
[Einc, Hinc] = PlaneWave_Excitation(r,k,omega_mu,Eo);

% -------------------------------------------------------------------------
%                 define EM vars and constants
% -------------------------------------------------------------------------

% properties at the given frequency
e_r = epsilon_r - 1j*sigma_e/(eo*omega);

% get domain dimensions
[L,M,N,~] = size(r);

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

tol = 1e-4;
form = 1;
solver = 'G';
precond = 'R';
[Jout,flag,relres,iter,resvec,solvert] = JVIE_Solver(Einc,e_r,r,freq,fN,tol,form,solver,precond);


% -------------------------------------------------------------------------
%                   Compute SAR
% -------------------------------------------------------------------------

Sout= Jout.*conj(Jout);
Sout = sum(Sout,4);
idxsar = find(abs(sigma_e(:))); % find position of elements where sigma is not zero
Sout(idxsar) = Sout(idxsar)./(sigma_e(idxsar));
idxsar = find(abs(rho(:))); % find position of elements where rho is not zero
Sout(idxsar) = Sout(idxsar)./(2*rho(idxsar));

% -------------------------------------------------------------------------
%                   Compute E field
% -------------------------------------------------------------------------

[Eout] = E_field_Nop(Jout,fN,Gram,freq,Einc);

% -------------------------------------------------------------------------
%                   Compute H field
% -------------------------------------------------------------------------

[Hout] = H_field_Kop(Jout,fK,Gram,Hinc);

% -------------------------------------------------------------------------
%                   Compute Global SAR and absorbed Power
% -------------------------------------------------------------------------

Pabs = Eout.*conj(Jout);
Pabs = sum(Pabs,4);
Pabs = 0.5*real(sum(Pabs(:)))*Gram;
Gsar = sum(Sout(:))*Gram;

% -------------------------------------------------------------------------
%                   Solve with MIE
% -------------------------------------------------------------------------

[E_Mie,H_Mie] = MIE_SERIES(r,Rad,er,se,freq);

% -------------------------------------------------------------------------
%              Compare results
% -------------------------------------------------------------------------

figidx = 1;

% mask fields outside sphere
Evie = zeros(L,M,N,3);
Hvie = zeros(L,M,N,3);
idx = find(abs(epsilon_r(:) -1));
Evie(idx) = Eout(idx);
Evie(L*M*N+idx) = Eout(L*M*N+idx);
Evie(2*L*M*N+idx) = Eout(2*L*M*N+idx);
Hvie(idx) = Hout(idx);
Hvie(L*M*N+idx) = Hout(L*M*N+idx);
Hvie(2*L*M*N+idx) = Hout(2*L*M*N+idx);

% -------------------------------------------------------------------------
%              Currents

RMS_J = sqrt(sum(Jout.*conj(Jout),4)/2);
scale = max(RMS_J(:));
scale = scale*1.05;
xcut = 0; ycut = 0; zcut = 0;
% plot views
[figidx,fighandle] = plot_3Dvector(RMS_J,r,xcut,ycut,zcut,figidx,scale,strcat(modelname, 'RMS(J) VIE'));

% -------------------------------------------------------------------------
%              E field

RMS_E = sqrt(sum(Evie.*conj(Evie),4)/2);
scale = max(RMS_E(:));
scale = scale*1.05;
xcut = 0; ycut = 0; zcut = 0;
% plot views
[figidx,fighandle] = plot_3Dvector(RMS_E,r,xcut,ycut,zcut,figidx,scale,strcat(modelname, 'RMS(E) VIE'));

RMS_E = abs(sqrt(sum(E_Mie.*conj(E_Mie),4)/2));
% plot views
[figidx,fighandle] = plot_3Dvector(RMS_E,r,xcut,ycut,zcut,figidx,scale,strcat(modelname, 'RMS(E) MIE'));

% -------------------------------------------------------------------------
%              H field

RMS_H = sqrt(sum(Hvie.*conj(Hvie),4)/2);
scale = max(RMS_H(:));
scale = scale*1.05;
xcut = 0; ycut = 0; zcut = 0;
% plot views
[figidx,fighandle] = plot_3Dvector(RMS_H,r,xcut,ycut,zcut,figidx,scale,strcat(modelname, 'RMS(H) VIE'));

RMS_H = abs(sqrt(sum(abs(H_Mie).*abs(H_Mie),4)/2));
% plot views
[figidx,fighandle] = plot_3Dvector(RMS_H,r,xcut,ycut,zcut,figidx,scale,strcat(modelname, 'RMS(H) MIE'));

