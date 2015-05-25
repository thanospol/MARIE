%%    Script to solve for circular coils in the presence of scattarer
% _________________________________________________________________________
%
%       Script to test the SIE+VIE and WIE+VIE solvers
%       2 port circular coil
%       Compares S-PARAMETERS w/ and wo/ homogeneous CUBE
%       SVIE vs WVIE results
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

clear all;
close all;
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



% % -------------------------------------------------------------------------
% % Define the COIL structure
% % -------------------------------------------------------------------------
% 
% % Empty coil structure
% COIL = struct('name', 'No Model Selected',...
%     'type', 'N',...
%     'Rhocoil', [],...
%     'Thickness', [],...
%     'index',  [],...
%     'etod', [], ...
%     'node', [], ...
%     'edge', [], ...
%     'elem', [], ...
%     'index_elem', [], ...
%     'Ct', [], ...
%     'Ln', [], ...
%     'Pn', [], ...
%     'port', [],...
%     'Pcoil', [],...
%     'Ncoil', [],...
%     'Dwire', []);


% -------------------------------------------------------------------------
% Generate the Surface Coil geometry for the msh
% -------------------------------------------------------------------------

if ispc
    % filename = '.\coils_geometry\Circularcoil_ref.msh'; % refined mesh
    SCOILfile = '.\coils_geometry\Circularcoil.msh'; % std mesh
    % filename = '.\coils_geometry\Circularcoil_coarse.msh'; % coarse mesh
else
    % filename = './coils_geometry/Circularcoil_ref.msh'; % refined mesh
    SCOILfile = './coils_geometry/Circularcoil.msh'; % std mesh
    % filename = './coils_geometry/Circularcoil_coarse.msh'; % coarse mesh
end

SCOIL = Import_SCOIL(SCOILfile);


% -------------------------------------------------------------------------
% Generate the Wire Coil geometry for the WIRE model
% -------------------------------------------------------------------------

% AWG 18 was used to model the wire
% diameter 1.024 mm
% area 0.823 mm2
% cooper:
%   sigma = 7.95e7 S/m
%   rho = 1.2579e-08 Ohm/m
%   delta = 3.26 um at 300MHz (skin depth)
%   R = rho/A = 0.0153 Ohm/m  (wiki Table of AWG18 0.02095 Ohm/m)
%   R(300MHz) = rho/(PI*(D-delta)*delta) = 1.2032 Ohm/m
%

WCOIL = struct('name', 'Circular loop wire',...
    'type', 'W',...
    'Rhocoil', [],...
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
    'Pcoil', [],...
    'Ncoil', [],...
    'Dwire', [],...
    'port', []);

WCOIL.Dwire = 0.500e-3; % diameter of the wire
WCOIL.Rhocoil = 1.2579e-08; % resistivity of copper
% deltaw = sqrt(2*rhow/(2*pi*f*mu)); % skin depth
% Rw = rhow/(pi*(Dw-deltaw)*deltaw); % resistance per meter of the wire

RAD = 0.040; % 40mm radius of the coil
NP = 100; % 100 segments
WCOIL.Pcoil = zeros(NP,3);
WCOIL.Ncoil = zeros(NP,3);
ALPHA = 2*pi/NP;
for ii = 1:NP-1
    
    xp = RAD*cos(ALPHA*ii);
    yp = RAD*sin(ALPHA*ii);
    zp = 0.12;
    
    WCOIL.Pcoil(ii,:) = [xp, yp, zp];
    WCOIL.Ncoil(ii+1,:) = [xp, yp, zp];
    
end

WCOIL.Ncoil(1,:) = [RAD,0,0.12];
WCOIL.Pcoil(NP,:) = [RAD,0,0.12];

nC = size(WCOIL.Pcoil,1);
WCOIL.port = [1, ceil((nC+1)/2)];


% -------------------------------------------------------------------------
% Generate the cube
% -------------------------------------------------------------------------

RHBMfile = 'CUBE_loop';

% size of the cube
Lside = 0.200; % 200 mm of side
Nvoxels = 20; % 50 voxels of discretization for the cube

% Domain is going to be at 2 extra voxel larger Lside
dx = Lside/Nvoxels;
x = -Nvoxels*dx/2:dx:(Nvoxels/2-1)*dx;
x = x + dx/2;
x = [x(1)-dx, x, x(end)+dx];
r = grid3d(x,x,x);

% Homogeneous CUBE
epsilon = 36;
sigma = 0.657;
density = 1000;
Center = [0,0,0];
Ldim = [Lside, Lside, Lside];

[idxS,epsilon_r,sigma_e,rho,~,~] = homogencube(r,Center,Ldim,epsilon,sigma,density);

RHBM = struct('name', 'Cube, side 20cm, eps_r=36, sig=0.657',...
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
% Compute the Frequency response
% -------------------------------------------------------------------------

NPORT = length(SCOIL.port);

fHz = linspace(2.5e8,3.5e8,5);
% fHz = 298.2e6;

Zsie_free = zeros(NPORT,NPORT,length(fHz));
Zmom_free = zeros(NPORT,NPORT,length(fHz));
Zsie_iter = zeros(NPORT,NPORT,length(fHz));
Zmom_iter = zeros(NPORT,NPORT,length(fHz));


for kk = 1:length(fHz)
    
    titerkk = tic;
    freq = fHz(kk);
    tol = 1e-3; coup = 5; gpu_flag = 1;
    
      
    % -------------------------------------------------------------------------
    % WIRE
    % -------------------------------------------------------------------------
  
    % generate the system matrix
    titer = tic;
    
    [ZPmom,~] = WIE_Solver(WCOIL,freq);
    
    Zmom_free(:,:,kk) = ZPmom;
    
    fprintf(1, '\n  -- Free space wire solver applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');

   
    % -------------------------------------------------------------------------
    % WIRE
    % -------------------------------------------------------------------------
    
    % generate the system matrix
    titer = tic;
    
    [ZPmomit,Jcmom,Jbmom] = VWIE_Solver(RHBM,WCOIL,freq,tol,coup,gpu_flag);
    
    Zmom_iter(:,:,kk) = ZPmomit;
    
    fprintf(1, '\n  -- Iterative wire solver applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');


    % -------------------------------------------------------------------------
    % Solve the free space case
    % -------------------------------------------------------------------------
    
    titer = tic;
    
    [ZPfree] = SIE_Solver(SCOIL,freq);
    
    Zsie_free(:,:,kk) = ZPfree;
    
    fprintf(1, '\n  -- Free space solver applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');

    % -------------------------------------------------------------------------
    % Solve using the iterative scheme
    % -------------------------------------------------------------------------
    
    titer = tic;
    
    [ZPiter,Jciter,Jbiter] = VSIE_Solver(RHBM,SCOIL,freq,tol,coup,gpu_flag);
    
    Zsie_iter(:,:,kk) = ZPiter;
    
    fprintf(1, '\n  -- Iterative method applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');
    
        
    % -------------------------------------------------------------------------
    % Done
    % -------------------------------------------------------------------------
    
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titerkk));
    fprintf(1, '\n');
    
    
end


Z0 = 50;
SPsie_free = z2s(Zsie_free, Z0);
SPmom_free = z2s(Zmom_free, Z0);
SPsie_iter = z2s(Zsie_iter, Z0);
SPmom_iter = z2s(Zmom_iter, Z0);

Tports = size(SPsie_free,1);
for ii = 1:Tports
    for jj = 1%:Tports
          
        figidx = 1000+ii*10+jj;
        figure(figidx)
        hold on
        plot(fHz,squeeze(real(SPsie_free(ii,jj,:))),'k-','LineWidth',2.0);
        plot(fHz,squeeze(real(SPsie_iter(ii,jj,:))),'k--','LineWidth',2.0);
        plot(fHz,squeeze(real(SPmom_free(ii,jj,:))),'b-','LineWidth',2.0);
        plot(fHz,squeeze(real(SPmom_iter(ii,jj,:))),'b--','LineWidth',2.0);
        set(gca, 'FontSize', 16);
        h_legend = legend('SIE Free', 'SVIE CUBE', 'WIE Free', 'WVIE CUBE');
        set(h_legend, 'FontSize', 16);
        h_xlabel = xlabel('Frequency'); set(h_xlabel, 'FontSize', 16);
        h_ylabel = ylabel(sprintf('Real(S_%d_%d)',ii,jj)); set(h_ylabel, 'FontSize', 16);
        grid on;
        hold off
        
        figidx = 2000+ii*10+jj;
        figure(figidx)
        hold on
        plot(fHz,squeeze(imag(SPsie_free(ii,jj,:))),'k-','LineWidth',2.0);
        plot(fHz,squeeze(imag(SPsie_iter(ii,jj,:))),'k--','LineWidth',2.0);
        plot(fHz,squeeze(imag(SPmom_free(ii,jj,:))),'b-','LineWidth',2.0);
        plot(fHz,squeeze(imag(SPmom_iter(ii,jj,:))),'b--','LineWidth',2.0);
        set(gca, 'FontSize', 16);
        h_legend = legend('SIE Free', 'SVIE CUBE', 'WIE Free', 'WVIE CUBE');
        set(h_legend, 'FontSize', 16);
        h_xlabel = xlabel('Frequency'); set(h_xlabel, 'FontSize', 16);
        h_ylabel = ylabel(sprintf('Imag(S_%d_%d)',ii,jj)); set(h_ylabel, 'FontSize', 16);
        grid on;
        hold off
        
    end
end



