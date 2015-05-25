%%    Script to solve for circular coils in the presence of scattarer
% _________________________________________________________________________
%
%       Script to test the SIE+VIE and WIE+VIE solvers
%       2 port circular coil
%       Compares FIELDS w/ and wo/ homogeneous CUBE
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

WCOIL = struct('Pcoil', [],...
               'Ncoil', [],...
               'Dwire', [],...
               'Rhocoil', [],...
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

freq = 298.2e6;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Solve with the WIRE approach (WVIE)
% -------------------------------------------------------------------------

% generate the system matrix
titer = tic;

% assign the coil structure
COIL = struct('name', 'Circular loop wire',...
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

tol = 1e-3;

[ZPw,Jcw,Jbw,Sbw,Ebw,Bbw,Gsarw,Pabsw] = MR_Solver(RHBM,COIL,freq,tol);

fprintf(1, '\n  -- VWIE solver applied');
fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
fprintf(1, '\n');


% -------------------------------------------------------------------------
% tune and stuff
% -------------------------------------------------------------------------

% with body
Cseries = 0.886*1e-12;
Cpar = 17.445*1e-12;

% % w/out body
% Cseries = 0.837*1e-12;
% Cpar = 64.045*1e-12;

Dport = 1;
Nports = size(ZPw,1);
VIN = zeros(Nports,1);
VIN(Dport,1) = -1;


% if you want tp change the tunning, matching, change the value of Cs
omega = 2*pi*freq;
Yt = 1j*omega*Cseries;
Yp = 1j*omega*Cpar;

% generate Ip for the given unitary Vp at port 5
Ycap = Yt*eye(Nports,Nports);
Ycap(Dport,Dport) = Yp; % the driving port

% get the loaded system with capacitors in parallel
YPw = inv(ZPw);
YL = YPw + Ycap;
ZL = inv(YL);

% port admittance, voltage and current of the loaded system
zporttun = ZL(Dport,Dport);
Yport = 1/zporttun;
Vport = (1 + 50*Yport)\VIN(Dport,1);
Iport = (VIN(Dport,1) - Vport)/50;

% excite the system with the current at the port
IIN = zeros(Nports,1);
IIN(Dport,1) = Iport;
Vptun = ZL*IIN; % get voltage at all ports
Iptun = YPw*Vptun; % generate the currents flowing into the coils

SPw = z2s(zporttun,50);
fprintf(1, '\n  S parameter after  tuning: %f (Z = %g + j %g)', 20*log10(abs(SPw)), real(zporttun), imag(zporttun));
fprintf(1, '\n  Port current after tuning:\n');
for ii = 1:length(Iptun)
    fprintf(1, '    Port %d: I = %g + j %g,   V = %g + j %g \n', ii, real(Iptun(ii)), imag(Iptun(ii)), real(Vptun(ii)), imag(Vptun(ii)));
end

[L,M,N,~] = size(Jbw);
Jbw = reshape(Jbw,3*L*M*N,Nports);
Jw = Jbw*Vptun; % currents in the body tunned
Jw = reshape(Jw,L,M,N,3);
Jbw = reshape(Jbw,L,M,N,3,Nports);

[L,M,N,~] = size(Ebw);
Ebw = reshape(Ebw,3*L*M*N,Nports);
Ew = Ebw*Vptun; % currents in the body tunned
Ew = reshape(Ew,L,M,N,3);
Ebw = reshape(Ebw,L,M,N,3,Nports);

Bbw = reshape(Bbw,3*L*M*N,Nports);
Bw = Bbw*Vptun; % currents in the body tunned
Bw = reshape(Bw,L,M,N,3);
Bbw = reshape(Bbw,L,M,N,3,Nports);

Sbw = reshape(Sbw,L*M*N,Nports);
Sw = Sbw*Vptun; % currents in the body tunned
Sw = reshape(Sw,L,M,N);
Sbw = reshape(Sbw,L,M,N,Nports);

Pw = Pabsw.'*Vptun;
Gw = Gsarw.'*Vptun;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Solve with SIE (VSIE)
% -------------------------------------------------------------------------

titer = tic;

% define the coil structure for the SIE model

COIL = struct('name', 'Circular Loop Surface model',...
    'type', 'S',...
    'Rhocoil', 0,...
    'Thickness', 1e-6,...
    'index',  SCOIL.index,...
    'etod', SCOIL.etod, ...
    'node', SCOIL.node, ...
    'edge', SCOIL.edge, ...
    'elem', SCOIL.elem, ...
    'index_elem', SCOIL.index_elem, ...
    'Ct', SCOIL.Ct, ...
    'Ln', SCOIL.Ln, ...
    'Pn', SCOIL.Pn, ...
    'port', SCOIL.port,...
    'Pcoil', [],...
    'Ncoil', [],...
    'Dwire', []);

tol = 1e-3;

[ZPs,Jcs,Jbs,Sbs,Ebs,Bbs,Gsars,Pabss] = MR_Solver(RHBM,COIL,freq,tol);

fprintf(1, '\n  -- VSIE solver applied');
fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
fprintf(1, '\n');

% -------------------------------------------------------------------------
% tune and stuff
% -------------------------------------------------------------------------

% with body
Cseries = 0.875*1e-12;
Cpar = 17.850*1e-12;

% % w/out body
% Cseries = 0.825*1e-12;
% Cpar = 72.555*1e-12;

Dport = 1;
VIN = zeros(Nports,1);
VIN(Dport,1) = 1;

% if you want tp change the tunning, matching, change the value of Cs
omega = 2*pi*freq;
Yt = 1j*omega*Cseries;
Yp = 1j*omega*Cpar;

% generate Ip for the given unitary Vp at port 5
Ycap = Yt*eye(Nports,Nports);
Ycap(Dport,Dport) = Yp; % the driving port

% get the loaded system with capacitors in parallel
YPs = inv(ZPs);
YL = YPs + Ycap;
ZL = inv(YL);

% port admittance, voltage and current of the loaded system
zporttun = ZL(Dport,Dport);
Yport = 1/zporttun;
Vport = (1 + 50*Yport)\VIN(Dport,1);
Iport = (VIN(Dport,1) - Vport)/50;

% excite the system with the current at the port
IIN = zeros(Nports,1);
IIN(Dport,1) = Iport;
Vptun = ZL*IIN; % get voltage at all ports
Iptun = YPs*Vptun; % generate the currents flowing into the coils


SPs = z2s(zporttun,50);
fprintf(1, '\n  S parameter after  tuning: %f (Z = %g + j %g)\n', 20*log10(abs(SPs)), real(zporttun), imag(zporttun));
fprintf(1, '\n  Port current after tuning:\n');
for ii = 1:length(Iptun)
    fprintf(1, '    Port %d: I = %g + j %g,   V = %g + j %g \n', ii, real(Iptun(ii)), imag(Iptun(ii)), real(Vptun(ii)), imag(Vptun(ii)));
end


[L,M,N,~] = size(Jbs);
Jbs = reshape(Jbs,3*L*M*N,Nports);
Js = Jbs*Iptun; % currents in the body tunned
Js = reshape(Js,L,M,N,3);
Jbs = reshape(Jbs,L,M,N,3,Nports);

[L,M,N,~] = size(Ebs);
Ebs = reshape(Ebs,3*L*M*N,Nports);
Es = Ebs*Iptun; % currents in the body tunned
Es = reshape(Es,L,M,N,3);
Ebs = reshape(Ebs,L,M,N,3,Nports);

Bbs = reshape(Bbs,3*L*M*N,Nports);
Bs = Bbs*Iptun; % currents in the body tunned
Bs = reshape(Bs,L,M,N,3);
Bbs = reshape(Bbs,L,M,N,3,Nports);

Sbs = reshape(Sbs,L*M*N,Nports);
Ss = Sbs*Iptun; % currents in the body tunned
Ss = reshape(Ss,L,M,N);
Sbs = reshape(Sbs,L,M,N,Nports);

Ps = Pabss.'*Iptun;
Gs = Gsars.'*Iptun;


% -------------------------------------------------------------------------
% Done
% -------------------------------------------------------------------------

maximumval = 0; % if 0: autoscale in the plot
minimumval = 0;
mapvec = maximumval*ones(L,M); % to set maximum
mapvec(L,M) = minimumval; % to set minimum

figure(1); vec = sum(conj(Jw).*Jw,4); imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(11,:,:)))); hold off; colorbar; % wire
figure(2); vec = sum(conj(Js).*Jw,4); imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(11,:,:)))); hold off; colorbar; % sie
figure(3); vec = sum(conj(Js-Jw).*(Js-Jw),4); imagesc(mapvec); hold on; imagesc(squeeze(abs(vec(11,:,:)))); hold off; colorbar; % error
