%%    Script for solving the S-parameters of a single 1-port coil 
% _________________________________________________________________________
%
%       Script to test the SIE solver with a single port planar coil
%       Compare with SEMCAD and SCUFF-RF results%       
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


% -------------------------------------------------------------------------
% Generate the Coil geometry to test
% -------------------------------------------------------------------------
if ispc
    SCOILfile = '.\coil_geometry\PLANAR_1coil_1port_SEMCAD.msh';
    % SCOILfile = '.\coil_geometry\PLANAR_1coil_1port_SEMCAD_Refined.msh';
else
    SCOILfile = './coil_geometry/PLANAR_1coil_1port_SEMCAD.msh';
    % SCOILfile = './coil_geometry/PLANAR_1coil_1port_SEMCAD_Refined.msh';
end

SCOIL = Import_SCOIL(SCOILfile);

% -------------------------------------------------------------------------
% Compute the Frequency response
% -------------------------------------------------------------------------

NPORT = length(SCOIL.port);

fHz = linspace(2.5e8,3.5e8,51);
ZP_sie = zeros(NPORT,NPORT,length(fHz));


for kk = 1:length(fHz)
    
    titerkk = tic;
    freq = fHz(kk);
    
    % -------------------------------------------------------------------------
    % Solve the free space case
    % -------------------------------------------------------------------------
    
    titer = tic;
    
    [ZP] = SIE_Solver(SCOIL,freq);
    
    ZP_sie(:,:,kk) = ZP;
    
    fprintf(1, '\n  -- Free space solver applied');
    fprintf(1, '\n          Elapsed time  = %g [sec]', toc(titer));
    fprintf(1, '\n');

    
end


Z0 = 50;
SP_sie = z2s(ZP_sie, Z0);

% -------------------------------------------------------------------------
% PLOT stuff
% -------------------------------------------------------------------------

if ispc
    load('.\semcad_data\PLANAR_1coil_1port_air_semcad.mat');
    load('.\scuff_data\PLANAR_1coil_1port_scuff.mat');
else
    load('./semcad_data/PLANAR_1coil_1port_air_semcad.mat');
    load('./scuff_data/PLANAR_1coil_1port_scuff.mat');
end
fHz = fHz/1e9;

ii = 1;
jj = 1;

figidx = 10000+ii*10+jj;
figure(figidx)
plot(ff_semcad/1e9,squeeze(real(SP_semcad(1,1,:))),'r-','LineWidth',3.0);
hold on
plot(ff_scuff/1e9,squeeze(real(SP_scuff(1,1,:))),'g-','LineWidth',3.0);
plot(fHz,squeeze(real(SP_sie(1,1,:))),'b-','LineWidth',3.0);
axis([0.25 0.35 0.9930 0.9985]);
set(gca, 'FontSize', 14);
h_legend = legend('SEMCAD', 'SCUFF RF', ...
    'SIE', 'Location', 'best');
set(h_legend, 'FontSize', 14);
h_xlabel = xlabel('Frequency (GHz)'); set(h_xlabel, 'FontSize', 16);
h_ylabel = ylabel(sprintf('S_%d_%d -- Real Part', 1, 1)); set(h_ylabel, 'FontSize', 16);
grid on;
hold off;

figidx = 20000+ii*10+jj;
figure(figidx)
plot(ff_semcad/1e9,squeeze(imag(SP_semcad(1,1,:))),'r-','LineWidth',3.0);
hold on
plot(ff_scuff/1e9,squeeze(imag(SP_scuff(1,1,:))),'g-','LineWidth',3.0);
plot(fHz,squeeze(imag(SP_sie(1,1,:))),'b-','LineWidth',3.0);
axis([0.25 0.35 -0.08 0.06]);
set(gca, 'FontSize', 14);
h_legend = legend('SEMCAD', 'SCUFF RF', ...
    'SIE', 'Location', 'best');
set(h_legend, 'FontSize', 14);
h_xlabel = xlabel('Frequency (GHz)'); set(h_xlabel, 'FontSize', 16);
h_ylabel = ylabel(sprintf('S_%d_%d -- Imaginary Part', 1, 1)); set(h_ylabel, 'FontSize', 16);
grid on;
hold off;


figidx = 30000+ii*10+jj;
figure(figidx)
plot(ff_semcad/1e9,squeeze(real(ZP_semcad(1,1,:))),'r-','LineWidth',3.0);
hold on
plot(ff_scuff/1e9,squeeze(real(ZP_scuff(1,1,:))),'g-','LineWidth',3.0);
plot(fHz,squeeze(real(ZP_sie(1,1,:))),'b-','LineWidth',3.0);
axis([0.25 0.35 0 5.5e4]);
set(gca, 'FontSize', 14);
h_legend = legend('SEMCAD', 'SCUFF RF', ...
    'SIE', 'Location', 'best');
set(h_legend, 'FontSize', 14);
h_xlabel = xlabel('Frequency (GHz)'); set(h_xlabel, 'FontSize', 16);
h_ylabel = ylabel(sprintf('Z_%d_%d -- Real Part', 1, 1)); set(h_ylabel, 'FontSize', 16);
grid on;
hold off;


figidx = 40000+ii*10+jj;
figure(figidx)
plot(ff_semcad/1e9,squeeze(imag(ZP_semcad(1,1,:))),'r-','LineWidth',3.0);
hold on
plot(ff_scuff/1e9,squeeze(imag(ZP_scuff(1,1,:))),'g-','LineWidth',3.0);
plot(fHz,squeeze(imag(ZP_sie(1,1,:))),'b-','LineWidth',3.0);
axis([0.25 0.35 -3e4 3e4]);
set(gca, 'FontSize', 14);
h_legend = legend('SEMCAD', 'SCUFF RF', ...
    'SIE', 'Location', 'best');
set(h_legend, 'FontSize', 14);
h_xlabel = xlabel('Frequency (GHz)'); set(h_xlabel, 'FontSize', 16);
h_ylabel = ylabel(sprintf('Z_%d_%d -- Imaginary Part', 1, 1)); set(h_ylabel, 'FontSize', 16);
grid on;
hold off;
