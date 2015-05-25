%%    MARIE - MAgnetic Resonance Integral Equation suite
% _________________________________________________________________________
%
%
%   Def.: Main script... to run MARIE just call this
%
% -------------------------------------------------------------------------
%
%     MARIE - Magnetic Resonance Integral Equation suite
%     Copyright (C)2014 Jorge Fernandez Villena / Athanasios G. Polimeridis 
%     Computational Prototyping Group
%     Research Laboratory of Electronics
%     Massachusetts Institute of Technology
%     contact: jvillena@mit.edu / thanos_p@mit.edu
% 
%     MARIE is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     MARIE is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% _________________________________________________________________________
%
%
%% DATA
%       RHBM  structure
%           name -  name of model
%           r - mapping of the internal edge number to dof number
%           epsilon_r - voxel dielectric
%           sigma_e  - voxel conductivity
%           rho  - voxel density
%           idxS - indexes of non-air material voxels
%           freqfN - freq of N operator circulant
%           freqfK - freq of K operator circulant
%           fN - N operator circulant
%           fK - K operator circulant
%           freqM - frequency of the MRGF
%           Dcoord - coordinates of DEIM points
%           Um, Sm, Vm - MRGFs for coil analysis
%           P - DEIM selection matrix
%           X - DEIM coupling weighting matrix
%           M - MRGF vectors for RHBM analysis
%       SCOIL structure
%           name - name of coil model
%           type - kind of coil model, wire of surface
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           ...
%       SOL structure
%           Jsol - body equivalent volumetric currents
%           Esol - electric fields
%           Bsol - magnetic fields
%           Ssol - local unweighted SAR
%           Gsar - global SAR
%           Pabs - absorved power
%           Zparam - port Z-parameters of coil
%           freq - frequency points of solutions
%       PUSE structure (same as solution, but with pulse applied)
%           Jcoil - coil current related basis functions
%           Jsol - body equivalent volumetric currents
%           Esol - electric fields
%           Bsol - magnetic fields
%           Ssol - local unweighted SAR
%           Gsar - global SAR
%           Pabs - absorved power
%           Ipulse - currents at the ports
%           Vpulse - voltages at the ports
%           freq - frequency points of solutions
%
%
% _________________________________________________________________________

clc;
clear
close all;
p = genpath(pwd);
addpath(p);


% initialize RHBM global variable
global RHBM;
RHBM.name = 'No Model Selected';
RHBM.r = [];
RHBM.epsilon_r = [];
RHBM.sigma_e = [];
RHBM.rho = [];
RHBM.idxS = [];

% initialize COIL global variable
global COIL;
COIL.name = 'No Model Selected';
COIL.type = 'N';

% initialize solution
global SOL;
SOL.Jsol = [];
SOL.Esol = [];
SOL.Bsol = [];
SOL.Ssol = [];
SOL.Gsar = [];
SOL.Pabs = [];
SOL.Zparam = [];
SOL.freq = [];

% initialize pulse solution
global PULSE;
PULSE.Jcoil = [];
PULSE.Jsol = [];
PULSE.Esol = [];
PULSE.Bsol = [];
PULSE.Ssol = [];
PULSE.Gsar = [];
PULSE.Pabs = [];
PULSE.Ipulse = [];
PULSE.Vpulse = [];
PULSE.freq = [];


% initialize figure counter
global FIGIDX;
FIGIDX = 0;

MARIE_MAIN;