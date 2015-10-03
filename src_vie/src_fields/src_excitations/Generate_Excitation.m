function [E,H] = Generate_Excitation(r,freq,Exctype,Pos,Comp)
%%    Function to Generate the fields due to a given excitation
% _________________________________________________________________________
%
%       Generates the E and H fields due to a given excitation
% _________________________________________________________________________
%
%% INPUT
%   r               4D (LxMxNx3) array with domain voxelized grid coordinates
%   freq            frequency
%   Exctype         'D' for dipole, 'P' for plane wave
%   Pos             position for dipole, direction for PW, .wsd file for Loop
%   Comp            dipole components, PW polarization, or loop current
%
%% OUTPUT
%   E               Electric field (LxMxNx3)
%   H               Magnetic field (LxMxNx3)
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
%

% -------------------------------------------------------------------------
% Obtain constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
omega = 2*pi*freq;
lambda = co/freq;
ko = 2*pi/lambda;
omega_mu = omega*mu;

% -------------------------------------------------------------------------
% Check the type
% -------------------------------------------------------------------------

switch Exctype
    
    case 'D'
        % Dipole
        [E,H] = Dipole_Excitation(r,ko,Comp,Pos);
        
    case 'P'
        % Plain wave
        k = ko*Pos;
        [E,H] = PlaneWave_Excitation(r,k,omega_mu,Comp);
        
    case 'L'
        % Constant Current Loop
        [E,H] = CCL_Excitation(r,ko,Comp,Pos);
        
    otherwise
        % zero fields
        E = 0*r;
        H = 0*r;
        
end
        

