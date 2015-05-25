function [E,H] = Dipole_Excitation(r,k,DipJ,DipCoord)
%%    Function to Generate a Hertzian Dipole excitation
% _________________________________________________________________________
%
%       Generates the E and H fields of an Hertzian dipole
%       with DipJ components in the position given by DipCoord
% _________________________________________________________________________
%
%% INPUT
%   r               4D (LxMxNx3) array with domain voxelized grid coordinates
%   k               vector with wavenumber in the medium
%   DipJ            3 Cartesian components of the Hertzian dipole
%   DipCoord        Coordinates of the position of the Hertzian dipole  
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
% Obtain coordinates of the domain
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);
Ocoord = reshape(r,L*M*N,3);

% -------------------------------------------------------------------------
% Generate E field
% -------------------------------------------------------------------------

[E] = E_field_DGF(DipJ,DipCoord,Ocoord,k);
E = reshape(E,L,M,N,3);

% -------------------------------------------------------------------------
% Generate H field
% -------------------------------------------------------------------------

[H] = H_field_DGF(DipJ,DipCoord,Ocoord,k);
H = reshape(H,L,M,N,3);


