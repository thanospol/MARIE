function [E,H] = CCL_Excitation(r,ko,Jc,cclfile)
%%    Function to Generate a Constant Current Loop excitation
% _________________________________________________________________________
%
%       Generates the E and H fields of an constant current loop
%       with Jc current amplitude in the loop given by .wsd file filename
% _________________________________________________________________________
%
%% INPUT
%   r               4D (LxMxNx3) array with domain voxelized grid coordinates
%   ko              wave number
%   Jc              current amplitude in the loop
%   cclfile         filename of the .wsd with the loop coordinates  
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
Dcoord = reshape(r,L*M*N,3);


% -------------------------------------------------------------------------
% Load loop coordinates as ASCii
% -------------------------------------------------------------------------

CCLdata = load(cclfile, '-ascii');

% -------------------------------------------------------------------------
%  prepare the coil current in the required format 
% -------------------------------------------------------------------------

Pcoil = CCLdata(:,1:3);
Ncoil = CCLdata(:,4:6);

clear CCLdata;

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

% [E] = E_field_DGF(Jcc,Ccoil,Dcoord,ko);
% [H] = H_field_DGF(Jcc,Ccoil,Dcoord,ko);
[E,H] = eval_DGF(Jcc,Ccoil,Dcoord,ko);

% in the whole domain, now reshape them to LMN3 size
E = reshape(E,L,M,N,3);
H = reshape(H,L,M,N,3);

