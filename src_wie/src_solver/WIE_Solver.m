function [ZP,Jc] = WIE_Solver(WCOIL,freq)
%%    Solver Routine for the WIE formulation
% _________________________________________________________________________
%
%   Solves the thin wire IE problem to obtain the Z parameters
%   Prepared for gapless structures: computes Y-parameters and inverts
%   based on the BEM approach on Wire Antennas (Harrington)
%       R.F. Harrington
%       Matrix Methods for Field Problems
%       Proc. IEEE 55(2): 136 - 149, Feb. 1967
%
% _________________________________________________________________________
%
%
%% INPUT
%   WCOIL structure
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           Dwire - diameter of wire
%           Rhocoil - resistivity of material
%           port - port definition
%   freq - frequency (Hz)
%
%
%% OUTPUT
%   ZP - Z-parameter matrix
%   Jc - Coil current coefficients
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________



% -------------------------------------------------------------------------
% Generate the WIE system

[Zo,F] = Assembly_WIE(WCOIL,freq);


% -------------------------------------------------------------------------
% form
%       YP = F.'*inv(Zo)*F;
%       [ZP] = inv(YP);

tic_i = tic;

Jc = Zo\F;

% Generate the Z parameter matrix
ZP = inv(F.'*Jc);
ZP = (ZP + ZP.')/2;

fprintf(1,'\n WIE system solve done,  Elapsed time  = %.2f [sec]' , toc(tic_i));
fprintf(1,'\n ----------------------------------------------------------\n');
