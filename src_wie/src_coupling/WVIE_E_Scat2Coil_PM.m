function [Eout] = WVIE_E_Scat2Coil_PM(Scoord,Ccoil,Dcoil,freq,Jb)
%%    Scatterer to Coil PM E field generation for the WIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING E field
%   due to VIE Jb current coefficients
%   (Jb values should already include the Gram scale, i.e. be integrated)
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       Ccoil - coordinates of the centers of coil segments 
%       Dcoil - Cartesian components of the segments
%       freq - frequency
%       Jb - VIE currents (values should be already integrated over vol)
%
%
%% Output
%       Eout - E field Nd vector with the contribution to all coil elements
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
%            Define EM constants
% -------------------------------------------------------------------------

co = 299792458;
lambda  = co/freq;
ko = 2*pi/lambda;

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

% compute field
[Eout] = E_field_DGF(Jb,Scoord,Ccoil,ko);

% scale by segment component and add
Eout = sum(Dcoil.*Eout,2);

% -------------------------------------------------------------------------
% done
% -------------------------------------------------------------------------


