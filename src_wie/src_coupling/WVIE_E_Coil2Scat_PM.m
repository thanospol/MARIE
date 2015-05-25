function [Eout] = WVIE_E_Coil2Scat_PM(Scoord,Ccoil,Dcoil,freq,Jc)
%%    Coil to Scatterer PM E field generation for the WIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING E field
%   due to WIE Jc current coefficients
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       Ccoil - coordinates of the centers of coil segments 
%       Dcoil - Cartesian components of the segments
%       freq - frequency
%       Jc - coil current coefficients
%
%
%% Output
%       Eout - E field 3xNo vector with the contribution of all segments
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
% Define variables and allocate space
% -------------------------------------------------------------------------

No = size(Scoord,1); % number of observation points

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

J = Dcoil; % scale the current by the length of each segment
for ii = 1:3
    J(:,ii) = J(:,ii).*Jc; % for each component
end

[Eout] = E_field_DGF(J,Ccoil,Scoord,ko); % compute the point-based fields

% -------------------------------------------------------------------------
% done
% -------------------------------------------------------------------------

Eout = reshape(Eout,3*No,1);


