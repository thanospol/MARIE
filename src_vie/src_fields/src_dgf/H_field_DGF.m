function [H] = H_field_DGF(J,Icoord,Ocoord,ko)
%%   Generation of H field via dyadic Green functions 
% _________________________________________________________________________
%
%       Evaluates the dyadic Green function
%       due to dipoles with current J in positions Icoord
%       to generate the H field in positions Ocoord
%
% _________________________________________________________________________
%
%% INPUT
%   J:      current dipoles (Ni x 3)
%   Icoord: vector with positions of current dipoles (Ni x 3)
%   Ocoord: vector with positions of observation points (No x 3)
%   ko:     wave number in the medium
%
%
%% OUTPUT
%   H:      Electric field (No x 3)
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
% Prepare data
% -------------------------------------------------------------------------

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);

% allocate space for output field
H = zeros(No,3); % x coomponent, y component, z component

% reshape J into 3 components just in case
J = reshape(J,Ni,3); % x coomponent, y component, z component

% -------------------------------------------------------------------------
% Get distance vectors
% -------------------------------------------------------------------------

X = zeros(No,Ni);
Y = zeros(No,Ni);
Z = zeros(No,Ni);

% loop to get the distances
for jj = 1:Ni
    X(:,jj) = Ocoord(:,1) - Icoord(jj,1);
    Y(:,jj) = Ocoord(:,2) - Icoord(jj,2);
    Z(:,jj) = Ocoord(:,3) - Icoord(jj,3);
end

% distance 3D
R2 = X.*X + Y.*Y + Z.*Z;
R = sqrt(R2);
R3 = R.*R2;
clear R2;

% -------------------------------------------------------------------------
% precompute multipliers
% -------------------------------------------------------------------------

% precompute value
const = 1j*ko*R;

% compute mult
mult = exp(-const)./(4*pi*R3);
mult = mult.*(const+1);

clear const; clear R; clear R2; clear R3;

X = mult.*X;
Y = mult.*Y;
Z = mult.*Z;

% -------------------------------------------------------------------------
% compute the Solution
% -------------------------------------------------------------------------

% x component
H(:,1) = Z*J(:,2) - Y*J(:,3);

% y component
H(:,2) = -Z*J(:,1) + X*J(:,3);

% z component
H(:,3) = Y*J(:,1) - X*J(:,2);


