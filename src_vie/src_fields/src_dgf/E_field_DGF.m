function [E] = E_field_DGF(J,Icoord,Ocoord,ko)
%%   Generation of E field via dyadic Green functions 
% _________________________________________________________________________
%
%       Evaluates the dyadic Green function
%       due to dipoles with current J in positions Icoord
%       to generate the E field in positions Ocoord
%
% _________________________________________________________________________
%
%% INPUT
%   J:      current dipoles (Ni x 3)
%   Icoord: vector with positions of current dipoles (Ni x 3)
%   Ocoord: vector with positions of observation points (No x 3)
%   ko:     wavelength number in the medium
%
%
%% OUTPUT
%   E:      Electric field (No x 3)
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

% EM constants
mu = 4*pi*1e-7;
co = 299792458;
lambda = 2*pi/ko;
f  = co/lambda;
omega = 2*pi*f;

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);

% allocate space for output field
E = zeros(No,3); % x coomponent, y component, z component

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

% -------------------------------------------------------------------------
% precompute multipliers
% -------------------------------------------------------------------------

% compute chi
chi = 1j*omega*mu*exp(-1j*ko*R)./(4*pi*ko*ko*R3);

% compute P and Q
P = 1j*ko*R + 1;
Q = ko*ko*R2 - P;
P = (Q - 2*P)./(R2);
Q = chi.*Q;
P = chi.*P;

clear chi; clear R; clear R2; clear R3;

% -------------------------------------------------------------------------
% compute the Solution
% -------------------------------------------------------------------------
 
% x component
E(:,1) = (P.*X.*X - Q)*J(:,1) + (P.*(X.*Y))*J(:,2) + (P.*(X.*Z))*J(:,3);

% y component
E(:,2) = (P.*(X.*Y))*J(:,1) + (P.*(Y.*Y) - Q)*J(:,2) + (P.*(Y.*Z))*J(:,3);

% z component
E(:,3) = (P.*(X.*Z))*J(:,1) + (P.*(Y.*Z))*J(:,2) + (P.*(Z.*Z) - Q)*J(:,3);

