function [Hout] = M_DGF_Triang(Ocoord,Ct,Ln,Pn,ko)
%%    H field DGF evaluation from RWG basis coefficients
% _________________________________________________________________________
%
%   Evaluate the H field Dyadic Green Function for each 
%   triangle element in a given set of observation points
%
% _________________________________________________________________________
%
%% Input
%       Ocoord - coordinates of the observation points (No x 3)
%       Ct - coordinates of the center of the triangle
%       Ln - values of the length of each side of the triangle
%       Pn - 3x3 matrix with coordinates of the rho vectors (Pn(:,1) == rho_1)
%       ko - wave number
%
%
%% Output
%       Hout - Tensor (No x 3 x 3) with the H field contribution of each edge of the element
%              Hout(:,3,2) is z component contribution of second edge 
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% obtain dimensions
No = size(Ocoord,1);

% allocate space for output field
Hout = zeros(No,3,3); % x coomponent, y component, z component, for each 


% Get distance vectors to the center of the triangle
X = Ocoord(:,1) - Ct(1);
Y = Ocoord(:,2) - Ct(2);
Z = Ocoord(:,3) - Ct(3);

% distance 3D
R2 = X.*X + Y.*Y + Z.*Z;
R = sqrt(R2);
R3 = R.*R2;

    
% For H DGF
% precompute value
const = 1j*ko*R;

% compute mult
mult = exp(-const)./(4*pi*R3);
mult = mult.*(const+1);

clear R3; clear const; clear R;

% multiplied distances for H field
Xh = mult.*X;
Yh = mult.*Y;
Zh = mult.*Z;


% compute the Solution

for ii = 1:3 % loop on the contribution of each edge 

    % get the function multiplying depending on the edge number
    Fn = Pn(:,ii)*Ln(ii)/2;

    % H field
    
    % x component ii element (multiplies by the 3 components of Fn of the ii edge)
    Hout(:,1,ii) = Zh*Fn(2) - Yh*Fn(3);
    
    % y component ii element (multiplies by the 3 components of Fn of the ii edge)
    Hout(:,2,ii) = -Zh*Fn(1) + Xh*Fn(3);
    
    % z component ii element (multiplies by the 3 components of Fn of the ii edge)
    Hout(:,3,ii) = Yh*Fn(1) - Xh*Fn(2);
    
    
end
    
    
% % % x component first element (multiplies by the 3 components of Pn of the first edge)
% % Gout(:,1,1) = (P.*X.*X - Q)*Pn(1,1)*Ln(1)/2 + (P.*(X.*Y))*Pn(2,1)*Ln(1)/2 + (P.*(X.*Z))*Pn(3,1)*Ln(1)/2;
% % % x component second element (multiplies by the 3 components of Pn of the second edge)
% % Gout(:,1,2) = (P.*X.*X - Q)*Pn(1,2)*Ln(2)/2 + (P.*(X.*Y))*Pn(2,2)*Ln(2)/2 + (P.*(X.*Z))*Pn(3,2)*Ln(2)/2;
% % % x component third element (multiplies by the 3 components of Pn of the third edge)
% % Gout(:,1,3) = (P.*X.*X - Q)*Pn(1,3)*Ln(3)/2 + (P.*(X.*Y))*Pn(2,3)*Ln(3)/2 + (P.*(X.*Z))*Pn(3,3)*Ln(3)/2;


