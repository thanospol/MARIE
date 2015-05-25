function [E,H] = eval_DGF(J,Icoord,Ocoord,ko)
%%   Function to evaluate the dyadic Green function due to dipole currents
% _________________________________________________________________________
%
%       Evaluates the dyadic Green function
%       due to dipoles with current J in positions Icoord
%       to generate the E and H fields in positions Ocoord
%       The function allows multiple (Ne) excitations
%
% _________________________________________________________________________
%
%% INPUT
%   J:      current dipoles (Ni x 3 x Ne)
%   Icoord: vector with positions of current dipoles (Ni x 3)
%   Ocoord: vector with positions of observation points (No x 3)
%   ko:     wavelength number in the medium
%
%
%% OUTPUT
%   E:      Electric field (No x 3 x Ne)
%   H:      Magnetic field (No x 3 x Ne)
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

% verify the dimensions of J
dim = size(J);
if (length(dim) == 2)
    % if Ne = 1, resize J to be a 3 dimensional vector
    % at the end E and H will be squeezed to 2D
    Ne = 1;
    J = reshape(J,dim(1),dim(2),1);
end

% obtain dimensions
Ni = size(Icoord,1);
No = size(Ocoord,1);

% allocate space for output field
Hx = zeros(No,Ne); 
Ex = zeros(No,Ne); 
Hy = zeros(No,Ne);
Ey = zeros(No,Ne);
Hz = zeros(No,Ne);
Ez = zeros(No,Ne);

J1 = zeros(1,Ne);
J2 = zeros(1,Ne);
J3 = zeros(1,Ne);

% -------------------------------------------------------------------------
% loop on the Icoord and apply the DGF
% -------------------------------------------------------------------------

for jj = 1:Ni
    
    % ---------------------------------------------------------------------
    % Get distance vectors
    
    X = Ocoord(:,1) - Icoord(jj,1);
    Y = Ocoord(:,2) - Icoord(jj,2);
    Z = Ocoord(:,3) - Icoord(jj,3);
         
    % distance 3D
    R2 = X.*X + Y.*Y + Z.*Z;
    R = sqrt(R2);
    R3 = R.*R2;
    
    % ---------------------------------------------------------------------
    % For E DGF
    % compute chi
    % chi = 1j*eta*dV*exp(-1j*ko*R)./(4*pi*ko*R3);
    chi = 1j*omega*mu*exp(-1j*ko*R)./(4*pi*ko*ko*R3);
    
    % compute P and Q
    P = 1j*ko*R + 1;
    Q = ko*ko*R2 - P;
    P = (Q - 2*P)./(R2);
    Q = chi.*Q;
    P = chi.*P;
    
    clear chi; clear R2;
    
    % ---------------------------------------------------------------------
    % For H DGF
    % precompute value
    const = 1j*ko*R;
    
    % compute mult
    mult = exp(-const)./(4*pi*R3);
    mult = mult.*(const+1);
    
    clear const; clear R; clear R3;
    
    % ---------------------------------------------------------------------
    % compute the Solution

    % get components of Jin
    J1(1,:) = J(jj,1,:);
    J2(1,:) = J(jj,2,:);
    J3(1,:) = J(jj,3,:);

    % ---------------------------------------------------------------------
    % E field    
    
    % x component
    Ex = Ex + (P.*X.*X - Q)*J1 + (P.*(X.*Y))*J2 + (P.*(X.*Z))*J3;
    % y component
    Ey = Ey + (P.*(X.*Y))*J1 + (P.*(Y.*Y) - Q)*J2 + (P.*(Y.*Z))*J3;
    % z component
    Ez = Ez + (P.*(X.*Z))*J1 + (P.*(Y.*Z))*J2 + (P.*(Z.*Z) - Q)*J3;
    
    % ---------------------------------------------------------------------
    % H field
    
    X = mult.*X;
    Y = mult.*Y;
    Z = mult.*Z;
    
    % x component
    Hx = Hx + Z*J2 - Y*J3;
    % y component
    Hy = Hy - Z*J1 + X*J3;
    % z component
    Hz = Hz + Y*J1 - X*J2;
    
end

% -------------------------------------------------------------------------
% Arrange output data
% -------------------------------------------------------------------------

% store in vector form with Ne columns
E = [Ex; Ey; Ez];
clear Ex; clear Ey; clear Ez;
H = [Hx; Hy; Hz];
clear Hx; clear Hy; clear Hz;

% reshape into Nox3xNe and squeeze in case Ne == 1
E = reshape(E,No,3,Ne);
E = squeeze(E);
H = reshape(H,No,3,Ne);
H = squeeze(H);


