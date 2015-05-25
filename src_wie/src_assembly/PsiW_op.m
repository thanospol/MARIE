function [Psi] = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k)
%%    Kernel operation for computing the MoM system of a thin wire
% _________________________________________________________________________
%
%   Def.: Compute the integral values (Psi) for the elements.
%   Formulas for the computation of the BEM value between two wire 
%   elements, based on the BEM approach on Wire Antennas (Harrington)
%       R.F. Harrington
%       Matrix Methods for Field Problems
%       Proc. IEEE 55(2): 136 - 149, Feb. 1967
%
% _________________________________________________________________________
%
%
%% INPUT
%   PPx is the x coord of positive end of the source segment
%   PPy is the y coord of positive end of the source segment
%   PPz is the z coord of positive end of the source segment
%   CPx is the x coord of center point of the source segment
%   CPy is the y coord of center point of the source segment
%   CPz is the z coord of center point of the source segment
%   OPx is the x coord of observation point
%   OPy is the y coord of observation point
%   OPz is the z coord of observation point
%   DM is a matrix with segment lengths
%   a is the wire radius
%   k is the wavelength number
%
%
%% OUTPUT
%   Psi is the value of the integral
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

lambda = 2*pi/k;

Psi = 0*PPx;

% % %%% --- MODIFICATION ----- 
% % Dval = mean(mean(DM)); % use mean segment length
% % %%% --- MODIFICATION ----- 
    
% ----------------------------------------
% Apply the rotation for the computation
% ----------------------------------------

% Assume Center in origin, with z axis crossing from Center to Positive
PPx = PPx - CPx; PPy = PPy - CPy; PPz = PPz - CPz;
OPx = OPx - CPx; OPy = OPy - CPy; OPz = OPz - CPz;


% apply the firt rotation in the xy plane
phi = -angle(PPx + 1j*PPy); % obtain the angle

PPxnew = cos(phi).*PPx - sin(phi).*PPy; % apply rotation to x coord of P
PPynew = sin(phi).*PPx + cos(phi).*PPy; % apply rotation to y coord of P

OPxnew = cos(phi).*OPx - sin(phi).*OPy; % apply rotation to x coord of O
OPynew = sin(phi).*OPx + cos(phi).*OPy; % apply rotation to y coord of O

% apply the second rotation in the xz plane
phi = -angle(PPz + 1j*PPxnew); % obtain the angle

PPx = cos(phi).*PPxnew + sin(phi).*PPz; % apply rotation to x coord of P
PPy = PPynew; % apply rotation to y coord of P
PPz = -sin(phi).*PPxnew + cos(phi).*PPz; % apply rotation to z coord of P

OPx = cos(phi).*OPxnew + sin(phi).*OPz; % apply rotation to x coord of O
OPy = OPynew; % apply rotation to y coord of O
OPz = -sin(phi).*OPxnew + cos(phi).*OPz; % apply rotation to z coord of O

% check the validity of the rotations

if ~isempty( find( (abs(PPx) + abs(PPy)) > 1e-10 ))
    fprintf(1, '\n Error in the rotation!\n');
end


% obtain values for the given new coordinates
Z = OPz;

Alpha = abs(PPz);


% % %%% --- MODIFICATION ----- 
% % Alpha = abs(PPz./PPz)*Dval/2; % use mean length of segment
% % %%% --- MODIFICATION ----- 


% ----------------------------------------
% compute matrix with distances
% ----------------------------------------

Dist = sqrt((OPx).^2 + (OPy).^2 + (OPz).^2);

% ----------------------------------------
% Same element
% ----------------------------------------
% 
% find elements that are the same
% [idx] = find(Dist <= 1e-8);
% 
% % apply operator
% % from Field Computation by Moment Methods,(book) Harrington, Sec.4.3
% % eq. 4-25
% 
% Psi(idx) = log(DM(idx)/a)./(2*pi*DM(idx))-1j*k/(4*pi);
% 
% % %%% --- MODIFICATION ----- 
% % Psi(idx) = log(Dval/a)./(2*pi*Dval)-1j*k/(4*pi); % use mean value
% % %%% --- MODIFICATION ----- 


% ----------------------------------------
% Close elements: near ineteraction
% ----------------------------------------

% find elements that are the same
% [idx] = find((Dist > 1e-8) & (Dist < lambda/5));
[idx] = find((Dist < lambda/5));

% apply operator
% from Matrix Methods for Field Problems, Harrington, Proc. IEEE
% 55(2): 136 - 149, Feb. 1967

if ~isempty(idx)
    
    Rho = sqrt(OPx.^2 + OPy.^2);
    [widx] = find( Rho < a );
    Rho(widx) = a; %Rho(widx) + a;
    R = sqrt(Rho.^2 + Z.^2);

    % compute the components I
    I1 = log( (Z(idx) + Alpha(idx) + sqrt(Rho(idx).^2 + (Z(idx)+Alpha(idx)).^2)) ./ (Z(idx) - Alpha(idx) + sqrt(Rho(idx).^2 + (Z(idx)-Alpha(idx)).^2)) );
    I2 = 2*Alpha(idx);
    I3 = (Alpha(idx) + Z(idx)) .* sqrt(Rho(idx).^2 + (Z(idx)+Alpha(idx)).^2) + (Alpha(idx) - Z(idx)) .* sqrt(Rho(idx).^2 + (Z(idx)-Alpha(idx)).^2) + I1.*(Rho(idx).^2);
    I3 = I3/2;
    I4 = 2*Alpha(idx).*(Rho(idx).^2) + (2*(Alpha(idx).^3) + 6*Alpha(idx).*(Z(idx).^2))/3;

    % compute the value of the integral Psi
    Psi(idx) = I1 - 1j*k*(I2 - R(idx).*I1) - (I3 - 2*R(idx).*I2 + I1.*(R(idx).^2)).* k^2 /2 + (I4 - 3*R(idx).*I3 + 3*I2.*(R(idx).^2) - I1.*(R(idx).^3))*1j*k^3 /6;
    Psi(idx) = Psi(idx).*exp(-1j*k.*R(idx))./(8*pi*Alpha(idx));

end
    
    
% ----------------------------------------
% Far elements: far ineteraction
% ----------------------------------------

% find elements that are in the far interaction
[idx] = find(Dist >= lambda/10);

% apply operator
% from Matrix Methods for Field Problems, Harrington, Proc. IEEE
% 55(2): 136 - 149, Feb. 1967

if ~isempty(idx)
    
    R = Dist;

    % compute the components A
    ZoverR = Z(idx)./R(idx);
    AlphaoverR = Alpha(idx)./R(idx);

    A0 = 1 + (AlphaoverR).^2 .* (3*(ZoverR).^2 - 1)/6 + (AlphaoverR).^4 .* ( 3 - 30*(ZoverR).^2 + 35*(ZoverR).^4 )./40;

    A1 = (AlphaoverR) .* (3*(ZoverR).^2 - 1)/6 + (AlphaoverR).^3 .* ( 3 - 30*(ZoverR).^2 + 35*(ZoverR).^4 )./40;

    A2 = -(ZoverR).^2 /6  - (AlphaoverR).^2 .*( 1 - 12*(ZoverR).^2 + 15*(ZoverR).^4 )./40;

    A3 = (AlphaoverR) .* ( 3*(ZoverR).^2 - 5*(ZoverR).^4 )./60;

    A4 = (ZoverR).^4 ./ 120;

    % compute the value of the integral Psi
    Psi(idx) = A0  + 1j*k*Alpha(idx).*A1 + ((k*Alpha(idx)).^2) .* A2 + 1j* ((k*Alpha(idx)).^3) .* A3 + ((k*Alpha(idx)).^4) .* A4;
    Psi(idx) = Psi(idx).*(exp(-1j*k*R(idx))./(4*pi*R(idx)));

end


