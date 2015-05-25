function [E,H] = WIE_Radiate(WCOIL,Jc,freq,r)
%%    Compute fields due to the solution to the VIE problem
% _________________________________________________________________________
%
%   Applies operators to generate the fields and other figures of merit
%
% _________________________________________________________________________
%
%% INPUT
%       WCOIL structure
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           Dwire - diameter of wire
%           Rhowire - resistivity of material
%           port - port definition
%       Jc - currents in the coil
%       freq - frequency
%       r - domain
%
%
%% OUTPUT
%       E        Solution electric field (LxMxNx3)
%       H        Solution magnetic field (LxMxNx3)
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

[L,M,N,~] = size(r);

% x,y,z coordinates of the domain
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

% %   Scoord are the coordinates of the RHBM points
% Scoord = [xd(idxS), yd(idxS), zd(idxS)];
%  Dcoord are the coordinates of the whole domain
Dcoord = [xd(:), yd(:), zd(:)];
clear xd; clear yd; clear zd;

% -------------------------------------------------------------------------
%         Compute incident fields due to Jc
% -------------------------------------------------------------------------

Ccoil = (WCOIL.Pcoil + WCOIL.Ncoil)/2; % center of coil segments
Dcoil = WCOIL.Pcoil - WCOIL.Ncoil; % length of each segment

J = Dcoil; % scale the current by the length of each segment
for ii = 1:3
    J(:,ii) = J(:,ii).*Jc; % for each component
end

% [E] = E_field_DGF(J,Ccoil,Dcoord,ko);
% [H] = H_field_DGF(J,Ccoil,Dcoord,ko);
[E,H] = eval_DGF(J,Ccoil,Dcoord,ko);

% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

% now reshape them to LMN3 size
E = reshape(E,L,M,N,3);
H = reshape(H,L,M,N,3);



