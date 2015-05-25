function [E,H] = SIE_Radiate(SCOIL,Jc,freq,r)
%%    Compute fields due to the solution to the SIE problem
% _________________________________________________________________________
%
%   Applies operators to generate the fields and other figures of merit
%
% _________________________________________________________________________
%
%% INPUT
%       SCOIL structure
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           ...
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
%                 define EM vars and constants
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);

% x,y,z coordinates of the domain
xd = r(:,:,:,1);
yd = r(:,:,:,2);
zd = r(:,:,:,3);

% Coordinates of all the domain
Dcoord = [xd(:), yd(:), zd(:)];
clear xd; clear yd; clear zd;


% -------------------------------------------------------------------------
%         Compute incident fields due to Jc
% -------------------------------------------------------------------------


LEVEL_DVrule = 5;

if (LEVEL_DVrule == 1)
    
    % Point matching for the E and H field
    [E,H] = SVIE_EM_Coil2Scat_PM_par(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.Ct,SCOIL.Ln,SCOIL.Pn,freq,Jc);
    
else
    
    % We use quadrature for the E field
    if (exist('ompQuadCoil2Scat', 'file') == 3)
        % we used the PARALLEL MEXED Quadrature based function...
        [E] = SVIE_E_Coil2Scat_QMEX_par(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,LEVEL_DVrule,Jc); % quadrature PARALLEL MEXED function
        % [E] = SVIE_E_Coil2Scat_QMEX(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,LEVEL_DVrule,Jc); % quadrature MEXED function
    else
        % we used the PARALLEL NOT MEXED Quadrature based function...
        [E] = SVIE_E_Coil2Scat_Q_par(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,LEVEL_DVrule,Jc); % quadrature PARALLEL NOT MEXED function
        % [E] = SVIE_E_Coil2Scat_Q(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.node,SCOIL.elem,freq,LEVEL_DVrule,Jc); % quadrature not mexed function
    end
    
    % Point matching for the H field
    [H] = SVIE_M_Coil2Scat_PM_par(Dcoord,SCOIL.index,SCOIL.etod,SCOIL.Ct,SCOIL.Ln,SCOIL.Pn,freq,Jc);
    
end


% -------------------------------------------------------------------------
%                 And it is done
% -------------------------------------------------------------------------

E = reshape(E,L,M,N,3);
H = reshape(H,L,M,N,3);

