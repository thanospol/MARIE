function [Eout] = SVIE_E_Scat2Coil_PM_par(Scoord,index,etod,Ct,Ln,Pn,freq,Jb)
%%    Scatterer to Coil PM E field generation for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING E field
%   due to VIE Jb current coefficients
%   (Jb values should already include the Gram scale, i.e. be integrated)
%   Parallel version (parfor)
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       Ct - coordinates of the center of the triangle
%       Ln - values of the length of each side of the triangle
%       Pn - 3x3 matrix with coordinates of the rho vectors (Pn(:,1) == rho_1)
%       freq - frequency
%       Jb - VIE currents (values should be already integrated over vol)
%
%
%% Output
%       Eout - E field Nd vector with the contribution to all edge elements
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

Ne = size(etod,2); % Number of elements
Nd = max(index); % number of dofs
No = size(Scoord,1); % number of observation points

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

Vout = zeros(3,Ne);

parfor ii = 1:Ne % loop on the elements
    
    % obtain the contribution for each element
    Gcoup = E_DGF_Triang(Scoord,Ct(:,ii),Ln(:,ii),Pn(:,:,ii),ko);
   
    % add each contribution multiplied by the sign given by etod
    % and multiplied by the current of the scatterer voxel
    
    Gcoup = reshape(Gcoup,3*No,3); % reshape to fit dimensions
    
    Vout(:,ii) = (Gcoup.'*Jb); % contribution of the 3 edges
    
end

% map the element contribution to the corresponding dof
Eout = zeros(Nd,1);

for ii = 1:Ne
    
    % get the sign and index for the contribution
    absnum = abs(etod(:,ii)); % internal index of the edge
    mult = etod(:,ii)./absnum; % +1 or -1
    
    idx1 = index(absnum(1));
    idx2 = index(absnum(2));
    idx3 = index(absnum(3));
    
    if (idx1)
        Eout(idx1,1) = Eout(idx1,1) + mult(1)*Vout(1,ii); % first edge
    end
    
    if (idx2)
        Eout(idx2,1) = Eout(idx2,1) + mult(2)*Vout(2,ii); % second edge
    end
    
    if (idx3)
        Eout(idx3,1) = Eout(idx3,1) + mult(3)*Vout(3,ii); % third edge
    end
    
end

% -------------------------------------------------------------------------
% done
% -------------------------------------------------------------------------


