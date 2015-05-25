function [Eout] = SVIE_E_Coil2Scat_PM(Scoord,index,etod,Ct,Ln,Pn,freq,Jc)
%%    Coil to Scatterer PM E field generation for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING E field
%   due to SIE Jc current coefficients
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
%       Jc - coil current coefficients
%
%
%% Output
%       Eout - E field 3xNo vector with the contribution of all edge elements
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
Ne = size(etod,2); % Number of elements

Eout = zeros(No,3);

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------


for ii = 1:Ne % loop on the elements
    
    % get the sign and index for the contribution
    absnum = abs(etod(:,ii)); % internal index of the edge
    mult = etod(:,ii)./absnum; % +1 or -1
    
    idx1 = index(absnum(1));
    idx2 = index(absnum(2));
    idx3 = index(absnum(3));
    
    % obtain the contribution for each element
    Ecoup = E_DGF_Triang(Scoord,Ct(:,ii),Ln(:,ii),Pn(:,:,ii),ko);
   
    % add each contribution multiplied by the sign given by etod
    % and multiplied by the current of the DOF
    
    if (idx1)
        Eout = Eout + Jc(idx1)*mult(1)*Ecoup(:,:,1); % first edge
    end
    
    if (idx2)
        Eout = Eout + Jc(idx2)*mult(2)*Ecoup(:,:,2); % second edge
    end
    
    if (idx3)
        Eout = Eout + Jc(idx3)*mult(3)*Ecoup(:,:,3); % third edge
    end
    
end


% -------------------------------------------------------------------------
% done
% -------------------------------------------------------------------------

Eout = reshape(Eout,3*No,1);


