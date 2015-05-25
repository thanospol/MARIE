function [Zbc] = Assembly_SCOUP_PM(Scoord,index,etod,Ct,Ln,Pn,freq) 
%%    PM coupling for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the POINT MATCHING Coupling SIE to VIE
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
%
%
%% Output
%       Zbc - Tensor (No x 3 x Nd) with the contribution of each edge of the element
%              Zbc(1000,3,200) is z contribution of 200-th edge to 1000-th element
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
Ne = size(etod,2); % number of elements
Nd = max(index); % number of dofs

Zbc = zeros(No,3,Nd);

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------


for ii = 1:Ne % loop on the elements
    
    %     r_1 = node(:,elem(1,ii)); % 3xNe with coordinates of the first node of all elements
    %     r_2 = node(:,elem(2,ii)); % 3xNe with coordinates of the first node of all elements
    %     r_3 = node(:,elem(3,ii)); % 3xNe with coordinates of the first node of all elements
    
    % get the sign and index for the contribution
    absnum = abs(etod(:,ii)); % internal index of the edge
    mult = etod(:,ii)./absnum; % +1 or -1
    
    idx1 = index(absnum(1));
    idx2 = index(absnum(2));
    idx3 = index(absnum(3));
       
    
    % obtain the contribution for each element
    Ecoup = E_DGF_Triang(Scoord,Ct(:,ii),Ln(:,ii),Pn(:,:,ii),ko);
    
    %     GJ1 = squeeze(Gcoup(:,:,1));
    %     GJ2 = squeeze(Gcoup(:,:,2));
    %     GJ3 = squeeze(Gcoup(:,:,3));
    %
    %     GGJ1(ii,:,:) = GJ1;
    %     GGJ2(ii,:,:) = GJ2;
    %     GGJ3(ii,:,:) = GJ3;
    
    
    % add each contribution multiplied by the sign given by etod
    
    if (idx1)
        Zbc(:,:,idx1) = Zbc(:,:,idx1) + mult(1)*Ecoup(:,:,1); % first edge
    end
    
    if (idx2)
        Zbc(:,:,idx2) = Zbc(:,:,idx2) + mult(2)*Ecoup(:,:,2); % second edge
    end
    
    if (idx3)
        Zbc(:,:,idx3) = Zbc(:,:,idx3) + mult(3)*Ecoup(:,:,3); % third edge
    end
    
    
end


% -------------------------------------------------------------------------
%             and it is done: report 
% -------------------------------------------------------------------------
    