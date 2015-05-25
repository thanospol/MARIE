function [Eout] = SVIE_E_Scat2Coil_Q_par(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,Jb)
%%    Scatterer to Coil Quadrature E field generation for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature-based E field
%   due to VIE Jb current coefficients
%   (Jb values should already include the Gram scale, i.e. be integrated)
%   Goes triangle by triangle, and calls the mexCoupling function for each
%   PARFOR version
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%       freq - frequency
%       LEVEL_DVrule - number of the quadrature rule (1 gives point matching)
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

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2*pi*freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
% eta = omega*mu/ko; %3.767303134617706e+002; 

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

No = size(Scoord,1); % number of observation points
Ne = size(elem,2); % number of elements
Nd = max(index); % number of dofs

Eout = zeros(Nd,1);


% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % enough
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule ( LEVEL_DVrule );

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

Jb = reshape(Jb,No,3); % reshape into the 3 components

for ii = 1:Ne % loop on the elements
    
    r1 = node(:,elem(1,ii)); % 3xNe with coordinates of the first node of all elements
    r2 = node(:,elem(2,ii)); % 3xNe with coordinates of the first node of all elements
    r3 = node(:,elem(3,ii)); % 3xNe with coordinates of the first node of all elements
    
    % get the sign and index for the contribution
    absnum = abs(etod(:,ii)); % internal index of the edge
    mult = etod(:,ii)./absnum; % +1 or -1

    idx1 = index(absnum(1));
    idx2 = index(absnum(2));
    idx3 = index(absnum(3));
    
    Eout1 = 0;
    Eout2 = 0;
    Eout3 = 0;

    parfor jj = 1:No % loop on the number of observation points
        
        % obtain the contribution for each element
        robs = Scoord(jj,:);
        [GJ1, GJ2, GJ3] = mexCoupling ( r1, r2, r3, robs, ko, Np_2D, Z1, Z2, Z3,  wp );

        %         GJ1 = scale*GJ1;
        %         GJ2 = scale*GJ2;
        %         GJ3 = scale*GJ3;
        %         GGJ1(ii,jj,:) = GJ1;
        %         GGJ2(ii,jj,:) = GJ2;
        %         GGJ3(ii,jj,:) = GJ3;

        
        % add each contribution multiplied by the sign given by etod
        
        Eout1 = Eout1 + mult(1)*(Jb(jj,:)*GJ1); % first edge
        Eout2 = Eout2 + mult(2)*(Jb(jj,:)*GJ2); % second edge
        Eout3 = Eout3 + mult(3)*(Jb(jj,:)*GJ3); % third edge
        
    end
    
    if (idx1)
        Eout(idx1,1) = Eout(idx1,1) + Eout1; % first edge
    end
    
    if (idx2)
        Eout(idx2,1) = Eout(idx2,1) + Eout2; % second edge
    end
    
    if (idx3)
        Eout(idx3,1) = Eout(idx3,1) + Eout3; % third edge
    end

end


% Zbc = Zbc/(1j*ko);
% scalefactor = (1j*eta)/(4*pi*ko);
% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  -e^scattered
% -------------------------------------------------------------------------
ce = 1i*omega*eo;
scalefactor = - 1 / ce / (4*pi);

Eout = scalefactor*Eout;
%  Zbc = (4*pi)/(1i*eta) * Zc2h;
 
    


