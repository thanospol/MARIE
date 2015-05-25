function  [Eout] = SVIE_E_Coil2Scat_Q(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,Jc)
%%    Coil to Scatterer Quadrature E field generation for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature-based E field
%   due to SIE Jc current coefficients
%   Goes triangle by triangle, and calls the mexCoupling function for each
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       freq - frequency
%       LEVEL_DVrule - level for the DUNAVANT rule
%       Jc - current coefficients for the coil
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

Eout = zeros(No,3);


% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % enough
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule ( LEVEL_DVrule );

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);

% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------


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

    for jj = 1:No % loop on the number of observation points
        
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
        
        if (idx1)
            Eout(jj,:) = Eout(jj,:) + Jc(idx1)*mult(1)*GJ1.'; % first edge
        end
        
        if (idx2)
            Eout(jj,:) = Eout(jj,:) + Jc(idx2)*mult(2)*GJ2.'; % second edge
        end
        
        if (idx3)
            Eout(jj,:) = Eout(jj,:) + Jc(idx3)*mult(3)*GJ3.'; % third edge
        end
        
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

Eout = reshape(Eout,3*No,1);
    