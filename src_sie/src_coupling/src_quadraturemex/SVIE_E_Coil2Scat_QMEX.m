function [Eout] = SVIE_E_Coil2Scat_QMEX(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,Jc) 
%%    Coil to Scatterer Quadrature E field generation for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature-based E field
%   due to SIE Jc current coefficients
%   Applies mexed function
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
omega = 2 * pi * freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
% eta = omega*mu/ko; %3.767303134617706e+002; 

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

NO = size(Scoord,1); % number of observation points
NE = size(elem,2); % number of elements
NC = max(index); % number of dofs

RO = Scoord.';

R1 = node(:,elem(1,:)); % 3xNe with coordinates of the first node of all elements
R2 = node(:,elem(2,:)); % 3xNe with coordinates of the first node of all elements
R3 = node(:,elem(3,:)); % 3xNe with coordinates of the first node of all elements
    
% get the sign and index for the contribution
ABSNUM = abs(etod(:,:)); % internal index of the edge
MULT = etod(:,:)./ABSNUM; % +1 or -1
IDX = index(ABSNUM);


% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------

% LEVEL_DVrule = 5; % 5 is more than enough for most cases
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule ( LEVEL_DVrule );

% [Np_2D,wp,Z1,Z2,Z3] = Gauss_2Dt(1);
% -------------------------------------------------------------------------
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

% call the mex fun
[Eout] = mexQuadCoil2Scat(R1(:),R2(:),R3(:),NE,RO(:),NO,IDX(:),MULT(:),NC,ko,Np_2D,Z1,Z2,Z3,wp,Jc(:));

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

Eout = scalefactor * Eout;
%  Zbc = (4*pi)/(1i*eta) * Zbc;   

