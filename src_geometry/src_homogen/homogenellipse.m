function [idx,epsilon_r,sigma_e,rho,mu_r,sigma_m] = homogenellipse(r,Cnt,Rad,Erat,e_r,s_e,dens,m_r,s_m)
%%    Defines a Ellipse, and assigns homogeneous properties
% _________________________________________________________________________
%
%       Defines a Ellipse in a given domain.
%       Optionally assigns material properties
%
% _________________________________________________________________________
%
%% INPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%   Cnt         Cartesian coordinates of the center of the ellipse 
%   Rad         Radius of the ellipse
%   Erat        Elliptic ratio (scale of x, y and z axis)
%
%
%% OPTIONAL INPUT
%   e_r         value of relative epsilon
%   s_e         value of electric conductivity
%   dens        value of the density
%   m_r         value of relative mu
%   s_m         value of magnetic conductivity
%
%
%% OUTPUT
%   idx         indexes of the positions
%   epsilon_r   3D (LxMxN) array with relative epsilon
%   sigma_e     3D (LxMxN) array with electric conductivity
%   rho         3D (LxMxN) array with the density
%   mu_r        3D (LxMxN) array with relative mu
%   sigma_m     3D (LxMxN) array with magnetic conductivity
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
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 5 || isempty(e_r))
   e_r = 1;
end
if(nargin < 6 || isempty(s_e))
   s_e = 0;
end
if(nargin < 7 || isempty(dens))
   dens = 0;
end
if(nargin < 8 || isempty(m_r))
   m_r = 1;
end
if(nargin < 9 || isempty(s_m))
   s_m = 0;
end



% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);
Erat = squeeze(Erat);
Cnt = squeeze(Cnt);

% define sphere
sphere = @(r)( ( (r(:,:,:,1) - Cnt(1))/Erat(1) ).^2 + ( (r(:,:,:,2) - Cnt(2))/Erat(2) ).^2 + ( (r(:,:,:,3) - Cnt(3))/Erat(3) ).^2 < Rad^2) ;
pointsphere = sphere(r); % ones for domain elements in the shell
idx = find(pointsphere(:)); % get indexes of elements

% -------------------------------------------------------------------------
% assign output
% -------------------------------------------------------------------------

epsilon_r = ones(L,M,N);
epsilon_r(idx) = e_r; 

mu_r = ones(L,M,N);
mu_r(idx) = m_r;

sigma_e = zeros(L,M,N);
sigma_e(idx) = s_e;

sigma_m = zeros(L,M,N);
sigma_m(idx) = s_m;

rho = zeros(L,M,N);
rho(idx) = dens;

