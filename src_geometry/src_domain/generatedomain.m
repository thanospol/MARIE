function [r] = generatedomain(res,x,y,z)
%%    Generates a 3D grid from the limits of a given domain
% _________________________________________________________________________
%
%       Generates a 3D cartesian grid of coordinates with given resolution
%       The resolution is fixed, and the final domain is the minimum for
%       the given resolution that encloses the specified dimensions
%
% _________________________________________________________________________
%
%% INPUT
%   res         resolution
%   x           minimum and maximum values of x
%   y           minimum and maximum values of y
%   z           minimum and maximum values of z
%
%
%% OUTPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
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

if(nargin < 2 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 3 || isempty(y))
   y = x;
end
if(nargin < 4 || isempty(z))
   z = x;
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

% just in case
x = squeeze(x);
y = squeeze(y);
z = squeeze(z);

% if the input is a single element, takes negative and positive
if length(x) == 1
    x = [-x x];
end
if length(y) == 1
    y = [-y y];
end
if length(z) == 1
    z = [-z z];
end

% -------------------------------------------------------------------------
% generate coordinate arrays
% -------------------------------------------------------------------------

% obtain length of each size
Dx = max(x) - min(x);
Dy = max(y) - min(y);
Dz = max(z) - min(z);

% obtain minimum number of cells in each direction
nx = ceil(Dx/2/res);
ny = ceil(Dy/2/res);
nz = ceil(Dz/2/res);

% compute centers of each array
Cx = (max(x) + min(x))/2;
Cy = (max(y) + min(y))/2;
Cz = (max(z) + min(z))/2;

% generate the arrays
x = Cx-nx*res:res:Cx+nx*res;
y = Cy-ny*res:res:Cy+ny*res;
z = Cz-nz*res:res:Cz+nz*res;

% -------------------------------------------------------------------------
% Generate grid
% -------------------------------------------------------------------------

r = grid3d(x,y,z);

