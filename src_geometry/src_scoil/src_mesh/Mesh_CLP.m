function [ct,ln,pn] = Mesh_CLP(node,elem)
%%  Pre-processing triangles to get info gor coupling
% _________________________________________________________________________
%
%   Pre-processing on the triangle elements to generate the information
%   for obtaining the coupling
%
% _________________________________________________________________________
%
%% Input
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%
%
%% Output
%       ct - coordinates of the baricenter of the triangle
%       ln - value of the length of each side of the triangle (3 values)
%       pn - 3 vectors with the components of each rho for the vertex
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
% Define variables and allocate space
% -------------------------------------------------------------------------

Ne = size(elem,2); % number of elements

ln = zeros(3,Ne); % lengths for n=1,2,3
pn = zeros(3,3,Ne); % the rows are 3 components for n=1,2,3 for each elem


% -------------------------------------------------------------------------
% (Avoid) loop on the elements and fill the information
% -------------------------------------------------------------------------

% get the coordinates of the nodes

r_1 = node(:,elem(1,:)); % 3xNe with coordinates of the first node of all elements
r_2 = node(:,elem(2,:)); % 3xNe with coordinates of the first node of all elements
r_3 = node(:,elem(3,:)); % 3xNe with coordinates of the first node of all elements


% get the values of the length of the edges

ledge = (r_2 - r_3); % vector with components of l_1
ledge = ledge.*ledge; % square the components
ln(1,:) = sum(ledge,1); % sum accross the rows to get the square of the length

ledge = (r_3 - r_1); % vector with components of l_2
ledge = ledge.*ledge; % square the components
ln(2,:) = sum(ledge,1); % sum accross the rows to get the square of the length

ledge = (r_1 - r_2); % vector with components of l_3
ledge = ledge.*ledge; % square the components
ln(3,:) = sum(ledge,1); % sum accross the rows to get the square of the length

ln = sqrt(ln); % get the element wise sq root with distances


% get the coordinates of the center of the elements

ct = (r_1 + r_2 + r_3)./3;


% get the components of the rho vectors

pn(:,1,:) = ct - r_1; % rho_1
pn(:,2,:) = ct - r_2; % rho_2
pn(:,3,:) = ct - r_3; % rho_3


    

    