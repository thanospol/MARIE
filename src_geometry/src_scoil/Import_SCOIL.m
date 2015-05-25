function [SCOIL] = Import_SCOIL(filename) 
%%    Loads a SCOIL structure from a mesh file
% _________________________________________________________________________
%
%   Reads the mesh file, and generates the data for the SCOIL
%
% _________________________________________________________________________
%
%
%% INPUT
%       filename - name of the file, with path, with .msh
%
%
%% OUTPUT
%       SCOIL structure with
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           Ct - coordinates of the center of the triangle
%           Ln - values of the length of each side of the triangle
%           Pn - 3x3 matrix with coordinates of the rho vectors (Pn(:,1) == rho_1)
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
%   Parse
% -------------------------------------------------------------------------

% call the mesh parser for surfaces
[node,~,e,elem,~] = Mesh_Parse(filename);

% -------------------------------------------------------------------------
% Pre-Processing
% -------------------------------------------------------------------------

[edge,etod,index,port,index_elem] = Mesh_PreProc(e,elem);

[Ct,Ln,Pn] = Mesh_CLP(node,elem);

% -------------------------------------------------------------------------
% Assign data to structure
% -------------------------------------------------------------------------

SCOIL = struct('index',  index,...
               'etod', etod, ...
               'node', node, ...
               'edge', edge, ...
               'elem', elem, ...
               'index_elem', index_elem, ...
               'Ct', Ct, ...
               'Ln', Ln, ...
               'Pn', Pn, ...
               'port', port);
           
           
