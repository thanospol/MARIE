function [Pcoil,Ncoil,port] = Import_WCOIL(filename) 
%%    Loads a WCOIL structure from a wire segment discretization (.wsd) file
% _________________________________________________________________________
%
%   Reads the file, and generates the data for the WCOIL
%
% _________________________________________________________________________
%
%
%% INPUT
%       filename - name of the file, with path, without .msh
%
%
%% OUTPUT
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           port - port definition
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
%   LOAD file as ASCII
% -------------------------------------------------------------------------

S = load(filename, '-ascii');

% first 3 values: x y and z coordinates of negative edge
Ncoil = S(:,1:3);

% second 3 values: x y and z coordinates of positive edge
Pcoil = S(:,4:6);

% last element, number of port (0 if no port)
[portindex] = find(squeeze(S(:,7)));
portvalues = squeeze(S(portindex,7));
idx = sort(portvalues, 'ascend');
port = portindex(idx);

