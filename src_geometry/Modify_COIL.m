function [COIL] = Modify_COIL(COIL,angle,trans) 
%%    Modify a COIL model
% _________________________________________________________________________
%
%   Gets the original COIL, modifies the position
%
% _________________________________________________________________________
%
%
%% INPUT
%       COIL - initial COIL
%       angle - rotation in degrees
%       trans = translation in x, y and z
%
%
%% OUTPUT
%       COIL -  COIL Model
%           name - name of the coil
%           type - 'S' if surface, 'W' if wire model
%           Rhocoil - resistivity of material
%           (for surface coil)
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition (either S or W, same field)
%           Ct - coordinates of the center of the triangle
%           Ln - values of the length of each side of the triangle
%           Pn - 3x3 matrix with coordinates of the rho vectors (Pn(:,1) == rho_1)
%           (for wire coil)
%           Pcoil - positive end of segment 
%           Ncoil - negative end of segment
%           Dwire - diameter of wire
%           port - port definition (either S or W, same field)
%
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2014
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________


if(nargin < 2)
   angle = [];
end
if(nargin < 3)
   trans = [];
end


% -------------------------------------------------------------------------
% check coil type
% -------------------------------------------------------------------------

switch COIL.type
    case 'S'
        
        % node has the coordinates of the triangle points: 3xNtriang
        
        % -----------------------------------------------------------------
        % Rotate the coordinates
        % -----------------------------------------------------------------
        
        if ~isempty(angle)
            
            phi = angle*pi/180; % change to radians
             
            node1 = COIL.node(1,:)*cos(phi) - COIL.node(2,:)*sin(phi);
            node2 = COIL.node(1,:)*sin(phi) + COIL.node(2,:)*cos(phi);
            COIL.node(1,:) = node1;
            COIL.node(2,:) = node2;
            
        end
                        
        % -----------------------------------------------------------------
        % Translate the coordinates
        % -----------------------------------------------------------------
        
        if ~isempty(trans)
            
            COIL.node(1,:) = COIL.node(1,:) + trans(1);
            COIL.node(2,:) = COIL.node(2,:) + trans(2);
            COIL.node(3,:) = COIL.node(3,:) + trans(3);
            
        end
        
        % -----------------------------------------------------------------
        % recompute the data of the triangles
        % -----------------------------------------------------------------
        [Ct,Ln,Pn] = Mesh_CLP(COIL.node,COIL.elem);
        COIL.Ct = Ct;
        COIL.Ln = Ln;
        COIL.Pn = Pn;
            
        
    case 'W'
                
        % COIL.Pcoil and COIL.Ncoil have the coordinates of segments
        
        % -----------------------------------------------------------------
        % Rotate the coordinates
        % -----------------------------------------------------------------
        
        if ~isempty(angle)
            
            phi = angle*pi/180; % change to radians
            
            pcoil1 = COIL.Pcoil(:,1)*cos(phi) - COIL.Pcoil(:,2)*sin(phi);
            pcoil2 = COIL.Pcoil(:,1)*sin(phi) + COIL.Pcoil(:,2)*cos(phi);
            COIL.Pcoil(:,1) = pcoil1;
            COIL.Pcoil(:,2) = pcoil2;
            
            ncoil1 = COIL.Ncoil(:,1)*cos(phi) - COIL.Ncoil(:,2)*sin(phi);
            ncoil2 = COIL.Ncoil(:,1)*sin(phi) + COIL.Ncoil(:,2)*cos(phi);
            COIL.Ncoil(:,1) = ncoil1;
            COIL.Ncoil(:,2) = ncoil2;
            
        end
                        
        % -----------------------------------------------------------------
        % Translate the coordinates
        % -----------------------------------------------------------------
        
        if ~isempty(trans)
            
            COIL.Pcoil(:,1) = COIL.Pcoil(:,1) + trans(1);
            COIL.Pcoil(:,2) = COIL.Pcoil(:,2) + trans(2);
            COIL.Pcoil(:,3) = COIL.Pcoil(:,3) + trans(3);
            
            COIL.Ncoil(:,1) = COIL.Ncoil(:,1) + trans(1);
            COIL.Ncoil(:,2) = COIL.Ncoil(:,2) + trans(2);
            COIL.Ncoil(:,3) = COIL.Ncoil(:,3) + trans(3);
            
        end
       
end

