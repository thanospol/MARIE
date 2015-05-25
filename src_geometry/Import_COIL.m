function [COIL] = Import_COIL(coilfile) 
%%    Loads a COIL structure from either .smm or .wmm file
% _________________________________________________________________________
%
%   Reads the surface marie model (.smm) or wire marie model (.wmm) file
%   and generates the data for the COIL structure
%
% _________________________________________________________________________
%
%
%% INPUT
%       coilfile - name of the file, with path
%
%% OUTPUT
%       COIL structure with
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



% -------------------------------------------------------------------------
% Open file and read subfiles
% -------------------------------------------------------------------------

last = length(coilfile);
extension = coilfile(last-3:last);
validext = 0;

% -------------------------------------------------------------------------
% Surface coil model
% -------------------------------------------------------------------------

if (strcmp(extension,'.smm'))
    validext = 1;
    
    fid = fopen(coilfile, 'r');
    if (fid < 0)
        fprintf(1, ' \n WARNING!: Unable to open file: %s\n', coilfile);
        COIL.name = 'No Model Selected';
        COIL.type = [];
        return;
    end
    
    
    modelname = fgetl(fid); % read first line with name
    
    line = fgetl(fid); % read second line with resistivity and thickness
    vec = sscanf(line,'%f'); % read node number and coordinates
    Rho = vec(1);
    Thickness = vec(2);
    
    line = fgetl(fid); % read line with file to parse
    nchar = length(line);
    extension = line(nchar-3:nchar);
    
    if strcmp(extension, '.msh') % mesh file
        
        [SCOIL] = Import_SCOIL(line);
        
    else
        
        if strcmp(extension, '.mat') % mat saved file
            
            load(line); % load the data
            % it should have a valid SCOIL
            
        else
            
            fprintf(1, ' \n WARNING!: Invalid file format\n');
            COIL.name = 'No Model Selected';
            COIL.type = [];
            return;
            
        end
        
    end
    
    % store the DATA
    COIL = struct('name', modelname,...
        'type', 'S',...
        'Rhocoil', Rho,...
        'Thickness', Thickness,...
        'index',  SCOIL.index,...
        'etod', SCOIL.etod, ...
        'node', SCOIL.node, ...
        'edge', SCOIL.edge, ...
        'elem', SCOIL.elem, ...
        'index_elem', SCOIL.index_elem, ...
        'Ct', SCOIL.Ct, ...
        'Ln', SCOIL.Ln, ...
        'Pn', SCOIL.Pn, ...
        'port', SCOIL.port,...
        'Pcoil', [],...
        'Ncoil', [],...
        'Dwire', []);
    
end






% -------------------------------------------------------------------------
% Wire coil model
% -------------------------------------------------------------------------

if (strcmp(extension,'.wmm'))
    validext = 1;
    
    fid = fopen(coilfile, 'r');
    if (fid < 0)
        fprintf(1, ' \n WARNING!: Unable to open file: %s\n', coilfile);
        COIL.name = 'No Model Selected';
        COIL.type = [];
        return;
    end
    
    
    modelname = fgetl(fid); % read first line with name
    
    line = fgetl(fid); % read second line with resistivity and thickness
    vec = sscanf(line,'%f'); % read node number and coordinates
    Rho = vec(1);
    Dwire = vec(2);
    
    line = fgetl(fid); % read line with coordinate file
    nchar = length(line);
    extension = line(nchar-3:nchar);
        
    if strcmp(extension, '.wsd') % wire segment discretization file
        
        [Pcoil,Ncoil,port] = Import_WCOIL(line);
        
    else
        
        if strcmp(extension, '.mat') % mat saved file
            
            load(line); % load the data
            % it should have a valid WCOIL
            
        else
            
            fprintf(1, ' \n WARNING!: Invalid file format\n');
            COIL.name = 'No Model Selected';
            COIL.type = [];
            return;
            
        end
        
    end
    
    % store the DATA
    COIL = struct('name', modelname,...
        'type', 'W',...
        'Rhocoil', Rho,...
        'Thickness', [],...
        'index',  [],...
        'etod', [], ...
        'node', [], ...
        'edge', [], ...
        'elem', [], ...
        'index_elem', [], ...
        'Ct', [], ...
        'Ln', [], ...
        'Pn', [], ...
        'Pcoil', Pcoil,...
        'Ncoil', Ncoil,...
        'Dwire', Dwire,...
        'port', port);
    
end


% -------------------------------------------------------------------------
% Invalid file
% -------------------------------------------------------------------------


if (validext == 0)
    fprintf(1, ' \n WARNING!: Wrong Format\n');
    COIL.name = 'No Model Selected';
    COIL.type = [];
    return;
end


           
           
