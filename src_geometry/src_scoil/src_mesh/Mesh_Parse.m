function [nodes,singlep,lines,triang,quad] = Mesh_Parse(meshfile)
%%    Parser for GMESH kind of discretization
% _________________________________________________________________________
%
%
%   Read basic surface mesh file from gmsh
%
% _________________________________________________________________________
%
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% voxelized
fid= fopen(meshfile, 'r');

while ~(feof(fid)) % it is not the end of the file, keep reading

       
        line = fgetl(fid); % read line
        
        % So far I do not give a S%*& about anything that is not a node or element 
        
        % NODDES ----------------------------------------------------------        
        if strcmp(line, '$Nodes') % start the nodes
            
            line = fgetl(fid); % read line
            nnodes = sscanf(line,'%d'); % read number of nodes
            
            % allocate space
            nodes = zeros(3,nnodes); % x, y and z coordinate of each node
            
            %IMPORTANT: in gmesh the node-number (first integer) is the one
            % given in the geo file and by gmesh, and it is  a postive (non-zero) integer,
            % but note that the node-numbers do not necessarily have to form a dense nor an ordered sequence. 
            nodeidx = zeros(nnodes,1); % index that translate node-number to the new number (position) in the 'nodes' matrix      
            
            for ii = 1:nnodes % read all nodes
                
                line = fgetl(fid); % read line
                vec = sscanf(line,'%f'); % read node number and coordinates
                nodes(1,ii) = vec(2);
                nodes(2,ii) = vec(3);
                nodes(3,ii) = vec(4);
                nodeidx(vec(1)) = ii; % in position node-number we store the new node number, which are now a sorted and compact sequence
                
            end
            
            line = fgetl(fid); % read line
            if ~(strcmp(line, '$EndNodes')) % read until the end of nodes
               fprintf(1, '\nWarning: number of nodes read do not match number of nodes given in mesh file! \n');
            end
            
        end
        
        % ELEMENTS --------------------------------------------------------
        if strcmp(line, '$Elements') % start the nodes
            
            line = fgetl(fid); % read line
            nelem = sscanf(line,'%d'); % read number of nodes
            
            % allocate space for nodes, edges and elements
            lines = zeros(3,nelem); % point 1, point 2, and physical flag of edge (line to which belongs)
            triang = zeros(4,nelem); % point 1, point 2, point 3 and physical flag of element (Surface to which belongs)
            quad = zeros(5,nelem); % point 1, point 2, point 3, point 4 and physical flag of element (Surface to which belongs)
            singlep = zeros(2,nelem);
            
            ecount = 0;
            tcount = 0;
            qcount = 0;
            pcount = 0;
            
            for ii = 1:nelem % read all nodes
                
                line = fgetl(fid); % read line
                vec = sscanf(line,'%d'); % read node number and values coordinates
                
                eltype = vec(2);
                
                switch eltype
                    
                    case 1 % edge
                 
                        ecount = ecount + 1;
                        % vec(3) give us the number of tags
                        ntags = vec(3);
                        lines(1,ecount) = nodeidx(vec(4+ntags)); % first point
                        lines(2,ecount) = nodeidx(vec(5+ntags)); % second point
                        lines(3,ecount) = vec(4); % phisical line
                    
                    case 2 % triangle
                        
                        tcount = tcount + 1;
                        % vec(3) give us the number of tags
                        ntags = vec(3);
                        triang(1,tcount) = nodeidx(vec(4+ntags)); % first point
                        triang(2,tcount) = nodeidx(vec(5+ntags)); % second point
                        triang(3,tcount) = nodeidx(vec(6+ntags)); % third point
                        triang(4,tcount) = vec(4); % phisical surface
                        
                    case 3 % quadrangle
                        
                        qcount = qcount + 1;
                        % vec(3) give us the number of tags
                        ntags = vec(3);
                        quad(1,qcount) = nodeidx(vec(4+ntags)); % first point
                        quad(2,qcount) = nodeidx(vec(5+ntags)); % second point
                        quad(3,qcount) = nodeidx(vec(6+ntags)); % third point
                        quad(4,qcount) = nodeidx(vec(7+ntags)); % fourth point
                        quad(5,qcount) = vec(4); % phisical surface
                        
                    case 15 % point
                        
                        pcount = pcount + 1;
                        % vec(3) give us the number of tags
                        ntags = vec(3);
                        singlep(1,pcount) = nodeidx(vec(4+ntags)); % first point
                        singlep(5,pcount) = vec(4); % phisical
                        
                    otherwise
                        
                        fprintf(1, '\nWarning: A non surface element was found!\n');
                        
                end % end switch
                        
                      
            end % end for reading elements
            
            % truncate the size of the matrices to the final number
            lines = lines(:,1:ecount); 
            triang = triang(:,1:tcount); 
            quad = quad(:,1:qcount); 
            singlep = singlep(:,1:pcount);
            
            line = fgetl(fid); % read line
            if ~(strcmp(line, '$EndElements')) % read until the end of elements
            
                fprintf(1, '\nWarning: number of elements read do not match number of elements given in mesh file! \n');
            
            else % we are done reading
                
                break; % get out of the loop                
            
            end
                       
        end
                
end

fclose(fid);

% report

fprintf(1, '\n\n -------------------------------------------------------');
fprintf(1, '\n  File:');
fprintf(1, '\n        %s', meshfile);
fprintf(1, '\n  Parse complete');
fprintf(1, '\n        %d nodes', nnodes);
fprintf(1, '\n        %d element points', pcount);
fprintf(1, '\n        %d lines', ecount);
fprintf(1, '\n        %d triangles', tcount);
fprintf(1, '\n        %d quadrangles', qcount);
fprintf(1, '\n -------------------------------------------------------\n\n');



            


