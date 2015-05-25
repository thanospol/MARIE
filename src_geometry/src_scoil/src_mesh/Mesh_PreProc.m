function [edge,etod,index,port,index_elem] = Mesh_PreProc(e,elem)
%%    Pre Processing for GMESH kind of discretization
% _________________________________________________________________________
%
%
%   Pre-process a mesh file from gmsh into SIE-friendly data
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
%

% Nn  = size(node,2); % Number of nodes
Nbd = size(e,2);    % Number of exterior (boundary) edges
Ne  = size(elem,2); % Number of elements


% -------------------------------------------------------------------------
%   Pre-processing: generate edge and etod: edge[2xNd], etod[3xNe]
%   Finds adjancencies
%   Assigns boundaries
% -------------------------------------------------------------------------

% t1 = clock;


first_node  = [3 1 2]; % define edge vectors clockwise
second_node = [2 3 1];

% for edge and etod generation
edge = zeros(2,3*Ne); % allocate maximun possible size
etod = zeros(3,Ne); % allocate size

% define array for boundary
% % kn=-1*ones(3*Ne,1); % allocate maximun possible size
kn=zeros(3*Ne,1); % allocate maximun possible size

% define arrays for adjacency
eparent = zeros(3*Ne,1); % allocate maximum possible size for adjacency check

NS = zeros(floor(Ne*Ne/2),2); nscount = 0;
VA = zeros(100*Ne,2); vacount = 0;
EA = zeros(3*Ne,2); eacount = 0;


% loop on elements
Nd = 0;
for ii=1:Ne
    
    count = Nd;
    flagpad = zeros(ii-1,1);
    flagead = zeros(ii-1,1);
    
 
    for in=1:3
        
        node1 = elem(first_node(in),ii);        
        node2 = elem(second_node(in),ii);
        flage  = 0;
        
        % check vertex adjacency
        R = elem(1,1:ii-1) - node1; flagpad(R == 0) = 1;
        R = elem(2,1:ii-1) - node1; flagpad(R == 0) = 1;
        R = elem(3,1:ii-1) - node1; flagpad(R == 0) = 1;
       
        for jj=1:Nd
            
            if (node1 == edge(1,jj)) && (node2 == edge(2,jj))
                etod(in,ii) = jj; % positive match found
                flage = 1; % edge adjacency
                flagead(eparent(jj)) = 2;
                kn(jj) = -1; % set boundary flag as internal node
            end
            
            if (node2 == edge(1,jj)) &&  (node1 == edge(2,jj))
                etod(in,ii) = -jj; % negative match found
                flage = 1; % edge adjacency
                flagead(eparent(jj)) = 2;
                kn(jj) = -1; % set boundary flag as internal node
            end
          
        end
        
        
        if (flage == 0) % new edge
            
            count = count+1;
            edge(1,count) = node1;
            edge(2,count) = node2;
            etod(in,ii) = count;
            
            eparent(count) = ii; % store the element to which the edge belongs
            
            
            for kk=1:Nbd % assign boundary
                bn1=e(1,kk);
                bn2=e(2,kk);
                if( (bn1 == node1) && (bn2 == node2) )
                    kn(count)=e(3,kk); % zero if external, positive integer if port
                end
                if( (bn2 == node1) && (bn1 == node2))
                    kn(count)=e(3,kk); % zero if external, positive integer if port
                end
            end
            
           
        end
        
    end
    
    
    % %     % check adjacency
       
    for kk = 1:length(flagead)
        
        adval = max(flagead(kk),flagpad(kk));
        switch adval
            case 1
                vacount = vacount+1;
                VA(vacount,1) = ii;
                VA(vacount,2) = kk;
            case 2
                eacount = eacount+1;
                EA(eacount,1) = ii;
                EA(eacount,2) = kk;
            case 0
                nscount = nscount+1;
                NS(nscount,1) = ii;
                NS(nscount,2) = kk;
        end
        
    end
    
    Nd = count;
    
end

edge = edge(:,1:Nd);
kn = kn(1:Nd);

% remove replicas in adjacency
ST = [1:Ne; 1:Ne].';
EA = EA(1:eacount,:);
VA = VA(1:vacount,:);
NS = NS(1:nscount,:);

% 
% t2 = clock;
% fprintf(1,'\n First loop done. elapsed time %f\n', etime(t2,t1));

% -------------------------------------------------------------------------
%   find the elements that are ports
% -------------------------------------------------------------------------

% for jj = 1:size(e,2)
% 
%     substractn1 = edge(1,:) - e(1,jj); % direct edge
%     substractn2 = edge(2,:) - e(2,jj);
%     sumnodes = abs(substractn1) + abs(substractn2);
%     idxport = find(sumnodes == 0);
%     if (idxport)
% %         fprintf(1,'\n (direct) port found it %d: port %d in edge %d', jj, e(3,jj), idxport);
%         kn(idxport) = e(3,jj); % zero if external, positive integer if port
%     end
%     
% 
%     substractn1 = edge(1,:) - e(2,jj); % reverse
%     substractn2 = edge(2,:) - e(1,jj);
%     sumnodes = abs(substractn1) + abs(substractn2);
%     idxport = find(sumnodes == 0);
%     if (idxport)
% %         fprintf(1,'\n (reverse) port found it %d: port %d in edge %d', jj, e(3,jj), idxport);
%         kn(idxport) = e(3,jj); % zero if external, positive integer if port
%     end
%     
% end


% -------------------------------------------------------------------------
%   Indexing & Number of unknown dofs (Nff)
% -------------------------------------------------------------------------

% t1 = clock;

% allocate space to index array
index=zeros(Nd,1);

% different types of edges, and sort them in ascending order
etype = sort(unique(kn)); 
idx = find(etype > 0);
etype = etype(idx);

portstruct.p = [];
portstruct.n = [];


% allocate space for the ports
if (~isempty(etype))
    nports = (length(etype))/2; % ports are pairs of edges, labeled with integer > 0
    port(nports) = portstruct;
end

% initialize counter and port number
ctr3 = 0;
pnum = 0;
for ii = 1:2:length(etype) % ports first! (ports do come in pairs)
    
    pnum = pnum+1; % increase port number
    
    % positive edge of the port 
    
    idx = find(kn == etype(ii)); % find indexes of elements with etype flag

    dofnum = ctr3+1:ctr3+length(idx); % create the dof number for each element
    
    index(idx) = dofnum; % assign the dof number to the index position
    
    port(pnum).p = dofnum; % store the dof number for the positive edge of port 
    
    ctr3 = ctr3+length(idx); % increase counter
    
    
    % negative edge of the port
    
    idx = find(kn == etype(ii+1)); % find indexes of elements with etype flag
    dofnum = ctr3+1:ctr3+length(idx); % create the dof number for each element
    
    index(idx) = dofnum; % assign the dof number to the index position
    
    port(pnum).n = dofnum; % store the dof number for the positive edge of port
    
    ctr3 = ctr3+length(idx); % increase counter
    
end

% now for the internal edges!
idx = find(kn == -1);
dofnum = ctr3+1:ctr3+length(idx); % create the dof number for each element
index(idx) = dofnum;

% and it is done

% t2 = clock;
% fprintf(1,'\n Third loop done. elapsed time %f\n', etime(t2,t1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Find adjacency                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index_elem.ST = ST;
index_elem.VA = VA;
index_elem.EA = EA;
index_elem.NS = NS;



