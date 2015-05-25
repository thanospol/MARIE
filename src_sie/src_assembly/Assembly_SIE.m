function [Z,F] = Assembly_SIE(SCOIL,freq) 
%%    Assembly the SIE system
% _________________________________________________________________________
%
%   Fucntion to generate the free-space Z matrix for the SIE
%   RWG functions are applied
%
% _________________________________________________________________________
%
%
%% INPUT
%       SCOIL structure with
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%       freq - frequency
%
%
%% OUTPUT
%       Z is the impedance matrix
%       F is the rhs for the port definition
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
%   

tini = clock;

% -------------------------------------------------------------------------
%            Define EM constants
% -------------------------------------------------------------------------

co = 299792458;
%
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
eta =  3.767303134617706e+002; 

% extract struct data
index = SCOIL.index;
etod = SCOIL.etod;
node = SCOIL.node;
edge = SCOIL.edge;
elem = SCOIL.elem;
Index_elem = SCOIL.index_elem;
port = SCOIL.port;

% -------------------------------------------------------------------------
%            Info 
% -------------------------------------------------------------------------

Nvars = max(index);
Nports = length(port);
ie_NS = Index_elem.NS(:,1); 
ie_ST = Index_elem.ST(:,1); 
ie_EA = Index_elem.EA(:,1); 
ie_VA = Index_elem.VA(:,1); 

fid = 1;
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n Assembling SIE System');
fprintf(fid, '\n');
fprintf(fid, '\n # DOFS:                %d', Nvars);
fprintf(fid, '\n # PORTS:               %d', Nports);
fprintf(fid, '\n # NS interactions:     %d',length(ie_NS));
fprintf(fid, '\n # EA interactions:     %d',length(ie_EA));
fprintf(fid, '\n # VA interactions:     %d',length(ie_VA));
fprintf(fid, '\n # ST interactions:     %d',length(ie_ST));
fprintf(fid, '\n Operating Frequency:   %.3f MHz',freq/1e6);
fprintf(fid, '\n');



% -------------------------------------------------------------------------
%            Order of Gauss Quadrature Integration 
% -------------------------------------------------------------------------

% Data
Eo     = 1;
B_o    = - Eo ;   % - sign because e^sca = -e^inc

% Quadrature
N_NS         = 4;  % NS
N_ST_psi     = 10; % ST
N_EA_theta   = 6;  % EA
N_EA_psi     = 6;
N_VA_theta_p = 6;  % VA
N_VA_theta_q = 6;
N_VA_psi     = 6;

% gather in a structure GL_order
GL_order.ST = N_ST_psi;
GL_order.EA = [N_EA_theta N_EA_psi];
GL_order.VA = [N_VA_theta_p N_VA_theta_q N_VA_psi];
GL_order.NS  = N_NS;

% -------------------------------------------------------------------------
%            Start the assembly procedure 
% -------------------------------------------------------------------------

%
Z = zeros(Nvars,Nvars);
F = zeros(Nvars,Nports);

% % isOpen = (matlabpool('size') > 0);
% % if (~isOpen)
% %     mycluster = parcluster;
% %     matlabpool('local',mycluster.NumWorkers);
% % end

% -------------------------------------------------------------------------
%             NS interactions
% -------------------------------------------------------------------------

tic 

[Z] = assembly_ns(index,etod,node,elem,Z,GL_order,Index_elem,ko);

TimeAssembly_NS = toc;

% -------------------------------------------------------------------------
%             EA interactions 
% -------------------------------------------------------------------------

tic 

[Z] = assembly_ea(index,etod,node,elem,Z,GL_order,Index_elem,ko);

TimeAssembly_EA = toc;

% -------------------------------------------------------------------------
%             VA interactions 
% -------------------------------------------------------------------------

tic 

[Z] = assembly_va(index,etod,node,elem,Z,GL_order,Index_elem,ko);

TimeAssembly_VA = toc;


% -------------------------------------------------------------------------
%             ST interactions 
% -------------------------------------------------------------------------

Z = Z + Z.'; 
%
tic 

[Z] = assembly_st(index,etod,node,elem,Z,GL_order,Index_elem,ko);

TimeAssembly_ST = toc;

% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  e^scattered
% -------------------------------------------------------------------------

Z = -(eta/(4*pi)) * Z;

% -------------------------------------------------------------------------
%             V assembly 
% -------------------------------------------------------------------------

tic

for pnum = 1:Nports
    [Vp] = excitation_coil(B_o,index,node,edge,port(pnum));
    F(:,pnum) = Vp;
end

TimeExcitation = toc;


% -------------------------------------------------------------------------
%             and it is done: report 
% -------------------------------------------------------------------------

tend = clock;

fprintf(fid, '\n Time in NS:            %.2f sec', TimeAssembly_NS);
fprintf(fid, '\n Time in EA:            %.2f sec', TimeAssembly_EA);
fprintf(fid, '\n Time in VA:            %.2f sec', TimeAssembly_VA);
fprintf(fid, '\n Time in ST:            %.2f sec', TimeAssembly_ST);
fprintf(fid, '\n Time in Port:          %.2f sec', TimeExcitation);
fprintf(fid, '\n -------------------------------------');
fprintf(fid, '\n Overall TIME:          %.3f sec', tend);
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n ');

