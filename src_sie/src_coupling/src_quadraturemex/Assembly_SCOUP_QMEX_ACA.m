function [U, V] = Assembly_SCOUP_QMEX_ACA(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,tol,order)
% -------------------------------------------------------------------------
%            Define EM constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;

omega = 2 * pi * freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% -------------------------------------------------------------------------
% Define variables and allocate space
% -------------------------------------------------------------------------

No = size(Scoord,1); % number of observation points
Ne = size(elem,2); % number of elements
Nd = max(index); % number of dofs

% -------------------------------------------------------------------------
% 1D cubature's number of points
% -------------------------------------------------------------------------
[ Np_2D, Z1, Z2, Z3, wp ] = dunavant_rule ( LEVEL_DVrule );

% -------------------------------------------------------------------------
% Instantiate variables for mex code
% -------------------------------------------------------------------------

RO = Scoord.';
R1 = node(:,elem(1,:)); % 3xNe with coordinates of the first node of all elements
R2 = node(:,elem(2,:)); % 3xNe with coordinates of the first node of all elements
R3 = node(:,elem(3,:)); % 3xNe with coordinates of the first node of all elements

idxNZ = zeros(Nd,3);
dtoe = zeros(3,3,Nd);
for n=1:Ne
	for k=1:3
		idx = index(abs(etod(k,n)));
		if idx
			idxNZ(idx,k) = idxNZ(idx,k)+1;
			dtoe(idxNZ(idx,k), k, idx) = sign(etod(k,n))*n;
		end
	end
end
dtoe = int64(dtoe(:));
etod = index(abs(etod(:,:))).*sign(etod(:,:));
etod = int64(etod(:));

M  = 3*No; % total number of rows
Mc = No;   % rows per component
N  = Nd;   % total number of columns
D = order;

% -------------------------------------------------------------------------
% Call mex function
% -------------------------------------------------------------------------
[U, V] = ompQuadCoil2ScatACA(R1(:),R2(:),R3(:),Ne,RO(:),No,etod,Nd,ko,Np_2D,Z1,Z2,Z3,wp,dtoe,tol,order);
U = cell2mat(U);
V = cell2mat(V);

% -------------------------------------------------------------------------
%             Final Z (with mult. constant) 
%
%      4*pi comes from omitting it in Green function!!
%       Z = discretization of  -e^scattered
% -------------------------------------------------------------------------
ce = 1i*omega*eo;
scalefactor = - 1 / ce / (4*pi);

if numel(U) < numel(V)
	U = scalefactor * U;
else
	V = conj(scalefactor) * V;
end
