function [U V SF RF] = Assembly_SCOUP_Q_aca(Scoord,index,etod,node,elem,freq,LEVEL_DVrule,tol,order) 
%%    Quadrature coupling for the SIE+VIE solver
% _________________________________________________________________________
%
%   Fucntion to generate the Quadrature Coupling SIE to VIE
%   Goes triangle by triangle, and calls the mexCoupling function for each
%   PARFOR version
%
% _________________________________________________________________________
%
%% Input
%       Scoord - coordinates of the observation points (No x 3)
%       node - coordinates of the nodes 
%       elem - 3 indexes of the nodes defining an element
%       etod - etod
%       index - mapping of the internal edge number to dof number
%       freq - frequency
%       LEVEL_DVrule - level for the DUNAVANT rule
%
%
%% Output
%       Zbc - Tensor (No x 3 x Nd) with the contribution of each edge of the element
%              Zbc(1000,3,200) is z contribution of 200-th edge to 1000-th element
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
%            Define EM constants
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2 * pi * freq;
lambda  = co/freq;
ko = 2*pi/lambda;

% Free-space impedance
% eta = omega*mu/ko; %3.767303134617706e+002; 

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
% loop on the elements and fill the matrix
% -------------------------------------------------------------------------

M  = 3*No; % total number of rows
Mc = No;   % rows per component
N  = Nd;   % total number of columns

D = order;

I = 1;
J = 0;

U = zeros(M, 0);
V = zeros(N, 0);

SF2  = zeros(1, 1);
RF2  = zeros(1, 1);
uvF2 = zeros(D, 1);

Ur = zeros(M, D);
Vr = zeros(N, D);

Vk = zeros(N, 3);
vk = zeros(N, 1);
uk = zeros(M, 1);

% dtoe
contributors = cell([Nd 3]);
for n = 1:Ne
	for k = 1:3
		idx = index(abs(etod(k,n)));
		if idx
			contributors{idx, k} = [contributors{idx, k} sign(etod(k,n))*n];
		end
	end
end

for k = 1:min(M,N)
    % --------------------------------------
	%         calculate v_k
	% --------------------------------------
	
	% calculate rows of A corresponding to I
	Ibase = mod(I-1, Mc)+1;
    parfor n = 1:Ne
		Vk_local = zeros(N,3);
		GJ = zeros(3, 3);
		r = node(:,elem(1:3,n));
        [GJ(:,1),GJ(:,2),GJ(:,3)] = mexCoupling(r(:,1), r(:,2), r(:,3), Scoord(Ibase,:), ko, Np_2D, Z1, Z2, Z3, wp);
		for d=1:3
			idx = index(abs(etod(d,n)));
			if idx
				for r=1:3
					Vk_local(idx, r) = Vk_local(idx, r) + (sign(etod(d,n))*GJ(r,d))';
				end
			end
		end
		Vk = Vk + Vk_local;
    end
	% calculate rows of R from rows of A
	for r=1:3
		Ir = Ibase+(r-1)*Mc;
		Vk(:,r) = Vk(:,r) - sum(V*U(Ir,:)', 2);
	end
	% find row with biggest element
	[Q J] = max(abs(Vk));
	[Q r] = max(abs([Vk(J(1),1) Vk(J(2),2), Vk(J(3), 3)]));
	% set new vk to row with largest element
	% and normalize the vector
	vk = Vk(:,r) / Vk(J(r), r);
	J = J(r);
	
	% --------------------------------------
	%           calculate u_k
	% --------------------------------------
	
	% calculate column of A
	for d=1:3
		for n=contributors{J,d}
			r = node(:, elem(1:3,abs(n)));
			parfor m=1:No
				GJ = zeros(3, 3);
				[GJ(:,1), GJ(:,2), GJ(:,3)] = mexCoupling(r(:,1), r(:,2), r(:,3), Scoord(m,:), ko, Np_2D, Z1, Z2, Z3, wp);
				uk_local = zeros(M,1);
				uk_local(m:Mc:m+2*Mc,1) = sign(n)*GJ(:,d);
				uk = uk + uk_local;
			end
		end
	end
	% calculate column of R from column of A
	uk = uk - sum(U*V(J,:)', 2);
	
	% --------------------------------------
	%      update norms & update basis
	% --------------------------------------
	
	% approximation norm
	uvF2(k+D) = uk'*uk .* vk'*vk;
	SF2(k+1) = SF2(k) + uvF2(k+D) + 2*sum(real((U'*uk).*(V'*vk)));
	% approximation norm
	Us = Ur(:, 2:end);
	Vs = Vr(:, 2:end);
	ukD = Ur(:,1);
	vkD = Vr(:,1);
	RF2(k+1) = RF2(k) + uvF2(k+D) - uvF2(k) + 2*sum(real((Us'*uk).*(Vs'*vk))) - 2*sum(real((Us'*ukD).*(Vs'*vkD)));
	% enlarge basis
	U = [U uk];
	V = [V vk];
	Vr = [Vs vk];
	Ur = [Us uk];
	
	% ---------------------------------------------
	%  terminate or get index of new row to expand
	% ---------------------------------------------
	
	if sqrt(RF2(k+1)) <= tol*sqrt(SF2(k+1))
		break;
	end
	
	uk = abs(uk);
	uk(I) = 0;
	[Q I] = max(uk);
	
	Vk(:) = 0;
	vk(:) = 0;
	uk(:) = 0;
end

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

SF = sqrt(SF2);
RF = sqrt(RF2);
end