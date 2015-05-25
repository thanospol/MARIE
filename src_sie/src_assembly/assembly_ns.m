function [Z_NS] = assembly_ns(index,etod,node,elem,ZR,GL_order,Index_elem,ko)
%%
%
A_o    = 1i*ko;
F_o    = 1/(1i*ko);
%
first_node  = [3 1 2];
second_node = [2 3 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Order of Gauss Quadrature Integration                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np_2D  = GL_order.NS;   
%%%%%%%%%%%%%%% Gauss quadrature rule for non-singular triangles %%%%%%%%%%
[Np,wt,Z,z1,z2,z3] = gauss_2d_sie(Np_2D);
W_t = wt'*wt;

Zi_j = zeros(Np,Np,3,3);
%
for ii=1:3
    for jj=1:3
        Zi_j(:,:,ii,jj) = Z(:,ii)*Z(:,jj)';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ie_NS = Index_elem.NS(:,1); 
je_NS = Index_elem.NS(:,2); 
n_NS_elem = length(ie_NS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Assembly                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%%%%%%%%%%%%%%%%%%%%%%%%% Main body of assembly %%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rmn_x = zeros(Np,Np);
Rmn_y = zeros(Np,Np);
Rmn_z = zeros(Np,Np);
%
Z_NS = ZR;
%
% mycluster = parcluster;
% matlabpool('local',mycluster.NumWorkers);
kko = ko;
%%%%%%%%%%%%%%%%%%%%%%%%% Non-Singular Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for index_NS = 1 : n_NS_elem
    
    ie = ie_NS(index_NS);
    je = je_NS(index_NS);
    
    % for ie=1:Ne
    %
    ao = abs(etod(:,ie));
    % coordinates of nodes of the observation triangle
    node_test_1 = elem(1,ie);
    node_test_2 = elem(2,ie);
    node_test_3 = elem(3,ie);
    %
    ro_1 = node(:,node_test_1);
    ro_2 = node(:,node_test_2);
    ro_3 = node(:,node_test_3);
    % vectors lo1, lo2, lo3
    l_obs = [ro_2-ro_3,ro_3-ro_1,ro_1-ro_2];
    %
    %     for je=(ie+1):Ne
    as = abs(etod(:,je));
    % coordinates of nodes of the source triangle
    node_basis_1 = elem(1,je);
    node_basis_2 = elem(2,je);
    node_basis_3 = elem(3,je);
    %
    rs_1 = node(:,node_basis_1);
    rs_2 = node(:,node_basis_2);
    rs_3 = node(:,node_basis_3);
    % ls1, ls2, ls3
    l_src = [rs_2-rs_3,rs_3-rs_1,rs_1-rs_2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ro_x = z1*ro_1(1)+z2*ro_2(1)+z3*ro_3(1);
    ro_y = z1*ro_1(2)+z2*ro_2(2)+z3*ro_3(2);
    ro_z = z1*ro_1(3)+z2*ro_2(3)+z3*ro_3(3);
    %
    rs_x = z1*rs_1(1)+z2*rs_2(1)+z3*rs_3(1);
    rs_y = z1*rs_1(2)+z2*rs_2(2)+z3*rs_3(2);
    rs_z = z1*rs_1(3)+z2*rs_2(3)+z3*rs_3(3);
    %
    for ii=1:Np
        Rmn_x(ii,:) = ro_x(ii)-rs_x;
        Rmn_y(ii,:) = ro_y(ii)-rs_y;
        Rmn_z(ii,:) = ro_z(ii)-rs_z;
    end
    %
    Rmn = sqrt(Rmn_x.^2 + Rmn_y.^2 + Rmn_z.^2);
    %
    GRmn = exp(-1i*kko*Rmn)./Rmn;
    W_GR = W_t .* GRmn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:3
        for j=1:3
            
            index_ao = index(ao(i));
            index_as = index(as(j));
            if (index_ao~=0)&&(index_as~=0)
                % sign of the dofs
                soi = sign(etod(i,ie));
                ssj = sign(etod(j,je));
                % i1, i2 and j1, j2
                i_1 = second_node(i);
                i_2 = first_node(i);
                j_1 = second_node(j);
                j_2 = first_node(j);
                % L_p,q = vec(Lp) * vec(Lq)
                Li1_j1 = l_obs(1,i_1)*l_src(1,j_1) + l_obs(2,i_1)*l_src(2,j_1) + l_obs(3,i_1)*l_src(3,j_1);
                Li1_j2 = l_obs(1,i_1)*l_src(1,j_2) + l_obs(2,i_1)*l_src(2,j_2) + l_obs(3,i_1)*l_src(3,j_2);
                Li2_j1 = l_obs(1,i_2)*l_src(1,j_1) + l_obs(2,i_2)*l_src(2,j_1) + l_obs(3,i_2)*l_src(3,j_1);
                Li2_j2 = l_obs(1,i_2)*l_src(1,j_2) + l_obs(2,i_2)*l_src(2,j_2) + l_obs(3,i_2)*l_src(3,j_2);
                % norm of Li, Lj
                L_i    = sqrt(l_obs(1,i)^2 + l_obs(2,i)^2 + l_obs(3,i)^2);
                L_j    = sqrt(l_src(1,j)^2 + l_src(2,j)^2 + l_src(3,j)^2);
                %
                Zi1_j1 = Zi_j(:,:,i_1,j_1);
                Zi1_j2 = Zi_j(:,:,i_1,j_2);
                Zi2_j1 = Zi_j(:,:,i_2,j_1);
                Zi2_j2 = Zi_j(:,:,i_2,j_2);
                %
                ZZ = Li1_j1*Zi2_j2-Li1_j2*Zi2_j1-Li2_j1*Zi1_j2+Li2_j2*Zi1_j1;
                ZZf = A_o*ZZ+4*F_o*ones(Np,Np);
                ZR_ie_je_ij = soi*ssj*L_i*L_j*(sum(sum(W_GR.*ZZf)));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Z_NS(index_ao,index_as) = Z_NS(index(ao(i)),index(as(j)))+ZR_ie_je_ij;
                %                     ZR(index_as,index_ao) = ZR(index(as(j)),index(ao(i)))+ZR_ie_je_ij;
            end %if ((index(ao(i))~=0)&&(index(as(j))~=0))
        end
    end %for i=1:3 for j=1:3
end
% ZR = ZR + ZR.';