function [Z_VA] = assembly_va(index,etod,node,elem,ZR,GL_order,Index_elem,ko)
%%
% %%%%%%%%%%%%%%% Gauss quadrature rule for singular triangles %%%%%%%%%%%%
N_theta_p = GL_order.VA(1);
N_theta_q = GL_order.VA(2);
N_psi     = GL_order.VA(3);
%
[ w_theta_p,z_theta_p ] = gauss_1d_sie ( N_theta_p );
[ w_theta_q,z_theta_q ] = gauss_1d_sie ( N_theta_q );
[ w_psi,z_psi ]         = gauss_1d_sie ( N_psi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vertex Adjacent Terms %%%%%%%%%%%%%%%%%%%%%%%%%%
ie_VA = Index_elem.VA(:,1); 
je_VA = Index_elem.VA(:,2); 
n_VA_elem = length(ie_VA);
%
Z_VA= ZR;

for index_VA = 1 : n_VA_elem
    ie = ie_VA(index_VA);
    je = je_VA(index_VA);
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
    %
    ro   = [ro_1 ro_2 ro_3];
    %
    p = [node_test_1  node_test_2  node_test_3];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    as = abs(etod(:,je));
    % coordinates of nodes of the source triangle
    node_basis_1 = elem(1,je);
    node_basis_2 = elem(2,je);
    node_basis_3 = elem(3,je);
    %
    rs_1 = node(:,node_basis_1);
    rs_2 = node(:,node_basis_2);
    rs_3 = node(:,node_basis_3);
    %
    rs   = [rs_1 rs_2 rs_3];
    %
    q = [node_basis_1  node_basis_2  node_basis_3];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Reorder the first 3 values of the elem ie and elem je
    %    so that the commom node appears in first position
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pinq] = ismember(p,q);            % pinq contains 1 in positions of p found in q
    [~,orderp] = sort(pinq,'descend'); % reorder so that the 1 appear first
    [qinp] = ismember(q,p);            % qinp contains 1 in positions of q found in p
    [~,orderq] = sort(qinp,'descend'); % reorder so that the 1 appear first
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           DEMCEM - DIRECT_WS_VA_RWG                   %
    [I_DE] = mexdirect_ws_va_rwg(ro(:,orderp(1)) , ro(:,orderp(2)) , ro(:,orderp(3)) , rs(:,orderq(2)) , rs(:,orderq(3)) , ko, N_theta_p, N_theta_q, N_psi, w_theta_p, z_theta_p, w_theta_q, z_theta_q, w_psi, z_psi);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:3
        i = orderp(ii);
        for jj=1:3
            j = orderq(jj);
            %
            index_ao = index(ao(i));
            index_as = index(as(j));
            if (index_ao~=0)&&(index_as~=0)
                % sign of the dofs
                soi=sign(etod(i,ie));
                ssj=sign(etod(j,je));
                %
                ZR_DE = soi*ssj* I_DE(jj + 3*(ii-1));
                %
                Z_VA(index_ao,index_as) = Z_VA(index(ao(i)),index(as(j)))+ZR_DE;
            end
        end
    end
end