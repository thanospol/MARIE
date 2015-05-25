function [ZR] = assembly_ea_par(index,etod,node,elem,ZR,GL_order,Index_elem,ko)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Order of Gauss Quadrature Integration                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_theta = GL_order.EA(1);
N_psi   = GL_order.EA(2);
%
[ w_theta , z_theta ] = gauss_1d_sie ( N_theta );
[ w_psi , z_psi ] = gauss_1d_sie ( N_psi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ie_EA = Index_elem.EA(:,1); 
je_EA = Index_elem.EA(:,2); 
n_EA_elem = length(ie_EA);
%

% semivectorized form
Z_EA_vector = zeros(9,n_EA_elem);
AO_index_vector = zeros(9,n_EA_elem);
AS_index_vector = zeros(9,n_EA_elem);


parfor index_EA = 1 : n_EA_elem
    ie = ie_EA(index_EA);
    je = je_EA(index_EA);
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
    %             Re-ordering for EA triangles
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lia,orderq] = ismember(p,q);               % lia contains 1 in positions of p found in q
                                                % orderq contains the positions of these
    [~,orderp] = sort(lia,'descend');           % reorder so that the 1 appear first
                                                % now orderp has the order in which the common elements
                                                % appear first in the p array
    orderq = orderq(orderp);                    % reorder orderq accordingly
    orderq(1:2) = orderq(2:-1:1);               % swapp the first and second entries of orderq
    orderq(3) = setdiff([1;2;3],orderq(1:2));   % set to 1 the third one (it was always zero, since only two nodes are shared)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                        DEMCEM                         %  
    [I_DE] = mexdirect_ws_ea_rwg(ro(:,orderp(1)) , ro(:,orderp(2)) , ro(:,orderp(3)), rs(:,orderq(3)) , ko, N_theta, N_psi , w_theta, z_theta, w_psi, z_psi);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Z_EA_local = zeros(9,1);
    AO_index_local = zeros(9,1);
    AS_index_local = zeros(9,1);
    localindex = 0;
    for ii=1:3
        i = orderp(ii);
        index_ao = index(ao(i));
        for jj=1:3
            j = orderq(jj);
            %
            index_as = index(as(j));
            if (index_ao && index_as)
                % sign of the dofs
                soi=sign(etod(i,ie));
                ssj=sign(etod(j,je));
                %
                ZR_DE = soi*ssj* I_DE(jj + 3*(ii-1));
                %
                % Z_EA(index_ao,index_as) = Z_EA(index(ao(i)),index(as(j)))+ZR_DE;
                
                localindex = localindex + 1;
                Z_EA_local(localindex,1) = ZR_DE;
                AO_index_local(localindex,1) = index_ao;
                AS_index_local(localindex,1) = index_as;
                
            end
        end
    end
    
    Z_EA_vector(:,index_EA) = Z_EA_local;
    AO_index_vector(:,index_EA) = AO_index_local;
    AS_index_vector(:,index_EA) = AS_index_local;
    
end

aoidx = nonzeros(AO_index_vector);
asidx = nonzeros(AS_index_vector);
zvals = nonzeros(Z_EA_vector);


for ii = 1:length(aoidx)
    ZR(aoidx(ii),asidx(ii)) = ZR(aoidx(ii),asidx(ii)) + zvals(ii);
end