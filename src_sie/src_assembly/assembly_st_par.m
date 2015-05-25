function [ZR] = assembly_st_par(index,etod,node,elem,ZR,GL_order,Index_elem,ko) 
%%
% %%%%%%%%%%%%%%% Gauss quadrature rule for singular triangles %%%%%%%%%%%%%%
Np_1D_WS_ST = GL_order.ST;
%%%%%%%%%%%%%%% Gauss quadrature rule for singular triangles %%%%%%%%%%%%%%
[w,z] = gauss_1d_sie(Np_1D_WS_ST);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edge Adjacent Terms %%%%%%%%%%%%%%%%%%%%%%%%%%
ie_ST = Index_elem.ST(:,1); 
% je_NS = Index_elem.NS(:,2); 
n_ST_elem = length(ie_ST);
%


% semivectorized form
Z_ST_vector = zeros(9,n_ST_elem);
AO_index_vector = zeros(9,n_ST_elem);
AS_index_vector = zeros(9,n_ST_elem);

parfor index_ST = 1 : n_ST_elem
    ie = ie_ST(index_ST);   %     je = je_ST(index_ST);
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
    %--------------------------------------------------------
    %                        DEMCEM                         %  
    [I_DE] = mexdirect_ws_st_rwg(ro_1,ro_2,ro_3,ko,Np_1D_WS_ST,w,z);
    %--------------------------------------------------------
    as = ao;
    
    Z_ST_local = zeros(9,1);
    AO_index_local = zeros(9,1);
    AS_index_local = zeros(9,1);
    localindex = 0;
    
    for i=1:3
        for j=1:3
            index_ao = index(ao(i));
            index_as = index(as(j));
            if (index_ao && index_as)
                % sign of the dofs
                soi=sign(etod(i,ie));
                ssj=sign(etod(j,ie));
                %
                ZR_DE = soi*ssj* I_DE(i + 3*(j-1));
                %
                % Z_ST(index_ao,index_as) = Z_ST(index_ao,index_as)+ZR_DE;
                localindex = localindex + 1;
                Z_ST_local(localindex,1) = ZR_DE;
                AO_index_local(localindex,1) = index_ao;
                AS_index_local(localindex,1) = index_as;
                                
            end
        end
    end
    
    Z_ST_vector(:,index_ST) = Z_ST_local;
    AO_index_vector(:,index_ST) = AO_index_local;
    AS_index_vector(:,index_ST) = AS_index_local;
    
end


aoidx = nonzeros(AO_index_vector);
asidx = nonzeros(AS_index_vector);
zvals = nonzeros(Z_ST_vector);

for ii = 1:length(aoidx)
    ZR(aoidx(ii),asidx(ii)) = ZR(aoidx(ii),asidx(ii)) + zvals(ii);
end