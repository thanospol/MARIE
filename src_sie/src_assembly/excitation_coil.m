function [V] = excitation_coil(B_o,index,node,edge,port)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Excitation Vector                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_node  = [3 1 2];
% second_node = [2 3 1];
%
Nff = max(index);
V  = zeros(Nff,1);
%
Ne_port1p = length(port(1).p);
%
for ind = 1:Ne_port1p
     index_edge = port(1).p(ind);
     real_edge  = find(index==index_edge);
     %
     point1 = edge(1,real_edge);
     point2 = edge(2,real_edge);
     %
     rpoint1 = node(:,point1);
     rpoint2 = node(:,point2);
     %
     L_edge = norm(rpoint2-rpoint1);
     %
     V(index_edge,1) = 0.5*L_edge;    
end
%
Ne_port1n = length(port(1).n);
%
for ind = 1:Ne_port1n
     index_edge = port(1).n(ind);
     real_edge  = find(index==index_edge);
     %
     point1 = edge(1,real_edge);
     point2 = edge(2,real_edge);
     %
     rpoint1 = node(:,point1);
     rpoint2 = node(:,point2);
     %
     L_edge = norm(rpoint2-rpoint1);
     %
     V(index_edge,1) = -0.5*L_edge;    
end

V = B_o*V;