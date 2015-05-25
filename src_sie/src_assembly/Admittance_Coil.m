function [Y_coil] = Admittance_Coil(I,index,node,edge,port)
%%
% first_node  = [3 1 2];
% second_node = [2 3 1];
%
% edge1_p1p = find(index==1);
% edge2_p1p = find(index==2);
%
Ne_port1p = length(port(1).p);

I_coilp = 0;

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
     I_coilp = I_coilp + L_edge * I(index_edge,1); 
     
end

Y_coil = I_coilp;

% %
% Ne_port1n = length(port(1).n);
% 
% I_coiln = 0;
% 
% for ind = 1:Ne_port1n
%      index_edge = port(1).n(ind);
%      real_edge  = find(index==index_edge);
%      %
%      point1 = edge(1,real_edge);
%      point2 = edge(2,real_edge);
%      %
%      rpoint1 = node(:,point1);
%      rpoint2 = node(:,point2);
%      %
%      L_edge = norm(rpoint2-rpoint1);
%      %
%      I_coiln = I_coiln + L_edge * I(index_edge,1); 
%      
% end
% 
% Y_coil = (I_coilp - I_coiln)/2;
