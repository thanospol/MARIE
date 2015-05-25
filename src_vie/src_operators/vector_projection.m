function [P] = vector_projection()
%%
x = [1,0,0];
y = [0,1,0];
z = [0,0,1];
%
KL = [-x;+x;-y;+y;-z;+z]';
%
P = zeros(6,6,6);
for kk=1:6
    k(1,:) = KL(:,kk);
    for ll=1:6
        l(1,:) = KL(:,ll);
        %    
        P(kk,ll,1) = dot( cross(k,x), cross(l,x) );
        P(kk,ll,2) = dot( cross(k,x), cross(l,y) );
        P(kk,ll,3) = dot( cross(k,x), cross(l,z) );
        %
        P(kk,ll,4) = dot( cross(k,y), cross(l,y) );
        P(kk,ll,5) = dot( cross(k,y), cross(l,z) );
        %
        P(kk,ll,6) = dot( cross(k,z), cross(l,z) );
                           
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%