
function [Npg,wt,Z1,Z2,Z3] = Gauss_2D_CIE(Np)


format long


[w,z] = Gauss_1D_CIE (Np);

Npg = Np * Np;

Z1 = zeros(Npg,1);
Z2 = zeros(Npg,1);
Z3 = zeros(Npg,1);
wt = zeros(Npg,1);
ctr=1;
for I=1:Np 
    for J=1:Np
        W = w(I)*w(J);
        x = z(I);
        y = z(J);
        zi = (1-y)/8;
%         Wi(I,J) = W*zi;
%         Z2_temp(I,J) = (1+y)/2;
%         Z3_temp(I,J) = (1-Z2_temp(I,J))*(1+x)/2;
%         Z1_temp(I,J) = 1-Z2_temp(I,J)-Z3_temp(I,J);
        
        wt(ctr,1) = W*zi;
        Z2(ctr,1) = (1+y)/2;
        Z3(ctr,1) = (1-Z2(ctr,1))*(1+x)/2;
        Z1(ctr,1) = 1-Z2(ctr,1)-Z3(ctr,1);
        ctr=ctr+1;
        
    end
end
% Ng=ctr-1;
% Npg=Ng;

% Z1=z1' ;
% Z2=z2' ;
% Z3=z3' ;
