function [Npg,wt,Z,z1,z2,z3] = gauss_2d_sie(Np)


format long


[w,z] = gauss_1d_sie (Np);
ctr=1;
for I=1:Np 
    for J=1:Np
        W = w(I)*w(J);
        x = z(I);
        y = z(J);
        zi = (1-y)/8;
        Wi(I,J) = W*zi;
        Z2(I,J) = (1+y)/2;
        Z3(I,J) = (1-Z2(I,J))*(1+x)/2;
        Z1(I,J) = 1-Z2(I,J)-Z3(I,J);
        
        wt(ctr)=Wi(I,J);
        z1(ctr)=Z1(I,J);
        z2(ctr)=Z2(I,J);
        z3(ctr)=Z3(I,J);
        ctr=ctr+1;
        
    end
end
Ng=ctr-1;
Npg=Ng;

Z=[z1' z2' z3'];
