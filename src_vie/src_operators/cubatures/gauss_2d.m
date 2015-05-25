function [Npg,w,u,v] = gauss_2d(Np)

[w1D,x1D] = gauss_1d (Np);
Npg = Np * Np;
w  = zeros(Npg,1);
u  = zeros(Npg,1);
v  = zeros(Npg,1);

ctr=1;

for I=1:Np
    for J=1:Np
        w (ctr,1)= w1D(I)*w1D(J);
        u(ctr,1) = (x1D(I) + 1) / 2;
        v(ctr,1) = (x1D(J) + 1) / 2;
        ctr=ctr+1;
    end
end
w = (1/4) * w ;
