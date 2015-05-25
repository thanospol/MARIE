function [bje, dbje, bji, dbji, bh, dbh] = calc_sphbessel(n, ka, k1a)
% calculate spherical bessel functions and the derivative of it's
% Riccati-bessel forms
% e stands for exterior, corresponding to ka
% i stands for interior, corresponding to k1a
% hankel function only calculates exterior with ka
% input: n : a vector of Bessel orders
%        ka: exterior scalar
%        k1a: interior scalar

ndegree = [n(1)-1, n];

bje = sqrt(pi/2/ka) * besselj(ndegree+0.5, ka);
bji = sqrt(pi/2/k1a) * besselj(ndegree+0.5, k1a);
bh = sqrt(pi/2/ka) * besselh(ndegree+0.5, 2, ka);

dbje = ka*bje(1:end-1) - n.*bje(2:end);
dbji = k1a*bji(1:end-1) - n.*bji(2:end);
dbh = ka*bh(1:end-1) - n.*bh(2:end);

bje = bje(2:end);
bji = bji(2:end);
bh = bh(2:end);



