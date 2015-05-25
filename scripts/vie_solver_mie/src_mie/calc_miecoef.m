function [an, bn, cn, dn] = calc_miecoef(n, ka, k1a, k1overk)
% calculate mie scattering coefficients
% input: n : a vector of Bessel orders
%        ka: exterior scalar
%        k1a: interior scalar
%        k1overk: k1/k


[bje, dbje, bji, dbji, bh, dbh] = calc_sphbessel(n, ka, k1a);

denorm1 = k1overk^2*bji.*dbh - bh.*dbji;
denorm2 = bji.*dbh - bh.*dbji;

an = (k1overk^2*bji.*dbje - bje.*dbji) ./ denorm1;
bn = (bji.*dbje - bje.*dbji) ./ denorm2;
cn = (bje.*dbh - bh.*dbje) ./ denorm2;
dn = k1overk*(bje.*dbh - bh.*dbje) ./ denorm1;