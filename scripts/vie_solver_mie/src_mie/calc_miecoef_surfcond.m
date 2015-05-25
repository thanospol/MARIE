function [an, bn, cn, dn] = calc_miecoef_surfcond(n, k, k1, a, sigma, omega, mu, mu1)
% calculate mie scattering coefficients
% input: n : a vector of Bessel orders
%        ka: exterior scalar
%        k1a: interior scalar
%        k1overk: k1/k


[bje, dbje, bji, dbji, bh, dbh] = calc_sphbessel(n, k*a, k1*a);

denorm1 = (k1/k)^2*mu*dbh.*bji - mu1*dbji.*bh - i*(k1/k)*sigma*omega*mu*mu1/(k1*k*a)*dbji.*dbh;
denorm2 = mu*dbji.*bh - mu1*bji.*dbh + i*sigma*omega*a*mu*mu1*bji.*bh;

an = ((k1/k)^2*mu*bji.*dbje - mu1*bje.*dbji - i*(k1/k)*sigma*omega*mu*mu1/(k1*k*a)*dbji.*dbje) ./ denorm1;
bn = (mu*dbji.*bje - mu1*dbje.*bji + i*sigma*omega*a*mu*mu1*bji.*bje) ./ denorm2;
cn = mu1*(dbje.*bh - bje.*dbh) ./ denorm2;
dn = mu1*k1/k*(bje.*dbh - dbje.*bh) ./ denorm1;