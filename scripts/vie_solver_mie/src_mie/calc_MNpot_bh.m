function [Mo, Me, No, Ne] = calc_MNpot_bh(theta, phi, r, k, n)
% calculate the vector M N potential for scattering while the Bessel
% function is the Hankel function
% --input: theta, phi, r are the target position, they are scalars
%          k is the wavenumber, n is a vector of degrees
% --output: Mo, Me, No, Ne: 3-by-length(n) array
%          the 3 rows are for r, theta, phi directions
rho = k*r;

[dummy, dummy, dummy, dummy, bh, dbh] = calc_sphbessel(n, rho, 0);

[Mo, Me, No, Ne] = calc_MNpot(theta, phi, rho, n, bh, dbh);


