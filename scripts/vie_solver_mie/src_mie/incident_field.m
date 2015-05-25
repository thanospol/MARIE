
function [Einc, Hinc] = incident_field(r, theta, phi, k, omega, mu, E0)
% to calculate incident field through the expansion


n = 1:40;
En = E0*(-1i).^n .* (2*n+1) ./(n.*(n+1));

[Mo, Me, No, Ne] = calc_MNpot_bj(theta, phi, r, k, n);

temp_Einc = Mo + 1i*Ne;
Einc = (En * temp_Einc.').';
    
temp_Hinc = -k/(omega*mu) * (Me - 1i*No);
Hinc = (En * temp_Hinc.').';