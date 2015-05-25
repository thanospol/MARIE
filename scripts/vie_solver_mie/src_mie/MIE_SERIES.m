function [E_Mie,H_Mie] = MIE_SERIES(r,Rad,Er,Se,freq)
%%   MIE series 
% _________________________________________________________________________
%
%   Generates E and H fields for a Sphere via MIE series
% _________________________________________________________________________
%
%
%
% -------------------------------------------------------------------------
%
%   A.G. Polimeridis -- thanos_p@mit.edu
%   J. Fernandez Villena -- jvillena@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*freq;
lambda = co/freq;
ko = 2*pi/lambda;
omega_mu = omega*mu;

% amplitude of incident field
E_inc = 1;

[L,M,N,~] = size(r);
eps = Er - 1j*Se/(eo*omega);
E_Mie = zeros(L,M,N,3);
H_Mie = zeros(L,M,N,3);
parfor ii = 1:L
    for jj = 1:M
        for kk = 1:N
            [phi,theta,rMie] = cart2sph(r(ii,jj,kk,1),r(ii,jj,kk,2),r(ii,jj,kk,3));
            theta = pi/2 - theta;
            [E_Mie_temp,H_Mie_temp] = polar_transmitted_field(theta, phi, rMie, ko, ko*sqrt(eps), omega_mu, Rad, E_inc);
            E_Mie(ii,jj,kk,:) =  Spher2Cartes_field(theta, phi)*E_Mie_temp;
            H_Mie(ii,jj,kk,:) =  Spher2Cartes_field(theta, phi)*H_Mie_temp;
        end
    end
end