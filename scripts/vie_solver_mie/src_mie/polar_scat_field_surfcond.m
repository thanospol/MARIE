function [Es, Hs] = polar_scat_field_surfcond(theta, phi, r, k, k1, omega, mu, mu1, a, sigma, E_inc)

%if(r<a); error('the target point should be outside the shpere');end

number_of_target = length(theta);

% scatter field
Es = zeros(3, number_of_target);
Hs = zeros(3, number_of_target);

% n = 1:40;
% 
% [an, bn, cn, dn] = calc_miecoef_surfcond(n, k, k1, a, sigma, omega, mu, mu1);



an=[];bn=[];
n = [];
temp_n = 1:40;

while 1
    [temp_an, temp_bn] = calc_miecoef_surfcond(temp_n, k, k1, a, sigma, omega, mu, mu1);
    an = [an temp_an]; 
    bn = [bn temp_bn];
    n = [n temp_n];
    
    nn=n(end);

    if (nn>100)
        break;
    end
    
    if(abs(an(1))<1e-6 || abs(bn(1))<1e-6)
        break;
    end
    
    if (abs(an(end)/an(1)*(2*nn+1)/(nn^2+nn))<1e-4 && abs(bn(end)/bn(1)*(2*nn+1)/(nn^2+nn))<1e-4)
        break;
    else
        temp_n=temp_n(end)+1:temp_n(end)+2*length(temp_n);
    end

end


an = ones(3, 1) * an;
bn = ones(3, 1) * bn;

En = E_inc*(-i).^n .* (2*n+1) ./(n.*(n+1));

for i_target = 1 : number_of_target

    [Mo, Me, No, Ne] = calc_MNpot_bh(theta(i_target), phi(i_target), r(i_target), k, n);

    temp_Es = -i*an.*Ne - bn.*Mo;
    Es(:, i_target) = (En * temp_Es.').';
    
    temp_Hs = k/(omega*mu) * (-i*bn.*No + an.*Me);
    Hs(:, i_target) = (En * temp_Hs.').';
    
end % for i_target

