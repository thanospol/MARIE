function [Et, Ht] = polar_transmitted_field_surfcond(theta, phi, r, k, k1, omega, mu, mu1, a, sigma, E_inc)

if(r>a); error('the target point should be inside the shpere');end

number_of_target = length(theta);

% scatter field
Et = zeros(3, number_of_target);
Ht = zeros(3, number_of_target);

% n = 1:40;
% 
% [an, bn, cn, dn] = calc_miecoef_surfcond(n, k, k1, a, sigma, omega, mu, mu1);


cn=[];dn=[];
n = [];
temp_n = 1:40;

while 1
    [an, bn, temp_cn, temp_dn] = calc_miecoef_surfcond(temp_n, k, k1, a, sigma, omega, mu, mu1);
    cn = [cn temp_cn]; 
    dn = [dn temp_dn];
    n = [n temp_n];
    
    nn=n(end);
    
    if (nn>100)
        break;
    end
    
    if(abs(cn(1))<1e-6 || abs(dn(1))<1e-6)
        break;
    end
    
    if (abs(cn(end)/cn(1)*(2*nn+1)/(nn^2+nn))<1e-4 && abs(dn(end)/dn(1)*(2*nn+1)/(nn^2+nn))<1e-4)
        break;
    else
        temp_n=temp_n(end)+1:temp_n(end)+2*length(temp_n);
    end

end





cn = ones(3, 1) * cn;
dn = ones(3, 1) * dn;

En = E_inc*(-i).^n .* (2*n+1) ./(n.*(n+1));

for i_target = 1 : number_of_target

    [Mo, Me, No, Ne] = calc_MNpot_bj(theta(i_target), phi(i_target), r(i_target), k1, n);

    temp_Et = cn.*Mo + i*dn.*Ne;
    Et(:, i_target) = (En * temp_Et.').';
    
    temp_Ht = -k1/(omega*mu1) * (dn.*Me - i*cn.*No);
    Ht(:, i_target) = (En * temp_Ht.').';
    
end % for i_target

