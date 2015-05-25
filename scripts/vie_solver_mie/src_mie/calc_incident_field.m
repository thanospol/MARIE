% calculate incident field to compare the MN expansion and the direct plane
% wave
clear all
clc

a = 1;
k = 5;
k1 = 20;
omega = 2*pi*1e10;
mu = 4*pi*1e-7;
mu1 = mu;
E0 = 1;
number_of_target = 100;

temp_theta = linspace(0, pi, number_of_target/2);
theta = [temp_theta, temp_theta(end:-1:1)];
phi = [0*ones(1, number_of_target/2), (pi)*ones(1, number_of_target/2)];
r = 5*ones(1, number_of_target);


theta_plot = [temp_theta, temp_theta+pi];
% direct incident field

temp_E = ones(3, 1) * (E0*exp(-1i*k*r.*cos(theta)));
Einc_direct = temp_E .* [sin(theta).*cos(phi); cos(theta).*cos(phi); -sin(phi)];
Hinc_direct = k/(omega*mu) * temp_E .* [sin(theta).*sin(phi); cos(theta).*sin(phi); cos(phi)];

% incident field by expansion

Einc_expansion = zeros(3, number_of_target);
Hinc_expansion = zeros(3, number_of_target);

n = 1:60;
En = E0*(-1i).^n .* (2*n+1) ./(n.*(n+1));

for i_target = 1 : number_of_target

    [Mo, Me, No, Ne] = calc_MNpot_bj(theta(i_target), phi(i_target), r(i_target), k, n);

    temp_Einc = Mo + 1i*Ne;
    Einc_expansion(:, i_target) = (En * temp_Einc.').';
    
    temp_Hinc = -k/(omega*mu) * (Me - 1i*No);
    Hinc_expansion(:, i_target) = (En * temp_Hinc.').';
    
end % for i_target


Esph_title = {'r', '$\theta$', '$\phi$'};
Hsph_title = {'r', '$\theta$', '$\phi$'};

for i_vector = 1:3
  
   figure(1)
   subplot(3, 2, (i_vector-1)*2+1);
   plot(theta_plot, real(Einc_direct(i_vector, :) - Einc_expansion(i_vector, :)), 'g.--');
   title([Esph_title{i_vector},' real']);
   
   
   subplot(3, 2, (i_vector-1)*2+2);
   plot(theta_plot, imag(Einc_direct(i_vector, :) - Einc_expansion(i_vector, :)), 'g.--'); 
   title([Esph_title{i_vector},' imag']);
   
   figure(2)
   subplot(3, 2, (i_vector-1)*2+1);
   plot(theta_plot, real(Hinc_direct(i_vector, :) - Hinc_expansion(i_vector, :)), 'g.--');
   title([Hsph_title{i_vector},' real']);
   
   
   subplot(3, 2, (i_vector-1)*2+2);
   plot(theta_plot, imag(Hinc_direct(i_vector, :) - Hinc_expansion(i_vector, :)), 'g.--'); 
   title([Hsph_title{i_vector},' imag']);
   
end