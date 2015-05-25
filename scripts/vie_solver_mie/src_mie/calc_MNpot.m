function [Mo, Me, No, Ne] = calc_MNpot(theta, phi, rho, n, bj, dbj)
% calculate the vector M N potential with the Bessel functions provided
% --input: theta, phi, rho=kr are the target position, they are scalars
%          k is the wavenumber, n is a vector of degrees
% --output: Mo, Me, No, Ne: 3-by-length(n) array
%          the 3 rows are for r, theta, phi directions.

pin = zeros(1, length(n));
taun = zeros(1, length(n));

% avoid singularity
if theta == 0; theta = theta+1e-6; end
if theta == pi; theta = theta-1e-6; end

for i = 1:length(n)
    
    % the definition of legendre function of Bohren's book is different as
    % Matlab's definition. Bohren didn't use Condon-Shortley phase (-1)^m
    temp = -legendre(n(i), cos(theta));
    
    pin(i) = temp(2)/sin(theta);

    taun(i) =-1/sin(theta)*(n(i)*(n(i)+1)*abs(sin(theta))*temp(1) + cos(theta)*temp(2));
    
end

% initialize
Mo = zeros(3, length(n));
Me = zeros(3, length(n));
No = zeros(3, length(n));
Ne = zeros(3, length(n));

Mo(1,:) = 0;
Mo(2,:) = cos(phi)*pin.*bj;
Mo(3,:) = -sin(phi)*taun.*bj;

Me(1,:) = 0;
Me(2,:) = -sin(phi)*pin.*bj;
Me(3,:) = -cos(phi)*taun.*bj;

No(1,:) = sin(phi)*n.*(n+1).*sin(theta).*pin.*bj/rho;
No(2,:) = sin(phi)*taun.*dbj/rho;
No(3,:) = cos(phi)*pin.*dbj/rho;

Ne(1,:) = cos(phi)*n.*(n+1).*sin(theta).*pin.*bj/rho;
Ne(2,:) = cos(phi)*taun.*dbj/rho;
Ne(3,:) = -sin(phi)*pin.*dbj/rho;
