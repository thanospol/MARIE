%%    Script to generate free space EM constants from frequency
% _________________________________________________________________________
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
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
eta =  3.767303134617706e+002; % Free-space impedance