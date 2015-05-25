function [Z,F] = Assembly_WIE(WCOIL,freq)
%%    Assembly the MoM system for a thin wire approximation
% _________________________________________________________________________
%
%   based on the BEM approach on Wire Antennas (Harrington)
%       R.F. Harrington
%       Matrix Methods for Field Problems
%       Proc. IEEE 55(2): 136 - 149, Feb. 1967
%
% _________________________________________________________________________
%
%
%% INPUT
%   WCOIL with data
%       Pcoil (Nseg x 3) with coordinates of positive points of segment
%       Ncoil (Nseg x 3) with coordinates of positive points of segment
%       Dwire - diameter of wire
%       Rhowire - resistivity of material
%       port - port definition
%   freq - frequency (Hz)
%
%
%% OUTPUT
%   Z is the impedance matrix
%   F is the rhs for the port definition
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________



% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------

mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
%
omega = 2*pi*freq;
lambda  = co/freq;
ko = 2*pi/lambda;


% -------------------------------------------------------------------------
% Assembly the matrices
% -------------------------------------------------------------------------
%

Pcoil = WCOIL.Pcoil;
Ncoil = WCOIL.Ncoil;
Dwire = WCOIL.Dwire;
Rhowire = WCOIL.Rhocoil;
Nports = length(WCOIL.port);
Nc = size(Pcoil,1); % number of coil segments

fid = 1;
fprintf(fid, '\n ----------------------------------------------------------');
fprintf(fid, '\n Assembling WIE System');
fprintf(fid, '\n');
fprintf(fid, '\n # DOFS:                %d', Nc);
fprintf(fid, '\n # PORTS:               %d', Nports);
fprintf(fid, '\n Operating Frequency:   %.3f MHz',freq/1e6);
fprintf(fid, '\n');


% -------------------------------------------------------------------------
% Generate the MoM Wire system

t1 = tic;

Dcoil = Pcoil - Ncoil; % obtain components of each segment (weight when appleid the current)


Z = MoMWire(Pcoil(:,1),Pcoil(:,2),Pcoil(:,3),Ncoil(:,1),Ncoil(:,2),Ncoil(:,3),Dwire/2,ko,mu,eo,omega);

% -------------------------------------------------------------------------
% with losses
if ~isempty(Rhowire)
    
    deltaw = sqrt(2*Rhowire/(2*pi*freq*mu)); % skin depth
    Rw = Rhowire/(pi*(Dwire-deltaw)*deltaw); % resistance per meter of the wire w/ skin effect
    Z = Z + Rw*diag(sqrt(diag(Dcoil*Dcoil.'))); % add contribution by length of each segment
    
end

% -------------------------------------------------------------------------
% assemble F matrix
F = zeros(Nc,Nports);
for ii = 1:Nports
    F(WCOIL.port(ii),ii) = 1;
end

% -------------------------------------------------------------------------
% And done
% -------------------------------------------------------------------------


fprintf(fid, '\n -------------------------------------');
fprintf(fid, '\n Overall TIME:          %.3f sec', toc(t1));
fprintf(fid, '\n');
fprintf(fid, '\n ----------------------------------------------------------\n ');


