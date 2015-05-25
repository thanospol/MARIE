function [Z] = MoMWire(Px,Py,Pz,Nx,Ny,Nz,a,k,relmu,releps,omega)
%%    Method of Moments for thin wire approximations
% _________________________________________________________________________
%
%   Def.: Formula for the computation of the BEM value between two wire 
%   elements, based on the BEM approach on Wire Antennas (Harrington)
%       R.F. Harrington
%       Matrix Methods for Field Problems
%       Proc. IEEE 55(2): 136 - 149, Feb. 1967
%
% _________________________________________________________________________
%
%
%% INPUT
%   Px is the x coord of positive end of the segments
%   Py is the y coord of positive end of the segments
%   Pz is the z coord of positive end of the segments
%   Nx is the x coord of negative end of the segments
%   Ny is the y coord of negative end of the segments
%   Nz is the z coord of negative end of the segments
%   a is the wire radius
%   k is the wavelength number
%   relmu relative permitivity
%   releps relative dielectric
%   omega is 2*pi*frequency of operation (Hz)
%
%
%% OUTPUT
%   Z is the impedance matrix
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
% Compute the coordinates in matrix form 
% -------------------------------------------------------------------------

% generate matrix coordinates

ns = length(Px);
PMx = zeros(ns,ns); PMy = zeros(ns,ns); PMz = zeros(ns,ns);
NMx = zeros(ns,ns); NMy = zeros(ns,ns); NMz = zeros(ns,ns);

for i = 1:ns
    
    PMx(:,i) = Px(i); % x coordinates of positive edge
    PMy(:,i) = Py(i); % y coordinates of positive edge
    PMz(:,i) = Pz(i); % z coordinates of positive edge
    NMx(:,i) = Nx(i); % x coordinates of negative edge
    NMy(:,i) = Ny(i); % y coordinates of negative edge
    NMz(:,i) = Nz(i); % z coordinates of negative edge
    
end

% compute length of the segments
DM = sqrt((PMx-NMx).^2 + (PMy-NMy).^2 + (PMz-NMz).^2);


% -------------------------------------------------------------------------
% Compute the integral values (Psi) for the elements
% -------------------------------------------------------------------------

% from Matrix Methods for Field Problems, Harrington, Proc. IEEE
% 55(2): 136 - 149, Feb. 1967



% -------------------------------------------------------------------------
% center to center Psi (n,m)

% center is the Center matrix
% positive is the Positive edge matrix
% observation point is the Center matrix

CPx = NMx + (PMx - NMx)./2; CPy = NMy + (PMy - NMy)./2; CPz = NMz + (PMz - NMz)./2;
PPx = PMx; PPy = PMy; PPz = PMz;
OPx = CPx.'; OPy = CPy.'; OPz = CPz.';

% Compute the integral values (Psi) for the elements
Psi_cc = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k);

% fid = fopen('Teststraightwire.txt', 'a');
% for i = 1:length(Psi_cc)
%     for j = 1:length(Psi_cc)
%         fprintf(fid, '\n Psicc(%d,%d) = %f', i,j, Psi_cc(i,j));
%     end
% end
% fprintf(fid, '\n\n');
% fclose(fid);



% Psi_cc

% -------------------------------------------------------------------------
% positive to positive Psi(n+,m+)

% center is the Positive matrix
% positive is the Positive edge matrix plus half of length
% observation point is the Positive matrix

CPx = PMx; CPy = PMy; CPz = PMz;
PPx = CPx + (PMx - NMx)./2; PPy = CPy + (PMy - NMy)./2; PPz = CPz + (PMz - NMz)./2;
OPx = PMx.'; OPy = PMy.'; OPz = PMz.';

% Compute the integral values (Psi) for the elements

Psi_pp = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k);

% fid = fopen('Teststraightwire.txt', 'a');
% for i = 1:length(Psi_pp)
%     for j = 1:length(Psi_pp)
%         fprintf(fid, '\n Psipp(%d,%d) = %f', i,j, Psi_pp(i,j));
%     end
% end
% fprintf(fid, '\n\n');
% fclose(fid);

% Psi_pp

% -------------------------------------------------------------------------
% positive to negative Psi(n+,m-)

% center is the Positive matrix
% positive is the Positive edge matrix plus half of length
% observation point is the Negative matrix

% CPx = PMx; CPy = PMy; CPz = PMz;
% PPx = CPx + (PMx - NMx)./2; PPy = CPy + (PMy - NMy)./2; PPz = NMz + (PMz - NMz)./2;
OPx = NMx.'; OPy = NMy.'; OPz = NMz.';

% Compute the integral values (Psi) for the elements

Psi_pn = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k);

% fid = fopen('Teststraightwire.txt', 'a');
% for i = 1:length(Psi_pn)
%     for j = 1:length(Psi_pn)
%         fprintf(fid, '\n Psipn(%d,%d) = %f', i,j, Psi_pn(i,j));
%     end
% end
% fprintf(fid, '\n\n');
% fclose(fid);

% Psi_pn

% -------------------------------------------------------------------------
% negative to positive Psi(n-,m+)


% center is the Negative matrix
% positive is the Center point
% observation point is the Positive matrix

CPx = NMx; CPy = NMy; CPz = NMz;
PPx = CPx + (PMx - NMx)./2; PPy = CPy + (PMy - NMy)./2; PPz = NMz + (PMz - NMz)./2;
OPx = PMx.'; OPy = PMy.'; OPz = PMz.';

% Compute the integral values (Psi) for the elements

Psi_np = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k);

% fid = fopen('Teststraightwire.txt', 'a');
% for i = 1:length(Psi_np)
%     for j = 1:length(Psi_np)
%         fprintf(fid, '\n Psinp(%d,%d) = %f', i,j, Psi_np(i,j));
%     end
% end
% fprintf(fid, '\n\n');
% fclose(fid);

% Psi_np


% -------------------------------------------------------------------------
% negative to negative Psi(n-,m-)


% center is the Negative matrix
% positive is the Center point
% observation point is the Negative matrix

% CPx = NMx; CPy = NMy; CPz = NMz;
% PPx = CPx + (PMx - NMx)./2; PPy = CPy + (PMy - NMy)./2; PPz = NMz + (PMz - NMz)./2;
OPx = NMx.'; OPy = NMy.'; OPz = NMz.';

% Compute the integral values (Psi) for the elements

Psi_nn = PsiW_op(PPx,PPy,PPz,CPx,CPy,CPz,OPx,OPy,OPz,DM,a,k);

% fid = fopen('Teststraightwire.txt', 'a');
% for i = 1:length(Psi_nn)
%     for j = 1:length(Psi_nn)
%         fprintf(fid, '\n Psinn(%d,%d) = %f', i,j, Psi_nn(i,j));
%     end
% end
% fprintf(fid, '\n\n');
% fclose(fid);

% Psi_nn


% -------------------------------------------------------------------------
% Compute the value of the impedance Zmn
% -------------------------------------------------------------------------

DLproduct = (PMx-NMx).*((PMx-NMx).') + (PMy-NMy).*((PMy-NMy).') + (PMz-NMz).*((PMz-NMz).');

Z = (1j * omega * relmu *(DLproduct)).* Psi_cc + (Psi_pp - Psi_np - Psi_pn + Psi_nn) ./ (1j * omega * releps);

