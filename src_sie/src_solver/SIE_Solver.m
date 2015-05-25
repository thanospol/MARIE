function [ZP,Jc] = SIE_Solver(SCOIL,freq)
%%    Solver Routine for the SIE formulation
% _________________________________________________________________________
%
%       Solves the SIE problem to obtain the Z parameters
%       Prepared for structures with gaps as ports
%       Generates the port Z parameters of the system
%       | Ztt    Zit.' | | Jt |   | F |
%       | Zit    Zii   | | Ji | = | 0 | Vin
%
%       Refs:
%           Eleftheriades and Mosig, IEEE TMTT 44(3):438-445 1996
%           Villena, Polimeridis et al, MRM 2014
% _________________________________________________________________________
%
%% INPUT
%       SCOIL structure with
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%           Ct - coordinates of the center of the triangle
%           Ln - values of the length of each side of the triangle
%           Pn - 3x3 matrix with coordinates of the rho vectors (Pn(:,1) == rho_1)
%       freq - frequency
%
%
%% OUTPUT
%       ZP - Z-parameters
%       Jc - current at the ports
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
% Generate the SIE system

[Zo,F] = Assembly_SIE_par(SCOIL,freq);

% -------------------------------------------------------------------------
% Solve for the Z parameters of the system
% -------------------------------------------------------------------------
%
%       The system is
%       | Ztt    Zit.' | | Jt |   | F |
%       | Zit    Zii   | | Ji | = | 0 | Vin
%
%       Needs to form
%       [ZP] = -Fpinv * Ztt * Fpinv.' + Fpinv * Zti * inv(Zii) * Zit * Fpinv.'
%
%       with
%       Fpinv = pinv(F)
%


% -------------------------------------------------------------------------
% split according to the port information

tic_i = tic;
idx = find(sum(abs(F),2));
Nports = size(F,2); % total number of external ports
Tports = length(idx); % number to terminals related to ports
Fp = F(1:Tports,1:Nports);
Fpinv = pinv(Fp); % pseudoinverse of the incidence matrix
Ztt = Zo(1:Tports,1:Tports); % terminal subblock
Zii = Zo(Tports+1:end,Tports+1:end); % internal subblock
Zit = Zo(Tports+1:end,1:Tports);

% Solve for the currents
Zit = Zit*Fpinv.';
xit = Zii\Zit;
Jc = [Fpinv.'; xit]; % terminal currents due to Ip and internal edges currents due to Ip

% Generate the Z parameter matrix
ZP = Fpinv*Ztt*Fpinv.' - Zit.'*xit;
ZP = (ZP + ZP.')/2;

fprintf(1,'\n SIE system solve done,  Elapsed time  = %.2f [sec]' , toc(tic_i));
fprintf(1,'\n ----------------------------------------------------------\n');
