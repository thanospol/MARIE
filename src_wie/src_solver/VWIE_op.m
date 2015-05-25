function [Vout] = VWIE_op(Jc,Zcc,fVIEsolver,ZbcL,ZbcR,freq,tau,Gram,Scoord,Ccoil,Dcoil,RHBM,tol,Jb,idxC)
%%    Operator for the VIE + WIE solver
% _________________________________________________________________________
%
%   Applies the MVP related to the VIE+WIE formulation:
%       f(Ic) = (Zcc - Zbc.'*inv(Zbb)*Zbc)*Ic
%
%       related to
%       | Zcc    Zbc.' | | Jc |   | F |
%       |              | |    | = |   | Vin
%       | Zbc    Zbb   | | Jb |   | 0 |
%
%
%   Uses a two step iterative scheme
%       Vrhs = tau*Gram*Zbc*Jc      this is the rhs for the VIE solver
%                                   Zbc*Jc gives point fields
%                                   need to integrate over the volume (Gram)
%       Jb = inv(Zbb)*Vrhs          this is solved iteratively
%       Zp = Zcc - Zbc.'*(Gram*Jb)  integrate the vol. currets, radiate to 
%                                   point fields and add contribution to Zcc
%
%
%       Refs:
%           Villena, Polimeridis et al, MRM 2014
%
% _________________________________________________________________________
%
%% INPUT
%       Jc - coil current 
%       Zcc - coil impedance matrix
%       fVIEsolver - function for solving the VIE solver
%       ZbcL - left factor of Zbc
%       ZbcR - right factor of Zbc
%       fG - circulant for the VIE
%       freq - frequency
%       tau - multiplier to convert fields to currents
%       Gram - gram matrix (volume to the voxel)
%       Scoord - coordinates of the scatterer (for on the fly coupling)
%       Ccoil - coordinates of the coil segments (for on the fly coupling)
%       Dcoil - length of coil segments (for on the fly coupling)
%       RHBM - body model struct
%       tol - tolerance
%       idxC - (optional) indexes of coil current coefficients relevant
%              (used for block cases, when only internal edges are solved)
%       Jb - initial guess
%
%
%% OUTPUT
%       Vout - solution of the MVP 
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________
   

% fast return if input is zero
if (nnz(Jc) == 0)
    Vout = Jc;
    return
end    


% -------------------------------------------------------------------------
%                 define EM vars
% -------------------------------------------------------------------------

tic_i = tic;

% -------------------------------------------------------------------------
% Transform the current Jc into equivalent fields in Scatterer  
% Transform the incident fields into eq. incident volumetric currents
% -------------------------------------------------------------------------

if ~(isempty(idxC)) % only some of the currents are used
    Nd = size(Ccoil,1); % number of dofs
    Jin = zeros(Nd,1);
    Jin(idxC) = Jc;
else
    Jin = Jc;
end

if isempty(ZbcR)
    % on the fly
    [J] = WVIE_E_Coil2Scat_PM(Scoord,Ccoil,Dcoil,freq,Jin); % compute the point-based fields
else
    if isempty(ZbcL)
        % single term
        J = ZbcR*Jin;
    else
        % compressed form (ACA or SVD)
        J = ZbcL*(ZbcR*Jin);
    end
end
clear Jin;

if ~(isempty(RHBM.M))
    % generate initial guess based on MRGFs
    Jb = RHBM.M*(RHBM.X*(RHBM.P.'*J)); % project Einc (J) on DEIM points (P and X), and scale basis
end

% transform into VIE input format: transform to current and integrate
J = tau.*(Gram*J);


% -------------------------------------------------------------------------
%  Solve to obtain equivalent VIE solution volumetric currents
% -------------------------------------------------------------------------

% BICGSTAB no precond, JVIE 2 formulation
J = fVIEsolver(J, tol, Jb);

% -------------------------------------------------------------------------
%  Transform into scattered fields and add direct contribution
% -------------------------------------------------------------------------

% to radiate we have to integrate the coefficients J over volume

if isempty(ZbcR)
    % on the fly
    [J] = WVIE_E_Scat2Coil_PM(Scoord,Ccoil,Dcoil,freq,Gram*J);
    if ~(isempty(idxC)) % only some of the currents are used
        J = J(idxC); % select the relevant indexes
    end
    Vout = Zcc*Jc - J; % scale by the length of the wire
else
    if isempty(ZbcL)
        % single term
        Vout = Zcc*Jc - ZbcR.'*(Gram*J);
    else
        % compressed form (ACA or SVD)
        Vout = Zcc*Jc - ZbcR.'*(ZbcL.'*(Gram*J));
    end
end


% -------------------------------------------------------------------------
%  Report
% -------------------------------------------------------------------------

Time= toc(tic_i);
fprintf(1,'\n Internal VIE+WIE operator applied, Elapsed time  = %g [sec]' , Time);





