function [Vout] = VSIE_op(Jc,Zcc,fVIEsolver,ZbcL,ZbcR,freq,tau,Gram,Scoord,RHBM,COIL,tol,Jb,idxC)
%%    Operator for the VIE + SIE solver
% _________________________________________________________________________
%
%   Applies the MVP related to the VIE+SIE formulation:
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
%       RHBM - body model struct
%       COIL - coil struct
%       tol - tolerance
%       Jb - initial guess
%       idxC - (optional) indexes of coil current coefficients relevant
%              (used for block cases, when only internal edges are solved)
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


tic_i = tic;

% some definitions:
LEVEL_DVrule = 3; % 5 is more than enough, 3 is ok

% -------------------------------------------------------------------------
% Transform the current Jc into equivalent fields in Scatterer  
% Transform the incident fields into eq. incident volumetric currents
% -------------------------------------------------------------------------

if ~(isempty(idxC)) % only some of the currents are used
    Nd = max(COIL.index); % number of dofs
    J = zeros(Nd,1);
    J(idxC) = Jc;
else
    J = Jc;
end

if isempty(ZbcR)
        % on the fly
    try
        % we used the PARALLEL MEXED Quadrature based function...
        J = SVIE_E_Coil2Scat_QMEX_par(Scoord,COIL.index,COIL.etod,COIL.node,COIL.elem,freq,LEVEL_DVrule,J); % compute the point-based fields
    catch
        % we used the PARALLEL NOT MEXED Quadrature based function...
        J = SVIE_E_Coil2Scat_Q_par(Scoord,COIL.index,COIL.etod,COIL.node,COIL.elem,freq,LEVEL_DVrule,J); % compute the point-based fields
    end
else
    if isempty(ZbcL)
        % single term
        J = ZbcR*J;
    else
        % compressed form (ACA or SVD)
        J = ZbcL*(ZbcR*J);
    end
end

if ~isempty(RHBM.M)
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
    try
        % we used the PARALLEL MEXED Quadrature based function...
        J = SVIE_E_Scat2Coil_QMEX_par(Scoord,COIL.index,COIL.etod,COIL.node,COIL.elem,freq,LEVEL_DVrule,Gram*J); % compute the point-based fields
    catch
        % we used the PARALLEL NOT MEXED Quadrature based function...
        J = SVIE_E_Scat2Coil_Q_par(Scoord,COIL.index,COIL.etod,COIL.node,COIL.elem,freq,LEVEL_DVrule,Gram*J); % compute the point-based fields
    end
    if ~(isempty(idxC)) % only some of the currents are used
         J = J(idxC);% select the relevant indexes
    end
    Vout = Zcc*Jc - J;
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
fprintf(1,'\n Internal VIE+SIE operator applied, Elapsed time  = %g [sec]' , Time);






